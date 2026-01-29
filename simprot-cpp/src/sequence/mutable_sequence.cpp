#include "simprot/sequence/mutable_sequence.hpp"

#include <utility>

namespace simprot {

MutableSequence::MutableSequence(const std::string& sequence,
                                 const std::vector<double>* rates,
                                 WichmannHillRNG& rng,
                                 double gamma_alpha) {
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        double rate = (rates && i < rates->size())
            ? (*rates)[i]
            : (gamma_alpha < 0.0 ? 1.0 : rng.gamma(gamma_alpha));
        append(sequence[i], rate);
    }
}

MutableSequence::MutableSequence(const std::string& sequence,
                                 const std::vector<double>& rates) {
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        double rate = (i < rates.size()) ? rates[i] : 1.0;
        append(sequence[i], rate);
    }
}

MutableSequence::MutableSequence(MutableSequence&& other) noexcept
    : head_(other.head_)
    , tail_(other.tail_)
    , size_(other.size_) {
    other.head_ = nullptr;
    other.tail_ = nullptr;
    other.size_ = 0;
}

MutableSequence& MutableSequence::operator=(MutableSequence&& other) noexcept {
    if (this != &other) {
        clear();
        head_ = other.head_;
        tail_ = other.tail_;
        size_ = other.size_;
        other.head_ = nullptr;
        other.tail_ = nullptr;
        other.size_ = 0;
    }
    return *this;
}

MutableSequence::~MutableSequence() {
    clear();
}

void MutableSequence::clear() {
    SequenceNode* current = head_;
    while (current) {
        SequenceNode* next = current->next;
        delete current;
        current = next;
    }
    head_ = nullptr;
    tail_ = nullptr;
    size_ = 0;
}

void MutableSequence::append(char residue, double rate) {
    auto* node = new SequenceNode(residue, rate);

    if (!head_) {
        head_ = tail_ = node;
    } else {
        node->prev = tail_;
        tail_->next = node;
        tail_ = node;
    }
    ++size_;
}

std::string MutableSequence::to_string() const {
    std::string result;
    result.reserve(size_);
    for (SequenceNode* node = head_; node; node = node->next) {
        result += node->residue;
    }
    return result;
}

std::vector<double> MutableSequence::rates() const {
    std::vector<double> result;
    result.reserve(size_);
    for (SequenceNode* node = head_; node; node = node->next) {
        result.push_back(node->rate);
    }
    return result;
}

void MutableSequence::normalize_rates() {
    if (size_ == 0) return;

    // Compute sum of rates
    double sum = 0.0;
    for (SequenceNode* node = head_; node; node = node->next) {
        sum += node->rate;
    }

    if (sum <= 0.0) return;

    // Normalize: rate_i = rate_i / sum * size
    // This makes the average rate = 1.0
    double factor = static_cast<double>(size_) / sum;
    for (SequenceNode* node = head_; node; node = node->next) {
        node->rate *= factor;
    }
}

SequenceNode* MutableSequence::insert_before(SequenceNode* before,
                                             std::size_t length,
                                             WichmannHillRNG& rng,
                                             const SubstitutionMatrix& matrix,
                                             double gamma_alpha) {
    if (length == 0) return before;

    // COMPATIBILITY: The original SIMPROT generates insertions in two phases:
    // 1. RandomSequence(size) - generates ALL amino acids first
    // 2. MakeSequenceList(seq, NULL) - generates ALL rates second
    // We must match this order to maintain RNG synchronization.

    // Phase 1: Generate all amino acids
    std::vector<char> residues;
    residues.reserve(length);
    for (std::size_t i = 0; i < length; ++i) {
        AminoAcidIndex aa = matrix.sample_from_frequencies(rng);
        residues.push_back(amino_acid_to_char(aa));
    }

    // Phase 2: Generate all rates and create nodes
    SequenceNode* insert_head = nullptr;
    SequenceNode* insert_tail = nullptr;

    for (std::size_t i = 0; i < length; ++i) {
        double rate = (gamma_alpha < 0.0) ? 1.0 : rng.gamma(gamma_alpha);

        auto* node = new SequenceNode(residues[i], rate);
        if (!insert_head) {
            insert_head = insert_tail = node;
        } else {
            node->prev = insert_tail;
            insert_tail->next = node;
            insert_tail = node;
        }
    }

    // Link the insertion into the main list
    if (before) {
        // Insert before the given node
        if (before->prev) {
            // Middle insertion
            before->prev->next = insert_head;
            insert_head->prev = before->prev;
        } else {
            // Prefix insertion (before is head)
            head_ = insert_head;
        }
        insert_tail->next = before;
        before->prev = insert_tail;
    } else {
        // Append at tail (before == nullptr)
        if (tail_) {
            tail_->next = insert_head;
            insert_head->prev = tail_;
            tail_ = insert_tail;
        } else {
            // Empty list
            head_ = insert_head;
            tail_ = insert_tail;
        }
    }

    size_ += length;
    return insert_head;
}

SequenceNode* MutableSequence::delete_at(SequenceNode* start, std::size_t count) {
    if (!start || count == 0) return start;

    bool deleting_head = (start == head_);
    SequenceNode* before_deletion = start->prev;

    // Delete nodes
    SequenceNode* current = start;
    std::size_t deleted = 0;
    while (current && deleted < count) {
        SequenceNode* next = current->next;
        if (next) {
            next->prev = current->prev;
        }
        delete current;
        current = next;
        ++deleted;
    }

    size_ -= deleted;

    // Update pointers
    if (size_ == 0) {
        head_ = tail_ = nullptr;
        return nullptr;
    }

    if (deleting_head) {
        head_ = current;
        if (head_) {
            head_->prev = nullptr;
        }
    } else if (before_deletion) {
        before_deletion->next = current;
    }

    // Update tail if we deleted to the end
    if (!current) {
        tail_ = before_deletion;
        if (tail_) {
            tail_->next = nullptr;
        }
    } else if (!current->next) {
        tail_ = current;
    }

    return current;
}

SequenceNode* MutableSequence::node_at(std::size_t index) const {
    if (index >= size_) return nullptr;

    SequenceNode* node = head_;
    for (std::size_t i = 0; i < index && node; ++i) {
        node = node->next;
    }
    return node;
}

void MutableSequence::clear_marks() {
    for (SequenceNode* node = head_; node; node = node->next) {
        node->mark = 0;
    }
}

void MutableSequence::for_each(const std::function<void(SequenceNode&)>& callback) {
    for (SequenceNode* node = head_; node; node = node->next) {
        callback(*node);
    }
}

void MutableSequence::for_each(const std::function<void(const SequenceNode&)>& callback) const {
    for (SequenceNode* node = head_; node; node = node->next) {
        callback(*node);
    }
}

} // namespace simprot
