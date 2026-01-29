#include "simprot/sequence/alignment.hpp"

#include <algorithm>
#include <unordered_set>

namespace simprot {

//==============================================================================
// Static factory method
//==============================================================================

Alignment Alignment::from_tree(const TreeNode& root) {
    return Alignment(root);
}

//==============================================================================
// Recursive constructor (GetAlignment equivalent)
//==============================================================================

Alignment::Alignment(const TreeNode& node) {
    // If the node is a leaf, create a single-sequence alignment
    if (node.is_leaf()) {
        if (!node.sequence.empty()) {
            sequences_.emplace_back(node.name, node.sequence);
        }
        return;
    }

    // Internal node: recursively build alignments for subtrees
    // and merge them using the profiles

    // Get alignment for left subtree and align with current node's profile
    Alignment left_aln(node.left ? Alignment(*node.left) : Alignment());
    Alignment left_profile(left_aln, node, ChildDirection::Left);

    // Get alignment for right subtree and align with current node's profile
    Alignment right_aln(node.right ? Alignment(*node.right) : Alignment());
    Alignment right_profile(right_aln, node, ChildDirection::Right);

    // Merge the two profile-aligned alignments
    Alignment merged(left_profile, right_profile);
    sequences_ = std::move(merged.sequences_);
}

//==============================================================================
// AlignProfile equivalent
//==============================================================================

Alignment::Alignment(const Alignment& child_aln,
                     const TreeNode& parent,
                     ChildDirection direction) {
    if (child_aln.empty()) {
        return;
    }

    // Get the appropriate profiles
    const std::string& profile_child = (direction == ChildDirection::Left)
        ? parent.left_profile_child
        : parent.right_profile_child;
    const std::string& profile_self = (direction == ChildDirection::Left)
        ? parent.left_profile_self
        : parent.right_profile_self;

    // If profiles are empty, just return the child alignment with parent prepended
    if (profile_child.empty() || profile_self.empty()) {
        sequences_.reserve(child_aln.sequences_.size() + 1);
        sequences_.emplace_back(parent.name, parent.sequence);
        for (const auto& seq : child_aln.sequences_) {
            sequences_.push_back(seq);
        }
        return;
    }

    // Create new alignment: parent + all child sequences
    const std::size_t num_child_seqs = child_aln.sequences_.size();
    sequences_.resize(num_child_seqs + 1);

    // First row is the parent
    sequences_[0].name = parent.name;

    // Copy child sequence names
    for (std::size_t i = 0; i < num_child_seqs; ++i) {
        sequences_[i + 1].name = child_aln.sequences_[i].name;
    }

    // Iterate over positions, aligning using profiles
    const std::size_t profile_len = profile_child.length();
    const std::size_t child_aln_len = child_aln.sequences_[0].sequence.length();

    std::size_t j = 0;  // Position in child alignment
    std::size_t k = 0;  // Position in profile

    while (j < child_aln_len && k < profile_len) {
        char child_char = child_aln.sequences_[0].sequence[j];
        char profile_char = profile_child[k];

        if (child_char == profile_char) {
            // Characters match: add profile_self to parent, child chars to children
            sequences_[0].sequence += profile_self[k];
            for (std::size_t i = 0; i < num_child_seqs; ++i) {
                sequences_[i + 1].sequence += child_aln.sequences_[i].sequence[j];
            }
            ++j;
            ++k;
        } else if (child_char == kGapChar) {
            // Child alignment has a gap: add gap to parent, propagate child gaps
            sequences_[0].sequence += kGapChar;
            for (std::size_t i = 0; i < num_child_seqs; ++i) {
                sequences_[i + 1].sequence += child_aln.sequences_[i].sequence[j];
            }
            ++j;
        } else {
            // Profile has extra character: add to parent, gaps to children
            sequences_[0].sequence += profile_self[k];
            for (std::size_t i = 0; i < num_child_seqs; ++i) {
                sequences_[i + 1].sequence += kGapChar;
            }
            ++k;
        }
    }

    // Handle remaining profile characters
    while (k < profile_len) {
        sequences_[0].sequence += profile_self[k];
        for (std::size_t i = 0; i < num_child_seqs; ++i) {
            sequences_[i + 1].sequence += kGapChar;
        }
        ++k;
    }

    // Handle remaining child alignment characters
    while (j < child_aln_len) {
        sequences_[0].sequence += kGapChar;
        for (std::size_t i = 0; i < num_child_seqs; ++i) {
            sequences_[i + 1].sequence += child_aln.sequences_[i].sequence[j];
        }
        ++j;
    }
}

//==============================================================================
// AlignAlignments equivalent
//==============================================================================

Alignment::Alignment(const Alignment& a, const Alignment& b) {
    if (a.empty()) {
        sequences_ = b.sequences_;
        return;
    }
    if (b.empty()) {
        sequences_ = a.sequences_;
        return;
    }

    // New alignment has all sequences from a, plus all but first from b
    // (first row of b is same as first row of a - the parent sequence)
    const std::size_t a_size = a.sequences_.size();
    const std::size_t b_size = b.sequences_.size();
    const std::size_t total_size = a_size + b_size - 1;

    sequences_.resize(total_size);

    // Copy names from a
    for (std::size_t i = 0; i < a_size; ++i) {
        sequences_[i].name = a.sequences_[i].name;
    }

    // Copy names from b (skipping first which is duplicate)
    for (std::size_t i = 1; i < b_size; ++i) {
        sequences_[a_size + i - 1].name = b.sequences_[i].name;
    }

    // Merge sequences using first rows as profiles
    const std::size_t a_len = a.sequences_[0].sequence.length();
    const std::size_t b_len = b.sequences_[0].sequence.length();

    std::size_t j = 0;  // Position in a
    std::size_t k = 0;  // Position in b

    while (j < a_len && k < b_len) {
        char a_char = a.sequences_[0].sequence[j];
        char b_char = b.sequences_[0].sequence[k];

        if (a_char != kGapChar && a_char == b_char) {
            // Same non-gap character: take from both
            for (std::size_t i = 0; i < a_size; ++i) {
                sequences_[i].sequence += a.sequences_[i].sequence[j];
            }
            for (std::size_t i = 1; i < b_size; ++i) {
                sequences_[a_size + i - 1].sequence += b.sequences_[i].sequence[k];
            }
            ++j;
            ++k;
        } else if (a_char == kGapChar) {
            // a has gap: take a's column, add gaps for b
            for (std::size_t i = 0; i < a_size; ++i) {
                sequences_[i].sequence += a.sequences_[i].sequence[j];
            }
            for (std::size_t i = 1; i < b_size; ++i) {
                sequences_[a_size + i - 1].sequence += kGapChar;
            }
            ++j;
        } else {
            // b differs or has gap: take b's column, add gaps for a
            for (std::size_t i = 0; i < a_size; ++i) {
                sequences_[i].sequence += kGapChar;
            }
            for (std::size_t i = 1; i < b_size; ++i) {
                sequences_[a_size + i - 1].sequence += b.sequences_[i].sequence[k];
            }
            ++k;
        }
    }

    // Handle remaining columns in a
    while (j < a_len) {
        for (std::size_t i = 0; i < a_size; ++i) {
            sequences_[i].sequence += a.sequences_[i].sequence[j];
        }
        for (std::size_t i = 1; i < b_size; ++i) {
            sequences_[a_size + i - 1].sequence += kGapChar;
        }
        ++j;
    }

    // Handle remaining columns in b
    while (k < b_len) {
        for (std::size_t i = 0; i < a_size; ++i) {
            sequences_[i].sequence += kGapChar;
        }
        for (std::size_t i = 1; i < b_size; ++i) {
            sequences_[a_size + i - 1].sequence += b.sequences_[i].sequence[k];
        }
        ++k;
    }
}

//==============================================================================
// Utility methods
//==============================================================================

std::vector<AlignedSequence> Alignment::leaf_sequences(const TreeNode& root) const {
    // Collect all leaf names
    std::unordered_set<std::string> leaf_names;

    std::function<void(const TreeNode&)> collect_leaves = [&](const TreeNode& node) {
        if (node.is_leaf()) {
            leaf_names.insert(node.name);
        } else {
            if (node.left) collect_leaves(*node.left);
            if (node.right) collect_leaves(*node.right);
        }
    };
    collect_leaves(root);

    // Filter sequences to only include leaves
    std::vector<AlignedSequence> result;
    result.reserve(leaf_names.size());

    for (const auto& seq : sequences_) {
        if (leaf_names.count(seq.name) > 0) {
            result.push_back(seq);
        }
    }

    return result;
}

Alignment Alignment::remove_gap_columns() const {
    if (sequences_.empty()) {
        return *this;
    }

    const std::size_t len = alignment_length();
    const std::size_t num_seqs = sequences_.size();

    // Find columns that are not all gaps
    std::vector<bool> keep_column(len, false);
    for (std::size_t col = 0; col < len; ++col) {
        for (std::size_t row = 0; row < num_seqs; ++row) {
            if (sequences_[row].sequence[col] != kGapChar) {
                keep_column[col] = true;
                break;
            }
        }
    }

    // Build new alignment without gap-only columns
    Alignment result;
    result.sequences_.resize(num_seqs);

    for (std::size_t row = 0; row < num_seqs; ++row) {
        result.sequences_[row].name = sequences_[row].name;
        result.sequences_[row].sequence.reserve(len);

        for (std::size_t col = 0; col < len; ++col) {
            if (keep_column[col]) {
                result.sequences_[row].sequence += sequences_[row].sequence[col];
            }
        }
    }

    return result;
}

std::vector<AlignedSequence> Alignment::ungapped_sequences() const {
    std::vector<AlignedSequence> result;
    result.reserve(sequences_.size());

    for (const auto& seq : sequences_) {
        std::string ungapped;
        ungapped.reserve(seq.sequence.length());

        for (char c : seq.sequence) {
            if (c != kGapChar) {
                ungapped += c;
            }
        }

        result.emplace_back(seq.name, std::move(ungapped));
    }

    return result;
}

} // namespace simprot
