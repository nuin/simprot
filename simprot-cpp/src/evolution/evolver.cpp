#include "simprot/evolution/evolver.hpp"

#include <algorithm>
#include <cstring>

namespace simprot {

SequenceEvolver::SequenceEvolver(const SimulationConfig& config, WichmannHillRNG& rng)
    : config_(config)
    , rng_(rng)
    , substitution_matrix_(create_substitution_matrix(config.substitution_model))
    , indel_model_(create_indel_model(
        config.indel_distribution,
        config.evol_scale,
        config.benner_k,
        config.indel_distribution_file.value_or(""))) {}

void SequenceEvolver::init_root_sequence(TreeNode& root) {
    std::size_t length = static_cast<std::size_t>(config_.root_sequence_length);
    root.sequence = generate_random_sequence(length);
    root.rates = generate_rates(length);
}

void SequenceEvolver::init_root_sequence(TreeNode& root, const std::string& sequence) {
    root.sequence = sequence;
    root.rates = generate_rates(sequence.size());
}

void SequenceEvolver::evolve(TreeNode& root) {
    // Make sure root has a sequence
    if (root.sequence.empty()) {
        throw std::runtime_error("Cannot evolve: root has no sequence");
    }

    // Evolve left subtree
    if (root.left) {
        mutate(root, *root.left, ChildDirection::Left);
        evolve(*root.left);
    }

    // Evolve right subtree
    if (root.right) {
        mutate(root, *root.right, ChildDirection::Right);
        evolve(*root.right);
    }
}

void SequenceEvolver::set_indel_callback(
    std::function<void(double, IndelType, int)> callback) {
    indel_callback_ = std::move(callback);
}

void SequenceEvolver::mutate(TreeNode& parent, TreeNode& child, ChildDirection direction) {
    // If branch length is zero, child is identical to parent
    if (child.distance <= 0.0) {
        child.sequence = parent.sequence;
        child.rates = parent.rates;

        // Set up trivial profiles (identical sequences)
        if (direction == ChildDirection::Left) {
            parent.left_profile_child = parent.sequence;
            parent.left_profile_self = parent.sequence;
        } else {
            parent.right_profile_child = parent.sequence;
            parent.right_profile_self = parent.sequence;
        }
        return;
    }

    // Convert parent sequence to mutable sequence
    MutableSequence seq(parent.sequence, parent.rates);

    // Normalize rates so average = 1.0
    seq.normalize_rates();

    // Initialize profiles
    std::string profile_child = parent.sequence;
    std::string profile_self = parent.sequence;

    // Get number of indels for this branch
    bool use_benner_scale = (config_.indel_distribution == IndelDistribution::Benner);
    int num_indels = compute_num_indels(
        child.distance,
        static_cast<int>(seq.size()),
        config_.indel_frequency,
        config_.evol_scale,
        use_benner_scale,
        rng_);

    // Perform each indel
    for (int i = 0; i < num_indels; ++i) {
        IndelType type = choose_indel_type(config_.indel_ratio, rng_);

        // Sample position and length
        int position = sample_indel_position(seq, type, rng_);
        double site_rate = get_indel_rate(seq, position, type);

        int indel_length = indel_model_->sample_length(
            static_cast<int>(seq.size()),
            child.distance,
            site_rate,
            rng_);

        if (indel_length > 0) {
            // Log the indel if callback is set
            if (indel_callback_) {
                indel_callback_(child.distance, type, indel_length);
            }

            // Mark positions and perform the indel
            mark_indel_positions(seq, position, indel_length, type);
            perform_indel(type, seq, profile_child, profile_self, indel_length);
        }
    }

    // Apply substitutions
    apply_substitutions(seq, profile_child, child.distance);

    // Store child sequence and rates
    child.sequence = seq.to_string();
    child.rates = seq.rates();

    // Store profiles in parent
    if (direction == ChildDirection::Left) {
        parent.left_profile_child = std::move(profile_child);
        parent.left_profile_self = std::move(profile_self);
    } else {
        parent.right_profile_child = std::move(profile_child);
        parent.right_profile_self = std::move(profile_self);
    }
}

void SequenceEvolver::perform_indel(
    IndelType type,
    MutableSequence& seq,
    std::string& profile_child,
    std::string& profile_self,
    int indel_size) {

    // Find the marked position
    SequenceNode* marked_node = nullptr;
    int marked_pos = -1;
    int idx = 0;

    seq.for_each([&](SequenceNode& node) {
        if (node.mark && !marked_node) {
            marked_node = &node;
            marked_pos = idx;
        }
        ++idx;
    });

    if (!marked_node) {
        seq.clear_marks();
        return;
    }

    // Map marked_pos to profile position (accounting for gaps)
    int profile_pos = 0;
    int seq_pos = 0;
    for (std::size_t i = 0; i < profile_child.size(); ++i) {
        if (profile_child[i] != kGapChar) {
            if (seq_pos == marked_pos) {
                profile_pos = static_cast<int>(i);
                break;
            }
            ++seq_pos;
        }
    }

    if (type == IndelType::Insertion) {
        // Perform insertion in sequence
        SequenceNode* insert_before = marked_node->next;
        SequenceNode* first_inserted = seq.insert_before(
            insert_before,
            static_cast<std::size_t>(indel_size),
            rng_,
            *substitution_matrix_,
            config_.gamma_alpha);

        // Update profiles
        // In the child profile, insert the new residues
        // In the self profile, insert gaps
        std::string inserted_chars;
        std::string gaps(indel_size, kGapChar);

        SequenceNode* node = first_inserted;
        for (int i = 0; i < indel_size && node; ++i) {
            inserted_chars += node->residue;
            node = node->next;
        }

        // Insert after the marked position in profiles
        int insert_pos = profile_pos + 1;
        if (insert_pos > static_cast<int>(profile_child.size())) {
            insert_pos = static_cast<int>(profile_child.size());
        }

        profile_child.insert(static_cast<std::size_t>(insert_pos), inserted_chars);
        profile_self.insert(static_cast<std::size_t>(insert_pos), gaps);

    } else {
        // Deletion
        // Find all marked nodes and delete them
        SequenceNode* delete_start = nullptr;
        int delete_count = 0;

        seq.for_each([&](SequenceNode& node) {
            if (node.mark) {
                if (!delete_start) {
                    delete_start = &node;
                }
                ++delete_count;
            }
        });

        if (delete_start && delete_count > 0) {
            // Update profiles before deleting
            // In child profile, replace residues with gaps
            // In self profile, keep the original residues
            int profile_idx = profile_pos;
            int deleted = 0;

            // Find the residues being deleted and replace with gaps
            for (std::size_t i = static_cast<std::size_t>(profile_idx);
                 i < profile_child.size() && deleted < delete_count;
                 ++i) {
                if (profile_child[i] != kGapChar) {
                    profile_child[i] = kGapChar;
                    ++deleted;
                }
            }

            // Perform deletion in sequence
            seq.delete_at(delete_start, static_cast<std::size_t>(delete_count));
        }
    }

    seq.clear_marks();
}

void SequenceEvolver::apply_substitutions(
    MutableSequence& seq,
    std::string& profile_child,
    double distance) {

    // Scale factor for PAM/JTT vs PMB
    double scale = (config_.substitution_model == SubstitutionModel::PMB)
        ? 1.0
        : 100.0;

    // Iterate through sequence and profile together
    std::size_t profile_idx = 0;

    seq.for_each([&](SequenceNode& node) {
        // Skip gaps in profile
        while (profile_idx < profile_child.size() &&
               profile_child[profile_idx] == kGapChar) {
            ++profile_idx;
        }

        if (profile_idx >= profile_child.size()) return;

        // Compute substitution
        double t = distance * node.rate * scale;
        AminoAcidIndex from = char_to_amino_acid(node.residue);
        AminoAcidIndex to = substitution_matrix_->sample_substitution(from, t, rng_);
        char new_residue = amino_acid_to_char(to);

        // Update sequence and profile
        node.residue = new_residue;
        profile_child[profile_idx] = new_residue;

        ++profile_idx;
    });
}

std::string SequenceEvolver::generate_random_sequence(std::size_t length) {
    std::string sequence;
    sequence.reserve(length);

    for (std::size_t i = 0; i < length; ++i) {
        AminoAcidIndex aa = substitution_matrix_->sample_from_frequencies(rng_);
        sequence += amino_acid_to_char(aa);
    }

    return sequence;
}

std::vector<double> SequenceEvolver::generate_rates(std::size_t length) {
    std::vector<double> rates;
    rates.reserve(length);

    // If gamma_alpha < 0, use equal rates (1.0)
    // Otherwise, sample from gamma distribution
    if (config_.gamma_alpha < 0.0) {
        rates.assign(length, 1.0);
    } else {
        for (std::size_t i = 0; i < length; ++i) {
            rates.push_back(rng_.gamma(config_.gamma_alpha));
        }
    }

    return rates;
}

} // namespace simprot
