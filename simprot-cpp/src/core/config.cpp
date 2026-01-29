#include "simprot/core/config.hpp"

#include <stdexcept>
#include <string>

namespace simprot {

void SimulationConfig::validate() const {
    // Root sequence length
    if (root_sequence_length <= 0) {
        throw std::invalid_argument(
            "root_sequence_length must be positive, got " +
            std::to_string(root_sequence_length));
    }

    // If root sequence provided, validate it
    if (root_sequence.has_value()) {
        const auto& seq = root_sequence.value();
        for (char c : seq) {
            try {
                (void)char_to_amino_acid(c);
            } catch (const std::invalid_argument&) {
                throw std::invalid_argument(
                    "Invalid amino acid '" + std::string(1, c) +
                    "' in root sequence");
            }
        }
    }

    // Gamma parameters
    // Note: gamma_alpha <= 0 means equal rates (valid)
    if (gamma_beta <= 0.0) {
        throw std::invalid_argument(
            "gamma_beta must be positive, got " + std::to_string(gamma_beta));
    }

    // Indel frequency
    if (indel_frequency < 0.0 || indel_frequency >= 1.0) {
        throw std::invalid_argument(
            "indel_frequency must be in [0, 1), got " +
            std::to_string(indel_frequency));
    }

    // Indel ratio
    if (indel_ratio < 0.0 || indel_ratio > 1.0) {
        throw std::invalid_argument(
            "indel_ratio must be in [0, 1], got " +
            std::to_string(indel_ratio));
    }

    // Max indel length
    if (max_indel_length <= 0) {
        throw std::invalid_argument(
            "max_indel_length must be positive, got " +
            std::to_string(max_indel_length));
    }

    // Evolution scale
    if (evol_scale <= 0.0) {
        throw std::invalid_argument(
            "evol_scale must be positive, got " + std::to_string(evol_scale));
    }

    // Tree branch scale
    if (tree_branch_scale <= 0.0) {
        throw std::invalid_argument(
            "tree_branch_scale must be positive, got " +
            std::to_string(tree_branch_scale));
    }

    // Variable branch gamma
    if (variable_branch_gamma < 0.0) {
        throw std::invalid_argument(
            "variable_branch_gamma must be non-negative, got " +
            std::to_string(variable_branch_gamma));
    }

    // Branch extinction probability
    if (branch_extinction_prob < 0.0 || branch_extinction_prob > 1.0) {
        throw std::invalid_argument(
            "branch_extinction_prob must be in [0, 1], got " +
            std::to_string(branch_extinction_prob));
    }

    // Extra terminal indels
    if (extra_terminal_indels < 0.0) {
        throw std::invalid_argument(
            "extra_terminal_indels must be non-negative, got " +
            std::to_string(extra_terminal_indels));
    }

    // Tree file is required
    if (tree_file.empty()) {
        throw std::invalid_argument("tree_file is required");
    }

    // Custom indel distribution requires frequencies
    if (indel_distribution == IndelDistribution::Custom &&
        custom_indel_frequencies.empty() &&
        !indel_distribution_file.has_value()) {
        throw std::invalid_argument(
            "Custom indel distribution requires either custom_indel_frequencies "
            "or indel_distribution_file");
    }
}

} // namespace simprot
