#ifndef SIMPROT_CORE_CONFIG_HPP
#define SIMPROT_CORE_CONFIG_HPP

/**
 * @file config.hpp
 * @brief Simulation configuration for SIMPROT.
 *
 * This header defines the SimulationConfig struct that holds all parameters
 * for a protein evolution simulation. This replaces the global variables
 * used in the original implementation.
 */

#include "simprot/core/types.hpp"

#include <optional>
#include <string>
#include <vector>

namespace simprot {

/**
 * @struct SimulationConfig
 * @brief Complete configuration for a SIMPROT simulation.
 *
 * All simulation parameters are encapsulated here, making it easy to
 * create reproducible simulations and avoiding global state.
 */
struct SimulationConfig {
    //==========================================================================
    // Random number generation
    //==========================================================================

    /** @brief Random seed for reproducibility. */
    int seed = 12345;

    //==========================================================================
    // Sequence parameters
    //==========================================================================

    /** @brief Length of the root sequence (amino acids). Default: 50. */
    int root_sequence_length = 50;

    /** @brief Optional user-provided root sequence (overrides random generation). */
    std::optional<std::string> root_sequence;

    //==========================================================================
    // Substitution model parameters
    //==========================================================================

    /** @brief Substitution matrix model (PAM, JTT, or PMB). Default: PMB. */
    SubstitutionModel substitution_model = SubstitutionModel::PMB;

    //==========================================================================
    // Rate variation parameters (Gamma distribution)
    //==========================================================================

    /**
     * @brief Alpha parameter for gamma-distributed site rates.
     *
     * Controls the amount of rate heterogeneity across sites.
     * - alpha < 1: High rate variation (some sites evolve much faster)
     * - alpha = 1: Exponential distribution
     * - alpha > 1: More uniform rates
     * - alpha = -1 (or <= 0): Equal rates (no variation)
     *
     * Default: 1.0
     */
    double gamma_alpha = 1.0;

    /** @brief Beta parameter for gamma distribution. Default: 1.0. */
    double gamma_beta = 1.0;

    //==========================================================================
    // Indel parameters
    //==========================================================================

    /**
     * @brief Indel frequency parameter.
     *
     * Controls the probability of indel events per site per unit branch length.
     * Higher values = more indels. Default: 0.03.
     */
    double indel_frequency = 0.03;

    /** @brief Indel length distribution model. Default: QianGoldstein. */
    IndelDistribution indel_distribution = IndelDistribution::QianGoldstein;

    /**
     * @brief Exponent for Benner (Zipfian) indel length distribution.
     *
     * Only used when indel_distribution = Benner.
     * Length probability ~ length^k. Default: -2.
     */
    int benner_k = -2;

    /**
     * @brief Ratio of insertions to deletions.
     *
     * Probability of insertion = indel_ratio.
     * Probability of deletion = 1 - indel_ratio.
     * Default: 0.5 (equal probability).
     */
    double indel_ratio = 0.5;

    /** @brief Maximum indel length. Default: 2048. */
    int max_indel_length = 2048;

    /**
     * @brief Evolution scale factor for indel length distribution.
     *
     * Used in Qian-Goldstein formula: normalDistance = rate * distance / evol_scale.
     * Default: 3.0.
     */
    double evol_scale = 3.0;

    /** @brief Custom indel length frequencies (for Custom distribution). */
    std::vector<double> custom_indel_frequencies;

    /**
     * @brief Extra terminal indels multiplier.
     *
     * If > 0, adds additional indels weighted toward sequence termini.
     * Default: 0.0 (disabled).
     */
    double extra_terminal_indels = 0.0;

    //==========================================================================
    // Tree scaling parameters
    //==========================================================================

    /**
     * @brief Scale factor for all branch lengths.
     *
     * All branch lengths read from the tree file are multiplied by this value.
     * Default: 1.0.
     */
    double tree_branch_scale = 1.0;

    /**
     * @brief Gamma shape for variable branch length scaling.
     *
     * If > 0, each branch length is multiplied by a gamma-distributed random
     * factor with this shape parameter. Default: 0.0 (disabled).
     */
    double variable_branch_gamma = 0.0;

    /**
     * @brief Probability of branch extinction.
     *
     * If > 0, each branch has this probability of being "extinct" (set to
     * zero length). Default: 0.0 (disabled).
     */
    double branch_extinction_prob = 0.0;

    //==========================================================================
    // Output options
    //==========================================================================

    /** @brief Output file for FASTA sequences (leaf nodes only). */
    std::optional<std::string> fasta_output_file;

    /** @brief Output file for true alignment (all sequences with gaps). */
    std::optional<std::string> alignment_output_file;

    /** @brief Output file for Phylip format alignment. */
    std::optional<std::string> phylip_output_file;

    /** @brief Output file for indel event log. */
    std::optional<std::string> indel_log_file;

    //==========================================================================
    // Input files
    //==========================================================================

    /** @brief Input tree file (Newick format). Required. */
    std::string tree_file;

    /** @brief Optional file with user-provided root sequence. */
    std::optional<std::string> root_sequence_file;

    /** @brief Optional file with custom indel length distribution. */
    std::optional<std::string> indel_distribution_file;

    //==========================================================================
    // Validation
    //==========================================================================

    /**
     * @brief Validate configuration parameters.
     * @throws std::invalid_argument if any parameter is invalid.
     */
    void validate() const;
};

} // namespace simprot

#endif // SIMPROT_CORE_CONFIG_HPP
