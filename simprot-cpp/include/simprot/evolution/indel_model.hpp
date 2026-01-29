#ifndef SIMPROT_EVOLUTION_INDEL_MODEL_HPP
#define SIMPROT_EVOLUTION_INDEL_MODEL_HPP

/**
 * @file indel_model.hpp
 * @brief Indel length distribution models for sequence evolution.
 *
 * This header defines the IndelModel interface and implementations for
 * different indel length distributions: Qian-Goldstein (4-term exponential),
 * Benner (Zipfian power-law), and custom distributions from file.
 */

#include "simprot/core/types.hpp"
#include "simprot/core/random.hpp"
#include "simprot/sequence/mutable_sequence.hpp"

#include <memory>
#include <vector>
#include <string>

namespace simprot {

/**
 * @class IndelModel
 * @brief Abstract base class for indel length distributions.
 *
 * Indel models define how the lengths of insertions and deletions are
 * sampled during sequence evolution. The model computes a cumulative
 * distribution function (CDF) based on sequence length, branch distance,
 * and site rate, then samples from this distribution.
 */
class IndelModel {
public:
    virtual ~IndelModel() = default;

    /**
     * @brief Get the name of this indel model.
     */
    [[nodiscard]] virtual std::string_view name() const = 0;

    /**
     * @brief Sample an indel length.
     * @param sequence_length Current sequence length.
     * @param branch_distance Branch length (evolutionary distance).
     * @param site_rate Gamma rate at the indel position.
     * @param rng Random number generator.
     * @return Sampled indel length (>= 1).
     */
    [[nodiscard]] virtual int sample_length(
        int sequence_length,
        double branch_distance,
        double site_rate,
        WichmannHillRNG& rng) const = 0;

protected:
    /**
     * @brief Maximum indel length (5% of sequence length, capped).
     */
    [[nodiscard]] static int max_indel_length(int sequence_length);

    /**
     * @brief Binary search in cumulative distribution to find sampled length.
     * @param cdf Cumulative distribution function.
     * @param x Random value in [0, 1].
     * @return Sampled length (1-indexed).
     */
    [[nodiscard]] static int binary_search_cdf(
        const std::vector<double>& cdf, double x);
};

/**
 * @class QianGoldsteinModel
 * @brief Qian-Goldstein 4-term exponential indel length model.
 *
 * The probability of indel length i is:
 *   P(i) = 1.027e-2 * exp(-i / (0.96 * d)) +
 *          3.031e-3 * exp(-i / (3.13 * d)) +
 *          6.141e-4 * exp(-i / (14.3 * d)) +
 *          2.090e-5 * exp(-i / (81.7 * d))
 *
 * where d = rate * distance / evol_scale
 *
 * Based on: Qian & Goldstein (2001) Proteins 45:102-104
 */
class QianGoldsteinModel final : public IndelModel {
public:
    /**
     * @brief Construct with evolution scale parameter.
     * @param evol_scale Scale factor for distance (default 3.0).
     */
    explicit QianGoldsteinModel(double evol_scale = 3.0);

    [[nodiscard]] std::string_view name() const override { return "Qian-Goldstein"; }

    [[nodiscard]] int sample_length(
        int sequence_length,
        double branch_distance,
        double site_rate,
        WichmannHillRNG& rng) const override;

private:
    double evol_scale_;
};

/**
 * @class BennerModel
 * @brief Benner/Zipfian power-law indel length model.
 *
 * The probability of indel length i is:
 *   P(i) = i^k
 *
 * where k is typically negative (default -2) to favor shorter indels.
 *
 * Based on: Benner et al. (1993) J. Mol. Biol. 229:1065-1082
 */
class BennerModel final : public IndelModel {
public:
    /**
     * @brief Construct with power-law exponent.
     * @param k Exponent (typically negative, default -2).
     */
    explicit BennerModel(int k = -2);

    [[nodiscard]] std::string_view name() const override { return "Benner"; }

    [[nodiscard]] int sample_length(
        int sequence_length,
        double branch_distance,
        double site_rate,
        WichmannHillRNG& rng) const override;

private:
    int k_;
};

/**
 * @class CustomIndelModel
 * @brief User-provided indel length distribution from file.
 *
 * Reads a file containing indel length frequencies (one per line).
 * Line 1 = frequency for length 1, Line 2 = frequency for length 2, etc.
 */
class CustomIndelModel final : public IndelModel {
public:
    /**
     * @brief Construct from a distribution file.
     * @param filename Path to the distribution file.
     */
    explicit CustomIndelModel(const std::string& filename);

    /**
     * @brief Construct from a frequency vector.
     * @param frequencies Frequencies for lengths 1, 2, 3, ...
     */
    explicit CustomIndelModel(std::vector<double> frequencies);

    [[nodiscard]] std::string_view name() const override { return "Custom"; }

    [[nodiscard]] int sample_length(
        int sequence_length,
        double branch_distance,
        double site_rate,
        WichmannHillRNG& rng) const override;

private:
    std::vector<double> cdf_;  // Pre-computed cumulative distribution
};

/**
 * @class NoIndelModel
 * @brief Model that produces no indels (substitution-only evolution).
 */
class NoIndelModel final : public IndelModel {
public:
    [[nodiscard]] std::string_view name() const override { return "None"; }

    [[nodiscard]] int sample_length(
        int sequence_length,
        double branch_distance,
        double site_rate,
        WichmannHillRNG& rng) const override {
        (void)sequence_length;
        (void)branch_distance;
        (void)site_rate;
        (void)rng;
        return 0;
    }
};

/**
 * @brief Factory function to create an indel model.
 * @param distribution The distribution type.
 * @param evol_scale Evolution scale for Qian-Goldstein.
 * @param benner_k Exponent for Benner model.
 * @param custom_file File path for custom distribution.
 * @return Unique pointer to the model instance.
 */
[[nodiscard]] std::unique_ptr<IndelModel> create_indel_model(
    IndelDistribution distribution,
    double evol_scale = 3.0,
    int benner_k = -2,
    const std::string& custom_file = "");

//==============================================================================
// Utility functions for indel operations
//==============================================================================

/**
 * @brief Compute the number of indels for a branch.
 *
 * Uses a Poisson-like process: for each site, the probability of an indel
 * event is f = 1 - exp(-indel_rate * distance / scale), where
 * indel_rate = -log(1 - indel_freq).
 *
 * @param distance Branch length.
 * @param sequence_length Current sequence length.
 * @param indel_freq Base indel frequency (e.g., 0.03).
 * @param evol_scale Evolution scale factor.
 * @param use_benner_scale If true, don't divide by evol_scale.
 * @param rng Random number generator.
 * @return Number of indels to perform.
 */
[[nodiscard]] int compute_num_indels(
    double distance,
    int sequence_length,
    double indel_freq,
    double evol_scale,
    bool use_benner_scale,
    WichmannHillRNG& rng);

/**
 * @brief Choose whether the next indel is an insertion or deletion.
 * @param indel_ratio Probability of insertion (0.5 = equal).
 * @param rng Random number generator.
 * @return IndelType::Insertion or IndelType::Deletion.
 */
[[nodiscard]] IndelType choose_indel_type(
    double indel_ratio,
    WichmannHillRNG& rng);

/**
 * @brief Sample an indel position in the sequence.
 *
 * Builds a cumulative distribution based on site rates and samples
 * from it using binary search.
 *
 * @param seq The sequence.
 * @param type Insertion or deletion.
 * @param rng Random number generator.
 * @return 0-based position index.
 */
[[nodiscard]] int sample_indel_position(
    const MutableSequence& seq,
    IndelType type,
    WichmannHillRNG& rng);

/**
 * @brief Get the rate at an indel position.
 * @param seq The sequence.
 * @param position 0-based position.
 * @param type Insertion or deletion (affects position interpretation).
 * @return Rate at the position.
 */
[[nodiscard]] double get_indel_rate(
    const MutableSequence& seq,
    int position,
    IndelType type);

/**
 * @brief Mark positions in the sequence for an indel operation.
 * @param seq The sequence.
 * @param position Starting position.
 * @param length Indel length.
 * @param type Insertion or deletion.
 */
void mark_indel_positions(
    MutableSequence& seq,
    int position,
    int length,
    IndelType type);

} // namespace simprot

#endif // SIMPROT_EVOLUTION_INDEL_MODEL_HPP
