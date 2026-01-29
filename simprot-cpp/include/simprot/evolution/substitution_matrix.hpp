#ifndef SIMPROT_EVOLUTION_SUBSTITUTION_MATRIX_HPP
#define SIMPROT_EVOLUTION_SUBSTITUTION_MATRIX_HPP

/**
 * @file substitution_matrix.hpp
 * @brief Amino acid substitution matrices for protein evolution.
 *
 * This header defines the SubstitutionMatrix interface and implementations
 * for PAM, JTT, and PMB matrices. The matrices are stored in eigendecomposed
 * form for efficient computation of transition probabilities.
 */

#include "simprot/core/types.hpp"
#include "simprot/core/random.hpp"

#include <array>
#include <cmath>
#include <memory>

namespace simprot {

/**
 * @class SubstitutionMatrix
 * @brief Abstract base class for substitution matrices.
 *
 * Substitution matrices define the probability of one amino acid substituting
 * for another over evolutionary time. The matrices are stored in eigendecomposed
 * form: eigenvalues (λ) and eigenvector matrix (V).
 *
 * The transition probability P(j|i,t) is computed as:
 *   P(j|i,t) = Σₖ V[k][j] · V[k][i] · exp(λ[k] · t) / π[j]
 *
 * where π[j] is the equilibrium frequency of amino acid j.
 */
class SubstitutionMatrix {
public:
    virtual ~SubstitutionMatrix() = default;

    /**
     * @brief Get the name of this substitution model.
     */
    [[nodiscard]] virtual std::string_view name() const = 0;

    /**
     * @brief Get the eigenvalues array.
     */
    [[nodiscard]] virtual const std::array<double, kNumAminoAcids>& eigenvalues() const = 0;

    /**
     * @brief Get the eigenvector matrix.
     *
     * probmat[k][i] is the k-th eigenvector's i-th component.
     */
    [[nodiscard]] virtual const std::array<std::array<double, kNumAminoAcids>, kNumAminoAcids>&
    eigenvectors() const = 0;

    /**
     * @brief Get equilibrium amino acid frequencies.
     *
     * These are derived from the eigenvector corresponding to eigenvalue 0.
     */
    [[nodiscard]] virtual const std::array<double, kNumAminoAcids>& frequencies() const = 0;

    /**
     * @brief Compute the substitution probability P(j|i,t).
     * @param from Source amino acid index (i).
     * @param to Target amino acid index (j).
     * @param time Evolutionary time/distance.
     * @return Probability of substituting 'from' with 'to'.
     */
    [[nodiscard]] double substitution_probability(
        AminoAcidIndex from, AminoAcidIndex to, double time) const;

    /**
     * @brief Sample a substitution given source amino acid and time.
     * @param from Source amino acid index.
     * @param time Evolutionary time/distance.
     * @param rng Random number generator.
     * @return Target amino acid index after substitution.
     *
     * This implements the GetSubstitution() logic from the original code.
     */
    [[nodiscard]] AminoAcidIndex sample_substitution(
        AminoAcidIndex from, double time, WichmannHillRNG& rng) const;

    /**
     * @brief Sample a random amino acid according to equilibrium frequencies.
     * @param rng Random number generator.
     * @return Sampled amino acid index.
     */
    [[nodiscard]] AminoAcidIndex sample_from_frequencies(WichmannHillRNG& rng) const;

protected:
    /**
     * @brief Initialize cumulative frequency distribution.
     *
     * Called by derived classes after setting up frequencies.
     */
    void init_cumulative_frequencies();

    std::array<double, kNumAminoAcids + 1> cumulative_frequencies_{};
};

/**
 * @brief Create a substitution matrix for the given model.
 * @param model The substitution model type.
 * @return Unique pointer to the matrix instance.
 */
[[nodiscard]] std::unique_ptr<SubstitutionMatrix>
create_substitution_matrix(SubstitutionModel model);

//==============================================================================
// Concrete matrix implementations
//==============================================================================

/**
 * @class PAMMatrix
 * @brief PAM (Point Accepted Mutation) substitution matrix.
 *
 * Based on Dayhoff et al. (1978).
 */
class PAMMatrix final : public SubstitutionMatrix {
public:
    PAMMatrix();

    [[nodiscard]] std::string_view name() const override { return "PAM"; }
    [[nodiscard]] const std::array<double, kNumAminoAcids>& eigenvalues() const override;
    [[nodiscard]] const std::array<std::array<double, kNumAminoAcids>, kNumAminoAcids>&
    eigenvectors() const override;
    [[nodiscard]] const std::array<double, kNumAminoAcids>& frequencies() const override;

private:
    std::array<double, kNumAminoAcids> frequencies_{};
};

/**
 * @class JTTMatrix
 * @brief JTT (Jones-Taylor-Thornton) substitution matrix.
 *
 * Based on Jones et al. (1992).
 */
class JTTMatrix final : public SubstitutionMatrix {
public:
    JTTMatrix();

    [[nodiscard]] std::string_view name() const override { return "JTT"; }
    [[nodiscard]] const std::array<double, kNumAminoAcids>& eigenvalues() const override;
    [[nodiscard]] const std::array<std::array<double, kNumAminoAcids>, kNumAminoAcids>&
    eigenvectors() const override;
    [[nodiscard]] const std::array<double, kNumAminoAcids>& frequencies() const override;

private:
    std::array<double, kNumAminoAcids> frequencies_{};
};

/**
 * @class PMBMatrix
 * @brief PMB (Probability Matrix from Blocks) substitution matrix.
 *
 * Based on Veerassamy et al. (2003).
 */
class PMBMatrix final : public SubstitutionMatrix {
public:
    PMBMatrix();

    [[nodiscard]] std::string_view name() const override { return "PMB"; }
    [[nodiscard]] const std::array<double, kNumAminoAcids>& eigenvalues() const override;
    [[nodiscard]] const std::array<std::array<double, kNumAminoAcids>, kNumAminoAcids>&
    eigenvectors() const override;
    [[nodiscard]] const std::array<double, kNumAminoAcids>& frequencies() const override;

private:
    std::array<double, kNumAminoAcids> frequencies_{};
};

} // namespace simprot

#endif // SIMPROT_EVOLUTION_SUBSTITUTION_MATRIX_HPP
