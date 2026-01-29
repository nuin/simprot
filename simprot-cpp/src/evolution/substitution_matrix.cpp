#include "simprot/evolution/substitution_matrix.hpp"
#include "simprot/evolution/matrix_data.hpp"

#include <algorithm>
#include <cmath>

namespace simprot {

//==============================================================================
// SubstitutionMatrix base class implementation
//==============================================================================

void SubstitutionMatrix::init_cumulative_frequencies() {
    const auto& freq = frequencies();
    cumulative_frequencies_[0] = 0.0;
    for (std::size_t i = 0; i < kNumAminoAcids; ++i) {
        cumulative_frequencies_[i + 1] = cumulative_frequencies_[i] + freq[i];
    }
}

double SubstitutionMatrix::substitution_probability(
    AminoAcidIndex from, AminoAcidIndex to, double time) const {

    const auto& eigvals = eigenvalues();
    const auto& eigvecs = eigenvectors();
    const auto& freq = frequencies();

    // P(j|i,t) = Σₖ V[k][j] · V[k][i] · exp(λ[k] · t) / π[j]
    double prob = 0.0;
    for (std::size_t k = 0; k < kNumAminoAcids; ++k) {
        prob += eigvecs[k][to] * eigvecs[k][from] * std::exp(eigvals[k] * time);
    }

    return prob / freq[to];
}

AminoAcidIndex SubstitutionMatrix::sample_substitution(
    AminoAcidIndex from, double time, WichmannHillRNG& rng) const {

    // This implements GetSubstitution() from the original code
    const auto& eigvals = eigenvalues();
    const auto& eigvecs = eigenvectors();
    const auto& freq = frequencies();

    // Pre-compute exp(λ[k] * t) for all eigenvalues
    std::array<double, kNumAminoAcids> exp_eigmat;
    for (std::size_t k = 0; k < kNumAminoAcids; ++k) {
        exp_eigmat[k] = std::exp(time * eigvals[k]);
    }

    // Get random threshold
    double x = rng.uniform();

    // Compute cumulative probability and find the amino acid
    double sum = 0.0;
    for (AminoAcidIndex j = 0; j < kNumAminoAcids && sum < x; ++j) {
        double prob = 0.0;
        for (std::size_t k = 0; k < kNumAminoAcids; ++k) {
            prob += eigvecs[k][j] * eigvecs[k][from] * exp_eigmat[k];
        }
        sum += prob / freq[j];

        if (sum >= x) {
            return j;
        }
    }

    // Should not reach here, but return last amino acid as fallback
    return kNumAminoAcids - 1;
}

AminoAcidIndex SubstitutionMatrix::sample_from_frequencies(WichmannHillRNG& rng) const {
    double x = rng.uniform();

    // Binary search in cumulative distribution
    auto it = std::upper_bound(
        cumulative_frequencies_.begin(),
        cumulative_frequencies_.end(),
        x);

    if (it == cumulative_frequencies_.begin()) {
        return 0;
    }

    auto idx = static_cast<AminoAcidIndex>(
        std::distance(cumulative_frequencies_.begin(), it) - 1);

    return std::min(idx, static_cast<AminoAcidIndex>(kNumAminoAcids - 1));
}

//==============================================================================
// PAMMatrix implementation
//==============================================================================

PAMMatrix::PAMMatrix() {
    // Compute equilibrium frequencies from the eigenvector with eigenvalue 0
    // (first row of eigenvector matrix)
    const auto& eigvecs = matrix_data::pam_eigenvectors;
    for (std::size_t i = 0; i < kNumAminoAcids; ++i) {
        frequencies_[i] = std::abs(eigvecs[0][i]);
    }
    init_cumulative_frequencies();
}

const std::array<double, kNumAminoAcids>& PAMMatrix::eigenvalues() const {
    return matrix_data::pam_eigenvalues;
}

const std::array<std::array<double, kNumAminoAcids>, kNumAminoAcids>&
PAMMatrix::eigenvectors() const {
    return matrix_data::pam_eigenvectors;
}

const std::array<double, kNumAminoAcids>& PAMMatrix::frequencies() const {
    return frequencies_;
}

//==============================================================================
// JTTMatrix implementation
//==============================================================================

JTTMatrix::JTTMatrix() {
    const auto& eigvecs = matrix_data::jtt_eigenvectors;
    for (std::size_t i = 0; i < kNumAminoAcids; ++i) {
        frequencies_[i] = std::abs(eigvecs[0][i]);
    }
    init_cumulative_frequencies();
}

const std::array<double, kNumAminoAcids>& JTTMatrix::eigenvalues() const {
    return matrix_data::jtt_eigenvalues;
}

const std::array<std::array<double, kNumAminoAcids>, kNumAminoAcids>&
JTTMatrix::eigenvectors() const {
    return matrix_data::jtt_eigenvectors;
}

const std::array<double, kNumAminoAcids>& JTTMatrix::frequencies() const {
    return frequencies_;
}

//==============================================================================
// PMBMatrix implementation
//==============================================================================

PMBMatrix::PMBMatrix() {
    const auto& eigvecs = matrix_data::pmb_eigenvectors;
    for (std::size_t i = 0; i < kNumAminoAcids; ++i) {
        frequencies_[i] = std::abs(eigvecs[0][i]);
    }
    init_cumulative_frequencies();
}

const std::array<double, kNumAminoAcids>& PMBMatrix::eigenvalues() const {
    return matrix_data::pmb_eigenvalues;
}

const std::array<std::array<double, kNumAminoAcids>, kNumAminoAcids>&
PMBMatrix::eigenvectors() const {
    return matrix_data::pmb_eigenvectors;
}

const std::array<double, kNumAminoAcids>& PMBMatrix::frequencies() const {
    return frequencies_;
}

//==============================================================================
// Factory function
//==============================================================================

std::unique_ptr<SubstitutionMatrix> create_substitution_matrix(SubstitutionModel model) {
    switch (model) {
        case SubstitutionModel::PAM:
            return std::make_unique<PAMMatrix>();
        case SubstitutionModel::JTT:
            return std::make_unique<JTTMatrix>();
        case SubstitutionModel::PMB:
            return std::make_unique<PMBMatrix>();
    }
    // Default to PMB
    return std::make_unique<PMBMatrix>();
}

} // namespace simprot
