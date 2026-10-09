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
    // InitSymbolCumulativeDensity in SIMPROT 1.04: running sums for 1..19
    // and a last entry forced to 1.0. For PAM the frequencies do not sum to
    // 1, so the array is not monotonic and the search below must follow the
    // original exactly.
    cumulative_frequencies_[0] = 0.0;
    for (std::size_t i = 1; i < kNumAminoAcids; ++i) {
        cumulative_frequencies_[i] = cumulative_frequencies_[i - 1] + freq[i - 1];
    }
    cumulative_frequencies_[kNumAminoAcids] = 1.0;
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

    // RandomSymbol in SIMPROT 1.04
    const auto& cdf = cumulative_frequencies_;
    int low = 0;
    int high = static_cast<int>(kNumAminoAcids);
    while (true) {
        int mid = (low + high) / 2;
        if (cdf[static_cast<std::size_t>(mid)] <= x) {
            if (cdf[static_cast<std::size_t>(mid) + 1] > x) {
                return static_cast<AminoAcidIndex>(mid);
            }
            // The original loops forever here when the array is not
            // monotonic around mid; stop instead (no reproducible run can
            // reach this).
            if (mid == low) return static_cast<AminoAcidIndex>(mid);
            low = mid;
        } else {
            high = mid;
        }
    }
}

namespace {

// MakeProtFreqs in SIMPROT 1.04: the frequencies are the absolute values of
// the eigenvector row whose eigenvalue is largest (first one on ties). For
// JTT and PMB that is row 0; the PAM data's largest eigenvalue is at index 10.
std::array<double, kNumAminoAcids> frequencies_from_eigenvectors(
    const std::array<double, kNumAminoAcids>& eigenvalues,
    const std::array<std::array<double, kNumAminoAcids>, kNumAminoAcids>& eigenvectors) {
    std::size_t maxeig = 0;
    for (std::size_t i = 0; i < kNumAminoAcids; ++i) {
        if (eigenvalues[i] > eigenvalues[maxeig]) maxeig = i;
    }
    std::array<double, kNumAminoAcids> freqs{};
    for (std::size_t i = 0; i < kNumAminoAcids; ++i) {
        freqs[i] = std::abs(eigenvectors[maxeig][i]);
    }
    return freqs;
}

}  // namespace

//==============================================================================
// PAMMatrix implementation
//==============================================================================

PAMMatrix::PAMMatrix() {
    frequencies_ = frequencies_from_eigenvectors(matrix_data::pam_eigenvalues,
                                                 matrix_data::pam_eigenvectors);
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
    frequencies_ = frequencies_from_eigenvectors(matrix_data::jtt_eigenvalues,
                                                 matrix_data::jtt_eigenvectors);
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
    frequencies_ = frequencies_from_eigenvectors(matrix_data::pmb_eigenvalues,
                                                 matrix_data::pmb_eigenvectors);
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
