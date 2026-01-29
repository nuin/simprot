#include "simprot/evolution/indel_model.hpp"

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cfloat>

namespace simprot {

namespace {
    constexpr int kMaxIndelLength = 2048;
    constexpr double kMaxIndelFraction = 0.05;  // 5% of sequence length
}

//==============================================================================
// IndelModel base class
//==============================================================================

int IndelModel::max_indel_length(int sequence_length) {
    int max_len = static_cast<int>(kMaxIndelFraction * sequence_length);
    return std::min(max_len, kMaxIndelLength);
}

int IndelModel::binary_search_cdf(const std::vector<double>& cdf, double x) {
    if (cdf.empty()) return 1;

    int low = 0;
    int high = static_cast<int>(cdf.size());

    while (true) {
        int mid = (low + high) / 2;
        if (mid == 0) return 1;
        if (mid == static_cast<int>(cdf.size())) return static_cast<int>(cdf.size());
        if (mid == low) return mid + 1;

        if (cdf[mid] <= x) {
            if (mid + 1 < static_cast<int>(cdf.size()) && cdf[mid + 1] > x) {
                return mid + 1;
            }
            low = mid;
        } else {
            high = mid;
        }
    }
}

//==============================================================================
// QianGoldsteinModel
//==============================================================================

QianGoldsteinModel::QianGoldsteinModel(double evol_scale)
    : evol_scale_(evol_scale > 0.0 ? evol_scale : 3.0) {}

int QianGoldsteinModel::sample_length(
    int sequence_length,
    double branch_distance,
    double site_rate,
    WichmannHillRNG& rng) const {

    // Compute normalized distance
    double normalized_dist = site_rate * branch_distance / evol_scale_;
    if (normalized_dist <= 0.0) {
        normalized_dist = 0.001;  // Avoid division by zero
    }

    int max_len = max_indel_length(sequence_length);

    // Build cumulative distribution using 4-term exponential formula
    std::vector<double> cdf(max_len + 1, 0.0);
    double sum = 0.0;

    for (int i = 1; i <= max_len; ++i) {
        double term1 = 1.027e-2 * std::exp(-i / (0.96 * normalized_dist));
        double term2 = 3.031e-3 * std::exp(-i / (3.13 * normalized_dist));
        double term3 = 6.141e-4 * std::exp(-i / (14.3 * normalized_dist));
        double term4 = 2.090e-5 * std::exp(-i / (81.7 * normalized_dist));

        double prob = term1 + term2 + term3 + term4;
        if (prob < DBL_EPSILON) break;

        cdf[i] = prob;
        sum += prob;
    }

    // Normalize and make cumulative
    if (sum > 0.0) {
        for (int i = 1; i <= max_len; ++i) {
            cdf[i] /= sum;
        }
        for (int i = 2; i <= max_len; ++i) {
            cdf[i] += cdf[i - 1];
            if (cdf[i - 1] >= 1.0 - DBL_EPSILON) break;
        }
    }

    // Sample from distribution
    double x = rng.uniform();
    return binary_search_cdf(cdf, x);
}

//==============================================================================
// BennerModel
//==============================================================================

BennerModel::BennerModel(int k) : k_(k) {}

int BennerModel::sample_length(
    int sequence_length,
    double branch_distance,
    double site_rate,
    WichmannHillRNG& rng) const {

    (void)branch_distance;  // Not used in Benner model
    (void)site_rate;

    int max_len = max_indel_length(sequence_length);

    // Build cumulative distribution using power-law: P(i) = i^k
    std::vector<double> cdf(max_len + 1, 0.0);
    double sum = 0.0;

    for (int i = 1; i <= max_len; ++i) {
        double prob = std::pow(static_cast<double>(i), static_cast<double>(k_));
        if (prob < DBL_EPSILON) break;

        cdf[i] = prob;
        sum += prob;
    }

    // Normalize and make cumulative
    if (sum > 0.0) {
        for (int i = 1; i <= max_len; ++i) {
            cdf[i] /= sum;
        }
        for (int i = 2; i <= max_len; ++i) {
            cdf[i] += cdf[i - 1];
            if (cdf[i - 1] >= 1.0 - DBL_EPSILON) break;
        }
    }

    // Sample from distribution
    double x = rng.uniform();
    return binary_search_cdf(cdf, x);
}

//==============================================================================
// CustomIndelModel
//==============================================================================

CustomIndelModel::CustomIndelModel(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open indel distribution file: " + filename);
    }

    std::vector<double> frequencies;
    double freq;
    while (file >> freq) {
        frequencies.push_back(freq);
    }

    if (frequencies.empty()) {
        throw std::runtime_error("Empty indel distribution file: " + filename);
    }

    // Build CDF
    double sum = 0.0;
    for (double f : frequencies) {
        sum += f;
    }

    cdf_.resize(frequencies.size() + 1, 0.0);
    if (sum > 0.0) {
        for (std::size_t i = 0; i < frequencies.size(); ++i) {
            cdf_[i + 1] = frequencies[i] / sum;
        }
        for (std::size_t i = 2; i < cdf_.size(); ++i) {
            cdf_[i] += cdf_[i - 1];
        }
    }
}

CustomIndelModel::CustomIndelModel(std::vector<double> frequencies) {
    if (frequencies.empty()) {
        throw std::runtime_error("Empty frequency vector for custom indel model");
    }

    double sum = 0.0;
    for (double f : frequencies) {
        sum += f;
    }

    cdf_.resize(frequencies.size() + 1, 0.0);
    if (sum > 0.0) {
        for (std::size_t i = 0; i < frequencies.size(); ++i) {
            cdf_[i + 1] = frequencies[i] / sum;
        }
        for (std::size_t i = 2; i < cdf_.size(); ++i) {
            cdf_[i] += cdf_[i - 1];
        }
    }
}

int CustomIndelModel::sample_length(
    int sequence_length,
    double branch_distance,
    double site_rate,
    WichmannHillRNG& rng) const {

    (void)sequence_length;
    (void)branch_distance;
    (void)site_rate;

    double x = rng.uniform();
    return binary_search_cdf(cdf_, x);
}

//==============================================================================
// Factory function
//==============================================================================

std::unique_ptr<IndelModel> create_indel_model(
    IndelDistribution distribution,
    double evol_scale,
    int benner_k,
    const std::string& custom_file) {

    switch (distribution) {
        case IndelDistribution::QianGoldstein:
            return std::make_unique<QianGoldsteinModel>(evol_scale);
        case IndelDistribution::Benner:
            return std::make_unique<BennerModel>(benner_k);
        case IndelDistribution::Custom:
            return std::make_unique<CustomIndelModel>(custom_file);
        case IndelDistribution::None:
            return std::make_unique<NoIndelModel>();
    }
    return std::make_unique<QianGoldsteinModel>(evol_scale);
}

//==============================================================================
// Utility functions
//==============================================================================

int compute_num_indels(
    double distance,
    int sequence_length,
    double indel_freq,
    double evol_scale,
    bool use_benner_scale,
    WichmannHillRNG& rng) {

    if (indel_freq <= 0.0 || indel_freq >= 1.0) {
        return 0;
    }

    // indel_rate = -log(1 - indel_freq)
    double indel_rate = -std::log(1.0 - indel_freq);

    // Probability of indel per site
    double f;
    if (use_benner_scale) {
        f = 1.0 - std::exp(-indel_rate * distance);
    } else {
        f = 1.0 - std::exp(-indel_rate * distance / evol_scale);
    }

    // Count number of indels (Poisson-like process)
    int num_indels = 0;
    for (int i = 0; i < sequence_length; ++i) {
        if (rng.uniform() < f) {
            ++num_indels;
        }
    }

    return num_indels;
}

IndelType choose_indel_type(double indel_ratio, WichmannHillRNG& rng) {
    return (rng.uniform() < indel_ratio)
        ? IndelType::Insertion
        : IndelType::Deletion;
}

int sample_indel_position(
    const MutableSequence& seq,
    IndelType type,
    WichmannHillRNG& rng) {

    int length = static_cast<int>(seq.size());
    if (length == 0) return 0;

    // COMPATIBILITY: Build CDF exactly as original SIMPROT does (InitCumulativeIndels)
    // For INSERTION: arraySize = length + 2, rate[0] is counted TWICE
    // For DELETION: arraySize = length + 1, each rate counted once

    int array_size;
    std::vector<double> cdf;
    double sum = 0.0;

    if (type == IndelType::Insertion) {
        // Original: cdf[0]=0, cdf[1]=rate[0], cdf[2]=rate[0], cdf[3]=rate[1], ...
        array_size = length + 2;
        cdf.resize(array_size, 0.0);

        // Get rates into a vector for easier indexing
        std::vector<double> rates;
        rates.reserve(length);
        seq.for_each([&](const SequenceNode& node) {
            rates.push_back(node.rate);
        });

        // Original behavior: rate[0] is added twice
        // cdf[1] = rate[0] (before loop in original)
        cdf[1] = rates[0];
        sum += rates[0];

        // Loop from i=2 with current=head in original
        // cdf[2] = rate[0], cdf[3] = rate[1], ..., cdf[length+1] = rate[length-1]
        for (int i = 2; i < array_size; ++i) {
            cdf[i] = rates[i - 2];
            sum += rates[i - 2];
        }
    } else {
        // Deletion: straightforward CDF
        array_size = length + 1;
        cdf.resize(array_size, 0.0);

        int idx = 1;
        seq.for_each([&](const SequenceNode& node) {
            cdf[idx] = node.rate;
            sum += node.rate;
            ++idx;
        });
    }

    // Normalize by sum
    if (sum > 0.0) {
        for (int i = 0; i < array_size; ++i) {
            cdf[i] /= sum;
        }
    }

    // Convert to cumulative
    for (int i = 1; i < array_size; ++i) {
        cdf[i] += cdf[i - 1];
    }

    // Sample position using binary search (GetIndelPosition)
    double x = rng.uniform();
    int search_length = array_size - 1;
    int low = 0, high = search_length;

    while (true) {
        int mid = (low + high) / 2;
        if (mid == 0) return 0;
        if (mid == search_length) return search_length - 1;
        if (mid == low) return mid;

        if (cdf[mid] <= x) {
            if (cdf[mid + 1] > x) return mid;
            low = mid;
        } else {
            high = mid;
        }
    }
}

double get_indel_rate(
    const MutableSequence& seq,
    int position,
    IndelType type) {

    // For insertion at position 0, use head's rate
    if (type == IndelType::Insertion && position == 0) {
        SequenceNode* head = seq.head();
        return head ? head->rate : 1.0;
    }

    // For insertion, the actual index is 1 smaller
    if (type == IndelType::Insertion && position > 0) {
        --position;
    }

    SequenceNode* node = seq.node_at(static_cast<std::size_t>(position));
    return node ? node->rate : 1.0;
}

void mark_indel_positions(
    MutableSequence& seq,
    int position,
    int length,
    IndelType type) {

    if (length <= 0) return;

    int seq_len = static_cast<int>(seq.size());
    if (position < 0 || position > seq_len) return;

    if (type == IndelType::Insertion) {
        // Mark the amino acid before the insertion position
        // If position == 0, mark only the first amino acid
        // If position == length, mark only the last amino acid
        int mark_pos = (position > 0) ? position - 1 : 0;
        SequenceNode* node = seq.node_at(static_cast<std::size_t>(mark_pos));
        if (node) {
            node->mark = 1;
        }
        // Also mark the position after (if exists and position > 0)
        if (position > 0 && position < seq_len) {
            SequenceNode* next_node = seq.node_at(static_cast<std::size_t>(position));
            if (next_node) {
                next_node->mark = 1;
            }
        }
    } else {
        // Deletion: mark all positions to be deleted
        for (int i = 0; i < length && (position + i) < seq_len; ++i) {
            SequenceNode* node = seq.node_at(static_cast<std::size_t>(position + i));
            if (node) {
                node->mark = 1;
            }
        }
    }
}

} // namespace simprot
