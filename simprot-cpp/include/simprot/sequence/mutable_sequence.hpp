#ifndef SIMPROT_SEQUENCE_MUTABLE_SEQUENCE_HPP
#define SIMPROT_SEQUENCE_MUTABLE_SEQUENCE_HPP

/**
 * @file mutable_sequence.hpp
 * @brief Doubly-linked list representation of sequences for efficient indels.
 *
 * This header defines MutableSequence, a doubly-linked list structure that
 * enables efficient insertion and deletion operations during sequence evolution.
 * Each node stores an amino acid, its evolutionary rate, and a marker for indel
 * operations.
 */

#include "simprot/core/types.hpp"
#include "simprot/core/random.hpp"
#include "simprot/evolution/substitution_matrix.hpp"

#include <memory>
#include <string>
#include <vector>
#include <functional>

namespace simprot {

/**
 * @class SequenceNode
 * @brief A node in the mutable sequence doubly-linked list.
 *
 * Each node represents a single amino acid site with its evolutionary rate
 * and a marker used during indel operations.
 */
class SequenceNode {
public:
    char residue{};          ///< Amino acid at this site
    double rate{1.0};        ///< Evolutionary rate at this site
    int mark{0};             ///< Marker for indel operations

    SequenceNode* prev{nullptr};
    SequenceNode* next{nullptr};

    SequenceNode() = default;
    explicit SequenceNode(char c, double r = 1.0)
        : residue(c), rate(r) {}
};

/**
 * @class MutableSequence
 * @brief Doubly-linked list for efficient sequence indel operations.
 *
 * This class provides a doubly-linked list representation of a protein sequence
 * that enables O(1) insertion and deletion at any position (given a node pointer).
 * It is used during the Mutate() operation to efficiently apply indels.
 *
 * Usage:
 * @code
 * MutableSequence seq("ACDEF", rates);
 * // Insert 3 residues before position
 * auto inserted = seq.insert_before(node_ptr, 3, rng, matrix);
 * // Delete 2 residues starting at position
 * auto next = seq.delete_at(node_ptr, 2);
 * @endcode
 */
class MutableSequence {
public:
    /**
     * @brief Construct an empty mutable sequence.
     */
    MutableSequence() = default;

    /**
     * @brief Construct from a sequence string with optional rates.
     * @param sequence Amino acid sequence string.
     * @param rates Evolutionary rates per site (nullptr = generate randomly).
     * @param rng Random number generator for rate generation.
     * @param gamma_alpha Alpha parameter for gamma distribution.
     */
    MutableSequence(const std::string& sequence,
                    const std::vector<double>* rates,
                    WichmannHillRNG& rng,
                    double gamma_alpha);

    /**
     * @brief Construct from sequence string with provided rates vector.
     * @param sequence Amino acid sequence string.
     * @param rates Evolutionary rates per site (must match sequence length).
     */
    MutableSequence(const std::string& sequence,
                    const std::vector<double>& rates);

    // Disable copy (nodes own raw pointers)
    MutableSequence(const MutableSequence&) = delete;
    MutableSequence& operator=(const MutableSequence&) = delete;

    // Enable move
    MutableSequence(MutableSequence&& other) noexcept;
    MutableSequence& operator=(MutableSequence&& other) noexcept;

    ~MutableSequence();

    /**
     * @brief Get the head node.
     */
    [[nodiscard]] SequenceNode* head() const noexcept { return head_; }

    /**
     * @brief Get the tail node.
     */
    [[nodiscard]] SequenceNode* tail() const noexcept { return tail_; }

    /**
     * @brief Get current sequence length.
     */
    [[nodiscard]] std::size_t size() const noexcept { return size_; }

    /**
     * @brief Check if sequence is empty.
     */
    [[nodiscard]] bool empty() const noexcept { return size_ == 0; }

    /**
     * @brief Convert to string representation.
     */
    [[nodiscard]] std::string to_string() const;

    /**
     * @brief Extract rates into a vector.
     */
    [[nodiscard]] std::vector<double> rates() const;

    /**
     * @brief Normalize rates so they sum to sequence length.
     *
     * After normalization: sum(rates) = size(), which means average rate = 1.0
     */
    void normalize_rates();

    /**
     * @brief Insert a random sequence before a given node.
     * @param before Node before which to insert (nullptr = append at tail).
     * @param length Number of residues to insert.
     * @param rng Random number generator.
     * @param matrix Substitution matrix for sampling residues.
     * @param gamma_alpha Gamma shape for rate sampling.
     * @return Pointer to the first inserted node.
     */
    SequenceNode* insert_before(SequenceNode* before,
                                std::size_t length,
                                WichmannHillRNG& rng,
                                const SubstitutionMatrix& matrix,
                                double gamma_alpha);

    /**
     * @brief Delete nodes starting at a given position.
     * @param start First node to delete.
     * @param count Number of nodes to delete.
     * @return Pointer to the node after the deleted segment (may be nullptr).
     */
    SequenceNode* delete_at(SequenceNode* start, std::size_t count);

    /**
     * @brief Get the node at a specific index.
     * @param index 0-based index.
     * @return Pointer to the node, or nullptr if out of range.
     */
    [[nodiscard]] SequenceNode* node_at(std::size_t index) const;

    /**
     * @brief Clear all marks on nodes.
     */
    void clear_marks();

    /**
     * @brief Iterate over all nodes with a callback.
     * @param callback Function called for each node.
     */
    void for_each(const std::function<void(SequenceNode&)>& callback);

    /**
     * @brief Iterate over all nodes with a const callback.
     * @param callback Function called for each node.
     */
    void for_each(const std::function<void(const SequenceNode&)>& callback) const;

private:
    SequenceNode* head_{nullptr};
    SequenceNode* tail_{nullptr};
    std::size_t size_{0};

    /**
     * @brief Append a new node at the tail.
     */
    void append(char residue, double rate);

    /**
     * @brief Clear all nodes and reset state.
     */
    void clear();
};

} // namespace simprot

#endif // SIMPROT_SEQUENCE_MUTABLE_SEQUENCE_HPP
