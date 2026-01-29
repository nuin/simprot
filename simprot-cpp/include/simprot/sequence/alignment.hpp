#ifndef SIMPROT_SEQUENCE_ALIGNMENT_HPP
#define SIMPROT_SEQUENCE_ALIGNMENT_HPP

/**
 * @file alignment.hpp
 * @brief True alignment reconstruction from evolved sequences.
 *
 * This header defines the Alignment class which builds true alignments
 * from the profile information stored during sequence evolution. The
 * alignment is reconstructed recursively by merging child alignments
 * with parent profiles.
 */

#include "simprot/core/types.hpp"
#include "simprot/tree/tree_node.hpp"

#include <string>
#include <vector>

namespace simprot {

/**
 * @struct AlignedSequence
 * @brief A single sequence in an alignment with its name.
 */
struct AlignedSequence {
    std::string name;      ///< Sequence identifier (taxon name)
    std::string sequence;  ///< Aligned sequence (may contain gaps)

    AlignedSequence() = default;
    AlignedSequence(std::string n, std::string s)
        : name(std::move(n)), sequence(std::move(s)) {}
};

/**
 * @class Alignment
 * @brief True alignment of sequences from a phylogenetic tree.
 *
 * This class reconstructs the true alignment from the profile information
 * stored at each tree node during evolution. The reconstruction uses a
 * recursive algorithm that:
 * 1. At leaves: creates a single-sequence alignment
 * 2. At internal nodes: merges left and right child alignments using
 *    the stored profiles to determine gap positions
 *
 * Usage:
 * @code
 * // After evolution is complete
 * Alignment aln = Alignment::from_tree(root);
 * for (const auto& seq : aln.sequences()) {
 *     std::cout << ">" << seq.name << "\n" << seq.sequence << "\n";
 * }
 * @endcode
 */
class Alignment {
public:
    /**
     * @brief Construct an empty alignment.
     */
    Alignment() = default;

    /**
     * @brief Construct alignment from a tree (recursive).
     * @param root Root node of the evolved tree.
     * @return The complete true alignment.
     *
     * This is the main entry point for alignment reconstruction.
     */
    [[nodiscard]] static Alignment from_tree(const TreeNode& root);

    /**
     * @brief Get all aligned sequences.
     */
    [[nodiscard]] const std::vector<AlignedSequence>& sequences() const {
        return sequences_;
    }

    /**
     * @brief Get number of sequences in alignment.
     */
    [[nodiscard]] std::size_t num_sequences() const {
        return sequences_.size();
    }

    /**
     * @brief Get alignment length (all sequences have same length).
     */
    [[nodiscard]] std::size_t alignment_length() const {
        return sequences_.empty() ? 0 : sequences_[0].sequence.length();
    }

    /**
     * @brief Check if alignment is empty.
     */
    [[nodiscard]] bool empty() const {
        return sequences_.empty();
    }

    /**
     * @brief Get sequences for leaf nodes only.
     * @param root The tree root (to determine which are leaves).
     * @return Vector of aligned sequences for leaves only.
     */
    [[nodiscard]] std::vector<AlignedSequence> leaf_sequences(const TreeNode& root) const;

    /**
     * @brief Remove columns that are all gaps.
     * @return New alignment with gap-only columns removed.
     */
    [[nodiscard]] Alignment remove_gap_columns() const;

    /**
     * @brief Get sequences without gaps (unaligned).
     * @return Vector of sequences with all gaps removed.
     */
    [[nodiscard]] std::vector<AlignedSequence> ungapped_sequences() const;

private:
    std::vector<AlignedSequence> sequences_;

    /**
     * @brief Recursive constructor for tree traversal.
     */
    explicit Alignment(const TreeNode& node);

    /**
     * @brief Align a child alignment with a parent node's profile.
     * @param child_aln The child alignment.
     * @param parent The parent node.
     * @param direction Left or right child.
     */
    Alignment(const Alignment& child_aln,
              const TreeNode& parent,
              ChildDirection direction);

    /**
     * @brief Merge two alignments using their first rows as profiles.
     * @param a First alignment.
     * @param b Second alignment.
     */
    Alignment(const Alignment& a, const Alignment& b);
};

} // namespace simprot

#endif // SIMPROT_SEQUENCE_ALIGNMENT_HPP
