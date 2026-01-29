#ifndef SIMPROT_TREE_TREE_NODE_HPP
#define SIMPROT_TREE_TREE_NODE_HPP

/**
 * @file tree_node.hpp
 * @brief Tree node structure for phylogenetic trees.
 *
 * This header defines the TreeNode class that represents nodes in the
 * evolutionary tree. Each node stores a sequence, evolutionary rates,
 * and alignment profiles used for constructing the true alignment.
 */

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace simprot {

/**
 * @class TreeNode
 * @brief A node in the phylogenetic tree.
 *
 * TreeNode represents both internal nodes and leaf nodes in a bifurcating
 * phylogenetic tree. Internal nodes have two children; leaf nodes have none.
 *
 * During simulation:
 * 1. The root node is initialized with a sequence
 * 2. Sequences evolve down the tree via Mutate()
 * 3. Profile strings track insertions/deletions for true alignment reconstruction
 *
 * Profile explanation:
 * - left_profile_child / left_profile_self: Track how this node's sequence
 *   aligns with its left child after mutation (child has the actual residues,
 *   self has gaps where insertions occurred)
 * - right_profile_child / right_profile_self: Same for right child
 */
class TreeNode {
public:
    //==========================================================================
    // Constructors
    //==========================================================================

    /** @brief Default constructor creates an empty node. */
    TreeNode() = default;

    /**
     * @brief Construct a named node.
     * @param name Node name (from Newick tree or auto-generated for internal nodes).
     */
    explicit TreeNode(std::string name);

    /**
     * @brief Construct a node with name and branch length.
     * @param name Node name.
     * @param distance Branch length to parent.
     */
    TreeNode(std::string name, double distance);

    // Move operations
    TreeNode(TreeNode&&) noexcept = default;
    TreeNode& operator=(TreeNode&&) noexcept = default;

    // No copying (tree structure uses unique_ptr)
    TreeNode(const TreeNode&) = delete;
    TreeNode& operator=(const TreeNode&) = delete;

    //==========================================================================
    // Tree structure
    //==========================================================================

    /** @brief Name of this node (taxon name for leaves, auto-generated for internal). */
    std::string name;

    /** @brief Branch length (evolutionary distance) to parent node. */
    double distance = 0.0;

    /** @brief Left child node (nullptr if leaf). */
    std::unique_ptr<TreeNode> left;

    /** @brief Right child node (nullptr if leaf). */
    std::unique_ptr<TreeNode> right;

    //==========================================================================
    // Sequence data
    //==========================================================================

    /**
     * @brief Amino acid sequence at this node (no gaps).
     *
     * For leaves, this is the final simulated sequence.
     * For internal nodes, this is the ancestral sequence.
     */
    std::string sequence;

    /**
     * @brief Site-specific evolutionary rates.
     *
     * Each position in the sequence has a rate drawn from a gamma distribution.
     * rates[i] corresponds to sequence[i].
     */
    std::vector<double> rates;

    //==========================================================================
    // Alignment profiles
    //==========================================================================

    /**
     * @brief Profile showing child sequence aligned to this node (left child).
     *
     * Contains the left child's sequence with gaps where the child had
     * deletions relative to this node.
     */
    std::string left_profile_child;

    /**
     * @brief Profile showing this node's sequence aligned to left child.
     *
     * Contains this node's sequence with gaps where the child had
     * insertions relative to this node.
     */
    std::string left_profile_self;

    /**
     * @brief Profile showing child sequence aligned to this node (right child).
     */
    std::string right_profile_child;

    /**
     * @brief Profile showing this node's sequence aligned to right child.
     */
    std::string right_profile_self;

    //==========================================================================
    // Status flags
    //==========================================================================

    /**
     * @brief Whether this branch is "extinct" (zero-length).
     *
     * Used for branch extinction simulation where some lineages don't evolve.
     */
    bool extinct = false;

    //==========================================================================
    // Query methods
    //==========================================================================

    /**
     * @brief Check if this node is a leaf (has no children).
     * @return true if this is a leaf node.
     */
    [[nodiscard]] bool is_leaf() const noexcept {
        return !left && !right;
    }

    /**
     * @brief Check if this node is internal (has children).
     * @return true if this is an internal node.
     */
    [[nodiscard]] bool is_internal() const noexcept {
        return left || right;
    }

    /**
     * @brief Check if this node has a sequence.
     * @return true if sequence is non-empty.
     */
    [[nodiscard]] bool has_sequence() const noexcept {
        return !sequence.empty();
    }

    /**
     * @brief Get sequence length.
     * @return Length of the sequence, or 0 if no sequence.
     */
    [[nodiscard]] std::size_t sequence_length() const noexcept {
        return sequence.size();
    }

    //==========================================================================
    // Tree operations
    //==========================================================================

    /**
     * @brief Count total nodes in subtree rooted at this node.
     * @return Number of nodes (including this one).
     */
    [[nodiscard]] std::size_t count_nodes() const;

    /**
     * @brief Count leaf nodes in subtree rooted at this node.
     * @return Number of leaf nodes.
     */
    [[nodiscard]] std::size_t count_leaves() const;

    /**
     * @brief Clear sequence and profile data (but keep tree structure).
     *
     * Used when running multiple simulations on the same tree.
     */
    void clear_sequences();

    /**
     * @brief Clear only the profile data.
     *
     * Profiles are only needed during alignment construction.
     */
    void clear_profiles();
};

} // namespace simprot

#endif // SIMPROT_TREE_TREE_NODE_HPP
