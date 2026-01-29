#ifndef SIMPROT_IO_NEWICK_READER_HPP
#define SIMPROT_IO_NEWICK_READER_HPP

/**
 * @file newick_reader.hpp
 * @brief Convenience wrapper for reading Newick format tree files.
 *
 * This header provides a simple interface for reading phylogenetic trees
 * from files or strings in Newick format.
 */

#include "simprot/tree/tree_node.hpp"
#include "simprot/tree/tree_parser.hpp"
#include "simprot/core/config.hpp"
#include "simprot/core/random.hpp"

#include <memory>
#include <string>

namespace simprot {

/**
 * @brief Read a tree from a Newick format file.
 * @param filename Path to the Newick file.
 * @return Root node of the parsed tree.
 * @throws ParseError if the format is invalid.
 * @throws std::runtime_error if the file cannot be read.
 */
[[nodiscard]] std::unique_ptr<TreeNode> read_newick_file(const std::string& filename);

/**
 * @brief Read a tree from a Newick format string.
 * @param newick The Newick format string.
 * @return Root node of the parsed tree.
 * @throws ParseError if the format is invalid.
 */
[[nodiscard]] std::unique_ptr<TreeNode> read_newick_string(const std::string& newick);

/**
 * @brief Read a tree with configuration options.
 * @param filename Path to the Newick file.
 * @param config Simulation configuration (for branch scaling, etc.).
 * @param rng Random number generator (for variable branch lengths).
 * @return Root node of the parsed tree.
 * @throws ParseError if the format is invalid.
 * @throws std::runtime_error if the file cannot be read.
 *
 * This function applies configuration options like:
 * - Branch length scaling (tree_branch_scale)
 * - Variable branch lengths (variable_branch_gamma)
 * - Branch extinction (branch_extinction_prob)
 */
[[nodiscard]] std::unique_ptr<TreeNode> read_newick_file(
    const std::string& filename,
    const SimulationConfig& config,
    WichmannHillRNG& rng);

/**
 * @brief Count the number of leaf nodes in a tree.
 * @param root Root of the tree.
 * @return Number of leaf nodes.
 */
[[nodiscard]] std::size_t count_leaves(const TreeNode& root);

/**
 * @brief Count the total number of nodes in a tree.
 * @param root Root of the tree.
 * @return Total number of nodes.
 */
[[nodiscard]] std::size_t count_nodes(const TreeNode& root);

/**
 * @brief Get the names of all leaf nodes.
 * @param root Root of the tree.
 * @return Vector of leaf node names.
 */
[[nodiscard]] std::vector<std::string> get_leaf_names(const TreeNode& root);

} // namespace simprot

#endif // SIMPROT_IO_NEWICK_READER_HPP
