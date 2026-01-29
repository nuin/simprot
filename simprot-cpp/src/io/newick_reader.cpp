#include "simprot/io/newick_reader.hpp"

#include <functional>

namespace simprot {

std::unique_ptr<TreeNode> read_newick_file(const std::string& filename) {
    NewickParser parser;
    return parser.parse_file(filename);
}

std::unique_ptr<TreeNode> read_newick_string(const std::string& newick) {
    NewickParser parser;
    return parser.parse(newick);
}

std::unique_ptr<TreeNode> read_newick_file(
    const std::string& filename,
    const SimulationConfig& config,
    WichmannHillRNG& rng) {

    NewickParser parser(&rng);

    // Apply configuration options
    parser.set_branch_scale(config.tree_branch_scale);
    parser.set_variable_branch_gamma(config.variable_branch_gamma);
    parser.set_extinction_probability(config.branch_extinction_prob);

    return parser.parse_file(filename);
}

std::size_t count_leaves(const TreeNode& root) {
    if (root.is_leaf()) {
        return 1;
    }

    std::size_t count = 0;
    if (root.left) count += count_leaves(*root.left);
    if (root.right) count += count_leaves(*root.right);
    return count;
}

std::size_t count_nodes(const TreeNode& root) {
    return root.count_nodes();
}

std::vector<std::string> get_leaf_names(const TreeNode& root) {
    std::vector<std::string> names;

    std::function<void(const TreeNode&)> collect = [&](const TreeNode& node) {
        if (node.is_leaf()) {
            names.push_back(node.name);
        } else {
            if (node.left) collect(*node.left);
            if (node.right) collect(*node.right);
        }
    };

    collect(root);
    return names;
}

} // namespace simprot
