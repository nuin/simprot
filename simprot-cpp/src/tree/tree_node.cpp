#include "simprot/tree/tree_node.hpp"

namespace simprot {

TreeNode::TreeNode(std::string node_name)
    : name(std::move(node_name)) {}

TreeNode::TreeNode(std::string node_name, double dist)
    : name(std::move(node_name)), distance(dist) {}

std::size_t TreeNode::count_nodes() const {
    std::size_t count = 1;  // This node

    if (left) {
        count += left->count_nodes();
    }
    if (right) {
        count += right->count_nodes();
    }

    return count;
}

std::size_t TreeNode::count_leaves() const {
    if (is_leaf()) {
        return 1;
    }

    std::size_t count = 0;
    if (left) {
        count += left->count_leaves();
    }
    if (right) {
        count += right->count_leaves();
    }

    return count;
}

void TreeNode::clear_sequences() {
    sequence.clear();
    rates.clear();
    clear_profiles();

    if (left) {
        left->clear_sequences();
    }
    if (right) {
        right->clear_sequences();
    }
}

void TreeNode::clear_profiles() {
    left_profile_child.clear();
    left_profile_self.clear();
    right_profile_child.clear();
    right_profile_self.clear();

    if (left) {
        left->clear_profiles();
    }
    if (right) {
        right->clear_profiles();
    }
}

} // namespace simprot
