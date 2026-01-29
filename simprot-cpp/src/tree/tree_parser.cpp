#include "simprot/tree/tree_parser.hpp"

#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace simprot {

std::unique_ptr<TreeNode> NewickParser::parse(std::string_view newick) {
    // Reset internal counter for each parse
    internal_count_ = 0;

    // Make a mutable copy for parsing
    std::string_view input = newick;
    skip_whitespace(input);

    if (input.empty()) {
        throw ParseError("Empty input");
    }

    auto root = parse_subtree(input);

    // Set root name
    root->name = "InternalRoot";

    return root;
}

std::unique_ptr<TreeNode> NewickParser::parse_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open tree file: " + filename);
    }

    // Read entire file into string, stripping newlines
    std::ostringstream ss;
    std::string line;
    while (std::getline(file, line)) {
        // Remove any whitespace and append
        for (char c : line) {
            if (!std::isspace(static_cast<unsigned char>(c))) {
                ss << c;
            }
        }
    }

    return parse(ss.str());
}

std::unique_ptr<TreeNode> NewickParser::parse_subtree(std::string_view& input) {
    auto node = std::make_unique<TreeNode>();
    bool expect_label = true;
    bool expect_distance = false;

    while (!input.empty()) {
        skip_whitespace(input);
        if (input.empty()) break;

        std::string token = get_next_token(input);
        if (token.empty()) break;

        char first_char = token[0];

        if (first_char == '(') {
            // Start of descendant list - parse left and right children
            node->left = parse_subtree(input);
            node->right = parse_subtree(input);
            expect_label = true;  // May have an internal node label after ')'

        } else if (first_char == ',' || first_char == ';' || first_char == ')') {
            // End of this subtree

            // If we haven't assigned a name yet, generate one for internal node
            if (node->name.empty()) {
                node->name = generate_internal_name();
            }

            return node;

        } else if (first_char == ':') {
            // Next token will be a distance
            expect_distance = true;
            expect_label = false;

        } else {
            // This is either a label or a distance
            if (expect_distance) {
                // Parse as distance
                try {
                    double base_length = std::stod(token);
                    node->distance = process_branch_length(base_length, *node);
                } catch (const std::exception&) {
                    throw ParseError("Invalid branch length: " + token);
                }
                expect_distance = false;

            } else if (expect_label) {
                // Parse as label
                node->name = token;
                expect_label = false;
            }
        }
    }

    // If we reach here without returning, something went wrong
    throw ParseError("Unexpected end of input");
}

std::string NewickParser::get_next_token(std::string_view& input) {
    skip_whitespace(input);

    if (input.empty()) {
        return "";
    }

    // Check if it's a single-character delimiter
    if (is_delimiter(input[0])) {
        char c = input[0];
        input.remove_prefix(1);
        return std::string(1, c);
    }

    // Otherwise, read until we hit a delimiter
    std::size_t end = 0;
    while (end < input.size() && !is_delimiter(input[end]) &&
           !std::isspace(static_cast<unsigned char>(input[end]))) {
        ++end;
    }

    std::string token(input.substr(0, end));
    input.remove_prefix(end);
    return token;
}

void NewickParser::skip_whitespace(std::string_view& input) {
    while (!input.empty() &&
           std::isspace(static_cast<unsigned char>(input[0]))) {
        input.remove_prefix(1);
    }
}

double NewickParser::process_branch_length(double base_length, TreeNode& node) {
    // IMPORTANT: Always consume RNG for extinction check to match original behavior.
    // The original simprot.cpp always calls rndu() at line 618 during tree parsing,
    // regardless of whether BranchExtinction is 0. This maintains RNG state sync.
    double extinction_roll = rng_ ? rng_->uniform() : 1.0;

    // Check for extinction
    if (extinction_prob_ > 0.0 && extinction_roll <= extinction_prob_) {
        node.extinct = true;
        node.name += " Neg";  // Mark as negative/extinct
        return 0.0;
    }

    // Apply branch scaling
    double length = base_length * branch_scale_;

    // Apply variable branch gamma if enabled
    if (rng_ && variable_branch_gamma_ > 0.0) {
        double gamma_factor = rng_->gamma(variable_branch_gamma_);
        length *= gamma_factor;
    }

    return length;
}

std::string NewickParser::generate_internal_name() {
    std::ostringstream ss;
    ss << "Internal_" << std::setw(2) << std::setfill('0') << internal_count_;
    ++internal_count_;
    return ss.str();
}

} // namespace simprot
