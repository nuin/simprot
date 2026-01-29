#ifndef SIMPROT_TREE_TREE_PARSER_HPP
#define SIMPROT_TREE_TREE_PARSER_HPP

/**
 * @file tree_parser.hpp
 * @brief Newick format tree parser.
 *
 * This header defines a recursive descent parser for Newick format phylogenetic
 * trees. The parser constructs a TreeNode structure from the string representation.
 *
 * Newick format: https://en.wikipedia.org/wiki/Newick_format
 * Example: ((A:0.1,B:0.2):0.3,(C:0.1,D:0.2):0.4);
 */

#include "simprot/tree/tree_node.hpp"
#include "simprot/core/random.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

namespace simprot {

/**
 * @class ParseError
 * @brief Exception thrown when Newick parsing fails.
 */
class ParseError : public std::runtime_error {
public:
    explicit ParseError(const std::string& message)
        : std::runtime_error("Parse error: " + message) {}
};

/**
 * @class NewickParser
 * @brief Recursive descent parser for Newick format trees.
 *
 * The parser reads a Newick-formatted string and constructs a binary tree
 * of TreeNode objects. Branch lengths are scaled according to configuration,
 * and optional random variation can be applied.
 *
 * Usage:
 * @code
 * NewickParser parser;
 * auto tree = parser.parse("((A:0.1,B:0.2):0.3,C:0.4);");
 * @endcode
 */
class NewickParser {
public:
    /**
     * @brief Construct a parser with default settings.
     */
    NewickParser() = default;

    /**
     * @brief Construct a parser with an RNG for variable branch lengths.
     * @param rng Pointer to RNG (parser does not take ownership).
     */
    explicit NewickParser(WichmannHillRNG* rng) : rng_(rng) {}

    /**
     * @brief Set the branch length scaling factor.
     * @param scale All branch lengths are multiplied by this value.
     */
    void set_branch_scale(double scale) { branch_scale_ = scale; }

    /**
     * @brief Enable variable branch lengths with gamma-distributed scaling.
     * @param gamma_shape Shape parameter for gamma distribution (0 = disabled).
     */
    void set_variable_branch_gamma(double gamma_shape) {
        variable_branch_gamma_ = gamma_shape;
    }

    /**
     * @brief Set probability of branch extinction.
     * @param prob Probability that a branch is set to zero length.
     */
    void set_extinction_probability(double prob) {
        extinction_prob_ = prob;
    }

    /**
     * @brief Parse a Newick format string into a tree.
     * @param newick The Newick format string.
     * @return Root node of the parsed tree.
     * @throws ParseError if the format is invalid.
     */
    [[nodiscard]] std::unique_ptr<TreeNode> parse(std::string_view newick);

    /**
     * @brief Parse a Newick format string from a file.
     * @param filename Path to the file containing the Newick tree.
     * @return Root node of the parsed tree.
     * @throws ParseError if the format is invalid.
     * @throws std::runtime_error if the file cannot be read.
     */
    [[nodiscard]] std::unique_ptr<TreeNode> parse_file(const std::string& filename);

private:
    /**
     * @brief Recursive parsing function.
     * @param pos Current position in the input string.
     * @return Parsed subtree.
     */
    std::unique_ptr<TreeNode> parse_subtree(std::string_view& input);

    /**
     * @brief Get the next token from the input.
     * @param input Input string view (modified to consume the token).
     * @return The token string.
     */
    std::string get_next_token(std::string_view& input);

    /**
     * @brief Skip whitespace in the input.
     * @param input Input string view (modified to skip whitespace).
     */
    void skip_whitespace(std::string_view& input);

    /**
     * @brief Check if a character is a Newick delimiter.
     */
    [[nodiscard]] static bool is_delimiter(char c) {
        return c == '(' || c == ')' || c == ',' || c == ':' || c == ';';
    }

    /**
     * @brief Process branch length with scaling and optional variation.
     * @param base_length The raw branch length from the Newick string.
     * @param node The node to potentially mark as extinct.
     * @return The processed branch length.
     */
    double process_branch_length(double base_length, TreeNode& node);

    /**
     * @brief Generate a unique name for an internal node.
     */
    std::string generate_internal_name();

    // Configuration
    double branch_scale_ = 1.0;
    double variable_branch_gamma_ = 0.0;
    double extinction_prob_ = 0.0;

    // RNG for variable branches (optional)
    WichmannHillRNG* rng_ = nullptr;

    // Counter for internal node names
    int internal_count_ = 0;
};

} // namespace simprot

#endif // SIMPROT_TREE_TREE_PARSER_HPP
