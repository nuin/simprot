#ifndef SIMPROT_EVOLUTION_EVOLVER_HPP
#define SIMPROT_EVOLUTION_EVOLVER_HPP

/**
 * @file evolver.hpp
 * @brief Main sequence evolution engine.
 *
 * This header defines the SequenceEvolver class, which orchestrates the
 * simulation of protein sequence evolution along a phylogenetic tree.
 * It implements the Evolve() recursive traversal and Mutate() operations
 * for substitutions and indels.
 */

#include "simprot/core/types.hpp"
#include "simprot/core/config.hpp"
#include "simprot/core/random.hpp"
#include "simprot/tree/tree_node.hpp"
#include "simprot/evolution/substitution_matrix.hpp"
#include "simprot/evolution/indel_model.hpp"
#include "simprot/sequence/mutable_sequence.hpp"

#include <memory>
#include <string>
#include <functional>

namespace simprot {

/**
 * @class SequenceEvolver
 * @brief Main engine for simulating protein sequence evolution.
 *
 * The evolver traverses a phylogenetic tree, applying mutations (substitutions
 * and indels) along each branch to simulate sequence evolution. It maintains
 * profile information at each node for true alignment reconstruction.
 *
 * Usage:
 * @code
 * SequenceEvolver evolver(config, rng);
 * evolver.init_root_sequence(tree_root);
 * evolver.evolve(tree_root);
 * @endcode
 */
class SequenceEvolver {
public:
    /**
     * @brief Construct an evolver with configuration and RNG.
     * @param config Simulation configuration.
     * @param rng Random number generator.
     */
    SequenceEvolver(const SimulationConfig& config, WichmannHillRNG& rng);

    /**
     * @brief Initialize the root sequence with random amino acids.
     * @param root Root node of the tree.
     *
     * Generates a random sequence of the configured length using the
     * substitution matrix's equilibrium frequencies. Also assigns
     * gamma-distributed rates to each site.
     */
    void init_root_sequence(TreeNode& root);

    /**
     * @brief Initialize the root sequence from a provided string.
     * @param root Root node of the tree.
     * @param sequence Initial sequence string.
     */
    void init_root_sequence(TreeNode& root, const std::string& sequence);

    /**
     * @brief Evolve sequences along the entire tree.
     * @param root Root node of the tree.
     *
     * Performs a preorder traversal, calling mutate() on each branch
     * to generate child sequences from parent sequences.
     */
    void evolve(TreeNode& root);

    /**
     * @brief Set a callback for indel logging.
     * @param callback Function called for each indel event.
     *
     * The callback receives: (distance, type, length)
     */
    void set_indel_callback(
        std::function<void(double, IndelType, int)> callback);

    /**
     * @brief Get the substitution matrix being used.
     */
    [[nodiscard]] const SubstitutionMatrix& substitution_matrix() const {
        return *substitution_matrix_;
    }

    /**
     * @brief Get the indel model being used.
     */
    [[nodiscard]] const IndelModel& indel_model() const {
        return *indel_model_;
    }

private:
    /**
     * @brief Mutate a parent sequence to produce a child sequence.
     * @param parent Parent tree node.
     * @param child Child tree node.
     * @param direction Whether child is left or right.
     *
     * This is the core mutation function that:
     * 1. Converts parent sequence to a MutableSequence
     * 2. Applies indels (updating profiles)
     * 3. Applies substitutions per site
     * 4. Stores the result in the child node
     */
    void mutate(TreeNode& parent, TreeNode& child, ChildDirection direction);

    /**
     * @brief Perform an indel operation and update profiles.
     * @param type Insertion or deletion.
     * @param seq The mutable sequence.
     * @param profile_child Profile for the child sequence.
     * @param profile_self Profile for the parent (self) sequence.
     * @param indel_size Size of the indel.
     */
    void perform_indel(
        IndelType type,
        MutableSequence& seq,
        std::string& profile_child,
        std::string& profile_self,
        int indel_size);

    /**
     * @brief Apply substitutions to all sites in a sequence.
     * @param seq The mutable sequence.
     * @param profile_child Profile being updated.
     * @param distance Branch length.
     */
    void apply_substitutions(
        MutableSequence& seq,
        std::string& profile_child,
        double distance);

    /**
     * @brief Generate a random sequence of given length.
     * @param length Sequence length.
     * @return Random sequence string.
     */
    [[nodiscard]] std::string generate_random_sequence(std::size_t length);

    /**
     * @brief Generate gamma-distributed rates for sites.
     * @param length Number of sites.
     * @return Vector of rates.
     */
    [[nodiscard]] std::vector<double> generate_rates(std::size_t length);

    // Configuration
    const SimulationConfig& config_;
    WichmannHillRNG& rng_;

    // Models
    std::unique_ptr<SubstitutionMatrix> substitution_matrix_;
    std::unique_ptr<IndelModel> indel_model_;

    // Indel callback
    std::function<void(double, IndelType, int)> indel_callback_;
};

} // namespace simprot

#endif // SIMPROT_EVOLUTION_EVOLVER_HPP
