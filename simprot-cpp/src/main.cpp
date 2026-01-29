/**
 * SIMPROT 2.0 - Modern C++20 Implementation
 *
 * A protein sequence evolution simulator that generates multiple sequence
 * alignments by simulating amino acid substitutions and insertion/deletion
 * events along phylogenetic trees.
 *
 * Original: Nuin PAS., Wang Z., and ERM Tillier. (2006)
 *           The accuracy of several multiple sequence alignments for proteins.
 *           BMC Bioinformatics 7:471.
 */

#include "simprot/core/types.hpp"
#include "simprot/core/config.hpp"
#include "simprot/core/random.hpp"
#include "simprot/tree/tree_node.hpp"
#include "simprot/io/newick_reader.hpp"
#include "simprot/io/fasta_writer.hpp"
#include "simprot/evolution/evolver.hpp"
#include "simprot/sequence/alignment.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>

namespace {

void print_usage(const char* program_name) {
    std::cerr << "SIMPROT 2.0 - Protein Sequence Evolution Simulator\n\n";
    std::cerr << "Usage: " << program_name << " [options]\n\n";
    std::cerr << "Required:\n";
    std::cerr << "  -f <file>    Input tree file (Newick format)\n\n";
    std::cerr << "Output options:\n";
    std::cerr << "  -a <file>    Output true alignment file (FASTA with gaps)\n";
    std::cerr << "  -s <file>    Output sequences file (FASTA without gaps)\n";
    std::cerr << "  -o <file>    Output indel log file\n\n";
    std::cerr << "Sequence options:\n";
    std::cerr << "  -r <int>     Root sequence length (default: 50)\n";
    std::cerr << "  -i <file>    Input root sequence file\n\n";
    std::cerr << "Model options:\n";
    std::cerr << "  -p <int>     Substitution model: 0=PAM, 1=JTT, 2=PMB (default: 2)\n";
    std::cerr << "  -x <float>   Gamma alpha for rate variation (-1 = equal rates, default: 1.0)\n";
    std::cerr << "  -g <float>   Indel frequency (default: 0.03)\n";
    std::cerr << "  -b <int>     Indel model: 0=Qian-Goldstein, 1=Benner, 3=None (default: 0)\n";
    std::cerr << "  -k <int>     Benner exponent k (default: -2)\n";
    std::cerr << "  -u <file>    Custom indel distribution file\n\n";
    std::cerr << "Tree options:\n";
    std::cerr << "  -t <float>   Branch length scale factor (default: 1.0)\n";
    std::cerr << "  -v <float>   Variable branch gamma (0 = disabled, default: 0)\n";
    std::cerr << "  -e <float>   Branch extinction probability (default: 0)\n\n";
    std::cerr << "Other options:\n";
    std::cerr << "  -S <int>     Random seed (default: 12345)\n";
    std::cerr << "  -h           Show this help message\n";
}

void print_version() {
    std::cout << "SIMPROT 2.0.0\n";
    std::cout << "Modern C++20 implementation of protein sequence evolution simulator\n";
    std::cout << "Original: Nuin, Wang & Tillier (2006) BMC Bioinformatics 7:471\n";
}

} // anonymous namespace

int main(int argc, char* argv[]) {
    simprot::SimulationConfig config;

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        }
        if (arg == "--version") {
            print_version();
            return 0;
        }

        // Options that take a value
        if (i + 1 >= argc && arg[0] == '-' && arg.length() == 2) {
            std::cerr << "Error: Option " << arg << " requires an argument\n";
            return 1;
        }

        if (arg == "-f") {
            config.tree_file = argv[++i];
        } else if (arg == "-a") {
            config.alignment_output_file = argv[++i];
        } else if (arg == "-s") {
            config.fasta_output_file = argv[++i];
        } else if (arg == "-o") {
            config.indel_log_file = argv[++i];
        } else if (arg == "-r") {
            config.root_sequence_length = std::atoi(argv[++i]);
        } else if (arg == "-i") {
            config.root_sequence_file = argv[++i];
        } else if (arg == "-p") {
            int model = std::atoi(argv[++i]);
            if (model < 0 || model > 2) {
                std::cerr << "Error: Invalid substitution model (must be 0, 1, or 2)\n";
                return 1;
            }
            config.substitution_model = static_cast<simprot::SubstitutionModel>(model);
        } else if (arg == "-x") {
            config.gamma_alpha = std::atof(argv[++i]);
        } else if (arg == "-g") {
            config.indel_frequency = std::atof(argv[++i]);
        } else if (arg == "-b") {
            int model = std::atoi(argv[++i]);
            if (model == 0) {
                config.indel_distribution = simprot::IndelDistribution::QianGoldstein;
            } else if (model == 1) {
                config.indel_distribution = simprot::IndelDistribution::Benner;
            } else if (model == 2) {
                config.indel_distribution = simprot::IndelDistribution::Custom;
            } else if (model == 3) {
                config.indel_distribution = simprot::IndelDistribution::None;
            } else {
                std::cerr << "Error: Invalid indel model (must be 0, 1, 2, or 3)\n";
                return 1;
            }
        } else if (arg == "-k") {
            config.benner_k = std::atoi(argv[++i]);
        } else if (arg == "-u") {
            config.indel_distribution_file = argv[++i];
            config.indel_distribution = simprot::IndelDistribution::Custom;
        } else if (arg == "-t") {
            config.tree_branch_scale = std::atof(argv[++i]);
        } else if (arg == "-v") {
            config.variable_branch_gamma = std::atof(argv[++i]);
        } else if (arg == "-e") {
            config.branch_extinction_prob = std::atof(argv[++i]);
        } else if (arg == "-S") {
            config.seed = std::atoi(argv[++i]);
        } else {
            std::cerr << "Error: Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    // Validate required arguments
    if (config.tree_file.empty()) {
        std::cerr << "Error: Tree file (-f) is required\n";
        print_usage(argv[0]);
        return 1;
    }

    // Validate configuration
    try {
        config.validate();
    } catch (const std::exception& e) {
        std::cerr << "Configuration error: " << e.what() << "\n";
        return 1;
    }

    // Initialize RNG
    simprot::WichmannHillRNG rng(config.seed);

    // Read tree
    std::unique_ptr<simprot::TreeNode> tree;
    try {
        tree = simprot::read_newick_file(config.tree_file, config, rng);
        std::cerr << "Loaded tree with " << simprot::count_nodes(*tree)
                  << " nodes (" << simprot::count_leaves(*tree) << " leaves)\n";
    } catch (const std::exception& e) {
        std::cerr << "Error reading tree: " << e.what() << "\n";
        return 1;
    }

    // Create evolver
    simprot::SequenceEvolver evolver(config, rng);

    // Set up indel logging if requested
    std::ofstream indel_log;
    if (config.indel_log_file) {
        indel_log.open(*config.indel_log_file);
        if (!indel_log) {
            std::cerr << "Error: Cannot open indel log file: " << *config.indel_log_file << "\n";
            return 1;
        }
        evolver.set_indel_callback([&](double distance, simprot::IndelType type, int length) {
            indel_log << ">distance " << distance << "\n";
            indel_log << (type == simprot::IndelType::Insertion ? "Ins " : "Del ")
                      << length << "\n";
        });
    }

    // Initialize root sequence
    if (config.root_sequence_file) {
        // Read root sequence from file
        std::ifstream seq_file(*config.root_sequence_file);
        if (!seq_file) {
            std::cerr << "Error: Cannot open root sequence file: " << *config.root_sequence_file << "\n";
            return 1;
        }
        std::string line, sequence;
        while (std::getline(seq_file, line)) {
            if (line.empty() || line[0] == '>') continue;
            sequence += line;
        }
        evolver.init_root_sequence(*tree, sequence);
        std::cerr << "Loaded root sequence of length " << sequence.length() << "\n";
    } else {
        evolver.init_root_sequence(*tree);
        std::cerr << "Generated random root sequence of length " << config.root_sequence_length << "\n";
    }

    // Evolve sequences
    std::cerr << "Evolving sequences...\n";
    evolver.evolve(*tree);
    std::cerr << "Evolution complete.\n";

    // Build true alignment
    std::cerr << "Building true alignment...\n";
    simprot::Alignment alignment = simprot::Alignment::from_tree(*tree);
    std::cerr << "Alignment: " << alignment.num_sequences() << " sequences, "
              << alignment.alignment_length() << " columns\n";

    // Write output files
    if (config.alignment_output_file) {
        try {
            // Get leaf sequences only for output
            auto leaf_seqs = alignment.leaf_sequences(*tree);
            simprot::write_fasta_file(*config.alignment_output_file, leaf_seqs);
            std::cerr << "Wrote alignment to: " << *config.alignment_output_file << "\n";
        } catch (const std::exception& e) {
            std::cerr << "Error writing alignment: " << e.what() << "\n";
            return 1;
        }
    }

    if (config.fasta_output_file) {
        try {
            std::ofstream fasta_out(*config.fasta_output_file);
            if (!fasta_out) {
                std::cerr << "Error: Cannot open FASTA output file: " << *config.fasta_output_file << "\n";
                return 1;
            }
            simprot::write_fasta_leaf_sequences(fasta_out, *tree);
            std::cerr << "Wrote sequences to: " << *config.fasta_output_file << "\n";
        } catch (const std::exception& e) {
            std::cerr << "Error writing sequences: " << e.what() << "\n";
            return 1;
        }
    }

    // Print summary if no output files specified
    if (!config.alignment_output_file && !config.fasta_output_file) {
        std::cerr << "\nNo output files specified. Use -a for alignment or -s for sequences.\n";
        std::cerr << "Leaf sequences:\n";
        simprot::write_fasta_leaf_sequences(std::cout, *tree);
    }

    std::cerr << "Done.\n";
    return 0;
}
