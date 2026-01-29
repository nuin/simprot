#include "simprot/io/fasta_writer.hpp"

#include <fstream>
#include <stdexcept>

namespace simprot {

void write_fasta_sequence(std::ostream& os,
                          const std::string& name,
                          const std::string& sequence,
                          std::size_t line_width) {
    // Write header line
    os << '>' << name << '\n';

    // Write sequence, optionally wrapping lines
    if (line_width == 0 || sequence.length() <= line_width) {
        os << sequence << '\n';
    } else {
        std::size_t pos = 0;
        while (pos < sequence.length()) {
            std::size_t len = std::min(line_width, sequence.length() - pos);
            os << sequence.substr(pos, len) << '\n';
            pos += len;
        }
    }
}

void write_fasta_sequences(std::ostream& os,
                           const std::vector<AlignedSequence>& sequences,
                           std::size_t line_width) {
    for (const auto& seq : sequences) {
        write_fasta_sequence(os, seq.name, seq.sequence, line_width);
    }
}

void write_fasta_alignment(std::ostream& os,
                           const Alignment& alignment,
                           std::size_t line_width) {
    write_fasta_sequences(os, alignment.sequences(), line_width);
}

void write_fasta_leaf_sequences(std::ostream& os,
                                const TreeNode& root,
                                std::size_t line_width) {
    // Recursive function to collect and write leaf sequences
    std::function<void(const TreeNode&)> write_leaves = [&](const TreeNode& node) {
        if (node.is_leaf()) {
            if (!node.sequence.empty()) {
                write_fasta_sequence(os, node.name, node.sequence, line_width);
            }
        } else {
            if (node.left) write_leaves(*node.left);
            if (node.right) write_leaves(*node.right);
        }
    };

    write_leaves(root);
}

void write_fasta_file(const std::string& filename,
                      const std::vector<AlignedSequence>& sequences,
                      std::size_t line_width) {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    write_fasta_sequences(file, sequences, line_width);
}

void write_fasta_file(const std::string& filename,
                      const Alignment& alignment,
                      std::size_t line_width) {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    write_fasta_alignment(file, alignment, line_width);
}

} // namespace simprot
