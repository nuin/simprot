#ifndef SIMPROT_IO_FASTA_WRITER_HPP
#define SIMPROT_IO_FASTA_WRITER_HPP

/**
 * @file fasta_writer.hpp
 * @brief FASTA format output for sequences and alignments.
 *
 * This header provides functions for writing sequences in FASTA format,
 * both for individual sequences and complete alignments.
 */

#include "simprot/sequence/alignment.hpp"
#include "simprot/tree/tree_node.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace simprot {

/**
 * @brief Write a single sequence in FASTA format.
 * @param os Output stream.
 * @param name Sequence name/identifier.
 * @param sequence The sequence string.
 * @param line_width Maximum characters per line (0 = no wrapping).
 */
void write_fasta_sequence(std::ostream& os,
                          const std::string& name,
                          const std::string& sequence,
                          std::size_t line_width = 60);

/**
 * @brief Write multiple sequences in FASTA format.
 * @param os Output stream.
 * @param sequences Vector of aligned sequences.
 * @param line_width Maximum characters per line (0 = no wrapping).
 */
void write_fasta_sequences(std::ostream& os,
                           const std::vector<AlignedSequence>& sequences,
                           std::size_t line_width = 60);

/**
 * @brief Write an alignment in FASTA format.
 * @param os Output stream.
 * @param alignment The alignment to write.
 * @param line_width Maximum characters per line (0 = no wrapping).
 */
void write_fasta_alignment(std::ostream& os,
                           const Alignment& alignment,
                           std::size_t line_width = 60);

/**
 * @brief Write leaf sequences from a tree in FASTA format (no gaps).
 * @param os Output stream.
 * @param root Root node of the tree.
 * @param line_width Maximum characters per line (0 = no wrapping).
 *
 * This writes only the sequences at leaf nodes, without any gap characters.
 */
void write_fasta_leaf_sequences(std::ostream& os,
                                const TreeNode& root,
                                std::size_t line_width = 60);

/**
 * @brief Write sequences to a FASTA file.
 * @param filename Output file path.
 * @param sequences Vector of aligned sequences.
 * @param line_width Maximum characters per line (0 = no wrapping).
 * @throws std::runtime_error if file cannot be opened.
 */
void write_fasta_file(const std::string& filename,
                      const std::vector<AlignedSequence>& sequences,
                      std::size_t line_width = 60);

/**
 * @brief Write alignment to a FASTA file.
 * @param filename Output file path.
 * @param alignment The alignment to write.
 * @param line_width Maximum characters per line (0 = no wrapping).
 * @throws std::runtime_error if file cannot be opened.
 */
void write_fasta_file(const std::string& filename,
                      const Alignment& alignment,
                      std::size_t line_width = 60);

} // namespace simprot

#endif // SIMPROT_IO_FASTA_WRITER_HPP
