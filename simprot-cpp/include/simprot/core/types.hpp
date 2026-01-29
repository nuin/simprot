#ifndef SIMPROT_CORE_TYPES_HPP
#define SIMPROT_CORE_TYPES_HPP

/**
 * @file types.hpp
 * @brief Core type definitions for SIMPROT.
 *
 * This header defines the fundamental types used throughout the SIMPROT
 * simulation, including amino acid representations, substitution models,
 * and indel types.
 */

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string_view>

namespace simprot {

/**
 * @brief Amino acid index type (0-19 for standard amino acids).
 */
using AminoAcidIndex = std::uint8_t;

/**
 * @brief Number of standard amino acids.
 */
inline constexpr std::size_t kNumAminoAcids = 20;

/**
 * @brief Standard amino acid one-letter codes in canonical order.
 *
 * This order matches the original SIMPROT/PAML convention:
 * A R N D C Q E G H I L K M F P S T W Y V
 */
inline constexpr std::array<char, kNumAminoAcids> kAminoAcidChars = {
    'A', 'R', 'N', 'D', 'C',
    'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P',
    'S', 'T', 'W', 'Y', 'V'
};

/**
 * @brief Convert amino acid index to one-letter code.
 * @param index Amino acid index (0-19).
 * @return One-letter amino acid code.
 * @throws std::out_of_range if index >= 20.
 */
[[nodiscard]] constexpr char amino_acid_to_char(AminoAcidIndex index) {
    if (index >= kNumAminoAcids) {
        throw std::out_of_range("Amino acid index out of range");
    }
    return kAminoAcidChars[index];
}

/**
 * @brief Convert one-letter code to amino acid index.
 * @param c One-letter amino acid code (case-insensitive).
 * @return Amino acid index (0-19).
 * @throws std::invalid_argument if not a valid amino acid.
 *
 * This function creates a lookup table at compile time for O(1) conversion.
 */
[[nodiscard]] constexpr AminoAcidIndex char_to_amino_acid(char c) {
    // Convert to uppercase
    if (c >= 'a' && c <= 'z') {
        c = static_cast<char>(c - 'a' + 'A');
    }

    // Lookup table: maps 'A'-'Z' to amino acid indices (255 = invalid)
    constexpr std::array<AminoAcidIndex, 26> lookup = {
        0,   // A -> 0
        255, // B -> invalid
        4,   // C -> 4
        3,   // D -> 3
        6,   // E -> 6
        13,  // F -> 13
        7,   // G -> 7
        8,   // H -> 8
        9,   // I -> 9
        255, // J -> invalid
        11,  // K -> 11
        10,  // L -> 10
        12,  // M -> 12
        2,   // N -> 2
        255, // O -> invalid
        14,  // P -> 14
        5,   // Q -> 5
        1,   // R -> 1
        15,  // S -> 15
        16,  // T -> 16
        255, // U -> invalid
        19,  // V -> 19
        17,  // W -> 17
        255, // X -> invalid
        18,  // Y -> 18
        255  // Z -> invalid
    };

    if (c < 'A' || c > 'Z') {
        throw std::invalid_argument("Invalid amino acid character");
    }

    AminoAcidIndex idx = lookup[static_cast<std::size_t>(c - 'A')];
    if (idx == 255) {
        throw std::invalid_argument("Invalid amino acid character");
    }
    return idx;
}

/**
 * @brief Substitution matrix model selection.
 *
 * These correspond to the `-p` command-line option in the original SIMPROT.
 */
enum class SubstitutionModel : std::uint8_t {
    PAM = 0,  ///< PAM (Dayhoff) matrix
    JTT = 1,  ///< JTT (Jones-Taylor-Thornton) matrix
    PMB = 2   ///< PMB (Probability Matrix from Blocks) matrix - default
};

/**
 * @brief Convert SubstitutionModel to string for display.
 */
[[nodiscard]] constexpr std::string_view to_string(SubstitutionModel model) {
    switch (model) {
        case SubstitutionModel::PAM: return "PAM";
        case SubstitutionModel::JTT: return "JTT";
        case SubstitutionModel::PMB: return "PMB";
    }
    return "Unknown";
}

/**
 * @brief Indel length distribution model.
 *
 * Controls how indel lengths are sampled during simulation.
 */
enum class IndelDistribution : std::uint8_t {
    QianGoldstein = 0,  ///< 4-term exponential (Qian & Goldstein) - default
    Benner = 1,         ///< Zipfian power-law (i^k)
    Custom = 2,         ///< User-provided distribution from file
    None = 3            ///< No indels (substitution only)
};

/**
 * @brief Convert IndelDistribution to string for display.
 */
[[nodiscard]] constexpr std::string_view to_string(IndelDistribution dist) {
    switch (dist) {
        case IndelDistribution::QianGoldstein: return "Qian-Goldstein";
        case IndelDistribution::Benner: return "Benner";
        case IndelDistribution::Custom: return "Custom";
        case IndelDistribution::None: return "None";
    }
    return "Unknown";
}

/**
 * @brief Type of insertion/deletion event.
 */
enum class IndelType : std::uint8_t {
    Deletion = 0,   ///< Deletion event (removes residues)
    Insertion = 1   ///< Insertion event (adds residues)
};

/**
 * @brief Gap character used in alignments.
 */
inline constexpr char kGapChar = '-';

/**
 * @brief Child direction in tree traversal.
 */
enum class ChildDirection : std::uint8_t {
    Left = 0,
    Right = 1
};

} // namespace simprot

#endif // SIMPROT_CORE_TYPES_HPP
