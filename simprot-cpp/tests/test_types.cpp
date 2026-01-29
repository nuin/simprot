#include <catch2/catch_test_macros.hpp>

#include "simprot/core/types.hpp"

using namespace simprot;

TEST_CASE("Amino acid character array", "[types]") {
    SECTION("Contains 20 standard amino acids") {
        REQUIRE(kAminoAcidChars.size() == 20);
    }

    SECTION("First amino acid is Alanine") {
        REQUIRE(kAminoAcidChars[0] == 'A');
    }

    SECTION("Array matches PAML order") {
        // A R N D C Q E G H I L K M F P S T W Y V
        REQUIRE(kAminoAcidChars[0] == 'A');
        REQUIRE(kAminoAcidChars[1] == 'R');
        REQUIRE(kAminoAcidChars[2] == 'N');
        REQUIRE(kAminoAcidChars[3] == 'D');
        REQUIRE(kAminoAcidChars[4] == 'C');
        REQUIRE(kAminoAcidChars[5] == 'Q');
        REQUIRE(kAminoAcidChars[6] == 'E');
        REQUIRE(kAminoAcidChars[7] == 'G');
        REQUIRE(kAminoAcidChars[8] == 'H');
        REQUIRE(kAminoAcidChars[9] == 'I');
        REQUIRE(kAminoAcidChars[10] == 'L');
        REQUIRE(kAminoAcidChars[11] == 'K');
        REQUIRE(kAminoAcidChars[12] == 'M');
        REQUIRE(kAminoAcidChars[13] == 'F');
        REQUIRE(kAminoAcidChars[14] == 'P');
        REQUIRE(kAminoAcidChars[15] == 'S');
        REQUIRE(kAminoAcidChars[16] == 'T');
        REQUIRE(kAminoAcidChars[17] == 'W');
        REQUIRE(kAminoAcidChars[18] == 'Y');
        REQUIRE(kAminoAcidChars[19] == 'V');
    }
}

TEST_CASE("amino_acid_to_char conversion", "[types]") {
    SECTION("Valid indices return correct characters") {
        REQUIRE(amino_acid_to_char(0) == 'A');
        REQUIRE(amino_acid_to_char(1) == 'R');
        REQUIRE(amino_acid_to_char(19) == 'V');
    }

    SECTION("Out of range index throws") {
        REQUIRE_THROWS_AS(amino_acid_to_char(20), std::out_of_range);
        REQUIRE_THROWS_AS(amino_acid_to_char(255), std::out_of_range);
    }
}

TEST_CASE("char_to_amino_acid conversion", "[types]") {
    SECTION("Valid uppercase characters convert correctly") {
        REQUIRE(char_to_amino_acid('A') == 0);
        REQUIRE(char_to_amino_acid('R') == 1);
        REQUIRE(char_to_amino_acid('N') == 2);
        REQUIRE(char_to_amino_acid('D') == 3);
        REQUIRE(char_to_amino_acid('C') == 4);
        REQUIRE(char_to_amino_acid('V') == 19);
    }

    SECTION("Lowercase characters convert correctly") {
        REQUIRE(char_to_amino_acid('a') == 0);
        REQUIRE(char_to_amino_acid('r') == 1);
        REQUIRE(char_to_amino_acid('v') == 19);
    }

    SECTION("Invalid characters throw") {
        REQUIRE_THROWS_AS(char_to_amino_acid('B'), std::invalid_argument);
        REQUIRE_THROWS_AS(char_to_amino_acid('J'), std::invalid_argument);
        REQUIRE_THROWS_AS(char_to_amino_acid('O'), std::invalid_argument);
        REQUIRE_THROWS_AS(char_to_amino_acid('U'), std::invalid_argument);
        REQUIRE_THROWS_AS(char_to_amino_acid('X'), std::invalid_argument);
        REQUIRE_THROWS_AS(char_to_amino_acid('Z'), std::invalid_argument);
        REQUIRE_THROWS_AS(char_to_amino_acid('1'), std::invalid_argument);
        REQUIRE_THROWS_AS(char_to_amino_acid('-'), std::invalid_argument);
    }

    SECTION("Round-trip conversion works") {
        for (AminoAcidIndex i = 0; i < kNumAminoAcids; ++i) {
            char c = amino_acid_to_char(i);
            REQUIRE(char_to_amino_acid(c) == i);
        }
    }
}

TEST_CASE("SubstitutionModel enum", "[types]") {
    SECTION("Enum values match CLI options") {
        REQUIRE(static_cast<int>(SubstitutionModel::PAM) == 0);
        REQUIRE(static_cast<int>(SubstitutionModel::JTT) == 1);
        REQUIRE(static_cast<int>(SubstitutionModel::PMB) == 2);
    }

    SECTION("to_string returns correct names") {
        REQUIRE(to_string(SubstitutionModel::PAM) == "PAM");
        REQUIRE(to_string(SubstitutionModel::JTT) == "JTT");
        REQUIRE(to_string(SubstitutionModel::PMB) == "PMB");
    }
}

TEST_CASE("IndelDistribution enum", "[types]") {
    SECTION("Enum values are sequential") {
        REQUIRE(static_cast<int>(IndelDistribution::QianGoldstein) == 0);
        REQUIRE(static_cast<int>(IndelDistribution::Benner) == 1);
        REQUIRE(static_cast<int>(IndelDistribution::Custom) == 2);
        REQUIRE(static_cast<int>(IndelDistribution::None) == 3);
    }

    SECTION("to_string returns correct names") {
        REQUIRE(to_string(IndelDistribution::QianGoldstein) == "Qian-Goldstein");
        REQUIRE(to_string(IndelDistribution::Benner) == "Benner");
        REQUIRE(to_string(IndelDistribution::Custom) == "Custom");
        REQUIRE(to_string(IndelDistribution::None) == "None");
    }
}

TEST_CASE("IndelType enum", "[types]") {
    SECTION("Deletion is 0, Insertion is 1") {
        REQUIRE(static_cast<int>(IndelType::Deletion) == 0);
        REQUIRE(static_cast<int>(IndelType::Insertion) == 1);
    }
}

TEST_CASE("Constants", "[types]") {
    SECTION("Gap character is dash") {
        REQUIRE(kGapChar == '-');
    }

    SECTION("Number of amino acids is 20") {
        REQUIRE(kNumAminoAcids == 20);
    }
}
