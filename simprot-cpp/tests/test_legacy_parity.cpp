// Regression tests for behaviour that must match SIMPROT 1.04 (legacy/) so
// that the same seed gives the same output. tests/legacy_parity/run.sh checks
// whole runs against the original; these pin the individual rules.

#include <catch2/catch_test_macros.hpp>

#include "simprot/core/random.hpp"
#include "simprot/evolution/indel_model.hpp"
#include "simprot/evolution/matrix_data.hpp"
#include "simprot/evolution/substitution_matrix.hpp"
#include "simprot/sequence/mutable_sequence.hpp"
#include "simprot/tree/tree_parser.hpp"

#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <vector>

using namespace simprot;

TEST_CASE("compute_num_indels draws one uniform per site even at frequency 0",
          "[legacy_parity]") {
    // GetNumIndels loops over every site with f = 0 instead of returning early
    WichmannHillRNG rng(42);
    WichmannHillRNG expected(42);

    REQUIRE(compute_num_indels(0.5, 25, 0.0, 3.0, false, rng) == 0);
    for (int i = 0; i < 25; ++i) {
        (void)expected.uniform();
    }
    REQUIRE(rng.uniform() == expected.uniform());
}

TEST_CASE("Indel lengths stop below the 5% cap", "[legacy_parity]") {
    // InitCumulativeIndelLength fills lengths 1 .. max-1, so a 100-residue
    // sequence (max 5) never gets an indel of length 5
    QianGoldsteinModel qg(3.0);
    BennerModel benner(-2);
    WichmannHillRNG rng(7);

    int qg_max = 0;
    int benner_max = 0;
    for (int i = 0; i < 20000; ++i) {
        int a = qg.sample_length(100, 2.0, 1.0, rng);
        int b = benner.sample_length(100, 2.0, 1.0, rng);
        REQUIRE(a >= 1);
        REQUIRE(b >= 1);
        qg_max = std::max(qg_max, a);
        benner_max = std::max(benner_max, b);
    }
    REQUIRE(qg_max == 4);
    REQUIRE(benner_max == 4);
}

TEST_CASE("Sequences shorter than 40 residues only get length-1 indels",
          "[legacy_parity]") {
    QianGoldsteinModel qg(3.0);
    WichmannHillRNG rng(11);
    for (int i = 0; i < 1000; ++i) {
        REQUIRE(qg.sample_length(39, 0.5, 1.0, rng) == 1);
    }
}

TEST_CASE("A deletion running past the end is moved back", "[legacy_parity]") {
    // MarkPositions: if indelPos + indelLength > length, start at length - indelLength
    MutableSequence seq("ACDEFGHIKL", std::vector<double>(10, 1.0));
    mark_indel_positions(seq, 8, 4, IndelType::Deletion);

    std::vector<int> marks;
    for (std::size_t i = 0; i < seq.size(); ++i) {
        marks.push_back(seq.node_at(i)->mark);
    }
    REQUIRE(marks == std::vector<int>{0, 0, 0, 0, 0, 0, 1, 1, 1, 1});
}

TEST_CASE("Amino acid frequencies use the largest-eigenvalue row",
          "[legacy_parity]") {
    // MakeProtFreqs picks the row of the largest eigenvalue: row 10 for PAM,
    // row 0 for JTT and PMB
    PAMMatrix pam;
    JTTMatrix jtt;
    for (std::size_t i = 0; i < kNumAminoAcids; ++i) {
        REQUIRE(pam.frequencies()[i] == std::abs(matrix_data::pam_eigenvectors[10][i]));
        REQUIRE(jtt.frequencies()[i] == std::abs(matrix_data::jtt_eigenvectors[0][i]));
    }
}

TEST_CASE("Rates are normalised by dividing, then multiplying",
          "[legacy_parity]") {
    const std::vector<double> rates = {0.3, 1.7, 0.24062775927899965, 2.9, 0.05};
    MutableSequence seq("ACDEF", rates);
    seq.normalize_rates();

    double sum = 0.0;
    for (double r : rates) sum += r;
    for (std::size_t i = 0; i < rates.size(); ++i) {
        double expected = rates[i];
        expected /= sum;
        expected *= static_cast<double>(rates.size());
        REQUIRE(seq.node_at(i)->rate == expected);
    }
}

TEST_CASE("Extinction leaves out the left child of an extinct node",
          "[legacy_parity]") {
    // SIMPROT 1.04 with seed 3 and branch extinction 0.3 on this tree makes
    // the internal branches above (A,B), (C,D) and ((A,B),(C,D)) extinct;
    // FlagTree then drops A and C (left children, one level only).
    WichmannHillRNG rng(3);
    NewickParser parser(&rng);
    parser.set_extinction_probability(0.3);
    auto tree = parser.parse(
        "(((A:0.1,B:0.2):0.15,(C:0.3,D:0.05):0.2):0.1,"
        "((E:0.2,F:0.1):0.25,(G:0.15,H:0.3):0.1):0.2);");

    std::map<std::string, bool> omitted;
    std::function<void(const TreeNode&)> walk = [&](const TreeNode& node) {
        if (node.is_leaf()) omitted[node.name] = node.omitted_from_output();
        if (node.left) walk(*node.left);
        if (node.right) walk(*node.right);
    };
    walk(*tree);

    REQUIRE(omitted == std::map<std::string, bool>{
        {"A", true}, {"B", false}, {"C", true}, {"D", false},
        {"E", false}, {"F", false}, {"G", false}, {"H", false}});
}
