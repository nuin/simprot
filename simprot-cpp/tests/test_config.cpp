#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "simprot/core/config.hpp"

using namespace simprot;
using Catch::Approx;

TEST_CASE("SimulationConfig defaults", "[config]") {
    SimulationConfig config;

    SECTION("Default seed") {
        REQUIRE(config.seed == 12345);
    }

    SECTION("Default sequence length") {
        REQUIRE(config.root_sequence_length == 50);
    }

    SECTION("Default substitution model is PMB") {
        REQUIRE(config.substitution_model == SubstitutionModel::PMB);
    }

    SECTION("Default gamma alpha") {
        REQUIRE(config.gamma_alpha == Approx(1.0));
    }

    SECTION("Default indel frequency") {
        REQUIRE(config.indel_frequency == Approx(0.03));
    }

    SECTION("Default indel distribution is Qian-Goldstein") {
        REQUIRE(config.indel_distribution == IndelDistribution::QianGoldstein);
    }

    SECTION("Default indel ratio is 0.5") {
        REQUIRE(config.indel_ratio == Approx(0.5));
    }

    SECTION("Default tree branch scale") {
        REQUIRE(config.tree_branch_scale == Approx(1.0));
    }

    SECTION("No root sequence by default") {
        REQUIRE_FALSE(config.root_sequence.has_value());
    }

    SECTION("No output files by default") {
        REQUIRE_FALSE(config.fasta_output_file.has_value());
        REQUIRE_FALSE(config.alignment_output_file.has_value());
        REQUIRE_FALSE(config.phylip_output_file.has_value());
    }
}

TEST_CASE("SimulationConfig modification", "[config]") {
    SimulationConfig config;

    SECTION("Can change seed") {
        config.seed = 42;
        REQUIRE(config.seed == 42);
    }

    SECTION("Can change substitution model") {
        config.substitution_model = SubstitutionModel::JTT;
        REQUIRE(config.substitution_model == SubstitutionModel::JTT);
    }

    SECTION("Can set root sequence") {
        config.root_sequence = "ACDEFGHIKLMNPQRSTVWY";
        REQUIRE(config.root_sequence.has_value());
        REQUIRE(config.root_sequence.value() == "ACDEFGHIKLMNPQRSTVWY");
    }

    SECTION("Can set output files") {
        config.fasta_output_file = "output.fasta";
        config.alignment_output_file = "alignment.fasta";
        REQUIRE(config.fasta_output_file.value() == "output.fasta");
        REQUIRE(config.alignment_output_file.value() == "alignment.fasta");
    }

    SECTION("Can set tree file") {
        config.tree_file = "tree.txt";
        REQUIRE(config.tree_file == "tree.txt");
    }
}

TEST_CASE("SimulationConfig validation", "[config]") {
    SimulationConfig config;
    config.tree_file = "test.tree";  // Required field

    SECTION("Valid default config passes") {
        REQUIRE_NOTHROW(config.validate());
    }

    SECTION("Missing tree file throws") {
        SimulationConfig bad_config;
        bad_config.tree_file = "";
        REQUIRE_THROWS_AS(bad_config.validate(), std::invalid_argument);
    }

    SECTION("Negative root sequence length throws") {
        config.root_sequence_length = -1;
        REQUIRE_THROWS_AS(config.validate(), std::invalid_argument);
    }

    SECTION("Zero root sequence length throws") {
        config.root_sequence_length = 0;
        REQUIRE_THROWS_AS(config.validate(), std::invalid_argument);
    }

    SECTION("Negative indel frequency throws") {
        config.indel_frequency = -0.01;
        REQUIRE_THROWS_AS(config.validate(), std::invalid_argument);
    }

    SECTION("Indel frequency > 1 throws") {
        config.indel_frequency = 1.5;
        REQUIRE_THROWS_AS(config.validate(), std::invalid_argument);
    }

    SECTION("Indel ratio < 0 throws") {
        config.indel_ratio = -0.1;
        REQUIRE_THROWS_AS(config.validate(), std::invalid_argument);
    }

    SECTION("Indel ratio > 1 throws") {
        config.indel_ratio = 1.1;
        REQUIRE_THROWS_AS(config.validate(), std::invalid_argument);
    }

    SECTION("Valid edge cases") {
        config.indel_frequency = 0.0;  // No indels is valid
        REQUIRE_NOTHROW(config.validate());

        config.indel_frequency = 0.99;  // Just under 1.0 is valid
        REQUIRE_NOTHROW(config.validate());

        config.indel_ratio = 0.0;  // All deletions
        REQUIRE_NOTHROW(config.validate());

        config.indel_ratio = 1.0;  // All insertions
        REQUIRE_NOTHROW(config.validate());
    }
}

TEST_CASE("SimulationConfig copy semantics", "[config]") {
    SimulationConfig config1;
    config1.seed = 42;
    config1.root_sequence_length = 100;
    config1.tree_file = "test.tree";
    config1.fasta_output_file = "output.fasta";

    SECTION("Copy constructor") {
        SimulationConfig config2 = config1;

        REQUIRE(config2.seed == 42);
        REQUIRE(config2.root_sequence_length == 100);
        REQUIRE(config2.tree_file == "test.tree");
        REQUIRE(config2.fasta_output_file.value() == "output.fasta");
    }

    SECTION("Copy assignment") {
        SimulationConfig config2;
        config2 = config1;

        REQUIRE(config2.seed == 42);
        REQUIRE(config2.root_sequence_length == 100);
    }

    SECTION("Copies are independent") {
        SimulationConfig config2 = config1;
        config2.seed = 999;
        config2.root_sequence_length = 200;

        REQUIRE(config1.seed == 42);
        REQUIRE(config1.root_sequence_length == 100);
    }
}
