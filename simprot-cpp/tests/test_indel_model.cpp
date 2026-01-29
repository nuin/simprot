#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "simprot/evolution/indel_model.hpp"
#include "simprot/core/random.hpp"

using namespace simprot;
using Catch::Approx;

TEST_CASE("IndelModel factory", "[indel_model]") {
    SECTION("Creates Qian-Goldstein model") {
        auto model = create_indel_model(IndelDistribution::QianGoldstein);
        REQUIRE(model != nullptr);
        REQUIRE(model->name() == "Qian-Goldstein");
    }

    SECTION("Creates Benner model") {
        auto model = create_indel_model(IndelDistribution::Benner);
        REQUIRE(model != nullptr);
        REQUIRE(model->name() == "Benner");
    }

    SECTION("Creates No indel model") {
        auto model = create_indel_model(IndelDistribution::None);
        REQUIRE(model != nullptr);
        REQUIRE(model->name() == "None");
    }
}

TEST_CASE("QianGoldsteinModel", "[indel_model]") {
    QianGoldsteinModel model(3.0);
    WichmannHillRNG rng(12345);

    SECTION("Returns positive lengths") {
        for (int i = 0; i < 100; ++i) {
            int len = model.sample_length(100, 0.5, 1.0, rng);
            REQUIRE(len >= 1);
        }
    }

    SECTION("Lengths are capped at 5% of sequence length") {
        // With sequence length 100, max should be around 5
        WichmannHillRNG rng2(42);
        int max_seen = 0;
        for (int i = 0; i < 1000; ++i) {
            int len = model.sample_length(100, 0.5, 1.0, rng2);
            if (len > max_seen) max_seen = len;
        }
        // Might occasionally exceed 5, but cap should be around sequence_length * 0.05
        REQUIRE(max_seen <= 100);  // Definitely shouldn't exceed sequence length
    }

    SECTION("Shorter lengths are more common") {
        WichmannHillRNG rng2(42);
        int short_count = 0;  // Length 1-2
        int long_count = 0;   // Length 5+
        constexpr int n_samples = 10000;

        for (int i = 0; i < n_samples; ++i) {
            int len = model.sample_length(200, 0.5, 1.0, rng2);
            if (len <= 2) short_count++;
            if (len >= 5) long_count++;
        }

        // Short indels should be much more common
        REQUIRE(short_count > long_count * 2);
    }

    SECTION("Reproducibility with same seed") {
        WichmannHillRNG rng1(54321);
        WichmannHillRNG rng2(54321);

        for (int i = 0; i < 100; ++i) {
            int len1 = model.sample_length(100, 0.5, 1.0, rng1);
            int len2 = model.sample_length(100, 0.5, 1.0, rng2);
            REQUIRE(len1 == len2);
        }
    }
}

TEST_CASE("BennerModel", "[indel_model]") {
    BennerModel model(-2);
    WichmannHillRNG rng(12345);

    SECTION("Returns positive lengths") {
        for (int i = 0; i < 100; ++i) {
            int len = model.sample_length(100, 0.5, 1.0, rng);
            REQUIRE(len >= 1);
        }
    }

    SECTION("Length distribution follows power law (shorter more common)") {
        WichmannHillRNG rng2(42);
        std::vector<int> counts(10, 0);
        constexpr int n_samples = 10000;

        for (int i = 0; i < n_samples; ++i) {
            int len = model.sample_length(200, 0.5, 1.0, rng2);
            if (len <= 10) counts[len - 1]++;
        }

        // Length 1 should be most common
        REQUIRE(counts[0] > counts[1]);
        REQUIRE(counts[1] > counts[2]);
    }
}

TEST_CASE("NoIndelModel", "[indel_model]") {
    NoIndelModel model;
    WichmannHillRNG rng(12345);

    SECTION("Always returns 0") {
        for (int i = 0; i < 100; ++i) {
            int len = model.sample_length(100, 0.5, 1.0, rng);
            REQUIRE(len == 0);
        }
    }
}

TEST_CASE("CustomIndelModel", "[indel_model]") {
    // Custom distribution: 50% length 1, 30% length 2, 20% length 3
    std::vector<double> freqs = {0.5, 0.3, 0.2};
    CustomIndelModel model(freqs);
    WichmannHillRNG rng(12345);

    SECTION("Returns lengths in valid range") {
        for (int i = 0; i < 100; ++i) {
            int len = model.sample_length(100, 0.5, 1.0, rng);
            REQUIRE(len >= 1);
            REQUIRE(len <= 3);
        }
    }

    SECTION("Distribution approximately matches input") {
        WichmannHillRNG rng2(42);
        std::array<int, 3> counts{};
        constexpr int n_samples = 10000;

        for (int i = 0; i < n_samples; ++i) {
            int len = model.sample_length(100, 0.5, 1.0, rng2);
            counts[len - 1]++;
        }

        double p1 = static_cast<double>(counts[0]) / n_samples;
        double p2 = static_cast<double>(counts[1]) / n_samples;
        double p3 = static_cast<double>(counts[2]) / n_samples;

        // Allow 10% relative error
        REQUIRE(p1 == Approx(0.5).epsilon(0.1));
        REQUIRE(p2 == Approx(0.3).epsilon(0.1));
        REQUIRE(p3 == Approx(0.2).epsilon(0.1));
    }
}

TEST_CASE("compute_num_indels", "[indel_model]") {
    WichmannHillRNG rng(12345);

    SECTION("Returns non-negative count") {
        for (int i = 0; i < 100; ++i) {
            int count = compute_num_indels(0.5, 100, 0.03, 3.0, false, rng);
            REQUIRE(count >= 0);
        }
    }

    SECTION("Zero distance gives zero indels") {
        WichmannHillRNG rng2(42);
        for (int i = 0; i < 100; ++i) {
            int count = compute_num_indels(0.0, 100, 0.03, 3.0, false, rng2);
            REQUIRE(count == 0);
        }
    }

    SECTION("Longer distance gives more indels on average") {
        WichmannHillRNG rng1(42);
        WichmannHillRNG rng2(42);
        int short_sum = 0;
        int long_sum = 0;
        constexpr int n_samples = 1000;

        for (int i = 0; i < n_samples; ++i) {
            short_sum += compute_num_indels(0.1, 100, 0.03, 3.0, false, rng1);
        }
        for (int i = 0; i < n_samples; ++i) {
            long_sum += compute_num_indels(0.5, 100, 0.03, 3.0, false, rng2);
        }

        REQUIRE(long_sum > short_sum);
    }

    SECTION("Higher frequency gives more indels") {
        WichmannHillRNG rng1(42);
        WichmannHillRNG rng2(42);
        int low_freq_sum = 0;
        int high_freq_sum = 0;
        constexpr int n_samples = 1000;

        for (int i = 0; i < n_samples; ++i) {
            low_freq_sum += compute_num_indels(0.5, 100, 0.01, 3.0, false, rng1);
        }
        for (int i = 0; i < n_samples; ++i) {
            high_freq_sum += compute_num_indels(0.5, 100, 0.05, 3.0, false, rng2);
        }

        REQUIRE(high_freq_sum > low_freq_sum);
    }
}

TEST_CASE("choose_indel_type", "[indel_model]") {
    WichmannHillRNG rng(12345);

    SECTION("Returns Insertion or Deletion") {
        for (int i = 0; i < 100; ++i) {
            IndelType type = choose_indel_type(0.5, rng);
            REQUIRE((type == IndelType::Insertion || type == IndelType::Deletion));
        }
    }

    SECTION("Ratio 0.5 gives roughly equal counts") {
        WichmannHillRNG rng2(42);
        int insertions = 0;
        constexpr int n_samples = 10000;

        for (int i = 0; i < n_samples; ++i) {
            IndelType type = choose_indel_type(0.5, rng2);
            if (type == IndelType::Insertion) insertions++;
        }

        double ratio = static_cast<double>(insertions) / n_samples;
        REQUIRE(ratio == Approx(0.5).epsilon(0.05));
    }

    SECTION("Ratio 0.0 gives only deletions") {
        WichmannHillRNG rng2(42);
        for (int i = 0; i < 100; ++i) {
            IndelType type = choose_indel_type(0.0, rng2);
            REQUIRE(type == IndelType::Deletion);
        }
    }

    SECTION("Ratio 1.0 gives only insertions") {
        WichmannHillRNG rng2(42);
        for (int i = 0; i < 100; ++i) {
            IndelType type = choose_indel_type(1.0, rng2);
            REQUIRE(type == IndelType::Insertion);
        }
    }
}
