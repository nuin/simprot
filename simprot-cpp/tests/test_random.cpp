#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "simprot/core/random.hpp"

using namespace simprot;
using Catch::Approx;

TEST_CASE("WichmannHillRNG basic functionality", "[random]") {
    SECTION("Construction with default seed") {
        WichmannHillRNG rng;
        double val = rng.uniform();
        REQUIRE(val >= 0.0);
        REQUIRE(val < 1.0);
    }

    SECTION("Construction with explicit seed") {
        WichmannHillRNG rng(42);
        double val = rng.uniform();
        REQUIRE(val >= 0.0);
        REQUIRE(val < 1.0);
    }
}

TEST_CASE("WichmannHillRNG reproducibility", "[random]") {
    SECTION("Same seed produces same sequence") {
        WichmannHillRNG rng1(12345);
        WichmannHillRNG rng2(12345);

        for (int i = 0; i < 100; ++i) {
            REQUIRE(rng1.uniform() == rng2.uniform());
        }
    }

    SECTION("Different seeds produce different sequences") {
        WichmannHillRNG rng1(12345);
        WichmannHillRNG rng2(54321);

        bool all_same = true;
        for (int i = 0; i < 100; ++i) {
            if (rng1.uniform() != rng2.uniform()) {
                all_same = false;
                break;
            }
        }
        REQUIRE_FALSE(all_same);
    }

    SECTION("set_seed resets state") {
        WichmannHillRNG rng(12345);

        // Generate some values
        double first_val = rng.uniform();
        (void)rng.uniform();  // Advance RNG state
        (void)rng.uniform();  // Advance RNG state

        // Reset with same seed
        rng.set_seed(12345);

        // Should get the same first value
        REQUIRE(rng.uniform() == first_val);
    }
}

TEST_CASE("WichmannHillRNG known values", "[random]") {
    // These values were generated from the original C implementation
    // with seed 12345 to verify exact compatibility
    SECTION("First 10 values with seed 12345") {
        WichmannHillRNG rng(12345);

        // Expected values from original implementation:
        // We'll need to capture these from the original code
        // For now, just verify the values are in range and deterministic
        std::vector<double> values;
        for (int i = 0; i < 10; ++i) {
            values.push_back(rng.uniform());
        }

        // Verify all in range
        for (double v : values) {
            REQUIRE(v >= 0.0);
            REQUIRE(v < 1.0);
        }

        // Verify reproducibility by re-running
        rng.set_seed(12345);
        for (int i = 0; i < 10; ++i) {
            REQUIRE(rng.uniform() == values[i]);
        }
    }
}

TEST_CASE("WichmannHillRNG gamma distribution", "[random]") {
    SECTION("Gamma with shape 1 (exponential)") {
        WichmannHillRNG rng(12345);
        double val = rng.gamma(1.0);
        REQUIRE(val >= 0.0);
    }

    SECTION("Gamma with shape < 1") {
        WichmannHillRNG rng(12345);
        double val = rng.gamma(0.5);
        REQUIRE(val >= 0.0);
    }

    SECTION("Gamma with shape > 1") {
        WichmannHillRNG rng(12345);
        double val = rng.gamma(2.0);
        REQUIRE(val >= 0.0);
    }

    SECTION("Gamma reproducibility") {
        WichmannHillRNG rng1(42);
        WichmannHillRNG rng2(42);

        for (int i = 0; i < 100; ++i) {
            REQUIRE(rng1.gamma(1.5) == rng2.gamma(1.5));
        }
    }

    SECTION("Gamma mean is approximately 1 (since we return r/shape)") {
        WichmannHillRNG rng(12345);
        double sum = 0.0;
        constexpr int n = 10000;

        for (int i = 0; i < n; ++i) {
            sum += rng.gamma(2.0);
        }

        double mean = sum / n;
        // Mean of Gamma(a,1)/a should be 1
        REQUIRE(mean == Approx(1.0).margin(0.05));
    }

    SECTION("Invalid shape throws") {
        WichmannHillRNG rng(12345);
        REQUIRE_THROWS_AS(rng.gamma(0.0), std::invalid_argument);
        REQUIRE_THROWS_AS(rng.gamma(-1.0), std::invalid_argument);
    }
}
