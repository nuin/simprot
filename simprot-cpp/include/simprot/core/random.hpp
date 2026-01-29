#ifndef SIMPROT_CORE_RANDOM_HPP
#define SIMPROT_CORE_RANDOM_HPP

/**
 * @file random.hpp
 * @brief Wichmann-Hill random number generator for reproducible simulations.
 *
 * This is an exact port of the PAML random number generator used in the original
 * SIMPROT. Using the same seed will produce identical sequences of random numbers,
 * which is essential for reproducibility.
 *
 * Reference:
 * Wichmann BA & Hill ID. 1982. An efficient and portable pseudo-random number
 * generator. Appl. Stat. 31:188-190
 */

namespace simprot {

/**
 * @class WichmannHillRNG
 * @brief Wichmann-Hill pseudo-random number generator.
 *
 * This class encapsulates the three-generator Wichmann-Hill PRNG, ported exactly
 * from the PAML package (via SIMPROT's random.c). The generator combines three
 * linear congruential generators to produce uniform random numbers in [0,1).
 *
 * The generator also includes methods for gamma-distributed random variates,
 * which are used for site-specific rate variation in protein evolution.
 */
class WichmannHillRNG {
public:
    /**
     * @brief Construct a new RNG with the given seed.
     * @param seed Initial seed value.
     */
    explicit WichmannHillRNG(int seed = 12345);

    /**
     * @brief Set the seed for the RNG.
     * @param seed New seed value.
     *
     * After calling this, the RNG state is reset and will produce a deterministic
     * sequence of values based on the seed.
     */
    void set_seed(int seed);

    /**
     * @brief Generate a uniform random number in [0,1).
     * @return A double uniformly distributed in [0,1).
     *
     * This is the core RNG method. The algorithm uses three linear congruential
     * generators whose outputs are combined.
     */
    [[nodiscard]] double uniform();

    /**
     * @brief Generate a gamma-distributed random variate.
     * @param shape The shape parameter (alpha) of the gamma distribution.
     * @return A gamma-distributed random value, scaled by 1/shape.
     *
     * Note: This returns r/s where r ~ Gamma(s,1), which gives the "rate"
     * interpretation used in SIMPROT for site-specific evolutionary rates.
     */
    [[nodiscard]] double gamma(double shape);

private:
    /**
     * @brief Gamma variate for shape < 1.
     * @param s Shape parameter (0 < s < 1).
     * @return Gamma-distributed random value (unscaled).
     */
    [[nodiscard]] double gamma_lt1(double s);

    /**
     * @brief Gamma variate for shape > 1.
     * @param s Shape parameter (s > 1).
     * @return Gamma-distributed random value (unscaled).
     */
    [[nodiscard]] double gamma_gt1(double s);

    // State variables for the three LCGs
    int x_rndu_{11};
    int y_rndu_{23};
    int z_rndu_{137};

    // Cached parameters for gamma_lt1 (shape < 1)
    // Original uses static with ss=10.0 initial value
    double gamma1_cached_shape_{10.0};
    double gamma1_a_{0.0};
    double gamma1_p_{0.0};
    double gamma1_uf_{0.0};
    double gamma1_d_{0.0};

    // Cached parameters for gamma_gt1 (shape > 1)
    // Original uses static with ss=0 initial value
    double gamma2_cached_shape_{0.0};
    double gamma2_b_{0.0};
    double gamma2_h_{0.0};
};

} // namespace simprot

#endif // SIMPROT_CORE_RANDOM_HPP
