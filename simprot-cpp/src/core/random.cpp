#include "simprot/core/random.hpp"

#include <cmath>
#include <stdexcept>

namespace simprot {

WichmannHillRNG::WichmannHillRNG(int seed) {
    set_seed(seed);
}

void WichmannHillRNG::set_seed(int seed) {
    // This exactly matches the SetSeed function from random.c
    z_rndu_ = 170 * (seed % 178) + 137;

    // Reset x and y to their initial values
    x_rndu_ = 11;
    y_rndu_ = 23;

    // Reset gamma caches to match original static initial values
    // gamma1 uses ss=10.0, gamma2 uses ss=0
    gamma1_cached_shape_ = 10.0;
    gamma2_cached_shape_ = 0.0;
}

double WichmannHillRNG::uniform() {
    // This exactly matches the rndu() function from random.c (non-FAST version)
    // The algorithm combines three linear congruential generators.

    x_rndu_ = 171 * (x_rndu_ % 177) -  2 * (x_rndu_ / 177);
    y_rndu_ = 172 * (y_rndu_ % 176) - 35 * (y_rndu_ / 176);
    z_rndu_ = 170 * (z_rndu_ % 178) - 63 * (z_rndu_ / 178);

    if (x_rndu_ < 0) x_rndu_ += 30269;
    if (y_rndu_ < 0) y_rndu_ += 30307;
    if (z_rndu_ < 0) z_rndu_ += 30323;

    double r = static_cast<double>(x_rndu_) / 30269.0
             + static_cast<double>(y_rndu_) / 30307.0
             + static_cast<double>(z_rndu_) / 30323.0;

    return r - static_cast<int>(r);
}

double WichmannHillRNG::gamma(double shape) {
    // This exactly matches rndgamma() from random.c
    double r = 0.0;

    if (shape <= 0.0) {
        throw std::invalid_argument("Gamma shape parameter must be positive");
    } else if (shape < 1.0) {
        r = gamma_lt1(shape);
    } else if (shape > 1.0) {
        r = gamma_gt1(shape);
    } else {
        // shape == 1.0: Gamma(1) is exponential
        r = -std::log(uniform());
    }

    // Return r/shape as in the original
    return r / shape;
}

double WichmannHillRNG::gamma_lt1(double s) {
    // This exactly matches rndgamma1() from random.c
    // Algorithm for shape parameter s < 1

    constexpr double small = 1e-37;
    double x = 0.0;

    // Recompute parameters only if shape changed (optimization from original)
    if (s != gamma1_cached_shape_) {
        gamma1_a_ = 1.0 - s;
        gamma1_p_ = gamma1_a_ / (gamma1_a_ + s * std::exp(-gamma1_a_));
        gamma1_uf_ = gamma1_p_ * std::pow(small / gamma1_a_, s);
        gamma1_d_ = gamma1_a_ * std::log(gamma1_a_);
        gamma1_cached_shape_ = s;
    }

    for (;;) {
        double r = uniform();
        double w;

        if (r > gamma1_p_) {
            x = gamma1_a_ - std::log((1.0 - r) / (1.0 - gamma1_p_));
            w = gamma1_a_ * std::log(x) - gamma1_d_;
        } else if (r > gamma1_uf_) {
            x = gamma1_a_ * std::pow(r / gamma1_p_, 1.0 / s);
            w = x;
        } else {
            return 0.0;
        }

        r = uniform();
        if (1.0 - r <= w && r > 0.0) {
            if (r * (w + 1.0) >= 1.0 || -std::log(r) <= w) {
                continue;
            }
        }
        break;
    }

    return x;
}

double WichmannHillRNG::gamma_gt1(double s) {
    // This exactly matches rndgamma2() from random.c
    // Algorithm for shape parameter s > 1

    double x;

    // Recompute parameters only if shape changed (optimization from original)
    if (s != gamma2_cached_shape_) {
        gamma2_b_ = s - 1.0;
        gamma2_h_ = std::sqrt(3.0 * s - 0.75);
        gamma2_cached_shape_ = s;
    }

    for (;;) {
        double r = uniform();
        double g = r - r * r;
        double f = (r - 0.5) * gamma2_h_ / std::sqrt(g);
        x = gamma2_b_ + f;

        if (x <= 0.0) {
            continue;
        }

        r = uniform();
        double d = 64.0 * r * r * g * g * g;

        if (d * x < x - 2.0 * f * f || std::log(d) < 2.0 * (gamma2_b_ * std::log(x / gamma2_b_) - f)) {
            break;
        }
    }

    return x;
}

} // namespace simprot
