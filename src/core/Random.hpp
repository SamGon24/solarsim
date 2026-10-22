// =============================================================================
//  SolarSim - Random.hpp
//  Deterministic pseudo random numbers: SplitMix64 for seeding and
//  xoshiro256** (Blackman & Vigna, 2018) as the main generator. Deterministic
//  output across compilers/platforms (unlike std::normal_distribution etc.).
// =============================================================================
#pragma once

#include "Constants.hpp"

#include <cmath>
#include <cstdint>

namespace solsim {

class SplitMix64 {
public:
    explicit SplitMix64(std::uint64_t seed) : m_state(seed) {}
    std::uint64_t next() {
        std::uint64_t z = (m_state += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

private:
    std::uint64_t m_state;
};

class Random {
public:
    explicit Random(std::uint64_t seed = 0x5EED5017A25ULL) { reseed(seed); }

    void reseed(std::uint64_t seed) {
        SplitMix64 sm(seed);
        for (auto& s : m_s) s = sm.next();
        m_hasSpare = false;
    }

    std::uint64_t nextU64() {
        const std::uint64_t result = rotl(m_s[1] * 5, 7) * 9;
        const std::uint64_t t = m_s[1] << 17;
        m_s[2] ^= m_s[0];
        m_s[3] ^= m_s[1];
        m_s[1] ^= m_s[2];
        m_s[0] ^= m_s[3];
        m_s[2] ^= t;
        m_s[3] = rotl(m_s[3], 45);
        return result;
    }

    // Uniform double in [0, 1) with 53 bits of randomness.
    double uniform() { return double(nextU64() >> 11) * (1.0 / 9007199254740992.0); }
    double uniform(double lo, double hi) { return lo + (hi - lo) * uniform(); }

    // Uniform integer in [lo, hi] (inclusive).
    int uniformInt(int lo, int hi) {
        const std::uint64_t span = std::uint64_t(hi - lo) + 1ULL;
        return lo + int(nextU64() % span);
    }

    // Standard normal via the Marsaglia polar method.
    double normal() {
        if (m_hasSpare) {
            m_hasSpare = false;
            return m_spare;
        }
        double u, v, s;
        do {
            u = uniform(-1.0, 1.0);
            v = uniform(-1.0, 1.0);
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);
        const double f = std::sqrt(-2.0 * std::log(s) / s);
        m_spare = v * f;
        m_hasSpare = true;
        return u * f;
    }
    double normal(double mean, double sigma) { return mean + sigma * normal(); }

    // Rayleigh distribution - the natural distribution of small orbital
    // eccentricities and inclinations in a dynamically relaxed population.
    double rayleigh(double sigma) {
        double u = uniform();
        if (u <= 0.0) u = 1e-300;
        return sigma * std::sqrt(-2.0 * std::log(u));
    }

    double angle() { return uniform(0.0, constants::TWO_PI); }

private:
    static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

    std::uint64_t m_s[4]{};
    double m_spare = 0.0;
    bool m_hasSpare = false;
};

} // namespace solsim
