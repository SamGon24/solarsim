// =============================================================================
//  SolarSim - Noise.hpp
//  Deterministic 3D gradient-free value noise and fractal Brownian motion,
//  evaluated on the unit sphere to texture planets without UV seams.
// =============================================================================
#pragma once

#include "../core/Vec3.hpp"

#include <cmath>
#include <cstdint>

namespace solsim::render {

class ValueNoise {
public:
    explicit ValueNoise(std::uint32_t seed = 1) : m_seed(seed) {}

    // Smooth noise in [0, 1].
    double noise(double x, double y, double z) const {
        const double fx = std::floor(x), fy = std::floor(y), fz = std::floor(z);
        const int ix = static_cast<int>(fx), iy = static_cast<int>(fy), iz = static_cast<int>(fz);
        const double tx = fade(x - fx), ty = fade(y - fy), tz = fade(z - fz);
        auto v = [&](int dx, int dy, int dz) { return hash(ix + dx, iy + dy, iz + dz); };
        const double x00 = lerp(v(0, 0, 0), v(1, 0, 0), tx);
        const double x10 = lerp(v(0, 1, 0), v(1, 1, 0), tx);
        const double x01 = lerp(v(0, 0, 1), v(1, 0, 1), tx);
        const double x11 = lerp(v(0, 1, 1), v(1, 1, 1), tx);
        return lerp(lerp(x00, x10, ty), lerp(x01, x11, ty), tz);
    }

    double noise(const Vec3& p) const { return noise(p.x, p.y, p.z); }

    // Fractal sum of `octaves` noise layers, normalised to [0, 1].
    double fbm(const Vec3& p, int octaves = 4, double lacunarity = 2.0, double gain = 0.5) const {
        double sum = 0.0, amp = 1.0, norm = 0.0;
        Vec3 q = p;
        for (int o = 0; o < octaves; ++o) {
            sum += amp * noise(q);
            norm += amp;
            amp *= gain;
            q = q * lacunarity + Vec3{17.31, -9.7, 4.13};
        }
        return sum / norm;
    }

    // Ridged multifractal: sharp crests (ice ridges, crater rims).
    double ridged(const Vec3& p, int octaves = 4) const {
        double sum = 0.0, amp = 0.5, norm = 0.0;
        Vec3 q = p;
        for (int o = 0; o < octaves; ++o) {
            const double n = 1.0 - std::fabs(2.0 * noise(q) - 1.0);
            sum += amp * n * n;
            norm += amp;
            amp *= 0.5;
            q = q * 2.1 + Vec3{3.7, 1.9, -6.2};
        }
        return sum / norm;
    }

private:
    static double fade(double t) { return t * t * t * (t * (t * 6.0 - 15.0) + 10.0); }
    static double lerp(double a, double b, double t) { return a + (b - a) * t; }

    double hash(int x, int y, int z) const {
        std::uint32_t h = m_seed;
        h ^= static_cast<std::uint32_t>(x) * 0x8DA6B343u;
        h ^= static_cast<std::uint32_t>(y) * 0xD8163841u;
        h ^= static_cast<std::uint32_t>(z) * 0xCB1AB31Fu;
        h ^= h >> 15;
        h *= 0x2C1B3C6Du;
        h ^= h >> 12;
        h *= 0x297A2D39u;
        h ^= h >> 15;
        return static_cast<double>(h & 0xFFFFFFu) / static_cast<double>(0xFFFFFFu);
    }

    std::uint32_t m_seed;
};

} // namespace solsim::render
