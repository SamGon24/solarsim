// =============================================================================
//  SolarSim - StarField.cpp
// =============================================================================
#include "StarField.hpp"

#include "../astro/Coordinates.hpp"
#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"
#include "../core/Random.hpp"

#include <cmath>

namespace solsim::render {

namespace {

struct CatalogStar {
    const char* name;
    double raDeg, decDeg, vmag, bv;
};

// Bright stars (J2000.0). Values rounded from the Yale Bright Star Catalogue.
const CatalogStar kBright[] = {
    {"Sirius", 101.287, -16.716, -1.46, 0.00},     {"Canopus", 95.988, -52.696, -0.74, 0.15},
    {"Arcturus", 213.915, 19.182, -0.05, 1.23},    {"Rigil Kentaurus", 219.902, -60.834, -0.27, 0.71},
    {"Vega", 279.234, 38.784, 0.03, 0.00},         {"Capella", 79.172, 45.998, 0.08, 0.80},
    {"Rigel", 78.634, -8.202, 0.13, -0.03},        {"Procyon", 114.825, 5.225, 0.34, 0.42},
    {"Achernar", 24.429, -57.237, 0.46, -0.16},    {"Betelgeuse", 88.793, 7.407, 0.50, 1.85},
    {"Hadar", 210.956, -60.373, 0.61, -0.23},      {"Altair", 297.696, 8.868, 0.76, 0.22},
    {"Acrux", 186.650, -63.099, 0.76, -0.24},      {"Aldebaran", 68.980, 16.509, 0.86, 1.54},
    {"Antares", 247.352, -26.432, 0.96, 1.83},     {"Spica", 201.298, -11.161, 0.97, -0.23},
    {"Pollux", 116.329, 28.026, 1.14, 1.00},       {"Fomalhaut", 344.413, -29.622, 1.16, 0.09},
    {"Deneb", 310.358, 45.280, 1.25, 0.09},        {"Mimosa", 191.930, -59.689, 1.25, -0.23},
    {"Regulus", 152.093, 11.967, 1.35, -0.11},     {"Adhara", 104.656, -28.972, 1.50, -0.21},
    {"Castor", 113.650, 31.888, 1.58, 0.03},       {"Shaula", 263.402, -37.104, 1.62, -0.22},
    {"Gacrux", 187.791, -57.113, 1.64, 1.59},      {"Bellatrix", 81.283, 6.350, 1.64, -0.22},
    {"Elnath", 81.573, 28.608, 1.65, -0.13},       {"Miaplacidus", 138.300, -69.717, 1.69, 0.07},
    {"Alnilam", 84.053, -1.202, 1.69, -0.18},      {"Alnair", 332.058, -46.961, 1.74, -0.13},
    {"Alnitak", 85.190, -1.943, 1.77, -0.21},      {"Alioth", 193.507, 55.960, 1.77, -0.02},
    {"Dubhe", 165.932, 61.751, 1.79, 1.07},        {"Mirfak", 51.081, 49.861, 1.79, 0.48},
    {"Wezen", 107.098, -26.393, 1.83, 0.68},       {"Sargas", 264.330, -42.998, 1.86, 0.40},
    {"Kaus Australis", 276.043, -34.385, 1.85, -0.03}, {"Avior", 125.628, -59.510, 1.86, 1.28},
    {"Alkaid", 206.885, 49.313, 1.86, -0.19},      {"Menkalinan", 89.882, 44.948, 1.90, 0.08},
    {"Atria", 252.166, -69.028, 1.91, 1.44},       {"Alhena", 99.428, 16.399, 1.93, 0.00},
    {"Peacock", 306.412, -56.735, 1.94, -0.20},    {"Polaris", 37.955, 89.264, 1.98, 0.60},
    {"Mirzam", 95.675, -17.956, 1.98, -0.23},      {"Alphard", 141.897, -8.659, 1.99, 1.44},
    {"Hamal", 31.793, 23.462, 2.01, 1.15},         {"Diphda", 10.897, -17.987, 2.04, 1.02},
    {"Nunki", 283.816, -26.297, 2.05, -0.22},      {"Menkent", 211.671, -36.370, 2.06, 1.01},
    {"Alpheratz", 2.097, 29.091, 2.06, -0.11},     {"Mirach", 17.433, 35.621, 2.05, 1.58},
    {"Saiph", 86.939, -9.670, 2.07, -0.18},        {"Kochab", 222.676, 74.156, 2.08, 1.47},
    {"Rasalhague", 263.734, 12.560, 2.08, 0.15},   {"Algol", 47.042, 40.956, 2.12, -0.05},
    {"Almach", 30.975, 42.330, 2.10, 1.37},        {"Denebola", 177.265, 14.572, 2.14, 0.09},
    {"Naos", 120.896, -40.003, 2.25, -0.27},       {"Mizar", 200.981, 54.925, 2.23, 0.02},
    {"Schedar", 10.127, 56.537, 2.24, 1.17},       {"Caph", 2.295, 59.150, 2.28, 0.34},
    {"Enif", 326.046, 9.875, 2.39, 1.52},          {"Scheat", 345.944, 28.083, 2.42, 1.67},
    {"Markab", 346.190, 15.205, 2.49, -0.04},
};

} // namespace

StarField::StarField(int proceduralCount, std::uint64_t seed) {
    using namespace solsim::constants;
    for (const auto& s : kBright) {
        Star st;
        st.direction = astro::poleFromRaDec(s.raDeg, s.decDeg);
        st.magnitude = static_cast<float>(s.vmag);
        st.color = starColor(s.bv);
        st.name = s.name;
        m_stars.push_back(st);
    }
    m_named = m_stars.size();

    // Galactic frame (J2000): north galactic pole and galactic centre.
    const Vec3 ngp = astro::poleFromRaDec(192.85948, 27.12825);
    const Vec3 gc = astro::poleFromRaDec(266.40500, -28.93617);
    const Vec3 gy = ngp.cross(gc).normalized();

    Random rng(seed);
    for (int k = 0; k < proceduralCount; ++k) {
        double b, l;
        if (rng.uniform() < 0.65) {
            // Disc population: Laplace-distributed latitude, longitude biased to the centre.
            const double u = rng.uniform(-1.0, 1.0);
            b = -math::sign(u) * 7.0 * std::log(1.0 - std::fabs(u) + 1e-12) * DEG2RAD;
            b = math::clamp(b, -HALF_PI, HALF_PI);
            l = rng.uniform() < 0.5 ? rng.normal(0.0, 55.0) * DEG2RAD : rng.angle();
        } else {
            b = std::asin(rng.uniform(-1.0, 1.0));
            l = rng.angle();
        }
        const Vec3 dir = gc * (std::cos(b) * std::cos(l)) + gy * (std::cos(b) * std::sin(l)) + ngp * std::sin(b);
        Star st;
        st.direction = dir.normalized();
        // Cumulative counts grow ~ 10^(0.35 m): inverse-transform sampling on [2.5, 9].
        const double u = rng.uniform(1e-6, 1.0);
        st.magnitude = static_cast<float>(math::clamp(9.0 + std::log10(u) / 0.35, 2.5, 9.0));
        st.color = starColor(math::clamp(rng.normal(0.65, 0.45), -0.3, 2.0));
        m_stars.push_back(st);
    }
}

} // namespace solsim::render
