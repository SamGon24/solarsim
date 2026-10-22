// =============================================================================
//  SolarSim - StarField.hpp
//  Background sky: the brightest real stars (J2000 positions, V magnitudes and
//  B-V colours) plus a procedural population concentrated toward the galactic
//  plane so the Milky Way band appears in the right place on the sky.
// =============================================================================
#pragma once

#include "Color.hpp"

#include "../core/Vec3.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace solsim::render {

struct Star {
    Vec3 direction;     // unit vector, ecliptic J2000
    float magnitude;    // visual magnitude
    Color color;        // normalised colour from B-V
    std::string name;   // empty for procedural stars
};

class StarField {
public:
    explicit StarField(int proceduralCount = 9000, std::uint64_t seed = 7);
    const std::vector<Star>& stars() const { return m_stars; }
    std::size_t namedCount() const { return m_named; }

private:
    std::vector<Star> m_stars;
    std::size_t m_named = 0;
};

} // namespace solsim::render
