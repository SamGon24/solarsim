// =============================================================================
//  SolarSim - LunarTheory.hpp
//  Geocentric position of the Moon from the truncated ELP-2000/82 series as
//  given in Meeus, "Astronomical Algorithms" (2nd ed.), chapter 47:
//  ~60 periodic terms in longitude & distance and ~60 in latitude, plus the
//  additive Venus/Jupiter/flattening terms. Accuracy ~10" in longitude and
//  ~4" in latitude - far better than Keplerian mean elements, because the
//  Moon's orbit is strongly perturbed by the Sun (evection 1.27 deg,
//  variation 0.66 deg, annual equation 0.19 deg ...).
//
//  Used to create the Moon's initial state vector; the N-body integration
//  then reproduces the same perturbations dynamically.
// =============================================================================
#pragma once

#include "../core/Vec3.hpp"

namespace solsim::astro {

struct LunarPosition {
    double longitudeDeg = 0.0;  // geocentric ecliptic longitude, mean equinox of date
    double latitudeDeg = 0.0;
    double distanceKm = 0.0;
};

LunarPosition lunarPositionOfDate(double jdTT);

// Geocentric position in the ecliptic & equinox of J2000 (AU).
Vec3 lunarPositionJ2000(double jdTT);
// Position and velocity (AU, AU/day) by central differences of the series.
void lunarStateJ2000(double jdTT, Vec3& r, Vec3& v);

} // namespace solsim::astro
