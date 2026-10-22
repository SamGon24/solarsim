// =============================================================================
//  SolarSim - Lambert.hpp
//  Lambert's problem: find the Keplerian arc connecting position r1 to r2 in a
//  given time of flight. Universal-variable formulation with bisection on
//  psi (Bate, Mueller & White 1971; Vallado, "Fundamentals of Astrodynamics",
//  algorithm 58) - slower than Izzo's method but extremely robust for the
//  single-revolution transfers used in porkchop plots.
// =============================================================================
#pragma once

#include "../core/Vec3.hpp"

namespace solsim::astro {

struct LambertSolution {
    Vec3 v1;           // departure velocity
    Vec3 v2;           // arrival velocity
    bool ok = false;
    int iterations = 0;
    double transferAngleDeg = 0.0;
};

LambertSolution solveLambert(const Vec3& r1, const Vec3& r2, double tof, double mu, bool prograde = true);

} // namespace solsim::astro
