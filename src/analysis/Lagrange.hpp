// =============================================================================
//  SolarSim - Lagrange.hpp
//  Circular restricted three-body problem (CR3BP) quantities for a primary /
//  secondary pair: the five Lagrange points (collinear L1-L3 from the quintic
//  equilibrium condition solved with Brent's method; triangular L4/L5), Hill
//  sphere, Laplace sphere of influence, Routh's stability criterion and a
//  numerical stability experiment with test particles.
// =============================================================================
#pragma once

#include "../physics/Simulation.hpp"

#include <array>
#include <string>

namespace solsim::analysis {

struct LagrangeSolution {
    std::string primary, secondary;
    double massRatio = 0.0;           // mu = m2 / (m1 + m2)
    double separation = 0.0;          // current distance (AU)
    std::array<double, 5> x{};        // rotating-frame coordinates (units of separation, barycentre origin)
    std::array<double, 5> y{};
    std::array<Vec3, 5> position{};   // barycentric positions now (AU)
    std::array<Vec3, 5> velocity{};   // co-rotating velocities (AU/day)
    std::array<double, 5> distanceFromSecondary{}; // AU
    double hillRadius = 0.0;          // AU, a (1-e) cbrt(m2 / 3 m1)
    double sphereOfInfluence = 0.0;   // AU, a (m2/m1)^(2/5)
    bool triangularStable = false;    // Routh: mu < 0.0385209
    double jacobiL1 = 0.0;            // Jacobi constant at L1 (normalised units)
};

// Collinear point solve for given mass ratio. which = 1, 2 or 3.
double collinearLagrangePoint(double mu, int which);

LagrangeSolution computeLagrangePoints(const physics::NBodySystem& sys, int primary, int secondary);

struct LibrationResult {
    int point = 0;
    double years = 0.0;
    double maxDeviationAU = 0.0;    // max distance from the instantaneous L point
    double finalDeviationAU = 0.0;
    double maxAngleDeg = 0.0;       // max angular excursion around the primary
    bool escaped = false;
};

// Places a massless particle at each requested point and integrates.
std::array<LibrationResult, 5> testLagrangeStability(const physics::Simulation& sim, int primary, int secondary,
                                                     double years, double perturbationAU = 0.0);

} // namespace solsim::analysis
