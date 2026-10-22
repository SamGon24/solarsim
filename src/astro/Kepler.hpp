// =============================================================================
//  SolarSim - Kepler.hpp
//  Solutions of the two-body problem:
//    * Kepler's equation for elliptic (M = E - e sin E), hyperbolic
//      (M = e sinh H - H) and parabolic (Barker's equation) orbits
//    * conversion between eccentric / hyperbolic / true anomalies
//    * Stumpff functions and the universal-variable propagator
//      (Battin / Vallado) used by the Wisdom-Holman integrator and the
//      mission planner to advance a body along a Keplerian arc exactly.
// =============================================================================
#pragma once

#include "../core/Vec3.hpp"

namespace solsim::astro {

struct KeplerSolveInfo {
    int iterations = 0;
    bool converged = false;
};

// Elliptic: returns eccentric anomaly E (rad) for mean anomaly M (rad), 0 <= e < 1.
double solveKeplerElliptic(double M, double e, KeplerSolveInfo* info = nullptr);
// Hyperbolic: returns hyperbolic anomaly H for mean anomaly M, e > 1.
double solveKeplerHyperbolic(double M, double e, KeplerSolveInfo* info = nullptr);
// Parabolic (Barker): W = 3 * sqrt(mu / (2 q^3)) * (t - tp); returns true anomaly.
double solveBarker(double W);

double trueFromEccentric(double E, double e);
double eccentricFromTrue(double nu, double e);
double meanFromEccentric(double E, double e);
double trueFromHyperbolic(double H, double e);
double hyperbolicFromTrue(double nu, double e);
double meanFromHyperbolic(double H, double e);

// Stumpff functions c2(psi) = (1 - cos sqrt psi)/psi, c3(psi) = (sqrt psi - sin sqrt psi)/psi^1.5
// with analytic continuation for psi < 0 and series near zero.
void stumpff(double psi, double& c2, double& c3);

struct UniversalResult {
    Vec3 r;
    Vec3 v;
    int iterations = 0;
    bool converged = false;
};

// Propagates the relative state (r0, v0) about a central body with
// gravitational parameter mu for time dt (any conic, any sign of dt).
UniversalResult propagateKepler(const Vec3& r0, const Vec3& v0, double mu, double dt);

} // namespace solsim::astro
