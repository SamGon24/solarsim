// =============================================================================
//  SolarSim - OrbitalElements.hpp
//  Classical (Keplerian) orbital elements and conversions to/from Cartesian
//  state vectors. Valid for all conic sections:
//
//     * ellipse   0 <= e < 1   (a > 0)
//     * parabola       e = 1   (a undefined -> perihelion distance q is used)
//     * hyperbola      e > 1   (a < 0)
//
//  To support every case with one representation the primary elements are the
//  perihelion distance q and the time of perihelion passage tp; the semi-major
//  axis and the mean anomaly are derived quantities.
// =============================================================================
#pragma once

#include "../core/Vec3.hpp"

#include <string>

namespace solsim::astro {

struct StateVector {
    Vec3 r; // position
    Vec3 v; // velocity
};

struct OrbitalElements {
    double q = 1.0;      // perihelion distance [AU]
    double e = 0.0;      // eccentricity
    double i = 0.0;      // inclination [rad]
    double node = 0.0;   // longitude of ascending node, Omega [rad]
    double argp = 0.0;   // argument of periapsis, omega [rad]
    double tp = 0.0;     // time of periapsis passage [JD]
    double epoch = 0.0;  // osculation epoch [JD] (informational)
    double mu = 0.0;     // gravitational parameter of the central body [AU^3/d^2]

    // ---- derived quantities ----
    bool isElliptic() const { return e < 1.0 - 1e-12; }
    bool isHyperbolic() const { return e > 1.0 + 1e-12; }
    bool isParabolic() const { return !isElliptic() && !isHyperbolic(); }

    double a() const;                 // semi-major axis (negative for hyperbolae, inf for parabola)
    double Q() const;                 // aphelion distance (inf if unbound)
    double p() const { return q * (1.0 + e); } // semi-latus rectum
    double meanMotion() const;        // rad/day
    double period() const;            // days (inf if unbound)
    double meanAnomalyAt(double jd) const;
    double trueAnomalyAt(double jd) const;
    double longitudeOfPeriapsis() const { return node + argp; } // varpi
    double meanLongitudeAt(double jd) const;  // node + argp + M (ellipse)
    double specificEnergy() const;            // -mu / 2a
    double specificAngularMomentum() const;   // sqrt(mu p)

    // Construct from (a, e, i, node, argp, M at epoch).
    static OrbitalElements fromMeanAnomaly(double a, double e, double i, double node, double argp,
                                           double meanAnomaly, double epoch, double mu);
    // Construct from (a, e, i, longitude of perihelion, mean longitude L at epoch).
    static OrbitalElements fromMeanLongitude(double a, double e, double i, double node, double varpi,
                                             double meanLongitude, double epoch, double mu);

    std::string describe() const;
};

// Cartesian state at Julian date `jd` from elements (reference frame of the elements).
StateVector elementsToState(const OrbitalElements& el, double jd);
// Position/velocity at a given true anomaly.
StateVector stateAtTrueAnomaly(const OrbitalElements& el, double nu);

// Osculating elements from a state vector at time `jd`.
OrbitalElements stateToElements(const Vec3& r, const Vec3& v, double mu, double jd);

// Perifocal -> inertial rotation matrix columns (P, Q, W unit vectors).
void perifocalBasis(double i, double node, double argp, Vec3& P, Vec3& Q, Vec3& W);

} // namespace solsim::astro
