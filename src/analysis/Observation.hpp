// =============================================================================
//  SolarSim - Observation.hpp
//  What an observer on Earth would see: geocentric (or topocentric) apparent
//  positions with light-time correction, RA/Dec (J2000 and of date),
//  elongation, phase angle, illuminated fraction, apparent magnitude,
//  angular diameter and horizontal coordinates.
// =============================================================================
#pragma once

#include "../physics/NBodySystem.hpp"

#include <string>

namespace solsim::analysis {

struct ObserverSite {
    bool enabled = false;
    double latitudeDeg = 0.0;
    double longitudeDeg = 0.0; // east positive
    double altitudeM = 0.0;
};

struct Observation {
    std::string name;
    Vec3 apparent;               // observer -> body, light-time corrected, ecliptic J2000 [AU]
    double distance = 0.0;       // AU
    double lightTimeMin = 0.0;   // minutes
    double helioDistance = 0.0;  // AU
    double raJ2000 = 0.0, decJ2000 = 0.0;   // degrees
    double raDate = 0.0, decDate = 0.0;     // degrees, mean equator & equinox of date
    double eclLon = 0.0, eclLat = 0.0;      // degrees, ecliptic of date (apparent longitude)
    double elongation = 0.0;     // Sun-observer-body angle, degrees
    bool eveningSky = false;     // body east of the Sun
    double phaseAngle = 0.0;     // Sun-body-observer angle, degrees
    double illuminated = 1.0;    // fraction of the disc lit
    bool hasMagnitude = false;
    double magnitude = 0.0;
    double angularDiameter = 0.0; // arcseconds
    bool hasHorizontal = false;
    double altitude = 0.0, azimuth = 0.0; // degrees
};

// Observes body `target` from body `observer` (normally Earth) at time jdTT.
// The system must be at time jdTT.
Observation observe(const physics::NBodySystem& sys, std::size_t target, std::size_t observer, double jdTT,
                    const ObserverSite& site = {});

// Light-time corrected vector from observer to target (barycentric linear back-propagation).
Vec3 apparentVector(const physics::NBodySystem& sys, std::size_t target, const Vec3& observerPos);

// Apparent geocentric ecliptic longitude of the Sun (of date, with aberration), degrees.
double apparentSolarLongitude(const physics::NBodySystem& sys, std::size_t earth, double jdTT);
// Apparent geocentric ecliptic longitude of a body (of date), degrees.
double apparentLongitude(const physics::NBodySystem& sys, std::size_t target, std::size_t earth, double jdTT);

// Visual magnitude (Meeus ch. 41 / Astronomical Almanac formulae).
bool apparentMagnitude(const std::string& name, double r, double delta, double phaseDeg, double ringTiltRad,
                       double& magOut);

} // namespace solsim::analysis
