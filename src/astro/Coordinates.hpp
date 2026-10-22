// =============================================================================
//  SolarSim - Coordinates.hpp
//  Reference frames and spherical astronomy:
//     * ecliptic J2000 <-> equatorial J2000 (ICRF-aligned) rotation
//     * IAU 1976 precession (J2000 -> mean equator/equinox of date)
//     * general precession in longitude and an abridged IAU 1980 nutation
//     * horizontal (alt/az) coordinates for an observer on Earth
//     * conversions to RA/Dec and ecliptic longitude/latitude
//
//  The simulation's inertial frame is the ECLIPTIC AND MEAN EQUINOX OF J2000
//  (x toward the vernal equinox, z toward the north ecliptic pole).
// =============================================================================
#pragma once

#include "../core/Mat3.hpp"
#include "../core/Vec3.hpp"

namespace solsim::astro {

struct Spherical {
    double lon = 0.0;  // rad, [0, 2pi)   (RA for equatorial frames)
    double lat = 0.0;  // rad, [-pi/2, pi/2] (Dec for equatorial frames)
    double r = 0.0;
};

Spherical toSpherical(const Vec3& v);
Vec3 fromSpherical(double lon, double lat, double r = 1.0);

// Rotation about the x-axis by the J2000 obliquity.
Vec3 eclipticToEquatorial(const Vec3& ecl);
Vec3 equatorialToEcliptic(const Vec3& equ);
const Mat3& eclipticToEquatorialMatrix();

// Unit vector (in the ecliptic J2000 frame) of a pole given by ICRF RA/Dec (deg).
Vec3 poleFromRaDec(double raDeg, double decDeg);

// Mean obliquity of the ecliptic of date (IAU 2006), radians.
double meanObliquity(double jdTT);

// IAU 1976 precession matrix: J2000 mean equatorial -> mean equatorial of date.
Mat3 precessionMatrix(double jdTT);

// General precession in ecliptic longitude since J2000 (radians).
double generalPrecessionLongitude(double jdTT);

// Abridged IAU 1980 nutation (~0.5" accuracy): returns dpsi and deps in radians.
void nutation(double jdTT, double& dpsi, double& deps);

// Converts a J2000 ecliptic vector to apparent ecliptic longitude of date
// (precession + nutation in longitude). Latitude change is negligible here.
double eclipticLongitudeOfDate(const Vec3& eclJ2000, double jdTT);

// Equatorial (of date) -> horizontal. Returns altitude & azimuth (rad), az
// measured from north through east. lstRad = local sidereal time.
void equatorialToHorizontal(double ra, double dec, double lstRad, double latRad, double& alt, double& az);

// Angular separation between two directions (radians), numerically robust.
double angularSeparation(const Vec3& a, const Vec3& b);

// Constellation-agnostic compass name for an azimuth: "N", "NNE", ...
const char* compassPoint(double azRad);

} // namespace solsim::astro
