// =============================================================================
//  SolarSim - Coordinates.cpp
// =============================================================================
#include "Coordinates.hpp"

#include "JulianDate.hpp"

#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"

#include <cmath>

namespace solsim::astro {

using namespace solsim::constants;

Spherical toSpherical(const Vec3& v) {
    Spherical s;
    s.r = v.norm();
    if (s.r == 0.0) return s;
    s.lon = math::wrapTwoPi(std::atan2(v.y, v.x));
    s.lat = std::atan2(v.z, std::sqrt(v.x * v.x + v.y * v.y));
    return s;
}

Vec3 fromSpherical(double lon, double lat, double r) {
    const double cl = std::cos(lat);
    return {r * cl * std::cos(lon), r * cl * std::sin(lon), r * std::sin(lat)};
}

const Mat3& eclipticToEquatorialMatrix() {
    static const Mat3 m = Mat3::rotX(OBLIQUITY_J2000);
    return m;
}

Vec3 eclipticToEquatorial(const Vec3& ecl) { return eclipticToEquatorialMatrix() * ecl; }

Vec3 equatorialToEcliptic(const Vec3& equ) { return eclipticToEquatorialMatrix().transposed() * equ; }

Vec3 poleFromRaDec(double raDeg, double decDeg) {
    const Vec3 equ = fromSpherical(raDeg * DEG2RAD, decDeg * DEG2RAD, 1.0);
    return equatorialToEcliptic(equ).normalized();
}

double meanObliquity(double jdTT) {
    const double T = julianCenturies(jdTT);
    const double arcsec = 84381.406 - 46.836769 * T - 0.0001831 * T * T + 0.00200340 * T * T * T;
    return arcsec * ARCSEC2RAD;
}

Mat3 precessionMatrix(double jdTT) {
    const double T = julianCenturies(jdTT);
    const double zeta = (2306.2181 * T + 0.30188 * T * T + 0.017998 * T * T * T) * ARCSEC2RAD;
    const double z = (2306.2181 * T + 1.09468 * T * T + 0.018203 * T * T * T) * ARCSEC2RAD;
    const double theta = (2004.3109 * T - 0.42665 * T * T - 0.041833 * T * T * T) * ARCSEC2RAD;
    // P = Rz(z) * Ry(-theta) * Rz(zeta)  (active-rotation convention)
    return Mat3::rotZ(z) * Mat3::rotY(-theta) * Mat3::rotZ(zeta);
}

double generalPrecessionLongitude(double jdTT) {
    const double T = julianCenturies(jdTT);
    return (5028.796195 * T + 1.1054348 * T * T) * ARCSEC2RAD;
}

void nutation(double jdTT, double& dpsi, double& deps) {
    const double T = julianCenturies(jdTT);
    // Fundamental arguments (degrees) - Meeus ch. 22.
    const double omega = math::wrap360(125.04452 - 1934.136261 * T) * DEG2RAD;  // Moon's node
    const double Ls = math::wrap360(280.4665 + 36000.7698 * T) * DEG2RAD;       // Sun mean longitude
    const double Lm = math::wrap360(218.3165 + 481267.8813 * T) * DEG2RAD;      // Moon mean longitude
    dpsi = (-17.20 * std::sin(omega) - 1.32 * std::sin(2 * Ls) - 0.23 * std::sin(2 * Lm) +
            0.21 * std::sin(2 * omega)) * ARCSEC2RAD;
    deps = (9.20 * std::cos(omega) + 0.57 * std::cos(2 * Ls) + 0.10 * std::cos(2 * Lm) -
            0.09 * std::cos(2 * omega)) * ARCSEC2RAD;
}

double eclipticLongitudeOfDate(const Vec3& eclJ2000, double jdTT) {
    double dpsi = 0.0, deps = 0.0;
    nutation(jdTT, dpsi, deps);
    const double lon = std::atan2(eclJ2000.y, eclJ2000.x);
    return math::wrapTwoPi(lon + generalPrecessionLongitude(jdTT) + dpsi);
}

void equatorialToHorizontal(double ra, double dec, double lst, double lat, double& alt, double& az) {
    const double H = lst - ra; // hour angle
    const double sinAlt = std::sin(lat) * std::sin(dec) + std::cos(lat) * std::cos(dec) * std::cos(H);
    alt = std::asin(math::clamp(sinAlt, -1.0, 1.0));
    // Azimuth from north, eastward.
    const double y = -std::sin(H) * std::cos(dec);
    const double x = std::sin(dec) * std::cos(lat) - std::cos(dec) * std::sin(lat) * std::cos(H);
    az = math::wrapTwoPi(std::atan2(y, x));
}

double angularSeparation(const Vec3& a, const Vec3& b) { return angleBetween(a, b); }

const char* compassPoint(double azRad) {
    static const char* pts[] = {"N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                                "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"};
    const double deg = math::wrap360(azRad * RAD2DEG);
    const int idx = static_cast<int>(std::floor((deg + 11.25) / 22.5)) % 16;
    return pts[idx];
}

} // namespace solsim::astro
