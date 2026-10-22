// =============================================================================
//  SolarSim - Observation.cpp
// =============================================================================
#include "Observation.hpp"

#include "../astro/Coordinates.hpp"
#include "../astro/JulianDate.hpp"
#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"
#include "../core/StringUtil.hpp"

#include <cmath>

namespace solsim::analysis {

using namespace solsim::constants;

Vec3 apparentVector(const physics::NBodySystem& sys, std::size_t target, const Vec3& observerPos) {
    // Iterate tau = |x(t - tau) - o| / c with x(t - tau) ~ x(t) - v tau.
    Vec3 d = sys.pos(target) - observerPos;
    double tau = 0.0;
    for (int k = 0; k < 3; ++k) {
        tau = d.norm() / C_AU_PER_DAY;
        d = (sys.pos(target) - sys.vel(target) * tau) - observerPos;
    }
    return d;
}

bool apparentMagnitude(const std::string& name, double r, double delta, double i, double ringTilt, double& m) {
    const double base = 5.0 * std::log10(r * delta);
    const std::string n = str::toLower(name);
    if (n == "mercury") m = -0.42 + base + 0.0380 * i - 0.000273 * i * i + 0.000002 * i * i * i;
    else if (n == "venus") m = -4.40 + base + 0.0009 * i + 0.000239 * i * i - 0.00000065 * i * i * i;
    else if (n == "mars") m = -1.52 + base + 0.016 * i;
    else if (n == "jupiter") m = -9.40 + base + 0.005 * i;
    else if (n == "saturn") {
        const double sb = std::sin(std::fabs(ringTilt));
        m = -8.88 + base + 0.044 * i - 2.60 * sb + 1.25 * sb * sb;
    } else if (n == "uranus") m = -7.19 + base;
    else if (n == "neptune") m = -6.87 + base;
    else if (n == "pluto") m = -1.00 + base;
    else if (n == "moon") m = 0.21 + base + 0.026 * std::fabs(i) + 4.0e-9 * std::pow(i, 4.0);
    else if (n == "ceres") m = 3.34 + base + 0.03 * i;
    else if (n == "vesta") m = 3.20 + base + 0.03 * i;
    else return false;
    return true;
}

Observation observe(const physics::NBodySystem& sys, std::size_t target, std::size_t observer, double jdTT,
                    const ObserverSite& site) {
    Observation ob;
    ob.name = sys.info(target).name;
    Vec3 obsPos = sys.pos(observer);

    // Topocentric correction (important for the Moon: parallax up to ~1 deg).
    double lst = 0.0;
    if (site.enabled) {
        const double jdUt = astro::ttToUTC(jdTT);
        lst = math::wrapTwoPi((astro::gmstDegrees(jdUt) + site.longitudeDeg) * DEG2RAD);
        const double phi = site.latitudeDeg * DEG2RAD;
        const double R = (EARTH_EQ_RADIUS_KM + site.altitudeM / 1000.0) / AU_KM;
        // Site vector in the equatorial frame of date -> J2000 -> ecliptic.
        const Vec3 siteEqDate{R * std::cos(phi) * std::cos(lst), R * std::cos(phi) * std::sin(lst), R * std::sin(phi)};
        const Vec3 siteEqJ2000 = astro::precessionMatrix(jdTT).transposed() * siteEqDate;
        obsPos += astro::equatorialToEcliptic(siteEqJ2000);
    }

    const Vec3 d = apparentVector(sys, target, obsPos);
    ob.apparent = d;
    ob.distance = d.norm();
    ob.lightTimeMin = ob.distance / C_AU_PER_DAY * 1440.0;
    const Vec3 sunVec = apparentVector(sys, 0, obsPos);
    const Vec3 helio = d - sunVec; // Sun -> body (as seen at emission time, approximately)
    ob.helioDistance = target == 0 ? 0.0 : helio.norm();

    const astro::Spherical eqJ = astro::toSpherical(astro::eclipticToEquatorial(d));
    ob.raJ2000 = eqJ.lon * RAD2DEG;
    ob.decJ2000 = eqJ.lat * RAD2DEG;
    const astro::Spherical eqD = astro::toSpherical(astro::precessionMatrix(jdTT) * astro::eclipticToEquatorial(d));
    ob.raDate = eqD.lon * RAD2DEG;
    ob.decDate = eqD.lat * RAD2DEG;
    ob.eclLon = astro::eclipticLongitudeOfDate(d, jdTT) * RAD2DEG;
    ob.eclLat = astro::toSpherical(d).lat * RAD2DEG;

    if (target != 0) {
        ob.elongation = astro::angularSeparation(d, sunVec) * RAD2DEG;
        const double dl = math::wrap180(std::atan2(d.y, d.x) * RAD2DEG - std::atan2(sunVec.y, sunVec.x) * RAD2DEG);
        ob.eveningSky = dl > 0.0;
        ob.phaseAngle = astro::angularSeparation(-helio, -d) * RAD2DEG;
        ob.illuminated = 0.5 * (1.0 + std::cos(ob.phaseAngle * DEG2RAD));
    }
    const double radius = sys.info(target).radius;
    if (ob.distance > radius) ob.angularDiameter = 2.0 * std::asin(radius / ob.distance) * RAD2ARCSEC;

    // Saturn's ring opening angle B as seen from the observer.
    const double ringTilt = std::asin(math::clamp(-d.normalized().dot(sys.info(target).pole.normalized()), -1.0, 1.0));
    if (target == 0) {
        ob.hasMagnitude = true;
        ob.magnitude = SUN_APP_MAG_V + 5.0 * std::log10(ob.distance);
    } else {
        ob.hasMagnitude = apparentMagnitude(ob.name, ob.helioDistance, ob.distance, ob.phaseAngle, ringTilt, ob.magnitude);
    }

    if (site.enabled) {
        double alt = 0.0, az = 0.0;
        astro::equatorialToHorizontal(eqD.lon, eqD.lat, lst, site.latitudeDeg * DEG2RAD, alt, az);
        // Standard atmospheric refraction (Bennett 1982), significant near the horizon.
        const double altDeg = alt * RAD2DEG;
        if (altDeg > -1.0) {
            const double R = 1.02 / std::tan((altDeg + 10.3 / (altDeg + 5.11)) * DEG2RAD) / 60.0;
            alt += R * DEG2RAD;
        }
        ob.hasHorizontal = true;
        ob.altitude = alt * RAD2DEG;
        ob.azimuth = az * RAD2DEG;
    }
    return ob;
}

double apparentSolarLongitude(const physics::NBodySystem& sys, std::size_t earth, double jdTT) {
    const Vec3 d = apparentVector(sys, 0, sys.pos(earth));
    double lon = astro::eclipticLongitudeOfDate(d, jdTT);
    // Annual aberration: -20.4898" / R (R in AU).
    lon -= 20.4898 * ARCSEC2RAD / d.norm();
    return math::wrap360(lon * RAD2DEG);
}

double apparentLongitude(const physics::NBodySystem& sys, std::size_t target, std::size_t earth, double jdTT) {
    const Vec3 d = apparentVector(sys, target, sys.pos(earth));
    return math::wrap360(astro::eclipticLongitudeOfDate(d, jdTT) * RAD2DEG);
}

} // namespace solsim::analysis
