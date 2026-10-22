// =============================================================================
//  SolarSim - test_astro.cpp : Kepler problem, elements, coordinates, theories
// =============================================================================
#include "TestFramework.hpp"

#include "astro/BodyCatalog.hpp"
#include "astro/Coordinates.hpp"
#include "astro/JulianDate.hpp"
#include "astro/Kepler.hpp"
#include "astro/Lambert.hpp"
#include "astro/LunarTheory.hpp"
#include "astro/OrbitalElements.hpp"
#include "core/Constants.hpp"

using namespace solsim;
using namespace solsim::constants;

TEST(Kepler, EllipticAllEccentricities) {
    for (double e : {0.0, 0.01, 0.2, 0.5, 0.9, 0.99, 0.999}) {
        for (double M = -10.0; M <= 10.0; M += 0.37) {
            astro::KeplerSolveInfo info;
            const double E = astro::solveKeplerElliptic(M, e, &info);
            CHECK(info.converged);
            CHECK_NEAR(E - e * std::sin(E), M, 1e-12);
        }
    }
}

TEST(Kepler, Hyperbolic) {
    for (double e : {1.01, 1.5, 3.0, 10.0})
        for (double M : {-50.0, -1.0, 0.0, 0.3, 5.0, 100.0}) {
            const double H = astro::solveKeplerHyperbolic(M, e);
            CHECK_NEAR(e * std::sinh(H) - H, M, 1e-9 * std::max(1.0, std::fabs(M)));
        }
}

TEST(Kepler, MeeusExample30a) {
    // Meeus example 30.a: e = 0.100, M = 5 deg -> E = 5.554589 deg
    const double E = astro::solveKeplerElliptic(5.0 * DEG2RAD, 0.1);
    CHECK_NEAR(E * RAD2DEG, 5.554589254, 1e-8);
}

TEST(Elements, RoundTripAllConics) {
    struct Case { double q, e, i, node, argp; };
    const Case cases[] = {
        {0.9, 0.0167, 0.1, 1.0, 2.0},   {0.3, 0.9671, 2.83, 1.02, 1.94}, // Halley-like retrograde
        {2.0, 0.0, 0.3, 0.5, 0.0},      {1.0, 0.2, 0.0, 0.0, 1.0},        // circular / equatorial
        {0.25, 1.2, 2.14, 0.43, 4.22},  {1.5, 1.0, 0.5, 2.0, 3.0},        // hyperbolic / parabolic
    };
    for (const auto& c : cases) {
        astro::OrbitalElements el;
        el.q = c.q; el.e = c.e; el.i = c.i; el.node = c.node; el.argp = c.argp;
        el.tp = J2000_JD; el.mu = GM_SUN;
        for (double dt : {-40.0, 0.0, 17.0, 250.0}) {
            const auto s = astro::elementsToState(el, J2000_JD + dt);
            const auto el2 = astro::stateToElements(s.r, s.v, GM_SUN, J2000_JD + dt);
            const auto s2 = astro::elementsToState(el2, J2000_JD + dt + 10.0);
            const auto s1 = astro::elementsToState(el, J2000_JD + dt + 10.0);
            CHECK_NEAR((s2.r - s1.r).norm(), 0.0, 1e-9);
            CHECK_NEAR(el2.e, el.e, 1e-9);
            CHECK_NEAR(el2.q, el.q, 1e-9);
        }
    }
}

TEST(Kepler, UniversalPropagatorMatchesElements) {
    astro::OrbitalElements el = astro::OrbitalElements::fromMeanAnomaly(1.7, 0.4, 0.2, 0.3, 0.9, 1.2, J2000_JD, GM_SUN);
    const auto s0 = astro::elementsToState(el, J2000_JD);
    for (double dt : {-900.0, -3.0, 0.5, 77.0, 1234.5, 20000.0}) {
        const auto u = astro::propagateKepler(s0.r, s0.v, GM_SUN, dt);
        const auto ref = astro::elementsToState(el, J2000_JD + dt);
        CHECK(u.converged);
        CHECK_NEAR((u.r - ref.r).norm(), 0.0, 1e-10);
        CHECK_NEAR((u.v - ref.v).norm(), 0.0, 1e-12);
    }
    // hyperbolic
    astro::OrbitalElements h;
    h.q = 0.3; h.e = 2.0; h.i = 0.5; h.tp = J2000_JD; h.mu = GM_SUN;
    const auto hs = astro::elementsToState(h, J2000_JD - 20.0);
    const auto hu = astro::propagateKepler(hs.r, hs.v, GM_SUN, 60.0);
    CHECK_NEAR((hu.r - astro::elementsToState(h, J2000_JD + 40.0).r).norm(), 0.0, 1e-9);
}

TEST(Kepler, LambertRecoversVelocity) {
    const Vec3 r1{1.0, 0.1, -0.02}, v1{-0.001, 0.0175, 0.0008};
    for (double tof : {30.0, 120.0, 250.0}) {
        const auto p = astro::propagateKepler(r1, v1, GM_SUN, tof);
        const auto sol = astro::solveLambert(r1, p.r, tof, GM_SUN, true);
        CHECK(sol.ok);
        CHECK_NEAR((sol.v1 - v1).norm(), 0.0, 1e-10);
        CHECK_NEAR((sol.v2 - p.v).norm(), 0.0, 1e-10);
    }
}

TEST(Coordinates, ObliquityAndFrames) {
    // The ecliptic north pole has Dec = 90 - obliquity, RA = 270 deg.
    const auto s = astro::toSpherical(astro::eclipticToEquatorial(Vec3::unitZ()));
    CHECK_NEAR(s.lat * RAD2DEG, 90.0 - OBLIQUITY_J2000_DEG, 1e-10);
    CHECK_NEAR(s.lon * RAD2DEG, 270.0, 1e-10);
    const Vec3 v{0.3, -0.7, 0.2};
    CHECK_NEAR((astro::equatorialToEcliptic(astro::eclipticToEquatorial(v)) - v).norm(), 0.0, 1e-15);
    // Precession: ~50.3"/yr in longitude.
    CHECK_NEAR(astro::generalPrecessionLongitude(J2000_JD + 36525.0) * RAD2ARCSEC, 5029.9, 1.0);
    // Meeus example 21.b: precession of theta Persei from J2000 to 2028 Nov 13.19 (approx check)
    const Mat3 P = astro::precessionMatrix(2462088.69);
    const Vec3 star = astro::fromSpherical(41.054063 * DEG2RAD, 49.227750 * DEG2RAD);
    const auto out = astro::toSpherical(P * star);
    CHECK_NEAR(out.lon * RAD2DEG, 41.547214, 2e-5);
    CHECK_NEAR(out.lat * RAD2DEG, 49.348483, 2e-5);
}

TEST(Catalog, StandishEarthAtJ2000) {
    const auto* emb = astro::BodyCatalog::find("Earth");
    CHECK(emb != nullptr);
    const auto el = astro::meanElementsAt(emb->mean, J2000_JD, GM_SUN);
    const auto s = astro::elementsToState(el, J2000_JD);
    // JPL Horizons, EMB heliocentric ecliptic J2000 at 2000-01-01.5 TDB: (-0.1771, 0.9672, 0.0000) AU
    CHECK_NEAR(s.r.x, -0.1771, 2e-4);
    CHECK_NEAR(s.r.y, 0.9672, 2e-4);
    CHECK_NEAR(s.r.z, 0.0, 1e-4);
    CHECK_EQ(astro::BodyCatalog::all().front().name, std::string("Sun"));
}

TEST(Catalog, LunarTheoryMeeus47a) {
    // Meeus example 47.a: 1992 April 12, 0h TD -> lambda = 133.162655 deg,
    // beta = -3.229126 deg, distance = 368409.7 km
    const auto p = astro::lunarPositionOfDate(2448724.5);
    CHECK_NEAR(p.longitudeDeg, 133.162655, 2e-3);
    CHECK_NEAR(p.latitudeDeg, -3.229126, 2e-3);
    CHECK_NEAR(p.distanceKm, 368409.7, 5.0);
}
