// =============================================================================
//  SolarSim - test_physics.cpp : integrators, conservation laws, builder
// =============================================================================
#include "TestFramework.hpp"

#include "analysis/Lagrange.hpp"
#include "analysis/Precession.hpp"
#include "astro/Kepler.hpp"
#include "core/Constants.hpp"
#include "physics/Simulation.hpp"
#include "physics/SystemBuilder.hpp"

using namespace solsim;
using namespace solsim::constants;

namespace {

physics::NBodySystem twoBody(double e) {
    physics::NBodySystem s;
    physics::BodyInfo sun;
    sun.name = "Sun";
    s.addBody(sun, GM_SUN, {0, 0, 0}, {0, 0, 0});
    physics::BodyInfo p;
    p.name = "Planet";
    p.parent = 0;
    const double r = 1.0 - e; // perihelion of a = 1 AU orbit
    const double v = std::sqrt(GM_SUN * (1 + e) / r);
    s.addBody(p, 3e-12, {r, 0, 0}, {0, v, 0});
    s.moveToBarycentricFrame();
    return s;
}

double positionError(const std::string& integrator, double dt, double days) {
    const auto sys = twoBody(0.3);
    physics::SimulationConfig cfg;
    cfg.integrator = integrator;
    cfg.dt = dt;
    physics::Simulation sim(sys, 0.0, cfg);
    const Vec3 r0 = sys.pos(1) - sys.pos(0), v0 = sys.vel(1) - sys.vel(0);
    sim.advance(days);
    const auto ref = astro::propagateKepler(r0, v0, sys.gm()[0] + sys.gm()[1], days);
    return ((sim.system().pos(1) - sim.system().pos(0)) - ref.r).norm();
}

} // namespace

TEST(Integrators, ConvergenceOrder) {
    // Halving the step must reduce the global error by ~2^order.
    struct C { const char* name; int order; double dt; };
    for (const C c : {C{"euler", 1, 0.01}, C{"symplectic-euler", 1, 0.01}, C{"leapfrog", 2, 0.5}, C{"rk4", 4, 2.0},
                      C{"yoshida4", 4, 2.0}, C{"wh", 2, 2.0}}) {
        const double days = c.order == 1 ? 50.0 : 400.0;
        const double e1 = positionError(c.name, c.dt, days);
        const double e2 = positionError(c.name, c.dt / 2.0, days);
        const double observed = std::log2(e1 / e2);
        CHECK_NEAR(observed, c.order, 0.35);
    }
}

TEST(Integrators, AdaptiveAndWisdomHolmanAccuracy) {
    CHECK(positionError("dopri5", 1.0, 1000.0) < 1e-8);
    // Wisdom-Holman solves the Kepler drift exactly: even 5-day steps through an
    // e = 0.3 orbit stay accurate (the residual comes from the small Sun-drift term).
    CHECK(positionError("wh", 5.0, 1000.0) < 1e-9);
    CHECK(positionError("wh", 5.0, 1000.0) < 1e-3 * positionError("leapfrog", 5.0, 1000.0));
}

TEST(Integrators, SymplecticEnergyIsBounded) {
    const auto sys = twoBody(0.5);
    for (const char* name : {"leapfrog", "yoshida4", "rk4"}) {
        physics::SimulationConfig cfg;
        cfg.integrator = name;
        cfg.dt = 3.0;
        physics::Simulation sim(sys, 0.0, cfg);
        double early = 0, late = 0;
        const int perYear = 37; // sample densely so every orbital phase is seen
        for (int k = 0; k < 200 * perYear; ++k) {
            sim.advance(365.25 / perYear);
            const double e = std::fabs(sim.invariantErrors().energyRel);
            if (k < 20 * perYear) early = std::max(early, e);
            if (k >= 180 * perYear) late = std::max(late, e);
        }
        if (std::string(name) == "rk4") CHECK(late > 5.0 * early);   // secular drift
        else CHECK(late < 1.5 * early + 1e-14);                      // bounded oscillation
    }
}

TEST(Integrators, TimeReversibility) {
    auto sys = twoBody(0.2);
    physics::SimulationConfig cfg;
    cfg.integrator = "leapfrog";
    cfg.dt = 1.0;
    physics::Simulation sim(sys, 0.0, cfg);
    const Vec3 start = sim.system().pos(1);
    sim.advance(1000.0);
    sim.advance(-1000.0);
    CHECK_NEAR((sim.system().pos(1) - start).norm(), 0.0, 1e-12);
}

TEST(Builder, DefaultSystemIsBarycentric) {
    physics::BuildOptions o;
    o.epochJD = J2000_JD;
    o.calibration = physics::CalibrationMode::None;
    const auto sys = physics::buildSystem(o);
    CHECK_EQ(sys.namedCount(), static_cast<std::size_t>(11));
    CHECK_EQ(sys.info(0).name, std::string("Sun"));
    CHECK_NEAR(sys.centerOfMass().norm(), 0.0, 1e-14);
    CHECK_NEAR(sys.centerOfMassVelocity().norm(), 0.0, 1e-18);
    // Earth-Moon distance must lie between perigee and apogee extremes.
    const double d = (sys.pos(static_cast<std::size_t>(sys.indexOf("Moon"))) - sys.pos(static_cast<std::size_t>(sys.indexOf("Earth")))).norm();
    CHECK(d * AU_KM > 356000.0 && d * AU_KM < 407000.0);
    // Sun's offset from the barycentre is dominated by Jupiter (~0.005-0.01 AU)
    CHECK(sys.pos(0).norm() > 0.002 && sys.pos(0).norm() < 0.012);
}

TEST(Builder, MoonsFoldMassWhenAbsent) {
    physics::BuildOptions a;
    a.calibration = physics::CalibrationMode::None;
    physics::BuildOptions b = a;
    b.moons = true;
    const auto sa = physics::buildSystem(a), sb = physics::buildSystem(b);
    const double jupA = sa.gm()[static_cast<std::size_t>(sa.indexOf("Jupiter"))];
    double jupB = sb.gm()[static_cast<std::size_t>(sb.indexOf("Jupiter"))];
    for (const char* m : {"Io", "Europa", "Ganymede", "Callisto"}) jupB += sb.gm()[static_cast<std::size_t>(sb.indexOf(m))];
    CHECK_NEAR(jupA / jupB, 1.0, 1e-12);
    CHECK_NEAR(sa.totalGM() / sb.totalGM(), 1.0, 1e-12);
}

TEST(Builder, PopulationsAreBound) {
    physics::BuildOptions o;
    o.calibration = physics::CalibrationMode::None;
    o.populations.asteroids = 300;
    o.populations.trojans = 100;
    o.populations.kuiper = 100;
    const auto sys = physics::buildSystem(o);
    CHECK_EQ(sys.particleCount(), static_cast<std::size_t>(500));
    for (std::size_t i = sys.namedCount(); i < sys.size(); ++i) {
        const Vec3 r = sys.pos(i) - sys.pos(0), v = sys.vel(i) - sys.vel(0);
        const double energy = 0.5 * v.norm2() - sys.gm()[0] / r.norm();
        CHECK(energy < 0.0);
    }
}

TEST(Physics, RelativisticPrecession) {
    physics::BuildOptions o;
    o.only = {"Sun", "Mercury"};
    o.calibration = physics::CalibrationMode::None;
    const auto sys = physics::buildSystem(o);
    analysis::PrecessionConfig pc;
    pc.years = 20;
    pc.dt = 0.05;
    pc.relativity = true;
    const double gr = analysis::measurePrecession(sys, pc).rateArcsecPerCentury;
    pc.relativity = false;
    const double newton = analysis::measurePrecession(sys, pc).rateArcsecPerCentury;
    CHECK_NEAR(gr - newton, 42.98, 0.3);
    CHECK_NEAR(newton, 0.0, 0.5);
}

TEST(Physics, LagrangePointsSunEarth) {
    physics::BuildOptions o;
    o.calibration = physics::CalibrationMode::None;
    const auto sys = physics::buildSystem(o);
    const auto L = analysis::computeLagrangePoints(sys, 0, sys.indexOf("Earth"));
    // L1/L2 lie ~1.5 million km from Earth.
    CHECK_NEAR(L.distanceFromSecondary[0] * AU_KM, 1.49e6, 0.04e6);
    CHECK_NEAR(L.distanceFromSecondary[1] * AU_KM, 1.50e6, 0.04e6);
    CHECK(L.triangularStable);
    CHECK_NEAR(analysis::collinearLagrangePoint(0.012150585, 1), 0.836915, 1e-5); // Earth-Moon L1
}
