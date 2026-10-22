// =============================================================================
//  SolarSim - test_render_analysis.cpp : camera, mapping, font, terminal
//  output, events and mission design.
// =============================================================================
#include "TestFramework.hpp"

#include "analysis/Events.hpp"
#include "analysis/Observation.hpp"
#include "analysis/Transfer.hpp"
#include "astro/JulianDate.hpp"
#include "core/Constants.hpp"
#include "physics/SystemBuilder.hpp"
#include "render/Camera.hpp"
#include "render/DisplayMapping.hpp"
#include "render/Font.hpp"
#include "render/TerminalRenderer.hpp"

using namespace solsim;
using namespace solsim::constants;

TEST(Render, CameraProjection) {
    render::Camera cam;
    cam.setViewport(800, 600);
    cam.setTarget({1, 2, 3});
    cam.setDistance(10);
    cam.setAngles(0.3, 0.4);
    render::Projected p;
    CHECK(cam.project({1, 2, 3}, p));
    CHECK_NEAR(p.x, 400.0, 1e-9);
    CHECK_NEAR(p.y, 300.0, 1e-9);
    CHECK_NEAR(p.depth, 10.0, 1e-9);
    // A ray through a projected point passes through the original point.
    const Vec3 w{1.5, 2.2, 2.7};
    CHECK(cam.project(w, p));
    const Vec3 dir = cam.rayDirection(p.x, p.y);
    const Vec3 toW = (w - cam.eye()).normalized();
    CHECK_NEAR((dir - toW).norm(), 0.0, 1e-9);
    // Points behind the camera are rejected.
    CHECK(!cam.project(cam.eye() - cam.forward(), p));
}

TEST(Render, DisplayMappingMonotonic) {
    render::DisplayMapping m;
    for (auto mode : {render::ScaleMode::Linear, render::ScaleMode::Sqrt, render::ScaleMode::Log}) {
        m.setMode(mode);
        CHECK_NEAR(m.radial(1.0), 1.0, 1e-12);
        double prev = 0;
        for (double r = 0.01; r < 100; r *= 1.3) {
            const double v = m.radial(r);
            CHECK(v > prev);
            prev = v;
        }
    }
}

TEST(Render, FontAndTerminal) {
    CHECK(render::BitmapFont::pixel('A', 2, 0));
    CHECK(!render::BitmapFont::pixel(' ', 2, 3));
    CHECK_EQ(render::BitmapFont::textWidth("abc", 2), (3 * 6 - 1) * 2);
    CHECK_EQ(render::TerminalRenderer::xterm256(0, 0, 0), 16);
    CHECK_EQ(render::TerminalRenderer::xterm256(255, 0, 0), 196);
    render::Image8 img;
    img.resize(10, 6);
    render::TerminalRenderer t;
    t.setMode(render::TerminalColorMode::Ascii);
    const std::string s = t.compose(img, {}, {"header"}, {"footer"}, 10, 5);
    CHECK(s.find("header") != std::string::npos);
    CHECK(s.find("footer") != std::string::npos);
}

TEST(Analysis, HohmannEarthMars) {
    const auto h = analysis::hohmann(GM_SUN, 1.0, 1.524);
    CHECK_NEAR(h.dv1, 2.94, 0.03);       // km/s
    CHECK_NEAR(h.dv2, 2.65, 0.03);
    CHECK_NEAR(h.timeOfFlight, 259.0, 2.0);
    CHECK_NEAR(h.phaseAngleDeg, 44.3, 0.5);
    CHECK_NEAR(h.synodicPeriod, 780.0, 5.0);
}

TEST(Analysis, JuneSolstice2024) {
    double jd0 = 0, ref = 0;
    std::string err;
    astro::parseDate("2024-06-01", jd0, err);
    astro::parseDate("2024-06-20 20:51", ref, err); // published UTC time
    physics::BuildOptions o;
    o.epochJD = jd0;
    o.calibration = physics::CalibrationMode::Moon;
    physics::Simulation sim(physics::buildSystem(o), jd0);
    analysis::EventSearchOptions eo;
    eo.startJD = jd0;
    eo.endJD = jd0 + 40;
    eo.planets = eo.moonPhases = eo.apsides = false;
    const auto ev = analysis::findEvents(sim, eo);
    bool found = false;
    for (const auto& e : ev)
        if (e.description == "June solstice") {
            found = true;
            CHECK_NEAR((e.jd - ref) * 1440.0, 0.0, 30.0); // within 30 minutes
        }
    CHECK(found);
}

TEST(Analysis, TotalSolarEclipseNewMoon) {
    // The new moon of the 2024-04-08 total solar eclipse: 18:21 UTC.
    double jd0 = 0, ref = 0;
    std::string err;
    astro::parseDate("2024-03-01", jd0, err);
    astro::parseDate("2024-04-08 18:21", ref, err);
    physics::BuildOptions o;
    o.epochJD = jd0;
    o.calibration = physics::CalibrationMode::Moon;
    physics::Simulation sim(physics::buildSystem(o), jd0);
    analysis::EventSearchOptions eo;
    eo.startJD = jd0;
    eo.endJD = jd0 + 45;
    eo.planets = eo.seasons = eo.apsides = false;
    const auto ev = analysis::findEvents(sim, eo);
    bool found = false;
    for (const auto& e : ev)
        if (e.description == "New Moon" && std::fabs(e.jd - ref) < 2.0) {
            found = true;
            CHECK_NEAR((e.jd - ref) * 1440.0, 0.0, 15.0);
        }
    CHECK(found);
}

TEST(Analysis, ObservationOfTheSun) {
    physics::BuildOptions o;
    o.epochJD = J2000_JD;
    o.calibration = physics::CalibrationMode::None;
    const auto sys = physics::buildSystem(o);
    const auto ob = analysis::observe(sys, 0, static_cast<std::size_t>(sys.indexOf("Earth")), J2000_JD);
    // 2000-01-01 12:00 TT: Sun RA ~ 18h45m, Dec ~ -23.0 deg, distance 0.9833 AU
    CHECK_NEAR(ob.raDate / 15.0, 18.75, 0.03);
    CHECK_NEAR(ob.decDate, -23.03, 0.05);
    CHECK_NEAR(ob.distance, 0.9833, 0.0005);
    CHECK_NEAR(ob.angularDiameter / 60.0, 32.5, 0.2);
}
