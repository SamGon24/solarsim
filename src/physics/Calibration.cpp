// =============================================================================
//  SolarSim - Calibration.cpp
// =============================================================================
#include "Calibration.hpp"

#include "Simulation.hpp"

#include "../astro/BodyCatalog.hpp"
#include "../astro/LunarTheory.hpp"
#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"
#include "../core/StringUtil.hpp"
#include "../core/Timer.hpp"

#include <cmath>

namespace solsim::physics {

using namespace solsim::constants;

bool parseCalibrationMode(const std::string& s, CalibrationMode& out) {
    const std::string k = str::toLower(s);
    if (k == "none" || k == "off" || k == "0") out = CalibrationMode::None;
    else if (k == "moon") out = CalibrationMode::Moon;
    else if (k == "all" || k == "full" || k == "planets") out = CalibrationMode::All;
    else return false;
    return true;
}

const char* calibrationModeName(CalibrationMode m) {
    switch (m) {
        case CalibrationMode::None: return "none";
        case CalibrationMode::Moon: return "moon";
        case CalibrationMode::All: return "all";
    }
    return "?";
}

namespace {

struct Target {
    std::size_t index = 0;
    bool isMoon = false;
    bool isBarycentric = false;      // Earth: compare the Earth-Moon barycentre
    const astro::BodyDefinition* def = nullptr;
    double meanMotion = 0.0;         // rad/day (reference)
    double scale = 1.0;
    double firstDrift = 0.0;
    double lastDrift = 0.0;
    std::vector<double> t, dl;
};

// Body (or EMB) position relative to the reference centre used by its theory.
Vec3 measured(const NBodySystem& s, const Target& tg, int moon) {
    if (tg.isMoon) return s.pos(tg.index) - s.pos(static_cast<std::size_t>(s.info(tg.index).parent));
    if (tg.isBarycentric && moon >= 0) {
        const std::size_t m = static_cast<std::size_t>(moon);
        const double gE = s.gm()[tg.index], gM = s.gm()[m];
        return (s.pos(tg.index) * gE + s.pos(m) * gM) / (gE + gM) - s.pos(0);
    }
    return s.pos(tg.index) - s.pos(0);
}

Vec3 reference(const NBodySystem& s, const Target& tg, double jd) {
    if (tg.isMoon) return astro::lunarPositionJ2000(jd);
    const auto el = astro::meanElementsAt(tg.def->mean, jd, s.gm()[0] + s.gm()[tg.index]);
    return astro::elementsToState(el, jd).r;
}

// Scales the body's velocity relative to its parent, preserving the momentum
// of its local subsystem (Earth-Moon pair or planet + moons).
void scaleVelocity(NBodySystem& sys, const Target& tg, double factor, int moon) {
    auto& v = sys.velocities();
    const auto& gm = sys.gm();
    if (tg.isMoon) {
        const std::size_t m = tg.index;
        const std::size_t e = static_cast<std::size_t>(sys.info(m).parent);
        const double gE = gm[e], gM = gm[m];
        const Vec3 vBary = (v[e] * gE + v[m] * gM) / (gE + gM);
        const Vec3 rel = (v[m] - v[e]) * factor;
        v[e] = vBary - rel * (gM / (gE + gM));
        v[m] = vBary + rel * (gE / (gE + gM));
        return;
    }
    // Planet: scale the heliocentric velocity of its system barycentre and move
    // all of its satellites with it.
    const std::size_t p = tg.index;
    Vec3 vSys = v[p] * gm[p];
    double mSys = gm[p];
    std::vector<std::size_t> members{p};
    for (std::size_t i = 0; i < sys.namedCount(); ++i)
        if (sys.info(i).parent == static_cast<int>(p)) {
            members.push_back(i);
            vSys += v[i] * gm[i];
            mSys += gm[i];
        }
    (void)moon;
    if (mSys <= 0.0) {
        v[p] = sys.vel(0) + (v[p] - sys.vel(0)) * factor;
        return;
    }
    const Vec3 vBary = vSys / mSys;
    const Vec3 newBary = sys.vel(0) + (vBary - sys.vel(0)) * factor;
    const Vec3 dv = newBary - vBary;
    for (std::size_t i : members) v[i] += dv;
}

} // namespace

namespace {

// Reduced copy for the calibration runs: the Sun, bodies orbiting the Sun
// directly, and Earth's Moon. Masses of omitted satellites are folded into
// their parents so the heliocentric dynamics are unchanged. Small moons with
// very short periods would otherwise force a tiny time step.
NBodySystem trimmedCopy(const NBodySystem& full) {
    NBodySystem out;
    for (std::size_t i = 0; i < full.namedCount(); ++i) {
        const auto& inf = full.info(i);
        const bool keep = i == 0 || inf.parent == 0 || str::iequals(inf.name, "Moon");
        if (!keep) continue;
        double gm = full.gm()[i];
        for (std::size_t j = 0; j < full.namedCount(); ++j) {
            const auto& c = full.info(j);
            if (c.parent == static_cast<int>(i) && !str::iequals(c.name, "Moon") && i != 0) gm += full.gm()[j];
        }
        BodyInfo copy = inf;
        copy.parent = inf.parent < 0 ? -1 : out.indexOf(full.info(static_cast<std::size_t>(inf.parent)).name);
        if (inf.parent >= 0 && copy.parent < 0) continue;
        out.addBody(copy, gm, full.pos(i), full.vel(i));
    }
    return out;
}

} // namespace

CalibrationReport calibrateInitialConditions(NBodySystem& sys, double epochJD, CalibrationMode mode, int iterations,
                                             const std::function<void(const std::string&)>& log) {
    CalibrationReport report;
    if (mode == CalibrationMode::None) return report;
    Stopwatch sw;
    // Standish's Table 1 is only valid for 1800-2050: outside it, calibrate the Moon only.
    const double year = 2000.0 + (epochJD - J2000_JD) / JULIAN_YEAR_D;
    if (mode == CalibrationMode::All && (year < 1800.0 || year > 2050.0)) {
        if (log) log("calibration: epoch outside 1800-2050, planets left uncalibrated");
        mode = CalibrationMode::Moon;
    }
    NBodySystem work = trimmedCopy(sys);
    const int moon = work.indexOf("Moon");
    const int earth = work.indexOf("Earth");

    std::vector<Target> targets;
    if (moon >= 0 && earth >= 0 && work.info(static_cast<std::size_t>(moon)).parent == earth) {
        const auto* d = astro::BodyCatalog::find("Moon");
        if (d && d->model == astro::OrbitModel::LunarTheory) {
            Target t;
            t.index = static_cast<std::size_t>(moon);
            t.isMoon = true;
            t.def = d;
            t.meanMotion = 481267.88123421 * DEG2RAD / JULIAN_CENTURY_D;
            targets.push_back(t);
        }
    }
    if (mode == CalibrationMode::All) {
        for (std::size_t i = 1; i < work.namedCount(); ++i) {
            if (work.info(i).parent != 0 || work.gm()[i] <= 0.0) continue;
            const auto* d = astro::BodyCatalog::find(work.info(i).name);
            if (!d || d->model != astro::OrbitModel::MeanElements) continue;
            Target t;
            t.index = i;
            t.def = d;
            t.isBarycentric = static_cast<int>(i) == earth;
            t.meanMotion = d->mean.dL * DEG2RAD / JULIAN_CENTURY_D;
            targets.push_back(t);
        }
    }
    if (targets.empty()) return report;

    // Span: one year resolves the Moon well; planets need longer baselines.
    const double span = mode == CalibrationMode::All ? 30.0 * JULIAN_YEAR_D : 365.25;
    const double step = mode == CalibrationMode::All ? 10.0 : 2.0;

    for (int it = 0; it < iterations; ++it) {
        SimulationConfig cfg;
        cfg.threads = 1;
        Simulation sim(work, epochJD, cfg);
        for (auto& tg : targets) {
            tg.t.clear();
            tg.dl.clear();
        }
        for (double t = 0.0; t <= span + 1e-9; t += step) {
            if (t > 0.0) sim.advance(step);
            for (auto& tg : targets) {
                const Vec3 a = measured(sim.system(), tg, moon);
                const Vec3 b = reference(sim.system(), tg, sim.time());
                tg.t.push_back(t);
                tg.dl.push_back(std::atan2(a.y, a.x) - std::atan2(b.y, b.x));
            }
        }
        for (auto& tg : targets) {
            math::unwrapAngles(tg.dl);
            const math::LinearFit fit = math::fitLine(tg.t, tg.dl);
            const double drift = fit.slope; // rad/day, N-body minus reference
            const double driftArcsecYr = drift * RAD2ARCSEC * JULIAN_YEAR_D;
            if (it == 0) tg.firstDrift = driftArcsecYr;
            tg.lastDrift = driftArcsecYr;
            // dn/n = -3 dv/v  ->  to remove a drift dn the velocity changes by dv/v = +dn / 3n
            const double factor = 1.0 + drift / (3.0 * tg.meanMotion);
            tg.scale *= factor;
            scaleVelocity(work, tg, factor, moon);
            if (log)
                log(str::format("calibration pass %d: %-8s drift %+9.2f\"/yr -> velocity x %.8f", it + 1,
                                work.info(tg.index).name.c_str(), driftArcsecYr, factor));
        }
        work.moveToBarycentricFrame();
    }

    // Apply the cumulative factors to the full system (targets matched by name).
    const int fullMoon = sys.indexOf("Moon");
    for (const auto& tg : targets) {
        Target full = tg;
        const int idx = sys.indexOf(work.info(tg.index).name);
        if (idx < 0) continue;
        full.index = static_cast<std::size_t>(idx);
        scaleVelocity(sys, full, tg.scale, fullMoon);
        CalibrationReport::Entry e;
        e.body = work.info(tg.index).name;
        e.driftBeforeArcsecPerYear = tg.firstDrift;
        e.driftAfterArcsecPerYear = tg.lastDrift;
        e.velocityScale = tg.scale;
        report.entries.push_back(e);
    }
    sys.moveToBarycentricFrame();
    report.wallSeconds = sw.seconds();
    return report;
}

} // namespace solsim::physics
