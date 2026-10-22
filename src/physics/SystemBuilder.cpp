// =============================================================================
//  SolarSim - SystemBuilder.cpp
// =============================================================================
#include "SystemBuilder.hpp"

#include "../astro/Coordinates.hpp"
#include "../astro/OrbitalElements.hpp"
#include "../core/Constants.hpp"
#include "../core/Log.hpp"
#include "../core/MathUtil.hpp"
#include "../core/Random.hpp"
#include "../core/StringUtil.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

namespace solsim::physics {

using namespace solsim::constants;
using astro::BodyDefinition;

BodyInfo makeBodyInfo(const BodyDefinition& def) {
    BodyInfo info;
    info.name = def.name;
    info.kind = def.kind;
    info.radius = def.radiusKm / AU_KM;
    info.rotationPeriod = def.rotationHours / 24.0;
    // Cosmetic initial rotation phase derived from the name (deterministic).
    std::uint64_t h = 1469598103934665603ULL;
    for (char c : def.name) h = (h ^ static_cast<unsigned char>(c)) * 1099511628211ULL;
    info.rotationPhase0 = double(h % 6283) / 1000.0;
    info.pole = astro::poleFromRaDec(def.poleRaDeg, def.poleDecDeg);
    info.albedo = def.albedo;
    info.typicalDistance = def.typicalDistanceAU();
    info.style = def.style;
    info.color = def.color;
    info.color2 = def.color2;
    info.rings = def.rings;
    info.approximate = def.approximate;
    info.notes = def.notes;
    return info;
}

namespace {

struct Node {
    BodyDefinition def;
    double gm = 0.0;         // effective GM (AU^3/d^2) after folding
    astro::StateVector rel;  // relative to parent (system barycentre -> parent)
    astro::StateVector abs;  // heliocentric
    bool hasState = false;
};

bool contains(const std::vector<std::string>& v, const std::string& name) {
    return std::any_of(v.begin(), v.end(), [&](const std::string& s) { return str::iequals(s, name); });
}

} // namespace

NBodySystem buildSystem(const BuildOptions& opt, BuildReport* report) {
    BuildReport localReport;
    BuildReport& rep = report ? *report : localReport;

    // ---- 1. candidate definitions -----------------------------------------------
    std::vector<BodyDefinition> pool = astro::BodyCatalog::all();
    for (const auto& ex : opt.extra) {
        auto it = std::find_if(pool.begin(), pool.end(), [&](const BodyDefinition& d) { return str::iequals(d.name, ex.name); });
        if (it != pool.end()) {
            *it = ex;
            rep.notes.push_back("custom definition overrides built-in '" + ex.name + "'");
        } else {
            pool.push_back(ex);
        }
    }
    auto findDef = [&](const std::string& name) -> const BodyDefinition* {
        for (const auto& d : pool)
            if (str::iequals(d.name, name)) return &d;
        return nullptr;
    };

    std::vector<std::string> wanted = astro::BodyCatalog::defaultSet();
    if (!opt.includeMoon) wanted.erase(std::remove(wanted.begin(), wanted.end(), "Moon"), wanted.end());
    if (opt.moons)
        for (const auto& n : astro::BodyCatalog::majorMoonSet()) wanted.push_back(n);
    if (opt.minor)
        for (const auto& n : astro::BodyCatalog::minorBodySet()) wanted.push_back(n);
    for (const auto& ex : opt.extra)
        if (!contains(wanted, ex.name)) wanted.push_back(ex.name);

    if (!opt.only.empty()) {
        for (const auto& n : opt.only)
            if (!findDef(n)) rep.warnings.push_back("--only: unknown body '" + n + "' ignored");
        std::vector<std::string> filtered;
        for (const auto& n : opt.only)
            if (findDef(n)) filtered.push_back(findDef(n)->name);
        wanted = filtered;
    }
    for (const auto& n : opt.exclude) {
        if (str::iequals(n, "Sun")) {
            rep.warnings.push_back("the Sun cannot be excluded");
            continue;
        }
        wanted.erase(std::remove_if(wanted.begin(), wanted.end(), [&](const std::string& s) { return str::iequals(s, n); }),
                     wanted.end());
    }
    if (!contains(wanted, "Sun")) wanted.insert(wanted.begin(), "Sun");

    // Add missing parents (transitively).
    for (bool changed = true; changed;) {
        changed = false;
        for (std::size_t k = 0; k < wanted.size(); ++k) {
            const BodyDefinition* d = findDef(wanted[k]);
            if (!d) continue;
            if (!d->parent.empty() && !contains(wanted, d->parent)) {
                if (!findDef(d->parent))
                    throw std::runtime_error("body '" + d->name + "' has unknown parent '" + d->parent + "'");
                wanted.push_back(findDef(d->parent)->name);
                rep.notes.push_back("added '" + d->parent + "' (parent of " + d->name + ")");
                changed = true;
            }
        }
    }

    // ---- 2. nodes in topological order (parents first) --------------------------
    std::vector<Node> nodes;
    std::set<std::string> placed;
    for (std::size_t guard = 0; nodes.size() < wanted.size() && guard < wanted.size() + 2; ++guard) {
        for (const auto& name : wanted) {
            const BodyDefinition* d = findDef(name);
            if (!d) {
                rep.warnings.push_back("unknown body '" + name + "' skipped");
                placed.insert(str::toLower(name));
                continue;
            }
            const std::string key = str::toLower(d->name);
            if (placed.count(key)) continue;
            if (!d->parent.empty() && !placed.count(str::toLower(d->parent))) continue;
            Node nd;
            nd.def = *d;
            nd.gm = d->gm();
            nodes.push_back(nd);
            placed.insert(key);
        }
    }
    if (nodes.empty() || !str::iequals(nodes.front().def.name, "Sun"))
        throw std::runtime_error("system must contain the Sun as its root");

    auto nodeIndex = [&](const std::string& name) -> int {
        for (std::size_t i = 0; i < nodes.size(); ++i)
            if (str::iequals(nodes[i].def.name, name)) return static_cast<int>(i);
        return -1;
    };

    // ---- 3. fold the mass of non-simulated catalog moons into their parent --------
    for (const auto& d : astro::BodyCatalog::all()) {
        if (d.kind != astro::BodyKind::Moon || nodeIndex(d.name) >= 0) continue;
        const int p = nodeIndex(d.parent);
        if (p < 0) continue;
        nodes[static_cast<std::size_t>(p)].gm += d.gm();
        rep.notes.push_back("GM of " + d.name + " folded into " + d.parent);
    }

    // ---- 4. relative states ---------------------------------------------------------
    for (std::size_t i = 1; i < nodes.size(); ++i) {
        Node& nd = nodes[i];
        const int p = nodeIndex(nd.def.parent);
        const Node& parent = nodes[static_cast<std::size_t>(p)];
        // For a planet orbiting the Sun the two-body mu uses the whole planetary system.
        double muSystem = nd.gm;
        for (const auto& c : nodes)
            if (str::iequals(c.def.parent, nd.def.name)) muSystem += c.gm;
        const double mu = parent.gm + muSystem;
        const BodyDefinition* parentDef = &parent.def;
        nd.rel = astro::relativeStateAt(nd.def, opt.epochJD, mu, parentDef);
        if (!nd.rel.r.isFinite() || !nd.rel.v.isFinite())
            throw std::runtime_error("could not compute initial state of '" + nd.def.name + "'");
    }

    // ---- 5. heliocentric states with system-barycentre corrections ------------------
    nodes[0].abs = {};
    for (std::size_t i = 1; i < nodes.size(); ++i) {
        Node& nd = nodes[i];
        const Node& parent = nodes[static_cast<std::size_t>(nodeIndex(nd.def.parent))];
        // The orbit model describes the barycentre of nd's own subsystem.
        double mSys = nd.gm;
        Vec3 shiftR{}, shiftV{};
        for (const auto& c : nodes) {
            if (!str::iequals(c.def.parent, nd.def.name) || c.gm <= 0.0) continue;
            mSys += c.gm;
            shiftR += c.rel.r * c.gm;
            shiftV += c.rel.v * c.gm;
        }
        nd.abs.r = parent.abs.r + nd.rel.r;
        nd.abs.v = parent.abs.v + nd.rel.v;
        if (mSys > nd.gm && mSys > 0.0) {
            nd.abs.r -= shiftR / mSys;
            nd.abs.v -= shiftV / mSys;
        }
    }

    // ---- 6. assemble the system -----------------------------------------------------
    NBodySystem sys;
    for (const auto& nd : nodes) {
        BodyInfo info = makeBodyInfo(nd.def);
        info.parent = nd.def.parent.empty() ? -1 : sys.indexOf(nd.def.parent);
        sys.addBody(info, nd.gm, nd.abs.r, nd.abs.v);
        if (nd.def.approximate) rep.notes.push_back(nd.def.name + ": approximate initial conditions");
    }

    sys.moveToBarycentricFrame();

    // Differential correction of the initial velocities (see Calibration.hpp).
    if (opt.calibration != CalibrationMode::None) {
        rep.calibration = calibrateInitialConditions(sys, opt.epochJD, opt.calibration, 2);
        for (const auto& e : rep.calibration.entries)
            rep.notes.push_back(str::format("calibrated %s: mean-motion drift %+.1f\"/yr corrected (v x %.8f)",
                                            e.body.c_str(), e.driftBeforeArcsecPerYear, e.velocityScale));
    }

    // Optional intruder star on a hyperbolic flyby.
    if (opt.intruder.enabled) {
        const auto& in = opt.intruder;
        const double gmStar = in.massSolar * GM_SUN;
        const double mu = sys.gm()[0] + gmStar;
        const double vinf = in.vInfKms * DAY_S / AU_KM;
        astro::OrbitalElements el;
        el.q = in.perihelionAU;
        el.e = 1.0 + in.perihelionAU * vinf * vinf / mu;
        el.i = in.inclinationDeg * DEG2RAD;
        el.node = 1.2;
        el.argp = 0.7;
        el.tp = opt.epochJD + in.daysToPerihelion;
        el.mu = mu;
        const astro::StateVector s = astro::elementsToState(el, opt.epochJD);
        BodyInfo info;
        info.name = "Intruder";
        info.kind = astro::BodyKind::Star;
        info.parent = 0;
        info.radius = 0.6 * SUN_RADIUS_KM / AU_KM;
        info.style = astro::SurfaceStyle::Star;
        info.color = 0xFFB070;
        info.color2 = 0xFF6020;
        info.typicalDistance = in.perihelionAU;
        info.notes = str::format("%.2f Msun star, q = %.1f AU, v_inf = %.1f km/s", in.massSolar, in.perihelionAU, in.vInfKms);
        sys.addBody(info, gmStar, sys.pos(0) + s.r, sys.vel(0) + s.v);
        rep.notes.push_back("intruder star: e = " + str::format("%.3f", el.e) + ", r0 = " + str::format("%.0f AU", s.r.norm()));
    }

    sys.moveToBarycentricFrame();

    // ---- 7. synthetic populations (test particles) ---------------------------------
    const auto& pop = opt.populations;
    if (pop.asteroids > 0) addMainBelt(sys, pop.asteroids, pop.seed, opt.epochJD);
    if (pop.trojans > 0) {
        if (sys.indexOf("Jupiter") < 0) rep.warnings.push_back("Trojans need Jupiter - skipped");
        else addJupiterTrojans(sys, pop.trojans, pop.seed + 1, opt.epochJD);
    }
    if (pop.kuiper > 0) addKuiperBelt(sys, pop.kuiper, pop.seed + 2, opt.epochJD);
    return sys;
}

// -----------------------------------------------------------------------------
// Populations
// -----------------------------------------------------------------------------
namespace {

void addFromElements(NBodySystem& sys, int group, const astro::OrbitalElements& el, double jd) {
    const astro::StateVector s = astro::elementsToState(el, jd);
    sys.addParticle(group, sys.pos(0) + s.r, sys.vel(0) + s.v);
}

// Mean longitude of a body's osculating heliocentric orbit.
bool heliocentricElements(const NBodySystem& sys, const std::string& name, double jd, astro::OrbitalElements& out) {
    const int idx = sys.indexOf(name);
    if (idx < 0) return false;
    const auto i = static_cast<std::size_t>(idx);
    out = astro::stateToElements(sys.pos(i) - sys.pos(0), sys.vel(i) - sys.vel(0), sys.gm()[0] + sys.gm()[i], jd);
    return true;
}

} // namespace

void addMainBelt(NBodySystem& sys, int count, std::uint64_t seed, double jd) {
    Random rng(seed);
    const int g = sys.addParticleGroup("main belt", 0x9A8F80);
    const double mu = sys.gm()[0];
    // Kirkwood gaps (mean-motion resonances with Jupiter): centre, half width [AU]
    const double gaps[][2] = {{2.065, 0.02}, {2.502, 0.03}, {2.825, 0.02}, {2.958, 0.015}, {3.279, 0.04}};
    int made = 0;
    while (made < count) {
        const double a = rng.uniform(2.1, 3.35);
        bool inGap = false;
        for (const auto& gp : gaps)
            if (std::fabs(a - gp[0]) < gp[1]) inGap = true;
        if (inGap) continue;
        const double e = std::min(0.35, rng.rayleigh(0.09));
        const double inc = std::min(35.0, rng.rayleigh(7.5)) * DEG2RAD;
        const auto el = astro::OrbitalElements::fromMeanAnomaly(a, e, inc, rng.angle(), rng.angle(), rng.angle(), jd, mu);
        addFromElements(sys, g, el, jd);
        ++made;
    }
}

void addJupiterTrojans(NBodySystem& sys, int count, std::uint64_t seed, double jd) {
    astro::OrbitalElements jup;
    if (!heliocentricElements(sys, "Jupiter", jd, jup)) return;
    Random rng(seed);
    const int g = sys.addParticleGroup("Jupiter Trojans", 0x7FB36B);
    const double mu = sys.gm()[0];
    const double lamJ = jup.meanLongitudeAt(jd);
    const double aJ = jup.a();
    for (int k = 0; k < count; ++k) {
        const double sign = (k % 2 == 0) ? 1.0 : -1.0; // L4 leads, L5 trails
        const double a = aJ * (1.0 + rng.normal(0.0, 0.008));
        const double e = std::min(0.2, rng.rayleigh(0.05));
        const double inc = std::min(40.0, rng.rayleigh(10.0)) * DEG2RAD;
        const double node = rng.angle();
        const double argp = rng.angle();
        const double lambda = lamJ + sign * 60.0 * DEG2RAD + rng.normal(0.0, 9.0) * DEG2RAD;
        const double M = lambda - node - argp;
        const auto el = astro::OrbitalElements::fromMeanAnomaly(a, e, inc, node, argp, math::wrapPi(M), jd, mu);
        addFromElements(sys, g, el, jd);
    }
}

void addKuiperBelt(NBodySystem& sys, int count, std::uint64_t seed, double jd) {
    Random rng(seed);
    const int g = sys.addParticleGroup("Kuiper belt", 0x6C8DC5);
    const double mu = sys.gm()[0];
    astro::OrbitalElements nep;
    const bool haveNeptune = heliocentricElements(sys, "Neptune", jd, nep);
    for (int k = 0; k < count; ++k) {
        const double node = rng.angle();
        const double argp = rng.angle();
        if (haveNeptune && rng.uniform() < 0.25) {
            // Plutino in the 3:2 mean-motion resonance: resonant angle
            // phi = 3 lambda - 2 lambda_N - varpi librates about 180 deg.
            const double a = nep.a() * std::pow(1.5, 2.0 / 3.0) + rng.normal(0.0, 0.1);
            const double e = rng.uniform(0.08, 0.3);
            const double inc = std::min(40.0, rng.rayleigh(10.0)) * DEG2RAD;
            const double varpi = node + argp;
            const double phi = constants::PI + rng.normal(0.0, 35.0) * DEG2RAD;
            const double lambda = (phi + 2.0 * nep.meanLongitudeAt(jd) + varpi) / 3.0 +
                                  rng.uniformInt(0, 2) * (constants::TWO_PI / 3.0);
            const auto el = astro::OrbitalElements::fromMeanAnomaly(a, e, inc, node, argp,
                                                                    math::wrapPi(lambda - varpi), jd, mu);
            addFromElements(sys, g, el, jd);
        } else {
            // Classical (cold + hot) population.
            const double a = rng.uniform(42.0, 47.8);
            const double e = std::min(0.2, rng.rayleigh(0.05));
            const double sigI = rng.uniform() < 0.5 ? 2.5 : 14.0;
            const double inc = std::min(45.0, rng.rayleigh(sigI)) * DEG2RAD;
            const auto el = astro::OrbitalElements::fromMeanAnomaly(a, e, inc, node, argp, rng.angle(), jd, mu);
            addFromElements(sys, g, el, jd);
        }
    }
}

} // namespace solsim::physics
