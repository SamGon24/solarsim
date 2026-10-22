// =============================================================================
//  SolarSim - Events.cpp
// =============================================================================
#include "Events.hpp"

#include "Observation.hpp"

#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"
#include "../core/StringUtil.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

namespace solsim::analysis {

using physics::Simulation;

namespace {

enum class TrackKind { RootRising, RootAny, Maximum, Minimum };

struct Tracker {
    TrackKind kind = TrackKind::RootAny;
    std::string category;
    std::string body;
    std::function<double(const Simulation&)> f;
    std::function<EventRecord(const Simulation&)> describe;
    double v[3] = {0, 0, 0}; // samples k-2, k-1, k
};

bool wanted(const EventSearchOptions& opt, const std::string& name) {
    if (opt.bodies.empty()) return true;
    return std::any_of(opt.bodies.begin(), opt.bodies.end(), [&](const std::string& b) { return str::iequals(b, name); });
}

bool isInner(const std::string& n) { return str::iequals(n, "Mercury") || str::iequals(n, "Venus"); }

} // namespace

std::vector<EventRecord> findEvents(const Simulation& simIn, const EventSearchOptions& opt,
                                    const std::function<void(double)>& progress) {
    std::vector<EventRecord> events;
    Simulation sim = simIn;
    const auto& sys0 = sim.system();
    const int earthI = sys0.indexOf("Earth");
    if (earthI < 0) return events;
    const std::size_t earth = static_cast<std::size_t>(earthI);
    const int moonI = sys0.indexOf("Moon");

    std::vector<Tracker> trackers;
    auto solarLon = [earth](const Simulation& s) { return apparentSolarLongitude(s.system(), earth, s.time()); };

    // ---- seasons -------------------------------------------------------------------
    if (opt.seasons) {
        static const char* names[] = {"March equinox", "June solstice", "September equinox", "December solstice"};
        for (int k = 0; k < 4; ++k) {
            Tracker t;
            t.kind = TrackKind::RootRising;
            t.category = "season";
            t.body = "Sun";
            t.f = [=](const Simulation& s) { return math::wrap180(solarLon(s) - 90.0 * k); };
            t.describe = [=](const Simulation& s) {
                EventRecord e;
                e.description = names[k];
                e.value = (s.system().pos(earth) - s.system().pos(0)).norm();
                return e;
            };
            trackers.push_back(t);
        }
    }

    // ---- lunar phases ----------------------------------------------------------------
    if (opt.moonPhases && moonI >= 0) {
        const std::size_t moon = static_cast<std::size_t>(moonI);
        static const char* names[] = {"New Moon", "First Quarter", "Full Moon", "Last Quarter"};
        for (int k = 0; k < 4; ++k) {
            Tracker t;
            t.kind = TrackKind::RootRising;
            t.category = "moon";
            t.body = "Moon";
            t.f = [=](const Simulation& s) {
                return math::wrap180(apparentLongitude(s.system(), moon, earth, s.time()) - solarLon(s) - 90.0 * k);
            };
            t.describe = [=](const Simulation& s) {
                EventRecord e;
                e.description = names[k];
                e.value = (s.system().pos(moon) - s.system().pos(earth)).norm() * constants::AU_KM;
                return e;
            };
            trackers.push_back(t);
        }
    }

    // ---- planets: conjunctions, oppositions, elongations ----------------------------
    if (opt.planets) {
        for (std::size_t i = 1; i < sys0.namedCount(); ++i) {
            const auto& inf = sys0.info(i);
            if (static_cast<int>(i) == earthI || static_cast<int>(i) == moonI) continue;
            if (inf.parent != 0) continue; // heliocentric bodies only
            if (inf.kind != astro::BodyKind::Planet && inf.kind != astro::BodyKind::DwarfPlanet) continue;
            if (!wanted(opt, inf.name)) continue;
            const std::string name = inf.name;
            auto lonDiff = [=](const Simulation& s) {
                return math::wrap180(apparentLongitude(s.system(), i, earth, s.time()) - solarLon(s));
            };
            {
                Tracker t;
                t.kind = TrackKind::RootAny;
                t.category = "conjunction";
                t.body = name;
                t.f = lonDiff;
                t.describe = [=](const Simulation& s) {
                    EventRecord e;
                    const Observation ob = observe(s.system(), i, earth, s.time());
                    const double sunDist = (s.system().pos(0) - s.system().pos(earth)).norm();
                    if (isInner(name)) e.description = name + (ob.distance < sunDist ? " inferior conjunction" : " superior conjunction");
                    else e.description = name + " in conjunction with the Sun";
                    e.value = ob.distance;
                    return e;
                };
                trackers.push_back(t);
            }
            if (!isInner(name)) {
                Tracker t;
                t.kind = TrackKind::RootAny;
                t.category = "opposition";
                t.body = name;
                t.f = [=](const Simulation& s) { return math::wrap180(lonDiff(s) - 180.0); };
                t.describe = [=](const Simulation& s) {
                    EventRecord e;
                    const Observation ob = observe(s.system(), i, earth, s.time());
                    e.description = name + " at opposition";
                    if (ob.hasMagnitude) e.description += str::format(" (mag %+.1f, %.3f AU)", ob.magnitude, ob.distance);
                    e.value = ob.distance;
                    return e;
                };
                trackers.push_back(t);
            } else {
                Tracker t;
                t.kind = TrackKind::Maximum;
                t.category = "elongation";
                t.body = name;
                t.f = [=](const Simulation& s) { return observe(s.system(), i, earth, s.time()).elongation; };
                t.describe = [=](const Simulation& s) {
                    EventRecord e;
                    const Observation ob = observe(s.system(), i, earth, s.time());
                    e.description = str::format("%s greatest %s elongation %.1f\xC2\xB0 (%s sky)", name.c_str(),
                                                ob.eveningSky ? "eastern" : "western", ob.elongation,
                                                ob.eveningSky ? "evening" : "morning");
                    e.value = ob.elongation;
                    return e;
                };
                trackers.push_back(t);
            }
        }
    }

    // ---- apsides -------------------------------------------------------------------------
    if (opt.apsides) {
        auto dist = [earth](const Simulation& s) { return (s.system().pos(earth) - s.system().pos(0)).norm(); };
        for (int k = 0; k < 2; ++k) {
            Tracker t;
            t.kind = k == 0 ? TrackKind::Minimum : TrackKind::Maximum;
            t.category = "apsis";
            t.body = "Earth";
            t.f = dist;
            t.describe = [=](const Simulation& s) {
                EventRecord e;
                e.value = dist(s);
                e.description = str::format("Earth at %s (%.6f AU)", k == 0 ? "perihelion" : "aphelion", e.value);
                return e;
            };
            trackers.push_back(t);
        }
    }
    if (trackers.empty()) return events;

    // ---- sampling loop ---------------------------------------------------------------------
    if (sim.time() < opt.startJD) sim.advanceTo(opt.startJD);
    const double h = std::max(0.01, opt.sampleDays);
    std::unique_ptr<Simulation> snap1, snap2; // snapshots at k-1 and k-2
    double t1 = 0.0, t2 = 0.0;
    int k = 0;

    auto evalAt = [](const Simulation& snap, double t, const Tracker& tr) {
        Simulation s = snap;
        s.advanceTo(t);
        return tr.f(s);
    };
    auto record = [&](const Simulation& snap, double t, const Tracker& tr) {
        Simulation s = snap;
        s.advanceTo(t);
        EventRecord e = tr.describe(s);
        e.jd = t;
        e.category = tr.category;
        e.body = tr.body;
        if (t >= opt.startJD && t <= opt.endJD) events.push_back(e);
    };

    for (double t = sim.time(); t <= opt.endJD + h; t += h, ++k) {
        if (k > 0) sim.advanceTo(t);
        if (progress && (k % 64) == 0) progress((t - opt.startJD) / std::max(1.0, opt.endJD - opt.startJD));
        for (auto& tr : trackers) {
            tr.v[0] = tr.v[1];
            tr.v[1] = tr.v[2];
            tr.v[2] = tr.f(sim);
            if (k >= 1 && (tr.kind == TrackKind::RootAny || tr.kind == TrackKind::RootRising)) {
                const double a = tr.v[1], b = tr.v[2];
                const bool crossing = (a < 0.0 && b >= 0.0) || (tr.kind == TrackKind::RootAny && a > 0.0 && b <= 0.0);
                if (crossing && std::fabs(a - b) < 90.0) {
                    auto g = [&](double tt) { return evalAt(*snap1, tt, tr); };
                    const math::RootResult r = math::brent(g, t1, t, 1e-6, 80);
                    record(*snap1, r.converged ? r.x : 0.5 * (t1 + t), tr);
                }
            }
            if (k >= 2 && (tr.kind == TrackKind::Maximum || tr.kind == TrackKind::Minimum)) {
                const double sgn = tr.kind == TrackKind::Maximum ? 1.0 : -1.0;
                const double a = sgn * tr.v[0], b = sgn * tr.v[1], c = sgn * tr.v[2];
                if (b > a && b >= c) {
                    auto g = [&](double tt) { return sgn * evalAt(*snap2, tt, tr); };
                    const double tm = math::goldenSectionMax(g, t2, t, 1e-5, 80);
                    record(*snap2, tm, tr);
                }
            }
        }
        snap2 = std::move(snap1);
        t2 = t1;
        snap1 = std::make_unique<Simulation>(sim);
        t1 = t;
    }
    std::sort(events.begin(), events.end(), [](const EventRecord& a, const EventRecord& b) { return a.jd < b.jd; });
    // Remove duplicates detected twice at sample boundaries.
    std::vector<EventRecord> unique;
    for (const auto& e : events) {
        if (!unique.empty() && unique.back().category == e.category && unique.back().body == e.body &&
            unique.back().description == e.description && std::fabs(unique.back().jd - e.jd) < 0.5)
            continue;
        unique.push_back(e);
    }
    return unique;
}

} // namespace solsim::analysis
