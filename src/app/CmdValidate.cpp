// =============================================================================
//  SolarSim - CmdValidate.cpp
//  `solarsim validate` - checks the simulator against the real world:
//    1. published dates/times of seasons, lunar phases (incl. eclipses),
//       oppositions and conjunctions
//    2. two-body motion vs the analytic Kepler solution
//    3. conservation of energy / angular momentum
//    4. Mercury's relativistic perihelion advance vs Einstein's formula
//    5. N-body vs Standish mean-element model (consistency, informative)
// =============================================================================
#include "AppContext.hpp"

#include "../analysis/Events.hpp"
#include "../analysis/Precession.hpp"
#include "../astro/BodyCatalog.hpp"
#include "../astro/JulianDate.hpp"
#include "../astro/Kepler.hpp"
#include "../core/Constants.hpp"
#include "../core/StringUtil.hpp"
#include "../core/Table.hpp"
#include "../core/Timer.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>

namespace solsim::app {

using namespace solsim::constants;

namespace {

struct Reference {
    const char* utc;          // published time (UTC); date-only entries use 12:00
    const char* description;  // prefix matched against the simulated event description
    double toleranceMin;      // PASS threshold
};

// Published values (USNO / IMCCE / NASA). Times in UTC.
const Reference kReferences[] = {
    {"2024-03-20 03:06", "March equinox", 60},      {"2024-06-20 20:51", "June solstice", 60},
    {"2024-09-22 12:44", "September equinox", 60},  {"2024-12-21 09:20", "December solstice", 60},
    {"2025-03-20 09:01", "March equinox", 60},      {"2025-06-21 02:42", "June solstice", 60},
    {"2025-09-22 18:19", "September equinox", 60},  {"2025-12-21 15:03", "December solstice", 60},
    {"2026-03-20 14:46", "March equinox", 60},      {"2026-06-21 08:24", "June solstice", 60},
    {"2026-09-23 00:05", "September equinox", 60},  {"2026-12-21 20:50", "December solstice", 60},
    {"2024-01-25 17:54", "Full Moon", 30},          {"2024-04-08 18:21", "New Moon", 30},   // total solar eclipse
    {"2024-10-17 11:26", "Full Moon", 30},          {"2025-03-14 06:55", "Full Moon", 30},  // total lunar eclipse
    {"2025-09-07 18:09", "Full Moon", 30},          {"2026-08-12 17:37", "New Moon", 30},   // total solar eclipse
    {"2023-11-03 12:00", "Jupiter at opposition", 1.5 * 1440}, {"2024-12-07 12:00", "Jupiter at opposition", 1.5 * 1440},
    {"2026-01-10 12:00", "Jupiter at opposition", 1.5 * 1440}, {"2023-08-27 12:00", "Saturn at opposition", 1.5 * 1440},
    {"2024-09-08 12:00", "Saturn at opposition", 1.5 * 1440},  {"2025-09-21 12:00", "Saturn at opposition", 1.5 * 1440},
    {"2025-01-16 12:00", "Mars at opposition", 1.5 * 1440},    {"2025-03-23 12:00", "Venus inferior conjunction", 1.5 * 1440},
};

struct Check {
    std::string name;
    bool pass = false;
    std::string detail;
};

} // namespace

int cmdValidate(AppContext& ctx) {
    std::vector<Check> checks;
    Stopwatch total;

    // ---- 1. real-world events ------------------------------------------------------------------
    {
        double from = 0, to = 0;
        std::string err;
        astro::parseDate("2023-06-01", from, err);
        astro::parseDate("2027-01-01", to, err);
        const double epochJD = ctx.args().has("epoch") ? ctx.epoch() : from;
        if (epochJD > from) throw UsageError("--epoch must be on or before 2023-06-01 for validation");
        auto sim = ctx.makeSimulation(epochJD);
        analysis::EventSearchOptions opt;
        opt.startJD = from;
        opt.endJD = to;
        opt.apsides = false;
        const auto events = analysis::findEvents(*sim, opt, [&](double f) {
            if (!ctx.quiet()) printProgress("searching events", f);
        });
        if (!ctx.quiet()) finishProgress();

        TextTable t({"Event", "Published (UTC)", "Simulated (UTC)", "Error", "Status"});
        int pass = 0, n = 0;
        double sumAbs[3] = {0, 0, 0};
        int cnt[3] = {0, 0, 0};
        for (const auto& ref : kReferences) {
            double jr = 0;
            astro::parseDate(ref.utc, jr, err);
            const analysis::EventRecord* best = nullptr;
            for (const auto& e : events)
                if (str::startsWith(e.description, ref.description) && std::fabs(e.jd - jr) < 5.0)
                    if (!best || std::fabs(e.jd - jr) < std::fabs(best->jd - jr)) best = &e;
            ++n;
            if (!best) {
                t.addRow({ref.description, ref.utc, "not found", "-", "FAIL"});
                continue;
            }
            const double errMin = (best->jd - jr) * 1440.0;
            const bool ok = std::fabs(errMin) <= ref.toleranceMin;
            pass += ok;
            const bool dateOnly = ref.toleranceMin >= 1440.0;
            const int cat = dateOnly ? 2 : (str::startsWith(ref.description, "Full") || str::startsWith(ref.description, "New") ? 1 : 0);
            sumAbs[cat] += std::fabs(errMin);
            ++cnt[cat];
            t.addRow({ref.description, dateOnly ? std::string(ref.utc).substr(0, 10) : ref.utc,
                      astro::formatUTC(best->jd, false).substr(0, 16),
                      dateOnly ? str::format("%+.2f d", errMin / 1440.0) : str::format("%+.1f min", errMin), ok ? "PASS" : "FAIL"});
        }
        t.setTitle("1) Published event times vs N-body simulation (epoch " + astro::formatDateOnly(epochJD) + ")");
        t.print(std::cout, TextTable::defaultStyle());
        std::printf("mean |error|: seasons %.1f min, lunar phases %.1f min, planetary events %.2f days\n\n",
                    cnt[0] ? sumAbs[0] / cnt[0] : 0.0, cnt[1] ? sumAbs[1] / cnt[1] : 0.0, cnt[2] ? sumAbs[2] / cnt[2] / 1440.0 : 0.0);
        checks.push_back({"published events", pass == n, str::format("%d/%d within tolerance", pass, n)});
    }

    // ---- 2. two-body motion vs analytic Kepler propagation ------------------------------------------
    {
        physics::BuildOptions bo = ctx.buildOptions(J2000_JD);
        bo.only = {"Sun", "Mercury"};
        bo.moons = bo.minor = false;
        bo.populations = {};
        bo.calibration = physics::CalibrationMode::None;
        const physics::NBodySystem sys = physics::buildSystem(bo);
        physics::SimulationConfig cfg;
        cfg.integrator = "yoshida4";
        cfg.dt = 0.05;
        physics::Simulation sim(sys, J2000_JD, cfg);
        const Vec3 r0 = sys.pos(1) - sys.pos(0), v0 = sys.vel(1) - sys.vel(0);
        const double mu = sys.gm()[0] + sys.gm()[1];
        sim.advance(10.0 * JULIAN_YEAR_D);
        const astro::UniversalResult k = astro::propagateKepler(r0, v0, mu, 10.0 * JULIAN_YEAR_D);
        const double errKm = ((sim.system().pos(1) - sim.system().pos(0)) - k.r).norm() * AU_KM;
        checks.push_back({"two-body vs Kepler (Mercury, 10 yr)", errKm < 10.0, str::format("position difference %.3f km", errKm)});
        std::printf("2) Two-body integration vs analytic universal-variable Kepler solution: %.3f km after 10 years\n", errKm);
    }

    // ---- 3. conservation laws --------------------------------------------------------------------------
    {
        const physics::NBodySystem sys = ctx.buildSystem(J2000_JD);
        physics::SimulationConfig cfg;
        cfg.integrator = "yoshida4";
        physics::Simulation sim(sys, J2000_JD, cfg);
        double maxE = 0;
        for (int k = 0; k < 50; ++k) {
            sim.advance(JULIAN_YEAR_D);
            maxE = std::max(maxE, std::fabs(sim.invariantErrors().energyRel));
        }
        const auto e = sim.invariantErrors();
        std::printf("3) Conservation over 50 years (yoshida4, dt %.3f d): max |dE/E| = %.2e, |dL|/L = %.2e, |dP| = %.1e\n",
                    sim.timeStep(), maxE, e.angularMomentumRel, e.momentumAbs);
        checks.push_back({"energy conservation (50 yr)", maxE < 1e-8, str::format("max |dE/E| = %.2e", maxE)});
        checks.push_back({"angular momentum conservation", e.angularMomentumRel < 1e-12, str::format("|dL|/L = %.2e", e.angularMomentumRel)});
    }

    // ---- 4. general relativity ----------------------------------------------------------------------------
    {
        physics::BuildOptions bo = ctx.buildOptions(J2000_JD);
        bo.only = {"Sun", "Mercury"};
        bo.populations = {};
        bo.calibration = physics::CalibrationMode::None;
        const physics::NBodySystem sys = physics::buildSystem(bo);
        analysis::PrecessionConfig pc;
        pc.years = 50;
        pc.dt = 0.04;
        const double newton = analysis::measurePrecession(sys, pc).rateArcsecPerCentury;
        pc.relativity = true;
        const double gr = analysis::measurePrecession(sys, pc).rateArcsecPerCentury;
        const double theory = analysis::grPrecessionArcsecPerCentury(sys.gm()[0], 0.38709927, 0.20563593);
        const double measured = gr - newton;
        std::printf("4) Mercury GR perihelion advance: measured %.3f\"/cy, Einstein %.3f\"/cy (%.2f%%)\n", measured, theory,
                    100.0 * (measured - theory) / theory);
        checks.push_back({"GR perihelion advance", std::fabs(measured - theory) < 0.01 * theory,
                          str::format("%.3f vs %.3f \"/cy", measured, theory)});
    }

    // ---- 5. consistency with the Standish mean-element model -----------------------------------------
    {
        const physics::NBodySystem sys = ctx.buildSystem(J2000_JD);
        physics::Simulation sim(sys, J2000_JD, ctx.simulationConfig());
        TextTable t({"Planet", "1 yr", "5 yr", "10 yr", "Standish accuracy"});
        const char* names[] = {"Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"};
        const char* acc[] = {"15\"", "20\"", "20\"", "40\"", "400\"", "600\"", "50\"", "10\""};
        std::vector<std::vector<double>> diff(8);
        double t0 = 0;
        for (double yr : {1.0, 5.0, 10.0}) {
            sim.advance((yr - t0) * JULIAN_YEAR_D);
            t0 = yr;
            const auto& s = sim.system();
            for (int k = 0; k < 8; ++k) {
                const int i = s.indexOf(names[k]);
                if (i < 0) { diff[static_cast<std::size_t>(k)].push_back(0.0); continue; }
                const auto* def = astro::BodyCatalog::find(names[k]);
                const auto el = astro::meanElementsAt(def->mean, sim.time(), s.gm()[0] + s.gm()[static_cast<std::size_t>(i)]);
                Vec3 nb = s.pos(static_cast<std::size_t>(i)) - s.pos(0);
                if (std::string(names[k]) == "Earth" && s.indexOf("Moon") >= 0) {
                    const std::size_t e = static_cast<std::size_t>(i), m = static_cast<std::size_t>(s.indexOf("Moon"));
                    nb = (s.pos(e) * s.gm()[e] + s.pos(m) * s.gm()[m]) / (s.gm()[e] + s.gm()[m]) - s.pos(0);
                }
                diff[static_cast<std::size_t>(k)].push_back(angleBetween(nb, astro::elementsToState(el, sim.time()).r) * RAD2ARCSEC);
            }
        }
        for (int k = 0; k < 8; ++k) {
            const auto& d = diff[static_cast<std::size_t>(k)];
            t.addRow({names[k], str::format("%.0f\"", d[0]), str::format("%.0f\"", d[1]), str::format("%.0f\"", d[2]), acc[k]});
        }
        t.setTitle("\n5) Heliocentric angular difference N-body vs Standish mean elements (informative)");
        t.print(std::cout, TextTable::defaultStyle());
        std::printf("Mean elements omit periodic perturbations (e.g. the Jupiter-Saturn 'great inequality'), so\n"
                    "differences of the order of the stated Standish accuracy are expected.\n");
    }

    // ---- summary ---------------------------------------------------------------------------------------------
    std::printf("\nSummary (%s):\n", str::formatDuration(total.seconds()).c_str());
    int failed = 0;
    for (const auto& c : checks) {
        std::printf("  [%s] %-38s %s\n", c.pass ? "PASS" : "FAIL", c.name.c_str(), c.detail.c_str());
        failed += !c.pass;
    }
    return failed == 0 ? 0 : 1;
}

} // namespace solsim::app
