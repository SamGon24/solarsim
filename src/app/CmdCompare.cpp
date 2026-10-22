// =============================================================================
//  SolarSim - CmdCompare.cpp
//  `solarsim compare` - integrator benchmark: every method integrates the same
//  initial conditions; errors are measured against a high-accuracy reference
//  (Dormand-Prince with rtol = 1e-13).
// =============================================================================
#include "AppContext.hpp"

#include "../core/Constants.hpp"
#include "../core/StringUtil.hpp"
#include "../core/Table.hpp"
#include "../core/Timer.hpp"
#include "../io/CsvWriter.hpp"
#include "../io/SvgWriter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>

namespace solsim::app {

using namespace solsim::constants;

namespace {

struct RunResult {
    std::string name, display;
    int order = 0;
    bool symplectic = false;
    std::uint64_t steps = 0, evals = 0;
    double wall = 0.0;
    double maxEnergy = 0.0, finalEnergy = 0.0, angMom = 0.0;
    std::vector<std::pair<double, double>> energyHistory; // (years, |dE/E|)
    std::vector<double> posErrKm;                         // per tracked body
    bool failed = false;
    std::string error;
};

} // namespace

int cmdCompare(AppContext& ctx) {
    const double years = ctx.args().getDouble("years", 10.0);
    if (!(years > 0.0 && years <= 2000.0)) throw UsageError("--years must be in (0, 2000]");
    const double jd = ctx.args().has("epoch") ? ctx.epoch() : J2000_JD;
    const physics::NBodySystem sys = ctx.buildSystem(jd);
    physics::SimulationConfig base = ctx.simulationConfig();
    const double dt = base.dt > 0.0 ? base.dt : physics::Simulation::autoTimeStep(sys);

    std::vector<std::string> names = physics::integratorNames();
    if (ctx.args().has("integrators") && str::toLower(ctx.args().getString("integrators")) != "all")
        names = ctx.args().getList("integrators");

    const std::vector<std::string> tracked = {"Mercury", "Earth", "Moon", "Jupiter"};
    std::vector<int> trackedIdx;
    for (const auto& n : tracked) trackedIdx.push_back(sys.indexOf(n));
    const double total = years * JULIAN_YEAR_D;
    const int samples = 200;

    if (!ctx.quiet()) {
        std::printf("Integrator comparison: %zu bodies, %.1f years, fixed dt = %.4f d (%s)\n", sys.size(), years, dt,
                    str::formatDuration(dt * DAY_S).c_str());
        std::printf("Reference: Dormand-Prince 5(4), rtol = 1e-13\n\n");
    }

    // Reference solution.
    physics::SimulationConfig refCfg = base;
    refCfg.integrator = "dopri5";
    refCfg.tolerance = 1e-13;
    refCfg.dt = dt;
    physics::Simulation ref(sys, jd, refCfg);
    {
        Stopwatch sw;
        for (int k = 1; k <= samples; ++k) {
            ref.advance(total / samples);
            if (!ctx.quiet()) printProgress("reference      ", double(k) / samples);
        }
        if (!ctx.quiet()) finishProgress();
    }

    auto relPos = [&](const physics::NBodySystem& s, int idx) {
        const auto i = static_cast<std::size_t>(idx);
        return s.info(i).parent >= 0 ? s.pos(i) - s.pos(static_cast<std::size_t>(s.info(i).parent)) : s.pos(i);
    };

    std::vector<RunResult> results;
    for (const auto& name : names) {
        RunResult rr;
        physics::SimulationConfig cfg = base;
        cfg.integrator = name;
        cfg.dt = dt;
        try {
            physics::Simulation sim(sys, jd, cfg);
            rr.name = sim.integrator().name();
            rr.display = sim.integrator().displayName();
            rr.order = sim.integrator().order();
            rr.symplectic = sim.integrator().isSymplectic();
            Stopwatch sw;
            for (int k = 1; k <= samples; ++k) {
                sim.advance(total / samples);
                const auto e = sim.invariantErrors();
                const double ae = std::fabs(e.energyRel);
                if (!std::isfinite(ae)) throw std::runtime_error("integration diverged");
                rr.maxEnergy = std::max(rr.maxEnergy, ae);
                rr.energyHistory.emplace_back(sim.elapsed() / JULIAN_YEAR_D, std::max(ae, 1e-17));
                if (!ctx.quiet()) printProgress(str::padRight(rr.name, 15), double(k) / samples);
            }
            if (!ctx.quiet()) finishProgress();
            rr.wall = sw.seconds();
            rr.steps = sim.totalSteps();
            rr.evals = sim.forceEvaluations();
            const auto e = sim.invariantErrors();
            rr.finalEnergy = e.energyRel;
            rr.angMom = e.angularMomentumRel;
            for (int idx : trackedIdx)
                rr.posErrKm.push_back(idx < 0 ? std::numeric_limits<double>::quiet_NaN() : (relPos(sim.system(), idx) - relPos(ref.system(), idx)).norm() * AU_KM);
        } catch (const std::exception& ex) {
            if (!ctx.quiet()) finishProgress();
            rr.failed = true;
            rr.error = ex.what();
            if (rr.name.empty()) rr.name = name;
        }
        results.push_back(rr);
    }

    std::vector<std::string> hdr = {"Integrator", "Order", "Sympl.", "Steps", "Force evals", "Time", "max|dE/E|", "dL/L"};
    for (const auto& n : tracked) hdr.push_back(n + " err");
    TextTable t(hdr);
    auto fmtKm = [](double km) {
        if (!std::isfinite(km)) return std::string("-");
        if (km < 1.0) return str::format("%.1f m", km * 1e3);
        if (km < 1e6) return str::format("%.1f km", km);
        return str::format("%.2e km", km);
    };
    for (const auto& r : results) {
        if (r.failed) {
            t.addRow({r.name, "-", "-", "-", "-", "-", "FAILED: " + r.error});
            continue;
        }
        std::vector<std::string> row = {r.name, std::to_string(r.order), r.symplectic ? "yes" : "no",
                                        std::to_string(r.steps), std::to_string(r.evals), str::formatDuration(r.wall),
                                        str::format("%.2e", r.maxEnergy), str::format("%.2e", r.angMom)};
        for (double e : r.posErrKm) row.push_back(fmtKm(e));
        t.addRow(row);
    }
    t.print(std::cout, TextTable::defaultStyle());
    if (!ctx.quiet())
        std::printf("Position errors after %.1f years vs the reference (Moon: geocentric, others: heliocentric).\n"
                    "Symplectic methods keep |dE/E| bounded; non-symplectic ones drift secularly.\n", years);

    // Optional SVG chart of the energy error history (log scale).
    if (ctx.args().has("plot")) {
        const double W = 900, H = 520, L = 80, R = 200, T = 40, B = 60;
        io::SvgWriter svg(W, H);
        svg.text(W / 2, 24, str::format("Relative energy error |dE/E| over %.0f years (dt = %.3f d)", years, dt), "#e8eef6", 16, "middle", "bold");
        const double lo = -17, hi = 0;
        auto X = [&](double y) { return L + (W - L - R) * (y / years); };
        auto Y = [&](double e) { return T + (H - T - B) * (1.0 - (std::log10(e) - lo) / (hi - lo)); };
        for (int p = static_cast<int>(lo); p <= static_cast<int>(hi); p += 2) {
            svg.line(L, Y(std::pow(10.0, p)), W - R, Y(std::pow(10.0, p)), "#2a3550", 1);
            svg.text(L - 8, Y(std::pow(10.0, p)) + 4, str::format("1e%d", p), "#8fa0bc", 11, "end");
        }
        for (int k = 0; k <= 5; ++k) {
            const double yy = years * k / 5.0;
            svg.text(X(yy), H - B + 18, str::format("%.0f", yy), "#8fa0bc", 11, "middle");
        }
        svg.text((L + W - R) / 2, H - 18, "years", "#8fa0bc", 12, "middle");
        const char* colors[] = {"#ff6b6b", "#ffa94d", "#ffd43b", "#69db7c", "#4dabf7", "#b197fc", "#f783ac"};
        int ci = 0;
        for (const auto& r : results) {
            if (r.failed) continue;
            std::vector<std::pair<double, double>> pts;
            for (const auto& h : r.energyHistory) pts.emplace_back(X(h.first), Y(std::max(1e-17, std::min(1.0, h.second))));
            const char* c = colors[ci++ % 7];
            svg.polyline(pts, c, 2);
            svg.text(W - R + 12, T + 20.0 * ci, r.name, c, 13);
        }
        std::string err;
        if (!svg.save(ctx.args().getString("plot"), err)) throw std::runtime_error(err);
        std::printf("wrote %s\n", ctx.args().getString("plot").c_str());
    }
    if (ctx.args().has("csv")) {
        io::CsvWriter csv(ctx.args().getString("csv"));
        csv.header({"integrator", "years", "abs_rel_energy_error"});
        for (const auto& r : results)
            for (const auto& h : r.energyHistory) csv.row({r.name, io::CsvWriter::num(h.first), io::CsvWriter::num(h.second)});
        std::printf("wrote %s\n", ctx.args().getString("csv").c_str());
    }
    return 0;
}

} // namespace solsim::app
