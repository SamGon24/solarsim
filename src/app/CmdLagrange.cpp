// =============================================================================
//  SolarSim - CmdLagrange.cpp
//  `solarsim lagrange` - Lagrange points of several two-body pairs and an
//  optional stability experiment with test particles.
// =============================================================================
#include "AppContext.hpp"

#include "../analysis/Lagrange.hpp"
#include "../astro/JulianDate.hpp"
#include "../core/Constants.hpp"
#include "../core/StringUtil.hpp"
#include "../core/Table.hpp"

#include <cstdio>
#include <iostream>
#include <utility>

namespace solsim::app {

using namespace solsim::constants;

namespace {
std::string fmtDist(double au) {
    const double km = au * AU_KM;
    if (km < 1e7) return str::format("%.0f km", km);
    return str::format("%.4f AU", au);
}
} // namespace

int cmdLagrange(AppContext& ctx) {
    const double jd = ctx.epoch();
    auto sim = ctx.makeSimulation(jd);
    const auto& sys = sim->system();

    std::vector<std::pair<std::string, std::string>> pairs;
    if (ctx.args().has("pair")) {
        const auto l = ctx.args().getList("pair");
        if (l.size() != 2) throw UsageError("--pair expects PRIMARY,SECONDARY (e.g. Sun,Jupiter)");
        pairs.emplace_back(l[0], l[1]);
    } else {
        pairs = {{"Sun", "Earth"}, {"Earth", "Moon"}, {"Sun", "Mars"}, {"Sun", "Jupiter"}, {"Sun", "Neptune"}};
    }

    for (const auto& pr : pairs) {
        const int a = sys.indexOf(pr.first), b = sys.indexOf(pr.second);
        if (a < 0 || b < 0) {
            std::printf("skipping %s-%s (not simulated)\n", pr.first.c_str(), pr.second.c_str());
            continue;
        }
        const analysis::LagrangeSolution s = analysis::computeLagrangePoints(sys, a, b);
        std::printf("\n%s-%s  (mu = m2/(m1+m2) = %.6e, separation %s)\n", s.primary.c_str(), s.secondary.c_str(),
                    s.massRatio, fmtDist(s.separation).c_str());
        TextTable t({"Point", "x (sep.)", "y (sep.)", "from " + s.secondary, "Notes"});
        t.setAlign(4, TextTable::Align::Left);
        const char* notes[5] = {"between the bodies - unstable (saddle)", "beyond the secondary - unstable",
                                "opposite side of the primary - unstable",
                                s.triangularStable ? "leads by 60 deg - stable (Trojans)" : "leads by 60 deg - UNSTABLE",
                                s.triangularStable ? "trails by 60 deg - stable (Trojans)" : "trails by 60 deg - UNSTABLE"};
        for (int k = 0; k < 5; ++k)
            t.addRow({"L" + std::to_string(k + 1), str::format("%+.8f", s.x[static_cast<std::size_t>(k)]),
                      str::format("%+.6f", s.y[static_cast<std::size_t>(k)]),
                      fmtDist(s.distanceFromSecondary[static_cast<std::size_t>(k)]), notes[k]});
        t.print(std::cout, TextTable::defaultStyle());
        std::printf("Hill sphere %s | Laplace sphere of influence %s | Jacobi constant C(L1) = %.6f\n",
                    fmtDist(s.hillRadius).c_str(), fmtDist(s.sphereOfInfluence).c_str(), s.jacobiL1);
        std::printf("Routh criterion: L4/L5 stable if mu < 0.0385209 -> %s\n", s.triangularStable ? "yes" : "no");

        const double years = ctx.args().getDouble("test-years", 0.0);
        if (years > 0.0) {
            std::printf("Stability experiment: massless probes released at L1..L5, integrated %.0f years...\n", years);
            const auto res = analysis::testLagrangeStability(*sim, a, b, years, 0.0);
            TextTable r({"Point", "max drift", "final drift", "max angle", "Result"});
            for (const auto& x : res)
                r.addRow({"L" + std::to_string(x.point), fmtDist(x.maxDeviationAU), fmtDist(x.finalDeviationAU),
                          str::format("%.2f deg", x.maxAngleDeg),
                          x.escaped ? "escaped (unstable)" : (x.maxAngleDeg < 30 ? "librates (stable)" : "wanders")});
            r.print(std::cout, TextTable::defaultStyle());
        }
    }
    if (!ctx.quiet() && !ctx.args().has("test-years"))
        std::printf("\nTip: add --pair Sun,Jupiter --test-years 200 to watch Trojan-like libration at L4/L5.\n");
    return 0;
}

} // namespace solsim::app
