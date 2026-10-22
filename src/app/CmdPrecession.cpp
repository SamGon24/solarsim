// =============================================================================
//  SolarSim - CmdPrecession.cpp
//  `solarsim precession` - reproduces the classic test of general relativity:
//  the anomalous 43"/century advance of Mercury's perihelion.
// =============================================================================
#include "AppContext.hpp"

#include "../analysis/Precession.hpp"
#include "../astro/OrbitalElements.hpp"
#include "../core/Constants.hpp"
#include "../core/StringUtil.hpp"
#include "../core/Table.hpp"

#include <cstdio>
#include <iostream>

namespace solsim::app {

using namespace solsim::constants;

int cmdPrecession(AppContext& ctx) {
    const std::string body = ctx.args().getString("body", "Mercury");
    const double years = ctx.args().getDouble("years", 100.0);
    if (!(years >= 5.0 && years <= 10000.0)) throw UsageError("--years must be between 5 and 10000");
    const double epochJD = ctx.args().has("epoch") ? ctx.epoch() : J2000_JD;

    physics::BuildOptions full = ctx.buildOptions(epochJD);
    physics::BuildOptions two = full;
    two.only = {"Sun", body};
    two.moons = two.minor = false;
    two.populations = {};
    two.calibration = physics::CalibrationMode::None;
    const physics::NBodySystem sysFull = physics::buildSystem(full);
    const physics::NBodySystem sysTwo = physics::buildSystem(two);
    if (sysFull.indexOf(body) < 0) throw UsageError("unknown body '" + body + "'");

    analysis::PrecessionConfig cfg;
    cfg.body = body;
    cfg.years = years;
    cfg.epochJD = epochJD;
    cfg.integrator = ctx.args().getString("integrator", "yoshida4");
    // Mercury's eccentric 88-day orbit needs a small step to keep the
    // integrator's own (numerical) precession far below 0.1"/cy.
    cfg.dt = ctx.args().getDouble("dt", 0.04);
    cfg.sampleDays = 2.0;

    const std::size_t bi = static_cast<std::size_t>(sysFull.indexOf(body));
    const auto el = astro::stateToElements(sysFull.pos(bi) - sysFull.pos(0), sysFull.vel(bi) - sysFull.vel(0),
                                           sysFull.gm()[0] + sysFull.gm()[bi], epochJD);
    const double theory = analysis::grPrecessionArcsecPerCentury(sysFull.gm()[0], el.a(), el.e);

    if (!ctx.quiet()) {
        std::printf("Perihelion precession of %s over %.0f years (%s, dt = %.3f d)\n", body.c_str(), years,
                    cfg.integrator.c_str(), cfg.dt);
        std::printf("Measured: slope of the osculating longitude of perihelion (Omega + omega), J2000 inertial frame.\n\n");
    }

    struct Run {
        const char* label;
        const physics::NBodySystem* sys;
        bool gr;
    };
    const Run runs[] = {{"Sun + body, Newtonian (control, expect ~0)", &sysTwo, false},
                        {"Sun + body, + general relativity", &sysTwo, true},
                        {"full Solar System, Newtonian", &sysFull, false},
                        {"full Solar System, + general relativity", &sysFull, true}};
    double rates[4] = {0, 0, 0, 0};
    TextTable t({"Model", "dvarpi/dt (\"/century)", "std. error", "wall time", "|dE/E|"});
    for (int k = 0; k < 4; ++k) {
        analysis::PrecessionConfig c = cfg;
        c.relativity = runs[k].gr;
        const auto r = analysis::measurePrecession(*runs[k].sys, c, [&](double f) {
            if (!ctx.quiet()) printProgress(str::format("run %d/4", k + 1), f);
        });
        if (!ctx.quiet()) {
            printProgress(str::format("run %d/4", k + 1), 1.0);
            finishProgress();
        }
        rates[k] = r.rateArcsecPerCentury;
        t.addRow({runs[k].label, str::format("%.3f", r.rateArcsecPerCentury), str::format("%.3f", r.stdErr),
                  str::formatDuration(r.wallSeconds), str::format("%.1e", std::fabs(r.finalEnergyError))});
    }
    t.print(std::cout, TextTable::defaultStyle());

    std::printf("\nRelativistic contribution (two-body):     %8.3f \"/cy\n", rates[1] - rates[0]);
    std::printf("Relativistic contribution (full system):  %8.3f \"/cy\n", rates[3] - rates[2]);
    std::printf("Einstein's prediction 6 pi GM/(c^2 a(1-e^2)): %5.3f \"/cy\n", theory);
    std::printf("Planetary (Newtonian) perturbations:      %8.3f \"/cy\n", rates[2] - rates[0]);
    if (str::iequals(body, "Mercury"))
        std::printf("\nReference values for Mercury (inertial frame): ~532\"/cy from the planets + 42.98\"/cy from GR\n"
                    "= ~575\"/cy. Le Verrier (1859) found the unexplained 43\"/cy that Einstein's 1915 theory resolved.\n");
    return 0;
}

} // namespace solsim::app
