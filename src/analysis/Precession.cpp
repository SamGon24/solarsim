// =============================================================================
//  SolarSim - Precession.cpp
// =============================================================================
#include "Precession.hpp"

#include "../astro/OrbitalElements.hpp"
#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"
#include "../core/Timer.hpp"
#include "../physics/Simulation.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace solsim::analysis {

using namespace solsim::constants;

double grPrecessionArcsecPerCentury(double gmSun, double a, double e) {
    const double c2 = C_AU_PER_DAY * C_AU_PER_DAY;
    const double perOrbit = 6.0 * PI * gmSun / (c2 * a * (1.0 - e * e));
    const double period = TWO_PI * std::sqrt(a * a * a / gmSun);
    return perOrbit * (JULIAN_CENTURY_D / period) * RAD2ARCSEC;
}

PrecessionResult measurePrecession(const physics::NBodySystem& system, const PrecessionConfig& cfg,
                                   const std::function<void(double)>& progress) {
    const int bi = system.indexOf(cfg.body);
    if (bi <= 0) throw std::invalid_argument("unknown body '" + cfg.body + "'");
    const std::size_t b = static_cast<std::size_t>(bi);

    physics::SimulationConfig sc;
    sc.integrator = cfg.integrator;
    sc.dt = cfg.dt;
    sc.relativity = cfg.relativity;
    sc.threads = 1;
    physics::Simulation sim(system, cfg.epochJD, sc);

    std::vector<double> t, varpi;
    const double total = cfg.years * JULIAN_YEAR_D;
    const std::size_t nSamples = static_cast<std::size_t>(total / cfg.sampleDays) + 1;
    t.reserve(nSamples);
    varpi.reserve(nSamples);
    Stopwatch sw;
    for (std::size_t k = 0; k < nSamples; ++k) {
        if (k > 0) sim.advance(cfg.sampleDays);
        const auto& s = sim.system();
        // Heliocentric osculating elements (two-body mu = GM_sun + GM_body).
        const auto el = astro::stateToElements(s.pos(b) - s.pos(0), s.vel(b) - s.vel(0), s.gm()[0] + s.gm()[b], sim.time());
        t.push_back((sim.time() - cfg.epochJD) / JULIAN_CENTURY_D);
        varpi.push_back(el.node + el.argp);
        if (progress && (k % 512) == 0) progress(double(k) / double(nSamples));
    }
    math::unwrapAngles(varpi);
    for (double& v : varpi) v *= RAD2ARCSEC;
    const math::LinearFit fit = math::fitLine(t, varpi);

    PrecessionResult r;
    r.rateArcsecPerCentury = fit.slope;
    r.stdErr = fit.slopeStdErr;
    r.r2 = fit.r2;
    r.samples = fit.n;
    r.wallSeconds = sw.seconds();
    r.finalEnergyError = sim.invariantErrors().energyRel;
    return r;
}

} // namespace solsim::analysis
