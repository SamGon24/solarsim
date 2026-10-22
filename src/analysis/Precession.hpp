// =============================================================================
//  SolarSim - Precession.hpp
//  Measures the secular advance of a planet's perihelion from the N-body
//  integration: samples the osculating heliocentric longitude of perihelion
//  varpi = Omega + omega, unwraps it and fits a straight line (least squares).
//
//  Classic result for Mercury (inertial J2000 frame, arcsec per century):
//      planetary perturbations  ~ 532
//      general relativity       ~  43   (Einstein 1915)
//      total                    ~ 575
//  GR prediction per orbit: dvarpi = 6 pi G M / (c^2 a (1 - e^2))
// =============================================================================
#pragma once

#include "../physics/NBodySystem.hpp"

#include <functional>
#include <string>

namespace solsim::analysis {

struct PrecessionConfig {
    std::string body = "Mercury";
    double years = 100.0;
    double sampleDays = 2.0;
    bool relativity = false;
    std::string integrator = "yoshida4";
    double dt = 0.0;          // 0 = automatic
    double epochJD = 2451545.0;
};

struct PrecessionResult {
    std::string label;
    double rateArcsecPerCentury = 0.0;
    double stdErr = 0.0;
    double r2 = 0.0;
    std::size_t samples = 0;
    double wallSeconds = 0.0;
    double finalEnergyError = 0.0;
};

PrecessionResult measurePrecession(const physics::NBodySystem& system, const PrecessionConfig& cfg,
                                   const std::function<void(double)>& progress = {});

// Analytic GR perihelion advance for an orbit (arcsec / century).
double grPrecessionArcsecPerCentury(double gmSun, double a, double e);

} // namespace solsim::analysis
