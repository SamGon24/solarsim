// =============================================================================
//  SolarSim - Transfer.hpp
//  Interplanetary mission design, from textbook to "flight dynamics":
//    1. analytic Hohmann transfer (circular, coplanar): delta-v, time of
//       flight, phase angle, synodic period
//    2. launch-window search with the real (N-body) planetary positions
//    3. PORKCHOP SCAN: Lambert arcs (Sun-only two-body) over a grid of launch
//       dates x flight times -> departure C3, arrival v_inf, total delta-v
//    4. N-body SINGLE SHOOTING: Newton iteration on the launch velocity with
//       a finite-difference Jacobian of the full N-body propagation (patched
//       conic: planetary gravity on the probe switched off)
//    5. B-PLANE TARGETING: secant iteration on the aim offset using flights
//       WITH the target's gravity until the requested periapsis is reached
//       (the role played by trajectory-correction manoeuvres in real missions)
//    6. departure/capture burns from parking orbits (Oberth effect), DLA
// =============================================================================
#pragma once

#include "../physics/Simulation.hpp"

#include <functional>
#include <string>
#include <vector>

namespace solsim::analysis {

struct HohmannTransfer {
    double r1 = 0, r2 = 0;            // AU
    double aTransfer = 0;             // AU
    double dv1 = 0, dv2 = 0;          // km/s (heliocentric burns)
    double timeOfFlight = 0;          // days
    double phaseAngleDeg = 0;         // target lead angle at launch
    double synodicPeriod = 0;         // days
};

HohmannTransfer hohmann(double mu, double r1, double r2);

struct MissionConfig {
    std::string from = "Earth";
    std::string to = "Mars";
    double searchStartJD = 0.0;
    double parkingAltitudeKm = 300.0;   // departure circular orbit
    double captureAltitudeKm = 400.0;   // desired periapsis at the target
    int maxIterations = 15;
    double toleranceKm = 25.0;          // N-body shooting tolerance
    double periapsisToleranceKm = 5.0;
    bool porkchop = true;               // optimise launch date & TOF with Lambert arcs
    double porkchopHalfWidthDays = 75.0;
    double porkchopStepDays = 2.0;
};

struct PorkchopGrid {
    std::vector<double> launchJD;       // columns
    std::vector<double> tofDays;        // rows
    std::vector<double> c3;             // row-major [tof][launch], km^2/s^2 (NaN if no solution)
    std::vector<double> totalDv;        // km/s (departure from parking orbit + capture)
    double bestLaunch = 0.0, bestTof = 0.0, bestC3 = 0.0, bestDv = 0.0;
    bool valid() const { return !launchJD.empty() && !tofDays.empty(); }
};

struct MissionResult {
    bool success = false;
    std::string message;
    HohmannTransfer analytic;
    PorkchopGrid porkchop;
    double windowJD = 0.0;              // Hohmann phase-angle window
    double launchJD = 0.0, arrivalJD = 0.0;
    double timeOfFlight = 0.0;          // days
    double transferAngleDeg = 0.0;
    double vInfDepart = 0.0;            // km/s
    double c3 = 0.0;                    // km^2/s^2
    double declinationDeg = 0.0;        // DLA of the departure asymptote
    double dvDepart = 0.0;              // km/s from the parking orbit
    double vInfArrive = 0.0;            // km/s
    double aimOffsetKm = 0.0;           // final B-plane offset
    double periapsisAltKm = 0.0;        // achieved at the target
    double closestApproachJD = 0.0;
    double dvCapture = 0.0;             // km/s into a circular orbit at periapsis
    bool impact = false;
    int iterations = 0;
    double finalMissKm = 0.0;
    std::vector<Vec3> trajectory;       // heliocentric samples (AU)
    std::vector<Vec3> departOrbit;      // heliocentric orbit samples of the departure body
    std::vector<Vec3> targetOrbit;      // heliocentric orbit samples of the target body
    Vec3 launchPosition, arrivalPosition;
};

using MissionLog = std::function<void(const std::string&)>;

MissionResult planMission(const physics::Simulation& baseSim, const MissionConfig& cfg, const MissionLog& log = {});

} // namespace solsim::analysis
