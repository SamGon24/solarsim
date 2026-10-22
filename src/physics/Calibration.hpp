// =============================================================================
//  SolarSim - Calibration.hpp
//  Differential correction of initial conditions (a simplified form of what
//  ephemeris makers call "orbit determination").
//
//  Initial state vectors derived from MEAN elements or truncated series are
//  not perfectly consistent with the N-body dynamics: a velocity error of a
//  few 1e-5 changes the osculating semi-major axis and therefore the mean
//  motion (dn/n = -3 dv/v), which makes the body drift linearly in longitude.
//
//  This module integrates the system over a calibration span, measures the
//  longitude drift of each body against its reference theory (ELP-2000/82 for
//  the Moon, Standish mean elements for the planets) with a least-squares
//  line fit, scales the initial velocity by (1 + dn / 3n) and repeats
//  (a body that runs ahead needs MORE orbital energy: larger a, slower n).
// =============================================================================
#pragma once

#include "NBodySystem.hpp"

#include <functional>
#include <string>
#include <vector>

namespace solsim::physics {

enum class CalibrationMode { None, Moon, All };
bool parseCalibrationMode(const std::string& s, CalibrationMode& out);
const char* calibrationModeName(CalibrationMode m);

struct CalibrationReport {
    struct Entry {
        std::string body;
        double driftBeforeArcsecPerYear = 0.0;
        double driftAfterArcsecPerYear = 0.0;
        double velocityScale = 1.0; // cumulative factor applied
    };
    std::vector<Entry> entries;
    double wallSeconds = 0.0;
};

// Modifies `sys` in place. `epochJD` is the time of the state vectors.
CalibrationReport calibrateInitialConditions(NBodySystem& sys, double epochJD, CalibrationMode mode,
                                             int iterations = 2,
                                             const std::function<void(const std::string&)>& log = {});

} // namespace solsim::physics
