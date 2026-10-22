// =============================================================================
//  SolarSim - SystemBuilder.hpp
//  Turns body definitions (built-in catalog + optional INI files) into an
//  N-body system with barycentric initial conditions at a chosen epoch.
//
//  Steps:
//   1. select bodies (default set, --moons, --minor, --only, --exclude, files)
//   2. fold the GM of non-simulated catalog moons into their parent planet
//   3. evaluate each orbit model relative to the parent (ecliptic J2000)
//   4. treat a planet's mean orbit as the orbit of its system barycentre and
//      displace planet & moons around it (e.g. Earth vs Earth-Moon barycentre)
//   5. shift everything to the Solar System barycentre (total momentum = 0)
//   6. add synthetic populations (asteroid belt, Trojans, Kuiper belt)
// =============================================================================
#pragma once

#include "Calibration.hpp"
#include "NBodySystem.hpp"

#include "../astro/BodyCatalog.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace solsim::physics {

struct PopulationOptions {
    int asteroids = 0;  // main belt test particles
    int trojans = 0;    // Jupiter L4 + L5 Trojans
    int kuiper = 0;     // classical KBOs + plutinos
    std::uint64_t seed = 20261021;
};

struct IntruderOptions {
    bool enabled = false;
    double massSolar = 0.5;           // solar masses
    double perihelionAU = 80.0;       // closest approach to the Sun
    double vInfKms = 3.0;             // hyperbolic excess speed
    double daysToPerihelion = 36525;  // from the epoch
    double inclinationDeg = 25.0;
};

struct BuildOptions {
    double epochJD = 2451545.0;
    bool moons = false;      // Galilean moons, Titan, Triton
    bool minor = false;      // Ceres, Vesta, Halley, Encke
    bool includeMoon = true; // Earth's Moon
    std::vector<std::string> only;     // restrict to these bodies (+ parents)
    std::vector<std::string> exclude;  // remove these bodies
    std::vector<astro::BodyDefinition> extra; // from --bodies files
    PopulationOptions populations;
    IntruderOptions intruder;
    CalibrationMode calibration = CalibrationMode::All;  // differential correction of initial velocities
};

struct BuildReport {
    std::vector<std::string> warnings;
    std::vector<std::string> notes;
    CalibrationReport calibration;
};

NBodySystem buildSystem(const BuildOptions& options, BuildReport* report = nullptr);

// Converts catalog data to runtime metadata.
BodyInfo makeBodyInfo(const astro::BodyDefinition& def);

// Population generators (append test particles; system must contain the Sun
// and - for Trojans/plutinos - Jupiter/Neptune).
void addMainBelt(NBodySystem& sys, int count, std::uint64_t seed, double jd);
void addJupiterTrojans(NBodySystem& sys, int count, std::uint64_t seed, double jd);
void addKuiperBelt(NBodySystem& sys, int count, std::uint64_t seed, double jd);

} // namespace solsim::physics
