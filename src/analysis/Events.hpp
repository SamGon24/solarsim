// =============================================================================
//  SolarSim - Events.hpp
//  Finds astronomical events by integrating the N-body system and locating
//  roots/extrema of geometric functions, refined with Brent's method and
//  golden-section search on short re-integrations from saved snapshots:
//
//    * equinoxes & solstices       apparent solar longitude = 0/90/180/270 deg
//    * oppositions & conjunctions  geocentric longitude difference body-Sun
//    * greatest elongations        maxima of the Sun-Earth-planet angle
//    * lunar phases                Moon-Sun longitude difference = 0/90/180/270
//    * Earth perihelion/aphelion   extrema of the Earth-Sun distance
//    * lunar perigee/apogee        extrema of the Earth-Moon distance
// =============================================================================
#pragma once

#include "../physics/Simulation.hpp"

#include <functional>
#include <string>
#include <vector>

namespace solsim::analysis {

struct EventRecord {
    double jd = 0.0;       // TT
    std::string category;  // "season", "opposition", "conjunction", ...
    std::string body;
    std::string description;
    double value = 0.0;    // context dependent (elongation deg, distance AU ...)
};

struct EventSearchOptions {
    double startJD = 0.0;
    double endJD = 0.0;
    double sampleDays = 0.5;
    bool seasons = true;
    bool planets = true;      // oppositions, conjunctions, elongations
    bool moonPhases = true;
    bool apsides = true;      // Earth perihelion/aphelion, lunar perigee/apogee
    std::vector<std::string> bodies; // restrict planetary events (empty = all)
};

// `sim` must be positioned at (or before) startJD; it is copied, not modified.
std::vector<EventRecord> findEvents(const physics::Simulation& sim, const EventSearchOptions& opt,
                                    const std::function<void(double)>& progress = {});

} // namespace solsim::analysis
