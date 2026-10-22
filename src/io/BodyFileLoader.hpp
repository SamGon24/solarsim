// =============================================================================
//  SolarSim - BodyFileLoader.hpp
//  Loads extra bodies from INI files (see data/*.ini for documented examples).
//
//  [body]
//  name      = Hale-Bopp          ; required
//  parent    = Sun                ; default Sun
//  kind      = comet              ; planet|dwarf|moon|asteroid|comet|spacecraft
//  ; mass: one of gm_km3s2, mass_kg, mass_earths, mass_suns (default: massless)
//  radius_km = 30
//  ; orbit (heliocentric/planetocentric ecliptic J2000):
//  q | a     = 0.914              ; AU   (a < 0 for hyperbolic orbits)
//  e         = 0.995086
//  i, node, argp (or peri = longitude of perihelion)       ; degrees
//  ; time reference - one of:
//  tp        = 1997-04-01T03:00   ; time of perihelion (date or JD)
//  M  + epoch                     ; mean anomaly at epoch (deg)
//  L  + epoch                     ; mean longitude at epoch (deg)
//  ; appearance (optional): color = #RRGGBB, color2, style, rings_inner_km ...
// =============================================================================
#pragma once

#include "../astro/BodyCatalog.hpp"

#include <string>
#include <vector>

namespace solsim::io {

struct LoadResult {
    std::vector<astro::BodyDefinition> bodies;
    std::vector<std::string> warnings;
    std::string error; // non-empty on failure
    bool ok() const { return error.empty(); }
};

LoadResult loadBodyFile(const std::string& path);
LoadResult loadBodyString(const std::string& text, const std::string& sourceName = "<string>");

} // namespace solsim::io
