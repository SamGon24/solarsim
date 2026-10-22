// =============================================================================
//  SolarSim - BodyCatalog.hpp
//  Built-in database of Solar System bodies: physical constants (GM, radius,
//  rotation, pole orientation, albedo), rendering hints and the orbital model
//  used to generate initial conditions.
//
//  Sources (values reproduced in BodyCatalog.cpp):
//    * GM values: JPL DE440 planetary & lunar ephemeris (Park et al. 2021)
//    * Radii / rotation / poles: IAU WGCCRE 2015 report (Archinal et al. 2018)
//    * Planetary mean elements: E.M. Standish, "Keplerian Elements for
//      Approximate Positions of the Major Planets" (JPL SSD), Table 1,
//      valid 1800 AD - 2050 AD
//    * Lunar mean elements: Meeus, "Astronomical Algorithms" (2nd ed.) ch. 47
//    * Minor bodies & moons: approximate osculating elements (flagged)
// =============================================================================
#pragma once

#include "OrbitalElements.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace solsim::astro {

enum class BodyKind : std::uint8_t { Star, Planet, DwarfPlanet, Moon, Asteroid, Comet, Spacecraft, Particle };
const char* kindName(BodyKind k);
bool parseKind(const std::string& s, BodyKind& out);

enum class SurfaceStyle : std::uint8_t {
    Star, Rocky, Cratered, Cloudy, Terran, Desert, GasBanded, IceGiant, Icy, Volcanic, Comet, Metallic
};
const char* surfaceStyleName(SurfaceStyle s);
bool parseSurfaceStyle(const std::string& s, SurfaceStyle& out);

struct RingSystem {
    bool present = false;
    double innerKm = 0.0;
    double outerKm = 0.0;
    double opacity = 1.0;  // overall visibility scale
};

// Mean elements with secular rates (per Julian century from J2000).
// Angles in degrees, a in AU.
struct MeanElements {
    double a = 0, e = 0, I = 0, L = 0, varpi = 0, node = 0;
    double da = 0, de = 0, dI = 0, dL = 0, dvarpi = 0, dnode = 0;
};

enum class OrbitModel : std::uint8_t {
    None,           // the central star
    MeanElements,   // Standish-style mean elements with rates (ecliptic J2000, relative to parent)
    Osculating,     // fixed conic: q, e, i, node, argp, tp  (ecliptic J2000, relative to parent)
    Equatorial,     // satellite in parent's equatorial plane: a, e, incl, lambda0 (approximate phase)
    LunarTheory     // Earth's Moon from the truncated ELP-2000/82 series (LunarTheory.hpp)
};

struct BodyDefinition {
    std::string name;
    std::string parent;          // empty for the Sun
    BodyKind kind = BodyKind::Planet;

    // ---- physical ----
    double gmKm3s2 = 0.0;        // gravitational parameter of the body alone
    double radiusKm = 0.0;       // mean radius
    double rotationHours = 0.0;  // sidereal rotation period (negative = retrograde)
    double poleRaDeg = 0.0;      // north pole, ICRF
    double poleDecDeg = 90.0;
    double albedo = 0.3;         // geometric albedo

    // ---- appearance ----
    SurfaceStyle style = SurfaceStyle::Rocky;
    std::uint32_t color = 0xAAAAAA;   // primary sRGB colour 0xRRGGBB
    std::uint32_t color2 = 0x666666;  // secondary (bands, continents...)
    RingSystem rings;

    // ---- orbit ----
    OrbitModel model = OrbitModel::None;
    MeanElements mean;                 // for MeanElements
    double q = 0, e = 0, i = 0, node = 0, argp = 0, tp = 0;   // Osculating (deg, JD)
    double aKm = 0, incl = 0, lambda0 = 0;                    // Equatorial (km, deg, deg at J2000)
    bool meanIsBarycenter = false;     // Standish 'EM Bary' row -> Earth-Moon barycentre
    bool approximate = false;          // flagged in `info` output
    std::string notes;

    // Typical orbital distance from parent (AU) - used for display scaling and dt selection.
    double typicalDistanceAU() const;
    double gm() const;                 // AU^3/day^2
    double radiusAU() const;
};

class BodyCatalog {
public:
    // The complete built-in database (stable order; the Sun is first).
    static const std::vector<BodyDefinition>& all();
    static const BodyDefinition* find(const std::string& name); // case-insensitive

    // Names of the default set (Sun, 8 planets, Moon, Pluto).
    static std::vector<std::string> defaultSet();
    static std::vector<std::string> majorMoonSet();   // Io, Europa, Ganymede, Callisto, Titan, Triton
    static std::vector<std::string> minorBodySet();   // Ceres, Vesta, Halley, Encke
};

// Evaluates Standish-style mean elements at a date -> osculating elements.
OrbitalElements meanElementsAt(const MeanElements& m, double jdTT, double mu);

// State of `def` relative to its parent at jdTT (ecliptic J2000 frame).
//   muCentral : GM(parent) + GM(body) - used to derive velocities
//   parentDef : required for OrbitModel::Equatorial (pole orientation)
StateVector relativeStateAt(const BodyDefinition& def, double jdTT, double muCentral,
                            const BodyDefinition* parentDef);

} // namespace solsim::astro
