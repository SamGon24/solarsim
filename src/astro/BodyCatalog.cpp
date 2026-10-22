// =============================================================================
//  SolarSim - BodyCatalog.cpp
//  Embedded data. See BodyCatalog.hpp for references. Every value marked
//  "approximate" is good enough for visualisation and qualitative dynamics but
//  should not be used for precise ephemerides.
// =============================================================================
#include "BodyCatalog.hpp"

#include "Coordinates.hpp"
#include "JulianDate.hpp"
#include "LunarTheory.hpp"

#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"
#include "../core/StringUtil.hpp"

#include <cmath>

namespace solsim::astro {

using namespace solsim::constants;

const char* kindName(BodyKind k) {
    switch (k) {
        case BodyKind::Star: return "star";
        case BodyKind::Planet: return "planet";
        case BodyKind::DwarfPlanet: return "dwarf planet";
        case BodyKind::Moon: return "moon";
        case BodyKind::Asteroid: return "asteroid";
        case BodyKind::Comet: return "comet";
        case BodyKind::Spacecraft: return "spacecraft";
        case BodyKind::Particle: return "particle";
    }
    return "?";
}

bool parseKind(const std::string& s, BodyKind& out) {
    const std::string k = str::toLower(str::trim(s));
    if (k == "star") out = BodyKind::Star;
    else if (k == "planet") out = BodyKind::Planet;
    else if (k == "dwarf" || k == "dwarf planet" || k == "dwarfplanet") out = BodyKind::DwarfPlanet;
    else if (k == "moon" || k == "satellite") out = BodyKind::Moon;
    else if (k == "asteroid" || k == "minor planet") out = BodyKind::Asteroid;
    else if (k == "comet") out = BodyKind::Comet;
    else if (k == "spacecraft" || k == "probe") out = BodyKind::Spacecraft;
    else if (k == "particle") out = BodyKind::Particle;
    else return false;
    return true;
}

const char* surfaceStyleName(SurfaceStyle s) {
    switch (s) {
        case SurfaceStyle::Star: return "star";
        case SurfaceStyle::Rocky: return "rocky";
        case SurfaceStyle::Cratered: return "cratered";
        case SurfaceStyle::Cloudy: return "cloudy";
        case SurfaceStyle::Terran: return "terran";
        case SurfaceStyle::Desert: return "desert";
        case SurfaceStyle::GasBanded: return "banded";
        case SurfaceStyle::IceGiant: return "icegiant";
        case SurfaceStyle::Icy: return "icy";
        case SurfaceStyle::Volcanic: return "volcanic";
        case SurfaceStyle::Comet: return "comet";
        case SurfaceStyle::Metallic: return "metallic";
    }
    return "rocky";
}

bool parseSurfaceStyle(const std::string& s, SurfaceStyle& out) {
    const std::string k = str::toLower(str::trim(s));
    static const SurfaceStyle all[] = {SurfaceStyle::Star,      SurfaceStyle::Rocky,     SurfaceStyle::Cratered,
                                       SurfaceStyle::Cloudy,    SurfaceStyle::Terran,    SurfaceStyle::Desert,
                                       SurfaceStyle::GasBanded, SurfaceStyle::IceGiant,  SurfaceStyle::Icy,
                                       SurfaceStyle::Volcanic,  SurfaceStyle::Comet,     SurfaceStyle::Metallic};
    for (SurfaceStyle st : all) {
        if (k == surfaceStyleName(st)) {
            out = st;
            return true;
        }
    }
    return false;
}

double BodyDefinition::gm() const { return gmKm3s2 * KM3S2_TO_AU3D2; }
double BodyDefinition::radiusAU() const { return radiusKm / AU_KM; }

double BodyDefinition::typicalDistanceAU() const {
    switch (model) {
        case OrbitModel::None: return 0.0;
        case OrbitModel::MeanElements: return mean.a;
        case OrbitModel::Osculating: return e < 1.0 ? q / (1.0 - e) : q;
        case OrbitModel::Equatorial: return aKm / AU_KM;
        case OrbitModel::LunarTheory: return 384399.0 / AU_KM;
    }
    return 0.0;
}

namespace {

BodyDefinition makeSun() {
    BodyDefinition b;
    b.name = "Sun";
    b.kind = BodyKind::Star;
    b.gmKm3s2 = GM_SUN_KM3S2;
    b.radiusKm = SUN_RADIUS_KM;
    b.rotationHours = 25.38 * 24.0;   // sidereal, equatorial (differential rotation!)
    b.poleRaDeg = 286.13;
    b.poleDecDeg = 63.87;
    b.style = SurfaceStyle::Star;
    b.color = 0xFFF1D6;
    b.color2 = 0xFF9A2E;
    b.model = OrbitModel::None;
    b.notes = "G2V main-sequence star, 99.86% of the system mass";
    return b;
}

BodyDefinition makePlanet(const char* name, BodyKind kind, double gm, double radius, double rotH, double poleRa,
                          double poleDec, double albedo, SurfaceStyle style, std::uint32_t c1, std::uint32_t c2,
                          const MeanElements& m, const char* notes) {
    BodyDefinition b;
    b.name = name;
    b.parent = "Sun";
    b.kind = kind;
    b.gmKm3s2 = gm;
    b.radiusKm = radius;
    b.rotationHours = rotH;
    b.poleRaDeg = poleRa;
    b.poleDecDeg = poleDec;
    b.albedo = albedo;
    b.style = style;
    b.color = c1;
    b.color2 = c2;
    b.model = OrbitModel::MeanElements;
    b.mean = m;
    b.notes = notes;
    return b;
}

BodyDefinition makeSatellite(const char* name, const char* parent, double gm, double radius, double rotH,
                             double aKm, double e, double incl, double lambda0, SurfaceStyle style,
                             std::uint32_t c1, std::uint32_t c2, const char* notes) {
    BodyDefinition b;
    b.name = name;
    b.parent = parent;
    b.kind = BodyKind::Moon;
    b.gmKm3s2 = gm;
    b.radiusKm = radius;
    b.rotationHours = rotH;
    b.albedo = 0.3;
    b.style = style;
    b.color = c1;
    b.color2 = c2;
    b.model = OrbitModel::Equatorial;
    b.aKm = aKm;
    b.e = e;
    b.incl = incl;
    b.lambda0 = lambda0;
    b.approximate = true;
    b.notes = notes;
    return b;
}

BodyDefinition makeMinor(const char* name, BodyKind kind, double gm, double radius, double rotH, double q,
                         double e, double i, double node, double argp, double tp, SurfaceStyle style,
                         std::uint32_t c1, std::uint32_t c2, const char* notes) {
    BodyDefinition b;
    b.name = name;
    b.parent = "Sun";
    b.kind = kind;
    b.gmKm3s2 = gm;
    b.radiusKm = radius;
    b.rotationHours = rotH;
    b.albedo = kind == BodyKind::Comet ? 0.04 : 0.1;
    b.style = style;
    b.color = c1;
    b.color2 = c2;
    b.model = OrbitModel::Osculating;
    b.q = q;
    b.e = e;
    b.i = i;
    b.node = node;
    b.argp = argp;
    b.tp = tp;
    b.approximate = true;
    b.notes = notes;
    return b;
}

std::vector<BodyDefinition> buildCatalog() {
    std::vector<BodyDefinition> v;
    v.push_back(makeSun());

    // ---------------------------------------------------------------------------------
    // Standish (JPL) Table 1: a[AU] e I[deg] L[deg] varpi[deg] node[deg] + rates per century
    // ---------------------------------------------------------------------------------
    v.push_back(makePlanet("Mercury", BodyKind::Planet, 22031.868551, 2439.7, 1407.5976, 281.0103, 61.4155, 0.142,
                           SurfaceStyle::Cratered, 0x9C9389, 0x6F6760,
                           {0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593,
                            0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081},
                           "Smallest planet; 3:2 spin-orbit resonance; GR perihelion advance 43\"/cy"));
    v.push_back(makePlanet("Venus", BodyKind::Planet, 324858.592, 6051.8, -5832.6, 272.76, 67.16, 0.689,
                           SurfaceStyle::Cloudy, 0xE8D7A6, 0xC9B07A,
                           {0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255,
                            0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418},
                           "Retrograde rotation; 92 bar CO2 atmosphere"));
    {
        BodyDefinition earth = makePlanet(
            "Earth", BodyKind::Planet, 398600.435507, 6371.0, 23.9344696, 0.0, 90.0, 0.434, SurfaceStyle::Terran,
            0x2F6FD0, 0x3F8F3A,
            {1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0, 0.00000562, -0.00004392,
             -0.01294668, 35999.37244981, 0.32327364, 0.0},
            "Mean elements describe the Earth-Moon barycentre");
        earth.meanIsBarycenter = true;
        v.push_back(earth);
    }
    {
        // Moon: mean elements relative to Earth (Meeus ch. 47 fundamental arguments).
        BodyDefinition moon;
        moon.name = "Moon";
        moon.parent = "Earth";
        moon.kind = BodyKind::Moon;
        moon.gmKm3s2 = 4902.800118;
        moon.radiusKm = 1737.4;
        moon.rotationHours = 655.7199;
        moon.poleRaDeg = 269.9949;
        moon.poleDecDeg = 66.5392;
        moon.albedo = 0.12;
        moon.style = SurfaceStyle::Cratered;
        moon.color = 0xB5B1AA;
        moon.color2 = 0x6E6B66;
        // Initial state from the ELP-2000/82 lunar theory; mean elements kept for reference.
        moon.model = OrbitModel::LunarTheory;
        moon.mean = {384399.0 / AU_KM, 0.0549, 5.145, 218.3164477, 83.3532465, 125.0445479,
                     0.0, 0.0, 0.0, 481267.88123421, 4069.0137287, -1934.1362891};
        moon.notes = "ELP-2000/82 (Meeus ch. 47); node regresses in 18.6 yr, perigee advances in 8.85 yr";
        v.push_back(moon);
    }
    v.push_back(makePlanet("Mars", BodyKind::Planet, 42828.375816, 3389.5, 24.6229, 317.68143, 52.88650, 0.170,
                           SurfaceStyle::Desert, 0xC1440E, 0xE8A77A,
                           {1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891,
                            0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343},
                           "Includes Phobos & Deimos in GM"));
    {
        BodyDefinition jup = makePlanet(
            "Jupiter", BodyKind::Planet, 126712764.1 - (5959.916 + 3202.739 + 9887.834 + 7179.289), 69911.0,
            9.9250, 268.056595, 64.495303, 0.538, SurfaceStyle::GasBanded, 0xD8C8A0, 0xA0704A,
            {5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909, -0.00011607, -0.00013253,
             -0.00183714, 3034.74612775, 0.21252668, 0.20469106},
            "GM excludes the Galilean moons (added back when they are not simulated)");
        jup.rings = {true, 92000.0, 226000.0, 0.015}; // faint dusty rings
        v.push_back(jup);
    }
    {
        BodyDefinition sat = makePlanet(
            "Saturn", BodyKind::Planet, 37940584.8418 - 8978.137, 58232.0, 10.656, 40.589, 83.537, 0.499,
            SurfaceStyle::GasBanded, 0xE3D3A3, 0xB89A62,
            {9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448, -0.00125060, -0.00050991,
             0.00193609, 1222.49362201, -0.41897216, -0.28867794},
            "Main rings C+B+A: 74,658 - 136,775 km");
        sat.rings = {true, 74658.0, 136775.0, 1.0};
        v.push_back(sat);
    }
    {
        BodyDefinition ura = makePlanet(
            "Uranus", BodyKind::Planet, 5794556.4, 25362.0, -17.24, 257.311, -15.175, 0.488, SurfaceStyle::IceGiant,
            0x9FD6E4, 0x7FB9CC,
            {19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503, -0.00196176,
             -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589},
            "Axial tilt 97.8 deg");
        ura.rings = {true, 41837.0, 51149.0, 0.18};
        v.push_back(ura);
    }
    v.push_back(makePlanet("Neptune", BodyKind::Planet, 6836527.10058 - 1428.495, 24622.0, 16.11, 299.36, 43.46,
                           0.442, SurfaceStyle::IceGiant, 0x4B70DD, 0x3354B0,
                           {30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574,
                            0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664},
                           "GM excludes Triton (added back when not simulated)"));
    v.push_back(makePlanet("Pluto", BodyKind::DwarfPlanet, 975.5, 1188.3, -153.2935, 132.993, -6.163, 0.52,
                           SurfaceStyle::Icy, 0xD9C0A0, 0x8A6D55,
                           {39.48211675, 0.24882730, 17.14001206, 238.92903833, 224.06891629, 110.30393684,
                            -0.00031596, 0.00005170, 0.00004818, 145.20780515, -0.04062942, -0.01183482},
                           "GM of the Pluto-Charon system; 3:2 resonance with Neptune"));

    // ---------------------------------------------------------------------------------
    // Major moons (approximate phases; orbits in the parent's equatorial plane).
    // Galilean mean longitudes from Meeus ch. 44 (preserves the Laplace resonance
    // lambda_Io - 3 lambda_Europa + 2 lambda_Ganymede = 180 deg).
    // ---------------------------------------------------------------------------------
    const double lamJ = 34.39644051;
    v.push_back(makeSatellite("Io", "Jupiter", 5959.916, 1821.6, 42.459, 421700.0, 0.0041, 0.05, lamJ + 163.8069,
                              SurfaceStyle::Volcanic, 0xE8D86A, 0xB5651D, "Most volcanically active body known"));
    v.push_back(makeSatellite("Europa", "Jupiter", 3202.739, 1560.8, 85.228, 671034.0, 0.0094, 0.47,
                              lamJ + 358.4140, SurfaceStyle::Icy, 0xDCCDB0, 0x9C7A5A,
                              "Subsurface ocean beneath an ice shell"));
    v.push_back(makeSatellite("Ganymede", "Jupiter", 9887.834, 2634.1, 171.71, 1070412.0, 0.0013, 0.20,
                              lamJ + 5.7176, SurfaceStyle::Cratered, 0xA59C8C, 0x6E6458,
                              "Largest moon in the Solar System"));
    v.push_back(makeSatellite("Callisto", "Jupiter", 7179.289, 2410.3, 400.54, 1882709.0, 0.0074, 0.19,
                              lamJ + 224.8092, SurfaceStyle::Cratered, 0x7A6E60, 0x3E3830,
                              "Most heavily cratered surface known"));
    v.push_back(makeSatellite("Titan", "Saturn", 8978.137, 2574.7, 382.69, 1221870.0, 0.0288, 0.35, 0.0,
                              SurfaceStyle::Cloudy, 0xD9A441, 0xB07A30, "Thick N2 atmosphere, methane lakes"));
    v.push_back(makeSatellite("Triton", "Neptune", 1428.495, 1353.4, -141.04, 354759.0, 0.000016, 156.885, 0.0,
                              SurfaceStyle::Icy, 0xD8C4C0, 0x9A8580,
                              "Retrograde orbit - likely a captured Kuiper belt object"));
    for (auto& b : v) {
        if (b.parent == "Jupiter") { b.poleRaDeg = 268.05; b.poleDecDeg = 64.49; }
        if (b.name == "Titan") { b.poleRaDeg = 39.48; b.poleDecDeg = 83.43; }
        if (b.name == "Triton") { b.poleRaDeg = 299.36; b.poleDecDeg = 41.17; }
    }

    // ---------------------------------------------------------------------------------
    // Minor bodies: q[AU] e i node argp[deg] tp[JD]
    // ---------------------------------------------------------------------------------
    {
        BodyDefinition ceres = makeMinor("Ceres", BodyKind::DwarfPlanet, 62.6284, 469.7, 9.074, 2.7675 * (1 - 0.0785),
                                         0.0785, 10.587, 80.255, 73.42, 2459920.5, SurfaceStyle::Cratered, 0x8C8A86,
                                         0x5E5C58, "Largest object in the asteroid belt");
        ceres.poleRaDeg = 291.418; ceres.poleDecDeg = 66.764;
        v.push_back(ceres);
        BodyDefinition vesta = makeMinor("Vesta", BodyKind::Asteroid, 17.2882, 262.7, 5.342, 2.3615 * (1 - 0.0887),
                                         0.0887, 7.142, 103.81, 151.66, 2459574.5, SurfaceStyle::Rocky, 0xA09A8F,
                                         0x6A655C, "Differentiated protoplanet; source of HED meteorites");
        vesta.poleRaDeg = 309.031; vesta.poleDecDeg = 42.235;
        v.push_back(vesta);
    }
    // Comet nuclei (~1e14 kg) are dynamically massless test particles here.
    v.push_back(makeMinor("Halley", BodyKind::Comet, 0.0, 5.5, 52.8, 0.58598, 0.96714, 162.26, 58.42,
                          111.33, 2446470.95, SurfaceStyle::Comet, 0xCFE3F5, 0x7FA9D0,
                          "1P/Halley - retrograde, period ~76 yr; perihelion 1986-02-09"));
    v.push_back(makeMinor("Encke", BodyKind::Comet, 0.0, 2.4, 11.0, 0.33593, 0.84833, 11.78, 334.57, 186.55,
                          2460239.5, SurfaceStyle::Comet, 0xCFE3F5, 0x7FA9D0,
                          "2P/Encke - shortest known period of a bright comet (3.3 yr)"));
    return v;
}

} // namespace

const std::vector<BodyDefinition>& BodyCatalog::all() {
    static const std::vector<BodyDefinition> catalog = buildCatalog();
    return catalog;
}

const BodyDefinition* BodyCatalog::find(const std::string& name) {
    for (const auto& b : all())
        if (str::iequals(b.name, name)) return &b;
    return nullptr;
}

std::vector<std::string> BodyCatalog::defaultSet() {
    return {"Sun", "Mercury", "Venus", "Earth", "Moon", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};
}

std::vector<std::string> BodyCatalog::majorMoonSet() { return {"Io", "Europa", "Ganymede", "Callisto", "Titan", "Triton"}; }

std::vector<std::string> BodyCatalog::minorBodySet() { return {"Ceres", "Vesta", "Halley", "Encke"}; }

OrbitalElements meanElementsAt(const MeanElements& m, double jdTT, double mu) {
    const double T = julianCenturies(jdTT);
    const double a = m.a + m.da * T;
    const double e = m.e + m.de * T;
    const double I = (m.I + m.dI * T) * DEG2RAD;
    const double L = (m.L + m.dL * T) * DEG2RAD;
    const double varpi = (m.varpi + m.dvarpi * T) * DEG2RAD;
    const double node = (m.node + m.dnode * T) * DEG2RAD;
    return OrbitalElements::fromMeanLongitude(a, e, I, node, varpi, math::wrapTwoPi(L), jdTT, mu);
}

StateVector relativeStateAt(const BodyDefinition& def, double jdTT, double muCentral,
                            const BodyDefinition* parentDef) {
    switch (def.model) {
        case OrbitModel::None: return {};
        case OrbitModel::MeanElements: {
            const OrbitalElements el = meanElementsAt(def.mean, jdTT, muCentral);
            return elementsToState(el, jdTT);
        }
        case OrbitModel::Osculating: {
            OrbitalElements el;
            el.q = def.q;
            el.e = def.e;
            el.i = def.i * DEG2RAD;
            el.node = def.node * DEG2RAD;
            el.argp = def.argp * DEG2RAD;
            el.tp = def.tp;
            el.mu = muCentral;
            el.epoch = jdTT;
            return elementsToState(el, jdTT);
        }
        case OrbitModel::LunarTheory: {
            StateVector s;
            lunarStateJ2000(jdTT, s.r, s.v);
            return s;
        }
        case OrbitModel::Equatorial: {
            // Build the parent's equatorial frame: n = ascending node of the equator on
            // the ecliptic, p = pole, m = p x n.
            Vec3 pole = Vec3::unitZ();
            if (parentDef) pole = poleFromRaDec(parentDef->poleRaDeg, parentDef->poleDecDeg);
            Vec3 nodeDir = Vec3::unitZ().cross(pole);
            if (nodeDir.norm() < 1e-9) nodeDir = Vec3::unitX();
            nodeDir = nodeDir.normalized();
            const Vec3 mDir = pole.cross(nodeDir).normalized();
            const double nodeLon = std::atan2(nodeDir.y, nodeDir.x);

            const double a = def.aKm / AU_KM;
            const double n = std::sqrt(muCentral / (a * a * a));
            const double lambda = def.lambda0 * DEG2RAD + n * (jdTT - J2000_JD);
            // Elements in the equatorial frame (node/argp = 0: angle measured from the equator node).
            OrbitalElements el = OrbitalElements::fromMeanAnomaly(a, def.e, def.incl * DEG2RAD, 0.0, 0.0,
                                                                  math::wrapPi(lambda - nodeLon), jdTT, muCentral);
            const StateVector local = elementsToState(el, jdTT);
            auto toEcl = [&](const Vec3& v) { return nodeDir * v.x + mDir * v.y + pole * v.z; };
            return {toEcl(local.r), toEcl(local.v)};
        }
    }
    return {};
}

} // namespace solsim::astro
