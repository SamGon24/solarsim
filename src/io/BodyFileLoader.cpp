// =============================================================================
//  SolarSim - BodyFileLoader.cpp
// =============================================================================
#include "BodyFileLoader.hpp"

#include "../astro/JulianDate.hpp"
#include "../astro/Kepler.hpp"
#include "../core/Constants.hpp"
#include "../core/IniParser.hpp"
#include "../core/MathUtil.hpp"
#include "../core/StringUtil.hpp"

#include <cmath>
#include <stdexcept>

namespace solsim::io {

using namespace solsim::constants;

namespace {

bool parseColor(const std::string& s, std::uint32_t& out) {
    std::string t = str::trim(s);
    if (!t.empty() && t[0] == '#') t = t.substr(1);
    if (str::startsWith(str::toLower(t), "0x")) t = t.substr(2);
    if (t.size() != 6) return false;
    std::uint32_t v = 0;
    for (char c : t) {
        v <<= 4;
        if (c >= '0' && c <= '9') v |= static_cast<std::uint32_t>(c - '0');
        else if (c >= 'a' && c <= 'f') v |= static_cast<std::uint32_t>(c - 'a' + 10);
        else if (c >= 'A' && c <= 'F') v |= static_cast<std::uint32_t>(c - 'A' + 10);
        else return false;
    }
    out = v;
    return true;
}

double parseDateValue(const IniSection& sec, const std::string& key) {
    double jd = 0.0;
    std::string err;
    if (!astro::parseDate(sec.get(key), jd, err))
        throw std::runtime_error("line " + std::to_string(sec.line) + ": '" + key + "': " + err);
    return jd;
}

astro::BodyDefinition parseBody(const IniSection& sec, std::vector<std::string>& warnings) {
    const std::string where = "[body] at line " + std::to_string(sec.line);
    astro::BodyDefinition b;
    b.name = sec.get("name");
    if (b.name.empty()) throw std::runtime_error(where + ": missing 'name'");
    b.parent = sec.get("parent", "Sun");
    if (str::iequals(b.name, "Sun")) throw std::runtime_error(where + ": the Sun cannot be redefined");
    if (sec.has("kind") && !astro::parseKind(sec.get("kind"), b.kind))
        throw std::runtime_error(where + ": unknown kind '" + sec.get("kind") + "'");
    if (!sec.has("kind")) b.kind = astro::BodyKind::Asteroid;

    // ---- mass ----
    if (sec.has("gm_km3s2")) b.gmKm3s2 = sec.getDouble("gm_km3s2", 0.0);
    else if (sec.has("mass_kg")) b.gmKm3s2 = sec.getDouble("mass_kg", 0.0) * G_SI * 1e-9;
    else if (sec.has("mass_earths")) b.gmKm3s2 = sec.getDouble("mass_earths", 0.0) * 398600.435507;
    else if (sec.has("mass_suns")) b.gmKm3s2 = sec.getDouble("mass_suns", 0.0) * GM_SUN_KM3S2;
    if (b.gmKm3s2 < 0.0) throw std::runtime_error(where + ": negative mass");

    // ---- physical / appearance ----
    b.radiusKm = sec.getDouble("radius_km", 10.0);
    b.rotationHours = sec.getDouble("rotation_hours", 24.0);
    b.poleRaDeg = sec.getDouble("pole_ra", 0.0);
    b.poleDecDeg = sec.getDouble("pole_dec", 90.0);
    b.albedo = sec.getDouble("albedo", 0.1);
    b.notes = sec.get("notes", "loaded from file");
    b.approximate = true;
    b.style = b.kind == astro::BodyKind::Comet ? astro::SurfaceStyle::Comet : astro::SurfaceStyle::Rocky;
    if (sec.has("style") && !astro::parseSurfaceStyle(sec.get("style"), b.style))
        warnings.push_back(where + ": unknown style '" + sec.get("style") + "' (using default)");
    b.color = b.kind == astro::BodyKind::Comet ? 0xCFE3F5 : 0xB0A090;
    b.color2 = 0x605850;
    if (sec.has("color") && !parseColor(sec.get("color"), b.color))
        warnings.push_back(where + ": bad color '" + sec.get("color") + "'");
    if (sec.has("color2") && !parseColor(sec.get("color2"), b.color2))
        warnings.push_back(where + ": bad color2 '" + sec.get("color2") + "'");
    if (sec.has("rings_outer_km")) {
        b.rings.present = true;
        b.rings.innerKm = sec.getDouble("rings_inner_km", b.radiusKm * 1.2);
        b.rings.outerKm = sec.getDouble("rings_outer_km", b.radiusKm * 2.0);
        b.rings.opacity = sec.getDouble("rings_opacity", 0.6);
    }

    // ---- orbit ----
    const std::string frame = str::toLower(sec.get("frame", "ecliptic"));
    if (frame == "equatorial") {
        b.model = astro::OrbitModel::Equatorial;
        b.aKm = sec.getDouble("a_km", 0.0);
        if (b.aKm <= 0.0) throw std::runtime_error(where + ": equatorial orbits need a_km > 0");
        b.e = sec.getDouble("e", 0.0);
        b.incl = sec.getDouble("i", 0.0);
        b.lambda0 = sec.getDouble("lambda0", 0.0);
        return b;
    }
    if (frame != "ecliptic") throw std::runtime_error(where + ": frame must be 'ecliptic' or 'equatorial'");

    b.model = astro::OrbitModel::Osculating;
    b.e = sec.getDouble("e", 0.0);
    if (b.e < 0.0) throw std::runtime_error(where + ": eccentricity must be >= 0");
    if (sec.has("q")) {
        b.q = sec.getDouble("q", 1.0);
    } else if (sec.has("a")) {
        const double a = sec.getDouble("a", 1.0);
        if (std::fabs(b.e - 1.0) < 1e-12) throw std::runtime_error(where + ": parabolic orbits need q, not a");
        b.q = a * (1.0 - b.e);
    } else if (sec.has("a_km")) {
        b.q = sec.getDouble("a_km", 0.0) / AU_KM * (1.0 - b.e);
    } else {
        throw std::runtime_error(where + ": orbit needs 'q' or 'a' (AU)");
    }
    if (!(b.q > 0.0)) throw std::runtime_error(where + ": perihelion distance must be positive (check a/e)");
    b.i = sec.getDouble("i", 0.0);
    b.node = sec.getDouble("node", 0.0);
    if (sec.has("argp")) b.argp = sec.getDouble("argp", 0.0);
    else if (sec.has("peri")) b.argp = sec.getDouble("peri", 0.0) - b.node; // longitude of perihelion
    else b.argp = 0.0;

    // Time reference -> time of perihelion passage.
    const double mu = GM_SUN; // used only to convert anomalies to times (parent mass refined later)
    double muCentral = mu;
    if (!str::iequals(b.parent, "Sun")) {
        if (const auto* p = astro::BodyCatalog::find(b.parent)) muCentral = p->gm();
    }
    if (sec.has("tp")) {
        b.tp = parseDateValue(sec, "tp");
    } else if (sec.has("m") || sec.has("l")) {
        const double epoch = sec.has("epoch") ? parseDateValue(sec, "epoch") : J2000_JD;
        if (std::fabs(b.e - 1.0) < 1e-9) throw std::runtime_error(where + ": use tp for parabolic orbits");
        const double a = b.q / (1.0 - b.e);
        const double n = std::sqrt(muCentral / std::fabs(a * a * a)); // rad/day
        double M = sec.has("m") ? sec.getDouble("m", 0.0) : sec.getDouble("l", 0.0) - (b.node + b.argp);
        M *= DEG2RAD;
        if (b.e < 1.0) M = math::wrapPi(M);
        b.tp = epoch - M / n;
    } else {
        warnings.push_back(where + ": no tp/M/L given - assuming perihelion at J2000");
        b.tp = J2000_JD;
    }
    return b;
}

} // namespace

namespace {
LoadResult parseDocument(const IniDocument& doc, const std::string& sourceName) {
    LoadResult res;
    try {
        for (const auto& sec : doc.sections()) {
            if (sec.name != "body") {
                res.warnings.push_back(sourceName + ": ignoring section [" + sec.name + "]");
                continue;
            }
            res.bodies.push_back(parseBody(sec, res.warnings));
        }
    } catch (const std::exception& e) {
        res.error = sourceName + ": " + e.what();
        res.bodies.clear();
    }
    return res;
}
} // namespace

LoadResult loadBodyString(const std::string& text, const std::string& sourceName) {
    IniDocument doc;
    std::string err;
    if (!doc.loadString(text, err, sourceName)) {
        LoadResult res;
        res.error = err;
        return res;
    }
    return parseDocument(doc, sourceName);
}

LoadResult loadBodyFile(const std::string& path) {
    IniDocument doc;
    std::string err;
    if (!doc.loadFile(path, err)) {
        LoadResult res;
        res.error = err;
        return res;
    }
    return parseDocument(doc, path);
}

} // namespace solsim::io
