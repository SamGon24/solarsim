// =============================================================================
//  SolarSim - DisplayMapping.cpp
// =============================================================================
#include "DisplayMapping.hpp"

#include "../core/Constants.hpp"
#include "../core/StringUtil.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace solsim::render {

const char* scaleModeName(ScaleMode m) {
    switch (m) {
        case ScaleMode::Linear: return "linear";
        case ScaleMode::Sqrt: return "sqrt";
        case ScaleMode::Log: return "log";
    }
    return "?";
}

bool parseScaleMode(const std::string& s, ScaleMode& out) {
    const std::string k = str::toLower(s);
    if (k == "linear" || k == "lin") out = ScaleMode::Linear;
    else if (k == "sqrt" || k == "root") out = ScaleMode::Sqrt;
    else if (k == "log" || k == "logarithmic") out = ScaleMode::Log;
    else return false;
    return true;
}

void DisplayMapping::prepare(const physics::NBodySystem& sys) {
    if (sys.revision() == m_preparedRevision && m_minChildDistance.size() == sys.namedCount()) return;
    m_preparedRevision = sys.revision();
    m_minChildDistance.assign(sys.namedCount(), std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < sys.namedCount(); ++i) {
        const int p = sys.info(i).parent;
        if (p <= 0) continue;
        double d = sys.info(i).typicalDistance;
        if (!(d > 0.0)) d = sys.relativePosition(i, p).norm();
        auto& slot = m_minChildDistance[static_cast<std::size_t>(p)];
        slot = std::min(slot, d);
    }
}

double DisplayMapping::radial(double r) const {
    if (m_trueScale && m_mode == ScaleMode::Linear) return r;
    switch (m_mode) {
        case ScaleMode::Linear: return r;
        case ScaleMode::Sqrt: return std::sqrt(r);
        case ScaleMode::Log: {
            constexpr double r0 = 0.5;
            static const double norm = std::log(1.0 + 1.0 / r0);
            return std::log(1.0 + r / r0) / norm;
        }
    }
    return r;
}

Vec3 DisplayMapping::mapHelio(const Vec3& rel) const {
    const double r = rel.norm();
    if (r <= 0.0) return rel;
    return rel * (radial(r) / r);
}

double DisplayMapping::displayRadius(const physics::NBodySystem& sys, std::size_t index) const {
    if (index >= sys.namedCount()) return 0.0;
    const double R = sys.info(index).radius;
    if (m_trueScale) {
        // True radius at the local radial compression of the body's position.
        return R;
    }
    const double rEarth = constants::EARTH_RADIUS_KM / constants::AU_KM;
    const double ratio = std::max(R / rEarth, 1e-4);
    return 0.022 * m_bodyScale * std::pow(ratio, 0.45);
}

Vec3 DisplayMapping::satelliteOffset(int parent, const Vec3& rel, const physics::NBodySystem& sys) const {
    if (m_trueScale) {
        // Keep true proportions relative to the planet: scale by local compression.
        return rel;
    }
    const double d = rel.norm();
    if (d <= 0.0) return rel;
    const std::size_t p = static_cast<std::size_t>(parent);
    double dmin = p < m_minChildDistance.size() ? m_minChildDistance[p] : d;
    if (!std::isfinite(dmin) || dmin <= 0.0) dmin = d;
    const double Rp = displayRadius(sys, p);
    // Logarithmic compression of the satellite system: innermost moon at ~2.6 R_p.
    const double D = Rp * std::max(1.3, 2.6 + 2.2 * std::log2(d / dmin));
    return rel * (D / d);
}

Vec3 DisplayMapping::relativeToParent(const physics::NBodySystem& sys, int parent, const Vec3& rel) const {
    if (parent <= 0) return mapHelio(rel);
    return position(sys, static_cast<std::size_t>(parent)) + satelliteOffset(parent, rel, sys);
}

Vec3 DisplayMapping::position(const physics::NBodySystem& sys, std::size_t index) const {
    if (sys.size() == 0) return {};
    if (index == 0) return {};
    const Vec3 sun = sys.pos(0);
    if (index >= sys.namedCount()) return mapHelio(sys.pos(index) - sun);
    const int p = sys.info(index).parent;
    if (p <= 0) return mapHelio(sys.pos(index) - sun);
    return relativeToParent(sys, p, sys.pos(index) - sys.pos(static_cast<std::size_t>(p)));
}

} // namespace solsim::render
