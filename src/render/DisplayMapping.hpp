// =============================================================================
//  SolarSim - DisplayMapping.hpp
//  The Solar System spans ~5 orders of magnitude (Moon orbit 0.0026 AU,
//  Neptune 30 AU) and bodies are tiny compared to their orbits, so a
//  physically linear picture shows almost nothing. This class maps physical
//  positions to "display space":
//
//    * heliocentric distances are compressed radially:  r' = f(r)
//        linear  f(r) = r
//        sqrt    f(r) = sqrt(r)                    (default)
//        log     f(r) = ln(1 + r/r0) / ln(1 + 1/r0)
//      every mode keeps f(1 AU) = 1
//    * satellite offsets are magnified logarithmically around their parent
//      so moons stay outside the (exaggerated) planet disc
//    * body radii are exaggerated with a power law (R/R_earth)^0.45
//
//  With --true-scale and the linear mode, the mapping is the identity.
// =============================================================================
#pragma once

#include "../core/Vec3.hpp"
#include "../physics/NBodySystem.hpp"

#include <string>
#include <vector>

namespace solsim::render {

enum class ScaleMode { Linear, Sqrt, Log };
const char* scaleModeName(ScaleMode m);
bool parseScaleMode(const std::string& s, ScaleMode& out);

class DisplayMapping {
public:
    void setMode(ScaleMode m) { m_mode = m; }
    ScaleMode mode() const { return m_mode; }
    void setTrueScale(bool on) { m_trueScale = on; }
    bool trueScale() const { return m_trueScale; }
    void setBodyScale(double s) { m_bodyScale = s; }
    double bodyScale() const { return m_bodyScale; }

    // Must be called when the system structure changes.
    void prepare(const physics::NBodySystem& sys);

    double radial(double r) const;           // f(r)
    Vec3 mapHelio(const Vec3& rel) const;    // heliocentric vector -> display
    // Display position of any point index (named bodies or particles).
    Vec3 position(const physics::NBodySystem& sys, std::size_t index) const;
    // Display position of a point given relative to `parent` (trails/orbits).
    Vec3 relativeToParent(const physics::NBodySystem& sys, int parent, const Vec3& rel) const;
    double displayRadius(const physics::NBodySystem& sys, std::size_t index) const;

private:
    Vec3 satelliteOffset(int parent, const Vec3& rel, const physics::NBodySystem& sys) const;

    ScaleMode m_mode = ScaleMode::Sqrt;
    bool m_trueScale = false;
    double m_bodyScale = 1.0;
    std::vector<double> m_minChildDistance; // per named body
    std::uint64_t m_preparedRevision = ~0ULL;
};

} // namespace solsim::render
