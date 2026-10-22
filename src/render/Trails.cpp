// =============================================================================
//  SolarSim - Trails.cpp
// =============================================================================
#include "Trails.hpp"

#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"

#include <cmath>

namespace solsim::render {

void TrailRecorder::setCapacity(std::size_t c) {
    m_capacity = c < 2 ? 2 : c;
    clear();
}

void TrailRecorder::clear() {
    m_trails.clear();
    m_revision = ~0ULL;
}

void TrailRecorder::sync(const physics::NBodySystem& sys) {
    if (m_revision == sys.revision() && m_trails.size() == sys.namedCount()) return;
    m_revision = sys.revision();
    m_trails.assign(sys.namedCount(), Trail{});
    for (std::size_t i = 0; i < sys.namedCount(); ++i) {
        Trail& t = m_trails[i];
        t.rel.assign(m_capacity, Vec3{});
        t.parent = sys.info(i).parent;
        if (t.parent < 0) continue;
        // Aim for ~1.25 orbits of history (at least one sample every 100 days).
        const double d = std::max(1e-6, sys.info(i).typicalDistance > 0 ? sys.info(i).typicalDistance
                                                                        : sys.relativePosition(i, t.parent).norm());
        const double mu = std::max(1e-30, sys.twoBodyMu(i));
        const double period = constants::TWO_PI * std::sqrt(d * d * d / mu);
        t.period = period;
        t.interval = math::clamp(1.25 * period / static_cast<double>(m_capacity), 1e-3, 100.0);
    }
}

void TrailRecorder::record(const physics::NBodySystem& sys, double time) {
    sync(sys);
    for (std::size_t i = 0; i < m_trails.size(); ++i) {
        Trail& t = m_trails[i];
        if (t.parent < 0) continue;
        const double gap = std::fabs(time - t.lastTime);
        if (t.started && gap < t.interval) continue;
        // When time advances by a large fraction of an orbit between samples the
        // polyline would be meaningless (random chords) - restart the trail.
        if (t.started && gap > t.period / 16.0) {
            t.count = 0;
            t.head = 0;
        }
        t.rel[t.head] = sys.pos(i) - sys.pos(static_cast<std::size_t>(t.parent));
        t.head = (t.head + 1) % m_capacity;
        if (t.count < m_capacity) ++t.count;
        t.lastTime = time;
        t.started = true;
    }
}

const Vec3& TrailRecorder::point(std::size_t body, std::size_t k) const {
    const Trail& t = m_trails[body];
    const std::size_t start = (t.head + m_capacity - t.count) % m_capacity;
    return t.rel[(start + k) % m_capacity];
}

} // namespace solsim::render
