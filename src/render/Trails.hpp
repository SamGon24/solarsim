// =============================================================================
//  SolarSim - Trails.hpp
//  Ring buffers of recent positions for every named body. Points are stored
//  RELATIVE TO THE PARENT at the time of recording, so a moon's trail is drawn
//  as a loop around its planet and a planet's trail as an arc around the Sun.
// =============================================================================
#pragma once

#include "../core/Vec3.hpp"
#include "../physics/NBodySystem.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace solsim::render {

class TrailRecorder {
public:
    explicit TrailRecorder(std::size_t capacity = 320) : m_capacity(capacity) {}

    void setCapacity(std::size_t c);
    void clear();
    // Records a sample for each body whose sampling interval has elapsed.
    void record(const physics::NBodySystem& sys, double time);

    std::size_t count(std::size_t body) const { return body < m_trails.size() ? m_trails[body].count : 0; }
    int parent(std::size_t body) const { return body < m_trails.size() ? m_trails[body].parent : -1; }
    // k = 0 is the oldest sample.
    const Vec3& point(std::size_t body, std::size_t k) const;

private:
    struct Trail {
        std::vector<Vec3> rel;
        std::size_t head = 0;   // next write position
        std::size_t count = 0;
        int parent = -1;
        double lastTime = 0.0;
        double interval = 1.0;
        double period = 1.0;    // orbital period estimate (days)
        bool started = false;
    };
    void sync(const physics::NBodySystem& sys);

    std::size_t m_capacity;
    std::vector<Trail> m_trails;
    std::uint64_t m_revision = ~0ULL;
};

} // namespace solsim::render
