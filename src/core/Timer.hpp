// =============================================================================
//  SolarSim - Timer.hpp
//  Wall-clock stopwatch based on std::chrono::steady_clock.
// =============================================================================
#pragma once

#include <chrono>

namespace solsim {

class Stopwatch {
public:
    using Clock = std::chrono::steady_clock;

    Stopwatch() : m_start(Clock::now()) {}
    void restart() { m_start = Clock::now(); }
    double seconds() const { return std::chrono::duration<double>(Clock::now() - m_start).count(); }
    double milliseconds() const { return seconds() * 1e3; }

private:
    Clock::time_point m_start;
};

// Exponential moving average - used for smoothed FPS / step-rate displays.
class RateMeter {
public:
    explicit RateMeter(double smoothing = 0.1) : m_alpha(smoothing) {}
    void sample(double value) {
        m_value = m_initialized ? (1.0 - m_alpha) * m_value + m_alpha * value : value;
        m_initialized = true;
    }
    double value() const { return m_value; }

private:
    double m_alpha;
    double m_value = 0.0;
    bool m_initialized = false;
};

} // namespace solsim
