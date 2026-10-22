// =============================================================================
//  SolarSim - Simulation.hpp
//  Owns the N-body system, the force models and an integrator, and advances
//  the simulated clock (JD, TDB). Copyable: copies are used by the event
//  finder to re-integrate short intervals while bracketing roots.
// =============================================================================
#pragma once

#include "Diagnostics.hpp"
#include "Forces.hpp"
#include "Integrators.hpp"
#include "NBodySystem.hpp"

#include <cstdint>
#include <functional>
#include <memory>
#include <string>

namespace solsim {
class ThreadPool;
}

namespace solsim::physics {

struct SimulationConfig {
    std::string integrator = "yoshida4";
    double dt = 0.0;          // days; 0 = automatic (from the shortest orbital period)
    double tolerance = 1e-10; // adaptive integrators only
    bool relativity = false;  // add the 1PN GR correction
    double softening = 0.0;   // AU (Plummer)
    unsigned threads = 0;     // 0 = hardware concurrency, 1 = single threaded
};

class Simulation {
public:
    Simulation(NBodySystem system, double epochJD, const SimulationConfig& config = {});
    ~Simulation();
    Simulation(const Simulation& other);
    Simulation& operator=(const Simulation& other);
    Simulation(Simulation&&) noexcept;
    Simulation& operator=(Simulation&&) noexcept;

    // Integrates by `days` (negative = backwards). Splits into integrator steps.
    void advance(double days);
    void advanceTo(double jd) { advance(jd - m_time); }

    double time() const { return m_time; }
    double epoch() const { return m_epoch; }
    double elapsed() const { return m_time - m_epoch; }

    double timeStep() const { return m_dt; }
    void setTimeStep(double dt);
    static double autoTimeStep(const NBodySystem& sys);

    NBodySystem& system() { return m_system; }
    const NBodySystem& system() const { return m_system; }
    ForceSet& forces() { return m_forces; }
    const ForceSet& forces() const { return m_forces; }
    Integrator& integrator() { return *m_integrator; }
    const Integrator& integrator() const { return *m_integrator; }
    // Throws std::invalid_argument for unknown names.
    void setIntegrator(const std::string& name);

    void setRelativity(bool on);
    bool relativity() const { return m_forces.has("gr1pn"); }
    NewtonianGravity& gravity();

    const Invariants& initialInvariants() const { return m_initial; }
    Invariants currentInvariants() const { return computeInvariants(m_system); }
    InvariantErrors invariantErrors() const;
    void rebaseInvariants(); // call after structural changes

    std::uint64_t totalSteps() const { return m_totalSteps; }
    std::uint64_t forceEvaluations() const { return m_forces.evaluations(); }
    double lastStepTaken() const { return m_lastStep; }

    // Called after every integrator step (optional; used by the precession tracker).
    using StepCallback = std::function<void(const Simulation&)>;
    void setStepCallback(StepCallback cb) { m_stepCallback = std::move(cb); }

    // Notify the simulation that the system changed externally.
    void systemModified();

    const SimulationConfig& config() const { return m_config; }

private:
    void configureForces();

    NBodySystem m_system;
    ForceSet m_forces;
    std::unique_ptr<Integrator> m_integrator;
    std::shared_ptr<ThreadPool> m_pool;
    SimulationConfig m_config;
    double m_epoch = 0.0;
    double m_time = 0.0;
    double m_dt = 0.1;
    double m_lastStep = 0.0;
    std::uint64_t m_totalSteps = 0;
    Invariants m_initial;
    StepCallback m_stepCallback;
};

} // namespace solsim::physics
