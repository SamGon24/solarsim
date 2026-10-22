// =============================================================================
//  SolarSim - Simulation.cpp
// =============================================================================
#include "Simulation.hpp"

#include "../astro/OrbitalElements.hpp"
#include "../core/Constants.hpp"
#include "../core/Log.hpp"
#include "../core/MathUtil.hpp"
#include "../core/ThreadPool.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace solsim::physics {

Simulation::Simulation(NBodySystem system, double epochJD, const SimulationConfig& config)
    : m_system(std::move(system)), m_config(config), m_epoch(epochJD), m_time(epochJD) {
    if (config.threads != 1 && m_system.particleCount() >= 256)
        m_pool = std::make_shared<ThreadPool>(config.threads);
    configureForces();
    setIntegrator(config.integrator);
    setTimeStep(config.dt > 0.0 ? config.dt : autoTimeStep(m_system));
    m_initial = computeInvariants(m_system);
}

Simulation::~Simulation() = default;

Simulation::Simulation(const Simulation& o)
    : m_system(o.m_system),
      m_forces(o.m_forces),
      m_integrator(o.m_integrator ? o.m_integrator->clone() : nullptr),
      m_pool(o.m_pool),
      m_config(o.m_config),
      m_epoch(o.m_epoch),
      m_time(o.m_time),
      m_dt(o.m_dt),
      m_lastStep(o.m_lastStep),
      m_totalSteps(o.m_totalSteps),
      m_initial(o.m_initial),
      m_stepCallback() {}

Simulation& Simulation::operator=(const Simulation& o) {
    if (this == &o) return *this;
    Simulation tmp(o);
    *this = std::move(tmp);
    return *this;
}

Simulation::Simulation(Simulation&&) noexcept = default;
Simulation& Simulation::operator=(Simulation&&) noexcept = default;

void Simulation::configureForces() {
    m_forces = ForceSet();
    m_forces.add(std::make_unique<NewtonianGravity>(m_config.softening, m_pool.get()));
    if (m_config.relativity) m_forces.add(std::make_unique<PostNewtonianGR>(0));
}

NewtonianGravity& Simulation::gravity() {
    NewtonianGravity* g = m_forces.gravity();
    if (!g) throw std::logic_error("simulation has no Newtonian gravity model");
    return *g;
}

void Simulation::setIntegrator(const std::string& name) {
    auto integ = createIntegrator(name);
    if (!integ) throw std::invalid_argument("unknown integrator '" + name + "'");
    integ->setTolerance(m_config.tolerance);
    m_integrator = std::move(integ);
    m_config.integrator = m_integrator->name();
    if (m_integrator->name() == "wh" && relativity())
        log::warn("the Wisdom-Holman integrator ignores velocity dependent forces (GR correction disabled in its steps)");
}

void Simulation::setRelativity(bool on) {
    if (on == relativity()) return;
    if (on) m_forces.add(std::make_unique<PostNewtonianGR>(0));
    else m_forces.remove("gr1pn");
    m_config.relativity = on;
}

void Simulation::setTimeStep(double dt) {
    if (!(dt > 0.0) || !std::isfinite(dt)) throw std::invalid_argument("time step must be positive");
    m_dt = dt;
}

double Simulation::autoTimeStep(const NBodySystem& sys) {
    // Resolve the shortest dynamical timescale: the circular period at the
    // periapsis distance of every named orbiting body.
    double best = 5.0 * 150.0; // cap: 5 days
    for (std::size_t i = 0; i < sys.namedCount(); ++i) {
        const int p = sys.info(i).parent;
        if (p < 0) continue;
        const Vec3 r = sys.relativePosition(i, p);
        const Vec3 v = sys.relativeVelocity(i, p);
        const double mu = sys.twoBodyMu(i);
        if (mu <= 0.0) continue;
        const astro::OrbitalElements el = astro::stateToElements(r, v, mu, 0.0);
        double d = std::min(r.norm(), el.q);
        d = std::max(d, 1e-6);
        const double period = constants::TWO_PI * std::sqrt(d * d * d / mu);
        best = std::min(best, period);
    }
    return math::clamp(best / 150.0, 1e-4, 5.0);
}

void Simulation::advance(double days) {
    if (days == 0.0 || !std::isfinite(days)) return;
    if (m_integrator->isAdaptive()) {
        double remaining = days;
        const double dir = days > 0 ? 1.0 : -1.0;
        int guard = 0;
        while (remaining * dir > 1e-12) {
            const double maxStep = dir * std::min(std::fabs(remaining), m_dt * 50.0);
            const double taken = m_integrator->step(m_system, m_forces, maxStep);
            remaining -= taken;
            m_time += taken;
            m_lastStep = taken;
            ++m_totalSteps;
            if (m_stepCallback) m_stepCallback(*this);
            if (++guard > 100000000) throw std::runtime_error("advance: step limit exceeded");
        }
        return;
    }
    const long long nSteps = std::max<long long>(1, static_cast<long long>(std::ceil(std::fabs(days) / m_dt - 1e-9)));
    const double h = days / static_cast<double>(nSteps);
    const double start = m_time;
    for (long long k = 0; k < nSteps; ++k) {
        m_integrator->step(m_system, m_forces, h);
        m_time = start + h * static_cast<double>(k + 1);
        m_lastStep = h;
        ++m_totalSteps;
        if (m_stepCallback) m_stepCallback(*this);
    }
}

InvariantErrors Simulation::invariantErrors() const { return compareInvariants(m_initial, computeInvariants(m_system)); }

void Simulation::rebaseInvariants() { m_initial = computeInvariants(m_system); }

void Simulation::systemModified() {
    m_integrator->reset();
    if (!m_pool && m_config.threads != 1 && m_system.particleCount() >= 256) {
        m_pool = std::make_shared<ThreadPool>(m_config.threads);
        if (auto* g = m_forces.gravity()) g->setThreadPool(m_pool.get());
    }
    rebaseInvariants();
}

} // namespace solsim::physics
