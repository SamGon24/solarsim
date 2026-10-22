// =============================================================================
//  SolarSim - Integrators.hpp
//  Numerical integrators for the second-order system  x'' = a(x, x').
//
//  name               order  symplectic  notes
//  ----------------   -----  ----------  ---------------------------------------
//  euler                1       no       explicit Euler - energy drifts quickly
//  symplectic-euler     1       yes      semi-implicit (kick then drift)
//  leapfrog             2       yes      Stoermer-Verlet drift-kick-drift
//  rk4                  4       no       classical Runge-Kutta, secular drift
//  yoshida4             4       yes      Yoshida (1990) triple-jump composition
//  dopri5             5(4)      no       Dormand-Prince with adaptive step size
//  wh                   2       yes      Wisdom-Holman map in democratic
//                                        heliocentric coordinates (Duncan,
//                                        Levison & Lee 1998) - Kepler drift is
//                                        solved exactly with universal variables
//
//  Symplectic integrators conserve a "shadow" Hamiltonian, so the energy error
//  stays bounded and oscillates instead of growing - the key property for
//  long-term Solar System integrations.
// =============================================================================
#pragma once

#include "Forces.hpp"
#include "NBodySystem.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace solsim::physics {

class Integrator {
public:
    virtual ~Integrator() = default;

    virtual std::string name() const = 0;
    virtual std::string displayName() const = 0;
    virtual int order() const = 0;
    virtual bool isSymplectic() const = 0;
    virtual bool isAdaptive() const { return false; }
    virtual int forceEvaluationsPerStep() const = 0;

    // Advances the system by at most `dt` days (negative = backwards in time).
    // Fixed-step methods always take exactly dt; adaptive methods may take a
    // smaller step and return the step actually taken.
    virtual double step(NBodySystem& sys, const ForceSet& forces, double dt) = 0;

    // Forget cached state (step-size history etc.).
    virtual void reset() {}
    virtual std::unique_ptr<Integrator> clone() const = 0;

    // Adaptive tolerance (relative); ignored by fixed-step methods.
    virtual void setTolerance(double /*rtol*/) {}
    virtual double tolerance() const { return 0.0; }
    // Step size suggested by an adaptive method (0 when not applicable).
    virtual double suggestedStep() const { return 0.0; }

    std::uint64_t steps() const { return m_steps; }
    std::uint64_t rejectedSteps() const { return m_rejected; }

protected:
    void ensureScratch(std::size_t n, std::size_t buffers);
    std::vector<std::vector<Vec3>> m_scratch;
    std::uint64_t m_steps = 0;
    std::uint64_t m_rejected = 0;
};

// ---- concrete integrators ---------------------------------------------------

class EulerIntegrator final : public Integrator {
public:
    std::string name() const override { return "euler"; }
    std::string displayName() const override { return "Explicit Euler"; }
    int order() const override { return 1; }
    bool isSymplectic() const override { return false; }
    int forceEvaluationsPerStep() const override { return 1; }
    double step(NBodySystem& sys, const ForceSet& forces, double dt) override;
    std::unique_ptr<Integrator> clone() const override { return std::make_unique<EulerIntegrator>(*this); }
};

class SymplecticEulerIntegrator final : public Integrator {
public:
    std::string name() const override { return "symplectic-euler"; }
    std::string displayName() const override { return "Symplectic Euler"; }
    int order() const override { return 1; }
    bool isSymplectic() const override { return true; }
    int forceEvaluationsPerStep() const override { return 1; }
    double step(NBodySystem& sys, const ForceSet& forces, double dt) override;
    std::unique_ptr<Integrator> clone() const override { return std::make_unique<SymplecticEulerIntegrator>(*this); }
};

class LeapfrogIntegrator final : public Integrator {
public:
    std::string name() const override { return "leapfrog"; }
    std::string displayName() const override { return "Leapfrog (Stoermer-Verlet DKD)"; }
    int order() const override { return 2; }
    bool isSymplectic() const override { return true; }
    int forceEvaluationsPerStep() const override { return 1; }
    double step(NBodySystem& sys, const ForceSet& forces, double dt) override;
    std::unique_ptr<Integrator> clone() const override { return std::make_unique<LeapfrogIntegrator>(*this); }
};

class RK4Integrator final : public Integrator {
public:
    std::string name() const override { return "rk4"; }
    std::string displayName() const override { return "Runge-Kutta 4"; }
    int order() const override { return 4; }
    bool isSymplectic() const override { return false; }
    int forceEvaluationsPerStep() const override { return 4; }
    double step(NBodySystem& sys, const ForceSet& forces, double dt) override;
    std::unique_ptr<Integrator> clone() const override { return std::make_unique<RK4Integrator>(*this); }
};

class Yoshida4Integrator final : public Integrator {
public:
    std::string name() const override { return "yoshida4"; }
    std::string displayName() const override { return "Yoshida 4th order"; }
    int order() const override { return 4; }
    bool isSymplectic() const override { return true; }
    int forceEvaluationsPerStep() const override { return 3; }
    double step(NBodySystem& sys, const ForceSet& forces, double dt) override;
    std::unique_ptr<Integrator> clone() const override { return std::make_unique<Yoshida4Integrator>(*this); }
};

class DormandPrinceIntegrator final : public Integrator {
public:
    explicit DormandPrinceIntegrator(double rtol = 1e-10) : m_rtol(rtol) {}
    std::string name() const override { return "dopri5"; }
    std::string displayName() const override { return "Dormand-Prince 5(4) adaptive"; }
    int order() const override { return 5; }
    bool isSymplectic() const override { return false; }
    bool isAdaptive() const override { return true; }
    int forceEvaluationsPerStep() const override { return 7; }
    double step(NBodySystem& sys, const ForceSet& forces, double dt) override;
    void reset() override { m_h = 0.0; }
    std::unique_ptr<Integrator> clone() const override { return std::make_unique<DormandPrinceIntegrator>(*this); }
    void setTolerance(double rtol) override { m_rtol = rtol; }
    double tolerance() const override { return m_rtol; }
    double suggestedStep() const override { return m_h; }

private:
    double m_rtol;
    double m_h = 0.0; // last successful |step| suggestion
};

class WisdomHolmanIntegrator final : public Integrator {
public:
    std::string name() const override { return "wh"; }
    std::string displayName() const override { return "Wisdom-Holman (democratic heliocentric)"; }
    int order() const override { return 2; }
    bool isSymplectic() const override { return true; }
    int forceEvaluationsPerStep() const override { return 2; }
    double step(NBodySystem& sys, const ForceSet& forces, double dt) override;
    std::unique_ptr<Integrator> clone() const override { return std::make_unique<WisdomHolmanIntegrator>(*this); }
    std::uint64_t keplerFailures() const { return m_keplerFailures; }

private:
    std::uint64_t m_keplerFailures = 0;
};

// ---- factory ------------------------------------------------------------------
std::unique_ptr<Integrator> createIntegrator(const std::string& name);
std::vector<std::string> integratorNames();

} // namespace solsim::physics
