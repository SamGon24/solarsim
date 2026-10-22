// =============================================================================
//  SolarSim - Integrators.cpp
// =============================================================================
#include "Integrators.hpp"

#include "../astro/Kepler.hpp"
#include "../core/StringUtil.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <stdexcept>

namespace solsim::physics {

void Integrator::ensureScratch(std::size_t n, std::size_t buffers) {
    if (m_scratch.size() < buffers) m_scratch.resize(buffers);
    for (std::size_t b = 0; b < buffers; ++b)
        if (m_scratch[b].size() != n) m_scratch[b].assign(n, Vec3{});
}

// -----------------------------------------------------------------------------
// Explicit Euler:  x1 = x0 + h v0 ;  v1 = v0 + h a(x0, v0)
// -----------------------------------------------------------------------------
double EulerIntegrator::step(NBodySystem& sys, const ForceSet& forces, double h) {
    auto& x = sys.positions();
    auto& v = sys.velocities();
    auto& a = sys.accelerations();
    forces.compute(sys, x.data(), v.data(), a.data());
    for (std::size_t i = 0; i < sys.size(); ++i) {
        x[i] += v[i] * h;
        v[i] += a[i] * h;
    }
    ++m_steps;
    return h;
}

// -----------------------------------------------------------------------------
// Symplectic Euler (kick-drift):  v1 = v0 + h a(x0) ;  x1 = x0 + h v1
// -----------------------------------------------------------------------------
double SymplecticEulerIntegrator::step(NBodySystem& sys, const ForceSet& forces, double h) {
    auto& x = sys.positions();
    auto& v = sys.velocities();
    auto& a = sys.accelerations();
    forces.compute(sys, x.data(), v.data(), a.data());
    for (std::size_t i = 0; i < sys.size(); ++i) {
        v[i] += a[i] * h;
        x[i] += v[i] * h;
    }
    ++m_steps;
    return h;
}

// -----------------------------------------------------------------------------
// Leapfrog, drift-kick-drift form (time reversible, symplectic, 2nd order)
// -----------------------------------------------------------------------------
double LeapfrogIntegrator::step(NBodySystem& sys, const ForceSet& forces, double h) {
    auto& x = sys.positions();
    auto& v = sys.velocities();
    auto& a = sys.accelerations();
    const std::size_t n = sys.size();
    const double hh = 0.5 * h;
    for (std::size_t i = 0; i < n; ++i) x[i] += v[i] * hh;
    forces.compute(sys, x.data(), v.data(), a.data());
    for (std::size_t i = 0; i < n; ++i) {
        v[i] += a[i] * h;
        x[i] += v[i] * hh;
    }
    ++m_steps;
    return h;
}

// -----------------------------------------------------------------------------
// Classical RK4 on y = (x, v), y' = (v, a(x, v))
// -----------------------------------------------------------------------------
double RK4Integrator::step(NBodySystem& sys, const ForceSet& forces, double h) {
    const std::size_t n = sys.size();
    ensureScratch(n, 10);
    auto& x = sys.positions();
    auto& v = sys.velocities();
    auto& xt = m_scratch[0];
    auto& vt = m_scratch[1];
    auto& k1x = m_scratch[2]; auto& k1v = m_scratch[3];
    auto& k2x = m_scratch[4]; auto& k2v = m_scratch[5];
    auto& k3x = m_scratch[6]; auto& k3v = m_scratch[7];
    auto& k4x = m_scratch[8]; auto& k4v = m_scratch[9];

    // k1
    forces.compute(sys, x.data(), v.data(), k1v.data());
    for (std::size_t i = 0; i < n; ++i) k1x[i] = v[i];
    // k2
    for (std::size_t i = 0; i < n; ++i) {
        xt[i] = x[i] + k1x[i] * (0.5 * h);
        vt[i] = v[i] + k1v[i] * (0.5 * h);
    }
    forces.compute(sys, xt.data(), vt.data(), k2v.data());
    for (std::size_t i = 0; i < n; ++i) k2x[i] = vt[i];
    // k3
    for (std::size_t i = 0; i < n; ++i) {
        xt[i] = x[i] + k2x[i] * (0.5 * h);
        vt[i] = v[i] + k2v[i] * (0.5 * h);
    }
    forces.compute(sys, xt.data(), vt.data(), k3v.data());
    for (std::size_t i = 0; i < n; ++i) k3x[i] = vt[i];
    // k4
    for (std::size_t i = 0; i < n; ++i) {
        xt[i] = x[i] + k3x[i] * h;
        vt[i] = v[i] + k3v[i] * h;
    }
    forces.compute(sys, xt.data(), vt.data(), k4v.data());
    for (std::size_t i = 0; i < n; ++i) k4x[i] = vt[i];

    const double s = h / 6.0;
    for (std::size_t i = 0; i < n; ++i) {
        x[i] += (k1x[i] + k2x[i] * 2.0 + k3x[i] * 2.0 + k4x[i]) * s;
        v[i] += (k1v[i] + k2v[i] * 2.0 + k3v[i] * 2.0 + k4v[i]) * s;
        sys.accelerations()[i] = k1v[i];
    }
    ++m_steps;
    return h;
}

// -----------------------------------------------------------------------------
// Yoshida 4th order symplectic composition (Yoshida 1990, Phys. Lett. A 150)
//   w1 = 1 / (2 - 2^(1/3)),  w0 = -2^(1/3) / (2 - 2^(1/3))
//   drift c1 kick d1 drift c2 kick d2 drift c3 kick d3 drift c4
// -----------------------------------------------------------------------------
double Yoshida4Integrator::step(NBodySystem& sys, const ForceSet& forces, double h) {
    static const double cbrt2 = std::cbrt(2.0);
    static const double w1 = 1.0 / (2.0 - cbrt2);
    static const double w0 = -cbrt2 / (2.0 - cbrt2);
    static const double c[4] = {0.5 * w1, 0.5 * (w0 + w1), 0.5 * (w0 + w1), 0.5 * w1};
    static const double d[3] = {w1, w0, w1};

    auto& x = sys.positions();
    auto& v = sys.velocities();
    auto& a = sys.accelerations();
    const std::size_t n = sys.size();
    for (int stage = 0; stage < 3; ++stage) {
        const double ch = c[stage] * h;
        for (std::size_t i = 0; i < n; ++i) x[i] += v[i] * ch;
        forces.compute(sys, x.data(), v.data(), a.data());
        const double dh = d[stage] * h;
        for (std::size_t i = 0; i < n; ++i) v[i] += a[i] * dh;
    }
    const double ch = c[3] * h;
    for (std::size_t i = 0; i < n; ++i) x[i] += v[i] * ch;
    ++m_steps;
    return h;
}

// -----------------------------------------------------------------------------
// Dormand-Prince 5(4) - embedded Runge-Kutta pair with local error control.
// -----------------------------------------------------------------------------
namespace dp {
constexpr double c2 = 1.0 / 5, c3 = 3.0 / 10, c4 = 4.0 / 5, c5 = 8.0 / 9;
constexpr double a21 = 1.0 / 5;
constexpr double a31 = 3.0 / 40, a32 = 9.0 / 40;
constexpr double a41 = 44.0 / 45, a42 = -56.0 / 15, a43 = 32.0 / 9;
constexpr double a51 = 19372.0 / 6561, a52 = -25360.0 / 2187, a53 = 64448.0 / 6561, a54 = -212.0 / 729;
constexpr double a61 = 9017.0 / 3168, a62 = -355.0 / 33, a63 = 46732.0 / 5247, a64 = 49.0 / 176,
                 a65 = -5103.0 / 18656;
constexpr double b1 = 35.0 / 384, b3 = 500.0 / 1113, b4 = 125.0 / 192, b5 = -2187.0 / 6784, b6 = 11.0 / 84;
// error coefficients e_i = b_i - b*_i
constexpr double e1 = 71.0 / 57600, e3 = -71.0 / 16695, e4 = 71.0 / 1920, e5 = -17253.0 / 339200,
                 e6 = 22.0 / 525, e7 = -1.0 / 40;
} // namespace dp

double DormandPrinceIntegrator::step(NBodySystem& sys, const ForceSet& forces, double dtMax) {
    using namespace dp;
    const std::size_t n = sys.size();
    ensureScratch(n, 18);
    auto& x = sys.positions();
    auto& v = sys.velocities();
    auto* kx = &m_scratch[0]; // kx[0..6]
    auto* kv = &m_scratch[7]; // kv[0..6]
    auto& xt = m_scratch[14];
    auto& vt = m_scratch[15];
    auto& xn = m_scratch[16];
    auto& vn = m_scratch[17];

    const double dir = dtMax < 0 ? -1.0 : 1.0;
    const double absMax = std::fabs(dtMax);
    if (absMax == 0.0) return 0.0;
    if (m_h <= 0.0) m_h = std::min(absMax, 0.05);

    const double atolX = m_rtol * 1e-2;  // AU
    const double atolV = m_rtol * 1e-4;  // AU/day

    // k1 at the start of the step (does not depend on h).
    forces.compute(sys, x.data(), v.data(), kv[0].data());
    for (std::size_t i = 0; i < n; ++i) kx[0][i] = v[i];

    for (int attempt = 0; attempt < 60; ++attempt) {
        const double absH = std::min(m_h, absMax);
        const bool clipped = absH < m_h;
        const double h = dir * absH;

        auto stage = [&](int s, std::initializer_list<double> coeffs) {
            for (std::size_t i = 0; i < n; ++i) {
                Vec3 dx{}, dv{};
                int j = 0;
                for (double cf : coeffs) {
                    if (cf != 0.0) {
                        dx += kx[j][i] * cf;
                        dv += kv[j][i] * cf;
                    }
                    ++j;
                }
                xt[i] = x[i] + dx * h;
                vt[i] = v[i] + dv * h;
            }
            forces.compute(sys, xt.data(), vt.data(), kv[s].data());
            for (std::size_t i = 0; i < n; ++i) kx[s][i] = vt[i];
        };
        (void)c2; (void)c3; (void)c4; (void)c5; // time-independent forces: c_i unused
        stage(1, {a21});
        stage(2, {a31, a32});
        stage(3, {a41, a42, a43});
        stage(4, {a51, a52, a53, a54});
        stage(5, {a61, a62, a63, a64, a65});
        // 5th-order solution (also the 7th stage abscissa: FSAL property)
        for (std::size_t i = 0; i < n; ++i) {
            xn[i] = x[i] + (kx[0][i] * b1 + kx[2][i] * b3 + kx[3][i] * b4 + kx[4][i] * b5 + kx[5][i] * b6) * h;
            vn[i] = v[i] + (kv[0][i] * b1 + kv[2][i] * b3 + kv[3][i] * b4 + kv[4][i] * b5 + kv[5][i] * b6) * h;
        }
        forces.compute(sys, xn.data(), vn.data(), kv[6].data());
        for (std::size_t i = 0; i < n; ++i) kx[6][i] = vn[i];

        // Error estimate (max-norm over bodies, scaled per body).
        double err = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const Vec3 ex = (kx[0][i] * e1 + kx[2][i] * e3 + kx[3][i] * e4 + kx[4][i] * e5 + kx[5][i] * e6 +
                             kx[6][i] * e7) * h;
            const Vec3 ev = (kv[0][i] * e1 + kv[2][i] * e3 + kv[3][i] * e4 + kv[4][i] * e5 + kv[5][i] * e6 +
                             kv[6][i] * e7) * h;
            const double sx = atolX + m_rtol * std::max(x[i].norm(), xn[i].norm());
            const double sv = atolV + m_rtol * std::max(v[i].norm(), vn[i].norm());
            err = std::max(err, std::max(ex.norm() / sx, ev.norm() / sv));
        }
        if (!std::isfinite(err)) err = 1e10;

        if (err <= 1.0) {
            for (std::size_t i = 0; i < n; ++i) {
                x[i] = xn[i];
                v[i] = vn[i];
                sys.accelerations()[i] = kv[6][i];
            }
            const double factor = std::min(5.0, std::max(0.2, 0.9 * std::pow(std::max(err, 1e-10), -0.2)));
            const double next = absH * factor;
            m_h = clipped ? std::max(m_h, next) : next;
            ++m_steps;
            return h;
        }
        ++m_rejected;
        m_h = absH * std::max(0.1, 0.9 * std::pow(err, -0.2));
        if (m_h < 1e-12) throw std::runtime_error("dopri5: step size underflow (close encounter?)");
    }
    throw std::runtime_error("dopri5: too many rejected steps");
}

// -----------------------------------------------------------------------------
// Wisdom-Holman, democratic heliocentric coordinates:
//   H = H_Kepler(Q, P) + H_interaction(Q) + H_sun(P)
//   step:  sun-drift(h/2)  kick(h/2)  Kepler(h)  kick(h/2)  sun-drift(h/2)
// Q_i = heliocentric positions, P_i/m_i = barycentric velocities u_i.
// -----------------------------------------------------------------------------
double WisdomHolmanIntegrator::step(NBodySystem& sys, const ForceSet& forces, double h) {
    const std::size_t n = sys.size();
    const std::size_t nM = sys.massiveCount();
    if (n < 2 || nM < 1) return h;
    ensureScratch(n, 3);
    auto& Q = m_scratch[0];
    auto& u = m_scratch[1];
    auto& aInt = m_scratch[2];
    auto& x = sys.positions();
    auto& v = sys.velocities();
    const auto& gm = sys.gm();
    const double m0 = gm[0];
    const double mTot = sys.totalGM();

    // Interaction gravity = Newtonian gravity without the central star as a source.
    NewtonianGravity interaction(0.0, nullptr);
    if (const NewtonianGravity* g = forces.gravity()) {
        interaction = *g;
    }
    interaction.setSkipSource(0);

    const Vec3 xcm = sys.centerOfMass();
    const Vec3 vcm = sys.centerOfMassVelocity();
    for (std::size_t i = 0; i < n; ++i) {
        Q[i] = x[i] - x[0];
        u[i] = v[i] - vcm;
    }
    Q[0] = Vec3{};

    auto momentumOverM0 = [&]() {
        Vec3 P{};
        for (std::size_t i = 1; i < nM; ++i) P += u[i] * gm[i];
        return P / m0;
    };
    auto sunDrift = [&](double tau) {
        const Vec3 shift = momentumOverM0() * tau;
        for (std::size_t i = 1; i < n; ++i) Q[i] += shift;
    };
    auto kick = [&](double tau) {
        for (std::size_t i = 0; i < n; ++i) aInt[i] = Vec3{};
        interaction.accumulate(sys, Q.data(), u.data(), aInt.data());
        forces.countEvaluation();
        for (std::size_t i = 1; i < n; ++i) u[i] += aInt[i] * tau;
    };

    sunDrift(0.5 * h);
    kick(0.5 * h);
    for (std::size_t i = 1; i < n; ++i) {
        const astro::UniversalResult r = astro::propagateKepler(Q[i], u[i], m0, h);
        if (!r.converged) ++m_keplerFailures;
        Q[i] = r.r;
        u[i] = r.v;
    }
    kick(0.5 * h);
    sunDrift(0.5 * h);

    // Back to barycentric coordinates.
    const Vec3 xcmNew = xcm + vcm * h;
    Vec3 weighted{};
    for (std::size_t i = 1; i < nM; ++i) weighted += Q[i] * gm[i];
    const Vec3 x0 = xcmNew - weighted / mTot;
    x[0] = x0;
    v[0] = vcm - momentumOverM0();
    for (std::size_t i = 1; i < n; ++i) {
        x[i] = Q[i] + x0;
        v[i] = u[i] + vcm;
    }
    ++m_steps;
    return h;
}

// -----------------------------------------------------------------------------
std::unique_ptr<Integrator> createIntegrator(const std::string& nameIn) {
    const std::string name = str::toLower(nameIn);
    if (name == "euler") return std::make_unique<EulerIntegrator>();
    if (name == "symplectic-euler" || name == "seuler" || name == "symplectic")
        return std::make_unique<SymplecticEulerIntegrator>();
    if (name == "leapfrog" || name == "verlet") return std::make_unique<LeapfrogIntegrator>();
    if (name == "rk4") return std::make_unique<RK4Integrator>();
    if (name == "yoshida4" || name == "yoshida") return std::make_unique<Yoshida4Integrator>();
    if (name == "dopri5" || name == "rk45" || name == "adaptive") return std::make_unique<DormandPrinceIntegrator>();
    if (name == "wh" || name == "wisdom-holman") return std::make_unique<WisdomHolmanIntegrator>();
    return nullptr;
}

std::vector<std::string> integratorNames() {
    return {"euler", "symplectic-euler", "leapfrog", "rk4", "yoshida4", "dopri5", "wh"};
}

} // namespace solsim::physics
