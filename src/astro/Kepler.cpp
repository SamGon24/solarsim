// =============================================================================
//  SolarSim - Kepler.cpp
// =============================================================================
#include "Kepler.hpp"

#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"

#include <cmath>

namespace solsim::astro {

using namespace solsim::constants;

double solveKeplerElliptic(double M, double e, KeplerSolveInfo* info) {
    // Reduce M to (-pi, pi] - Kepler's equation is periodic in M.
    const double Mr = math::wrapPi(M);
    const double offset = M - Mr;
    // Starting value (Danby 1988): E0 = M + 0.85 e sign(sin M) is robust for all e < 1.
    double E = (e < 0.8) ? Mr : Mr + 0.85 * e * math::sign(std::sin(Mr));
    KeplerSolveInfo local;
    for (int i = 0; i < 64; ++i) {
        // Halley's method (cubic convergence).
        const double sE = std::sin(E), cE = std::cos(E);
        const double f = E - e * sE - Mr;
        const double fp = 1.0 - e * cE;
        const double fpp = e * sE;
        const double delta = f / (fp - 0.5 * f * fpp / fp);
        E -= delta;
        local.iterations = i + 1;
        if (std::fabs(delta) < 1e-15) {
            local.converged = true;
            break;
        }
    }
    if (info) *info = local;
    return E + offset;
}

double solveKeplerHyperbolic(double M, double e, KeplerSolveInfo* info) {
    // Initial guess (Vallado): H0 = asinh(M/e) is good for large |M|.
    double H = std::asinh(M / e);
    if (std::fabs(M) < 1e-3) H = M / (e - 1.0);
    KeplerSolveInfo local;
    for (int i = 0; i < 100; ++i) {
        const double sH = std::sinh(H), cH = std::cosh(H);
        const double f = e * sH - H - M;
        const double fp = e * cH - 1.0;
        const double fpp = e * sH;
        double delta = f / (fp - 0.5 * f * fpp / fp);
        // Damp huge steps (can happen for |M| very large with poor H0).
        if (std::fabs(delta) > 1.0) delta = math::sign(delta);
        H -= delta;
        local.iterations = i + 1;
        if (std::fabs(delta) < 1e-14 * std::max(1.0, std::fabs(H))) {
            local.converged = true;
            break;
        }
    }
    if (info) *info = local;
    return H;
}

double solveBarker(double W) {
    // Closed-form solution of Barker's equation:  tan(nu/2) = Y - 1/Y,
    // Y = cbrt(W/2 + sqrt(W^2/4 + 1))   (Meeus ch. 34)
    const double Y = std::cbrt(0.5 * W + std::sqrt(0.25 * W * W + 1.0));
    const double s = Y - 1.0 / Y;
    return 2.0 * std::atan(s);
}

double trueFromEccentric(double E, double e) {
    return 2.0 * std::atan2(std::sqrt(1.0 + e) * std::sin(0.5 * E), std::sqrt(1.0 - e) * std::cos(0.5 * E));
}

double eccentricFromTrue(double nu, double e) {
    return 2.0 * std::atan2(std::sqrt(1.0 - e) * std::sin(0.5 * nu), std::sqrt(1.0 + e) * std::cos(0.5 * nu));
}

double meanFromEccentric(double E, double e) { return E - e * std::sin(E); }

double trueFromHyperbolic(double H, double e) {
    return 2.0 * std::atan(std::sqrt((e + 1.0) / (e - 1.0)) * std::tanh(0.5 * H));
}

double hyperbolicFromTrue(double nu, double e) {
    return 2.0 * std::atanh(std::sqrt((e - 1.0) / (e + 1.0)) * std::tan(0.5 * nu));
}

double meanFromHyperbolic(double H, double e) { return e * std::sinh(H) - H; }

void stumpff(double psi, double& c2, double& c3) {
    if (psi > 1e-6) {
        const double sp = std::sqrt(psi);
        c2 = (1.0 - std::cos(sp)) / psi;
        c3 = (sp - std::sin(sp)) / (psi * sp);
    } else if (psi < -1e-6) {
        const double sp = std::sqrt(-psi);
        c2 = (1.0 - std::cosh(sp)) / psi;
        c3 = (std::sinh(sp) - sp) / (-psi * sp);
    } else {
        // Taylor series: c2 = 1/2 - psi/24 + psi^2/720 ; c3 = 1/6 - psi/120 + psi^2/5040
        c2 = 0.5 - psi / 24.0 + psi * psi / 720.0;
        c3 = 1.0 / 6.0 - psi / 120.0 + psi * psi / 5040.0;
    }
}

UniversalResult propagateKepler(const Vec3& r0v, const Vec3& v0v, double mu, double dtIn) {
    UniversalResult out;
    out.r = r0v;
    out.v = v0v;
    if (dtIn == 0.0) {
        out.converged = true;
        return out;
    }
    const double sqrtMu = std::sqrt(mu);
    const double r0 = r0v.norm();
    const double v0sq = v0v.norm2();
    const double rdotv = r0v.dot(v0v);
    const double alpha = 2.0 / r0 - v0sq / mu; // = 1/a

    double dt = dtIn;
    // For ellipses remove whole periods: improves conditioning for long arcs.
    if (alpha > 1e-12) {
        const double period = TWO_PI / (sqrtMu * alpha * std::sqrt(alpha));
        if (std::fabs(dt) > period) dt = std::fmod(dt, period);
    }

    // Initial guess for the universal anomaly chi.
    double chi;
    if (alpha > 1e-12) {
        chi = sqrtMu * dt * alpha;
    } else if (alpha < -1e-12) {
        const double a = 1.0 / alpha;
        const double s = math::sign(dt);
        const double num = -2.0 * mu * alpha * dt;
        const double den = rdotv + s * std::sqrt(-mu * a) * (1.0 - r0 * alpha);
        chi = s * std::sqrt(-a) * std::log(std::fabs(num / den));
        if (!std::isfinite(chi)) chi = sqrtMu * dt / r0;
    } else {
        // Near-parabolic: use the cubic (Barker-like) estimate.
        const Vec3 h = r0v.cross(v0v);
        const double p = h.norm2() / mu;
        const double s = 0.5 * std::atan(1.0 / (3.0 * std::sqrt(mu / (p * p * p)) * dt));
        const double w = std::atan(std::cbrt(std::tan(s)));
        chi = std::sqrt(p) * 2.0 / std::tan(2.0 * w);
    }

    double c2 = 0.5, c3 = 1.0 / 6.0, psi = 0.0, r = r0;
    bool converged = false;
    int it = 0;
    for (; it < 100; ++it) {
        psi = chi * chi * alpha;
        stumpff(psi, c2, c3);
        r = chi * chi * c2 + rdotv / sqrtMu * chi * (1.0 - psi * c3) + r0 * (1.0 - psi * c2);
        const double f = sqrtMu * dt - chi * chi * chi * c3 - rdotv / sqrtMu * chi * chi * c2 - r0 * chi * (1.0 - psi * c3);
        double delta = f / r;
        // Guard against wild steps in the first iterations.
        const double limit = 0.5 * std::max(1.0, std::fabs(chi));
        if (std::fabs(delta) > limit) delta = math::sign(delta) * limit;
        chi += delta;
        if (std::fabs(delta) < 1e-13 * std::max(1.0, std::fabs(chi))) {
            converged = true;
            ++it;
            break;
        }
    }
    psi = chi * chi * alpha;
    stumpff(psi, c2, c3);
    r = chi * chi * c2 + rdotv / sqrtMu * chi * (1.0 - psi * c3) + r0 * (1.0 - psi * c2);

    const double f = 1.0 - chi * chi / r0 * c2;
    const double g = dt - chi * chi * chi / sqrtMu * c3;
    const double gdot = 1.0 - chi * chi / r * c2;
    const double fdot = sqrtMu / (r * r0) * chi * (psi * c3 - 1.0);

    out.r = r0v * f + v0v * g;
    out.v = r0v * fdot + v0v * gdot;
    out.iterations = it;
    out.converged = converged;
    return out;
}

} // namespace solsim::astro
