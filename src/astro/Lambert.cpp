// =============================================================================
//  SolarSim - Lambert.cpp
// =============================================================================
#include "Lambert.hpp"

#include "Kepler.hpp"

#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"

#include <cmath>

namespace solsim::astro {

using namespace solsim::constants;

LambertSolution solveLambert(const Vec3& r1, const Vec3& r2, double tof, double mu, bool prograde) {
    LambertSolution sol;
    const double r1n = r1.norm(), r2n = r2.norm();
    if (r1n <= 0 || r2n <= 0 || tof <= 0 || mu <= 0) return sol;

    double cosDnu = math::clamp(r1.dot(r2) / (r1n * r2n), -1.0, 1.0);
    double dnu = std::acos(cosDnu);
    const double crossZ = r1.cross(r2).z;
    // Short way for prograde motion when the transfer sweeps counter-clockwise.
    if ((prograde && crossZ < 0.0) || (!prograde && crossZ >= 0.0)) dnu = TWO_PI - dnu;
    sol.transferAngleDeg = dnu * RAD2DEG;

    const double A = std::sin(dnu) * std::sqrt(r1n * r2n / (1.0 - std::cos(dnu)));
    if (std::fabs(A) < 1e-14) return sol; // 0 or 180 deg: plane undefined

    double psi = 0.0, c2 = 0.5, c3 = 1.0 / 6.0;
    double psiUp = 4.0 * PI * PI, psiLow = -4.0 * PI;
    const double sqrtMu = std::sqrt(mu);
    double y = 0.0;
    for (int it = 0; it < 200; ++it) {
        y = r1n + r2n + A * (psi * c3 - 1.0) / std::sqrt(c2);
        if (A > 0.0 && y < 0.0) {
            // Raise psi_low until y becomes positive.
            int guard = 0;
            while (y < 0.0 && guard++ < 100) {
                psi += 0.1;
                stumpff(psi, c2, c3);
                y = r1n + r2n + A * (psi * c3 - 1.0) / std::sqrt(c2);
            }
            psiLow = psi;
        }
        const double chi = std::sqrt(y / c2);
        const double dt = (chi * chi * chi * c3 + A * std::sqrt(y)) / sqrtMu;
        sol.iterations = it + 1;
        if (std::fabs(dt - tof) < 1e-9 * std::max(1.0, tof)) {
            sol.ok = true;
            break;
        }
        if (dt <= tof) psiLow = psi;
        else psiUp = psi;
        psi = 0.5 * (psiUp + psiLow);
        stumpff(psi, c2, c3);
        if (psiUp - psiLow < 1e-14) {
            sol.ok = true;
            break;
        }
    }
    if (!sol.ok || y <= 0.0) {
        sol.ok = false;
        return sol;
    }
    const double f = 1.0 - y / r1n;
    const double g = A * std::sqrt(y / mu);
    const double gdot = 1.0 - y / r2n;
    sol.v1 = (r2 - r1 * f) / g;
    sol.v2 = (r2 * gdot - r1) / g;
    return sol;
}

} // namespace solsim::astro
