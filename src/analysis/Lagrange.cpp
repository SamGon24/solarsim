// =============================================================================
//  SolarSim - Lagrange.cpp
// =============================================================================
#include "Lagrange.hpp"

#include "../astro/OrbitalElements.hpp"
#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"

#include <cmath>
#include <stdexcept>

namespace solsim::analysis {

double collinearLagrangePoint(double mu, int which) {
    // Equilibrium along the x axis of the rotating frame (barycentre at 0,
    // primary at -mu, secondary at 1 - mu):
    //   f(x) = x - (1-mu)(x+mu)/|x+mu|^3 - mu(x-1+mu)/|x-1+mu|^3 = 0
    auto f = [mu](double x) {
        const double d1 = x + mu, d2 = x - 1.0 + mu;
        return x - (1.0 - mu) * d1 / std::pow(std::fabs(d1), 3.0) - mu * d2 / std::pow(std::fabs(d2), 3.0);
    };
    const double eps = 1e-12;
    math::RootResult r;
    switch (which) {
        case 1: r = math::brent(f, -mu + 1e-6, 1.0 - mu - eps * 1e3, 1e-15, 300); break;
        case 2: r = math::brent(f, 1.0 - mu + eps * 1e3, 2.5, 1e-15, 300); break;
        case 3: r = math::brent(f, -2.5, -mu - 1e-6, 1e-15, 300); break;
        default: throw std::invalid_argument("collinear point must be 1, 2 or 3");
    }
    if (!r.converged && !std::isfinite(r.x)) throw std::runtime_error("Lagrange point solve failed");
    return r.x;
}

LagrangeSolution computeLagrangePoints(const physics::NBodySystem& sys, int primary, int secondary) {
    if (primary < 0 || secondary < 0 || primary == secondary) throw std::invalid_argument("invalid body pair");
    const std::size_t p = static_cast<std::size_t>(primary), s = static_cast<std::size_t>(secondary);
    const double m1 = sys.gm()[p], m2 = sys.gm()[s];
    if (m1 <= 0.0 || m2 <= 0.0) throw std::invalid_argument("both bodies need mass for Lagrange points");

    LagrangeSolution sol;
    sol.primary = sys.info(p).name;
    sol.secondary = sys.info(s).name;
    sol.massRatio = m2 / (m1 + m2);
    const double mu = sol.massRatio;

    const Vec3 r = sys.pos(s) - sys.pos(p);
    const Vec3 v = sys.vel(s) - sys.vel(p);
    const double a = r.norm();
    sol.separation = a;
    const Vec3 xhat = r / a;
    const Vec3 zhat = r.cross(v).normalized();
    const Vec3 yhat = zhat.cross(xhat);
    const Vec3 bary = (sys.pos(p) * m1 + sys.pos(s) * m2) / (m1 + m2);
    const Vec3 vbary = (sys.vel(p) * m1 + sys.vel(s) * m2) / (m1 + m2);
    const Vec3 omega = r.cross(v) / (a * a); // instantaneous angular velocity

    sol.x[0] = collinearLagrangePoint(mu, 1);
    sol.x[1] = collinearLagrangePoint(mu, 2);
    sol.x[2] = collinearLagrangePoint(mu, 3);
    sol.x[3] = 0.5 - mu;
    sol.y[3] = std::sqrt(3.0) / 2.0;   // L4 leads the secondary
    sol.x[4] = 0.5 - mu;
    sol.y[4] = -std::sqrt(3.0) / 2.0;  // L5 trails

    for (int k = 0; k < 5; ++k) {
        const Vec3 off = (xhat * sol.x[static_cast<std::size_t>(k)] + yhat * sol.y[static_cast<std::size_t>(k)]) * a;
        sol.position[static_cast<std::size_t>(k)] = bary + off;
        sol.velocity[static_cast<std::size_t>(k)] = vbary + omega.cross(off);
        const double dx = sol.x[static_cast<std::size_t>(k)] - (1.0 - mu), dy = sol.y[static_cast<std::size_t>(k)];
        sol.distanceFromSecondary[static_cast<std::size_t>(k)] = std::sqrt(dx * dx + dy * dy) * a;
    }

    const astro::OrbitalElements el = astro::stateToElements(r, v, m1 + m2, 0.0);
    const double semi = el.isElliptic() ? el.a() : a;
    sol.hillRadius = semi * (1.0 - std::min(el.e, 0.99)) * std::cbrt(m2 / (3.0 * m1));
    sol.sphereOfInfluence = semi * std::pow(m2 / m1, 0.4);
    sol.triangularStable = mu < 0.5 * (1.0 - std::sqrt(23.0 / 27.0)); // Routh / Gascheau: 0.0385209

    // Jacobi constant at L1: C = x^2 + y^2 + 2(1-mu)/r1 + 2 mu/r2 (zero velocity)
    const double xl1 = sol.x[0];
    sol.jacobiL1 = xl1 * xl1 + 2.0 * (1.0 - mu) / std::fabs(xl1 + mu) + 2.0 * mu / std::fabs(xl1 - 1.0 + mu);
    return sol;
}

std::array<LibrationResult, 5> testLagrangeStability(const physics::Simulation& simIn, int primary, int secondary,
                                                     double years, double perturbationAU) {
    physics::Simulation sim = simIn;
    auto& sys = sim.system();
    const LagrangeSolution sol0 = computeLagrangePoints(sys, primary, secondary);
    std::array<int, 5> idx{};
    for (int k = 0; k < 5; ++k) {
        physics::BodyInfo info;
        info.name = "L" + std::to_string(k + 1) + "-probe";
        info.kind = astro::BodyKind::Particle;
        info.parent = primary;
        info.radius = 1e-9;
        const Vec3 kick{perturbationAU, 0.0, 0.0};
        idx[static_cast<std::size_t>(k)] =
            sys.addBody(info, 0.0, sol0.position[static_cast<std::size_t>(k)] + kick, sol0.velocity[static_cast<std::size_t>(k)]);
    }
    sim.systemModified();

    std::array<LibrationResult, 5> res{};
    for (int k = 0; k < 5; ++k) {
        res[static_cast<std::size_t>(k)].point = k + 1;
        res[static_cast<std::size_t>(k)].years = years;
    }
    const double total = years * constants::JULIAN_YEAR_D;
    const LagrangeSolution s0 = computeLagrangePoints(sys, primary, secondary);
    const double sampleDays = std::max(1.0, total / 4000.0);
    for (double t = 0.0; t < total; t += sampleDays) {
        sim.advance(sampleDays);
        const LagrangeSolution s = computeLagrangePoints(sys, primary, secondary);
        const Vec3 pp = sys.pos(static_cast<std::size_t>(primary));
        for (int k = 0; k < 5; ++k) {
            auto& r = res[static_cast<std::size_t>(k)];
            const Vec3 probe = sys.pos(static_cast<std::size_t>(idx[static_cast<std::size_t>(k)]));
            const double dev = (probe - s.position[static_cast<std::size_t>(k)]).norm();
            r.maxDeviationAU = std::max(r.maxDeviationAU, dev);
            r.finalDeviationAU = dev;
            const double ang = angleBetween(probe - pp, s.position[static_cast<std::size_t>(k)] - pp) * constants::RAD2DEG;
            r.maxAngleDeg = std::max(r.maxAngleDeg, ang);
            if (dev > 0.5 * s0.separation) r.escaped = true;
        }
    }
    return res;
}

} // namespace solsim::analysis
