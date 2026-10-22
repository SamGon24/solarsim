// =============================================================================
//  SolarSim - Transfer.cpp
// =============================================================================
#include "Transfer.hpp"

#include "../astro/Coordinates.hpp"
#include "../astro/JulianDate.hpp"
#include "../astro/Lambert.hpp"
#include "../astro/OrbitalElements.hpp"
#include "../core/Constants.hpp"
#include "../core/Mat3.hpp"
#include "../core/MathUtil.hpp"
#include "../core/StringUtil.hpp"

#include <cmath>
#include <limits>

namespace solsim::analysis {

using namespace solsim::constants;
using physics::Simulation;

namespace {

constexpr double KMS_PER_AUD = AU_KM / DAY_S; // 1 AU/day in km/s (~1731.46)

Vec3 helioPos(const physics::NBodySystem& s, std::size_t i) { return s.pos(i) - s.pos(0); }
Vec3 helioVel(const physics::NBodySystem& s, std::size_t i) { return s.vel(i) - s.vel(0); }

void sampleOrbit(const physics::NBodySystem& s, std::size_t i, std::vector<Vec3>& out) {
    const astro::OrbitalElements el = astro::stateToElements(helioPos(s, i), helioVel(s, i), s.twoBodyMu(i), 0.0);
    if (!el.isElliptic()) return;
    Vec3 P, Q, W;
    astro::perifocalBasis(el.i, el.node, el.argp, P, Q, W);
    const double a = el.a(), b = a * std::sqrt(1.0 - el.e * el.e);
    for (int k = 0; k <= 360; ++k) {
        const double E = TWO_PI * k / 360.0;
        out.push_back(P * (a * (std::cos(E) - el.e)) + Q * (b * std::sin(E)));
    }
}

// Daily table of heliocentric states with cubic Hermite interpolation.
class StateTable {
public:
    void build(Simulation sim, std::size_t a, std::size_t b, double start, double end) {
        m_start = start;
        sim.advanceTo(start);
        for (double t = start; t <= end + 1.0; t += 1.0) {
            if (t > start) sim.advanceTo(t);
            const auto& s = sim.system();
            m_ra.push_back(helioPos(s, a));
            m_va.push_back(helioVel(s, a));
            m_rb.push_back(helioPos(s, b));
            m_vb.push_back(helioVel(s, b));
        }
    }
    bool covers(double t) const { return t >= m_start && t <= m_start + double(m_ra.size()) - 2.0; }
    void stateA(double t, Vec3& r, Vec3& v) const { interp(m_ra, m_va, t, r, v); }
    void stateB(double t, Vec3& r, Vec3& v) const { interp(m_rb, m_vb, t, r, v); }

private:
    void interp(const std::vector<Vec3>& R, const std::vector<Vec3>& V, double t, Vec3& r, Vec3& v) const {
        const double x = t - m_start;
        std::size_t i = static_cast<std::size_t>(std::floor(x));
        if (i + 1 >= R.size()) i = R.size() - 2;
        const double s = x - double(i);
        const double h00 = 2 * s * s * s - 3 * s * s + 1, h10 = s * s * s - 2 * s * s + s;
        const double h01 = -2 * s * s * s + 3 * s * s, h11 = s * s * s - s * s;
        r = R[i] * h00 + V[i] * h10 + R[i + 1] * h01 + V[i + 1] * h11;
        const double d00 = 6 * s * s - 6 * s, d10 = 3 * s * s - 4 * s + 1, d01 = -6 * s * s + 6 * s, d11 = 3 * s * s - 2 * s;
        v = R[i] * d00 + V[i] * d10 + R[i + 1] * d01 + V[i + 1] * d11;
    }
    double m_start = 0.0;
    std::vector<Vec3> m_ra, m_va, m_rb, m_vb;
};

// Oberth-effect burn from a circular orbit of radius r (km) to hyperbolic excess vinf (km/s).
double burnFromCircular(double muKm, double rKm, double vinf) {
    return std::sqrt(vinf * vinf + 2.0 * muKm / rKm) - std::sqrt(muKm / rKm);
}

} // namespace

HohmannTransfer hohmann(double mu, double r1, double r2) {
    HohmannTransfer h;
    h.r1 = r1;
    h.r2 = r2;
    h.aTransfer = 0.5 * (r1 + r2);
    const double v1 = std::sqrt(mu / r1), v2 = std::sqrt(mu / r2);
    const double vp = std::sqrt(mu * (2.0 / r1 - 1.0 / h.aTransfer));
    const double va = std::sqrt(mu * (2.0 / r2 - 1.0 / h.aTransfer));
    h.dv1 = std::fabs(vp - v1) * KMS_PER_AUD;
    h.dv2 = std::fabs(v2 - va) * KMS_PER_AUD;
    h.timeOfFlight = PI * std::sqrt(h.aTransfer * h.aTransfer * h.aTransfer / mu);
    const double n2 = std::sqrt(mu / (r2 * r2 * r2)), n1 = std::sqrt(mu / (r1 * r1 * r1));
    h.phaseAngleDeg = math::wrapPi(PI - n2 * h.timeOfFlight) * RAD2DEG;
    h.synodicPeriod = TWO_PI / std::fabs(n1 - n2);
    return h;
}

MissionResult planMission(const Simulation& baseSim, const MissionConfig& cfg, const MissionLog& logIn) {
    MissionLog log = logIn ? logIn : [](const std::string&) {};
    MissionResult res;
    Simulation sim = baseSim;
    const auto& sys0 = sim.system();
    const int depI = sys0.indexOf(cfg.from), tgtI = sys0.indexOf(cfg.to);
    if (depI <= 0 || tgtI <= 0 || depI == tgtI) {
        res.message = "departure and target must be two different bodies orbiting the Sun";
        return res;
    }
    const std::size_t dep = static_cast<std::size_t>(depI), tgt = static_cast<std::size_t>(tgtI);
    if (sys0.info(dep).parent != 0 || sys0.info(tgt).parent != 0 || sys0.gm()[dep] <= 0 || sys0.gm()[tgt] <= 0) {
        res.message = "both bodies must be massive and orbit the Sun directly";
        return res;
    }
    const double muSun = sys0.gm()[0];
    const double muD = sys0.gm()[dep] * AU3D2_TO_KM3S2, muT = sys0.gm()[tgt] * AU3D2_TO_KM3S2;
    const double RdKm = sys0.info(dep).radius * AU_KM, RtKm = sys0.info(tgt).radius * AU_KM;
    const double rPark = RdKm + cfg.parkingAltitudeKm, rCapture = RtKm + cfg.captureAltitudeKm;

    // ---- 1. analytic Hohmann numbers from the osculating semi-major axes -------------------
    if (sim.time() < cfg.searchStartJD) sim.advanceTo(cfg.searchStartJD);
    const auto elDep = astro::stateToElements(helioPos(sim.system(), dep), helioVel(sim.system(), dep), sim.system().twoBodyMu(dep), sim.time());
    const auto elTgt = astro::stateToElements(helioPos(sim.system(), tgt), helioVel(sim.system(), tgt), sim.system().twoBodyMu(tgt), sim.time());
    res.analytic = hohmann(muSun, elDep.a(), elTgt.a());
    log(str::format("Hohmann (circular, coplanar): a_t = %.4f AU, dv1 = %.3f km/s, dv2 = %.3f km/s, TOF = %.1f d, "
                    "phase angle = %.1f deg, synodic period = %.0f d",
                    res.analytic.aTransfer, res.analytic.dv1, res.analytic.dv2, res.analytic.timeOfFlight,
                    res.analytic.phaseAngleDeg, res.analytic.synodicPeriod));

    // ---- 2. Hohmann launch window: required phase angle with the real positions ------------
    const double phaseReq = res.analytic.phaseAngleDeg * DEG2RAD;
    auto phaseError = [&](const physics::NBodySystem& s) {
        const Vec3 a = helioPos(s, dep), b = helioPos(s, tgt);
        return math::wrapPi(std::atan2(b.y, b.x) - std::atan2(a.y, a.x) - phaseReq);
    };
    const double searchSpan = std::min(6000.0, 1.2 * res.analytic.synodicPeriod + 30.0);
    double prevErr = phaseError(sim.system()), prevT = sim.time();
    Simulation prevSnap = sim;
    double window = -1.0;
    for (double t = sim.time() + 1.0; t <= cfg.searchStartJD + searchSpan; t += 1.0) {
        sim.advanceTo(t);
        const double err = phaseError(sim.system());
        if (((prevErr < 0 && err >= 0) || (prevErr > 0 && err <= 0)) && std::fabs(err - prevErr) < HALF_PI) {
            auto g = [&](double tt) {
                Simulation s = prevSnap;
                s.advanceTo(tt);
                return phaseError(s.system());
            };
            const math::RootResult r = math::brent(g, prevT, t, 1e-5, 60);
            window = r.converged ? r.x : t;
            break;
        }
        prevErr = err;
        prevT = t;
        prevSnap = sim;
    }
    if (window < 0) {
        res.message = "no launch window found in the search interval";
        return res;
    }
    res.windowJD = window;
    log("Hohmann phase-angle window: " + astro::formatUTC(window, false));

    // ---- 3. porkchop scan with Lambert arcs ------------------------------------------------------
    double launch = window, tof = res.analytic.timeOfFlight;
    Vec3 lambertV1;
    bool haveLambert = false;
    {
        const double half = cfg.porkchopHalfWidthDays;
        const double l0 = std::max(cfg.searchStartJD, window - half), l1 = window + half;
        const double tMin = 0.5 * res.analytic.timeOfFlight, tMax = 1.35 * res.analytic.timeOfFlight;
        StateTable table;
        table.build(baseSim, dep, tgt, l0, l1 + tMax + 2.0);
        PorkchopGrid& pc = res.porkchop;
        for (double L = l0; L <= l1 + 1e-9; L += cfg.porkchopStepDays) pc.launchJD.push_back(L);
        const int nT = 60;
        for (int k = 0; k < nT; ++k) pc.tofDays.push_back(tMin + (tMax - tMin) * k / (nT - 1));
        pc.c3.assign(pc.launchJD.size() * pc.tofDays.size(), std::numeric_limits<double>::quiet_NaN());
        pc.totalDv = pc.c3;
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t r = 0; r < pc.tofDays.size(); ++r) {
            for (std::size_t c = 0; c < pc.launchJD.size(); ++c) {
                const double L = pc.launchJD[c], T = pc.tofDays[r];
                if (!table.covers(L) || !table.covers(L + T)) continue;
                Vec3 r1, v1p, r2, v2p;
                table.stateA(L, r1, v1p);
                table.stateB(L + T, r2, v2p);
                const astro::LambertSolution ls = astro::solveLambert(r1, r2, T, muSun, true);
                if (!ls.ok) continue;
                const double vinfD = (ls.v1 - v1p).norm() * KMS_PER_AUD;
                const double vinfA = (ls.v2 - v2p).norm() * KMS_PER_AUD;
                const double c3 = vinfD * vinfD;
                const double dv = burnFromCircular(muD, rPark, vinfD) + burnFromCircular(muT, rCapture, vinfA);
                pc.c3[r * pc.launchJD.size() + c] = c3;
                pc.totalDv[r * pc.launchJD.size() + c] = dv;
                if (dv < best) {
                    best = dv;
                    pc.bestLaunch = L;
                    pc.bestTof = T;
                    pc.bestC3 = c3;
                    pc.bestDv = dv;
                    lambertV1 = ls.v1;
                    haveLambert = true;
                }
            }
        }
        if (cfg.porkchop && haveLambert) {
            launch = pc.bestLaunch;
            tof = pc.bestTof;
            log(str::format("porkchop optimum (%zux%zu Lambert arcs): launch %s, TOF %.1f d, C3 = %.2f km^2/s^2, total dv = %.3f km/s",
                            pc.launchJD.size(), pc.tofDays.size(), astro::formatUTC(launch, false).c_str(), tof, pc.bestC3, pc.bestDv));
        }
    }

    // ---- 4. state at launch ----------------------------------------------------------------------
    Simulation launchSnap = baseSim;
    launchSnap.advanceTo(launch);
    const auto& ls = launchSnap.system();
    res.launchJD = launch;
    res.timeOfFlight = tof;
    const double arrival = launch + tof;
    res.arrivalJD = arrival;
    const Vec3 rDep = helioPos(ls, dep), vDep = helioVel(ls, dep);
    res.launchPosition = rDep;

    Vec3 v0;
    if (cfg.porkchop && haveLambert) {
        v0 = lambertV1;
    } else {
        Simulation probe = launchSnap;
        probe.advance(tof);
        const double r1 = rDep.norm(), r2 = helioPos(probe.system(), tgt).norm();
        const double at = 0.5 * (r1 + r2);
        v0 = rDep.cross(vDep).cross(rDep).normalized() * std::sqrt(muSun * (2.0 / r1 - 1.0 / at));
    }

    // Spacecraft flight helper (patched conic: departure body's gravity is never applied).
    auto makeFlight = [&](const Vec3& vHelio, bool targetGravity, int& scIndex) {
        Simulation s = launchSnap;
        auto& sys = s.system();
        physics::BodyInfo info;
        info.name = "Spacecraft";
        info.kind = astro::BodyKind::Spacecraft;
        info.parent = 0;
        info.radius = 1e-9;
        info.color = 0x66FFCC;
        scIndex = sys.addBody(info, 0.0, sys.pos(dep), sys.vel(0) + vHelio);
        s.gravity().excludePair(static_cast<int>(dep), scIndex);
        if (!targetGravity) s.gravity().excludePair(static_cast<int>(tgt), scIndex);
        s.systemModified();
        return s;
    };
    auto missAt = [&](const Vec3& v, const Vec3& offset, Vec3* vRel) {
        int sc = 0;
        Simulation s = makeFlight(v, false, sc);
        s.advanceTo(arrival);
        const auto& sys = s.system();
        const std::size_t i = static_cast<std::size_t>(sc);
        if (vRel) *vRel = sys.vel(i) - sys.vel(tgt);
        return (sys.pos(i) - sys.pos(tgt)) - offset;
    };

    // ---- 5. Newton-Raphson single shooting ----------------------------------------------------
    auto shoot = [&](Vec3& v, const Vec3& offset) {
        double missKm = std::numeric_limits<double>::infinity();
        for (int it = 0; it < cfg.maxIterations; ++it) {
            const Vec3 F0 = missAt(v, offset, nullptr);
            missKm = F0.norm() * AU_KM;
            ++res.iterations;
            log(str::format("  shooting iteration %2d: miss %.1f km", res.iterations, missKm));
            if (missKm < cfg.toleranceKm) break;
            const double delta = 1e-7; // AU/day (~0.17 m/s)
            Mat3 J;
            for (int c = 0; c < 3; ++c) {
                Vec3 dv{};
                dv[static_cast<std::size_t>(c)] = delta;
                const Vec3 col = (missAt(v + dv, offset, nullptr) - F0) / delta;
                J.m[0][c] = col.x;
                J.m[1][c] = col.y;
                J.m[2][c] = col.z;
            }
            Mat3 Jinv;
            if (!J.inverse(Jinv)) break;
            Vec3 step = -(Jinv * F0);
            const double maxStep = 0.15 * v.norm();
            if (step.norm() > maxStep) step = step * (maxStep / step.norm());
            v += step;
        }
        return missKm;
    };
    log("N-body shooting to the centre of " + cfg.to + " (planetary gravity on the probe off):");
    shoot(v0, Vec3{});

    // ---- 6. B-plane targeting with the target's gravity ----------------------------------------
    Vec3 vRel;
    missAt(v0, Vec3{}, &vRel);
    const double vinfA = vRel.norm() * KMS_PER_AUD;
    Vec3 bDir = vRel.cross(Vec3::unitZ());
    if (bDir.norm() < 1e-12) bDir = Vec3::unitX();
    bDir = bDir.normalized();

    struct Encounter {
        double periKm = 0.0, timeJD = 0.0;
    };
    // Full-gravity flight around the arrival date; returns the closest approach.
    auto encounter = [&](const Vec3& v) {
        int sc = 0;
        Simulation s = makeFlight(v, true, sc);
        s.setIntegrator("dopri5");
        s.integrator().setTolerance(1e-12);
        const std::size_t i = static_cast<std::size_t>(sc);
        s.advanceTo(arrival - 8.0);
        Encounter e;
        e.periKm = std::numeric_limits<double>::infinity();
        double d0 = 0, d1 = 0;
        const double h = 0.002;
        for (int n = 0; s.time() < arrival + 8.0; ++n) {
            s.advanceTo(s.time() + h);
            const double d = (s.system().pos(i) - s.system().pos(tgt)).norm() * AU_KM;
            if (d < e.periKm) {
                e.periKm = d;
                e.timeJD = s.time();
            }
            if (n >= 2 && d1 < d0 && d1 <= d) {
                const double o = math::parabolicVertexOffset(d0, d1, d);
                const double dm = d1 - 0.25 * (d0 - d) * o;
                if (dm < e.periKm && dm > 0) {
                    e.periKm = dm;
                    e.timeJD = s.time() - h + o * h;
                }
            }
            d0 = d1;
            d1 = d;
            if (d < 0.2 * RtKm) break; // impact
        }
        return e;
    };

    // Patched-conic first guess: b = r_p sqrt(1 + 2 mu / (r_p v_inf^2))
    double b = rCapture * std::sqrt(1.0 + 2.0 * muT / (rCapture * vinfA * vinfA));
    log(str::format("B-plane targeting for a %.0f km periapsis (v_inf = %.3f km/s, patched-conic b = %.0f km):",
                    cfg.captureAltitudeKm, vinfA, b));
    double bPrev = 0.0, fPrev = 0.0;
    Encounter enc;
    Vec3 vBest = v0;
    for (int it = 0; it < 10; ++it) {
        Vec3 v = v0;
        res.finalMissKm = shoot(v, bDir * (b / AU_KM));
        enc = encounter(v);
        vBest = v;
        const double f = enc.periKm - rCapture;
        log(str::format("  aim offset %8.1f km -> periapsis altitude %8.1f km", b, enc.periKm - RtKm));
        if (std::fabs(f) < cfg.periapsisToleranceKm) break;
        double bNext;
        if (it == 0) bNext = b * (f < 0 ? 1.15 : 0.9);
        else bNext = b - f * (b - bPrev) / (f - fPrev);
        if (!std::isfinite(bNext) || bNext <= 0.0) bNext = 0.5 * b;
        bPrev = b;
        fPrev = f;
        b = bNext;
        v0 = v;
    }
    v0 = vBest;
    res.aimOffsetKm = b;

    // ---- 7. results --------------------------------------------------------------------------------
    const Vec3 vInfDep = v0 - vDep;
    res.vInfDepart = vInfDep.norm() * KMS_PER_AUD;
    res.c3 = res.vInfDepart * res.vInfDepart;
    const Vec3 eq = astro::eclipticToEquatorial(vInfDep.normalized());
    res.declinationDeg = std::asin(math::clamp(eq.z, -1.0, 1.0)) * RAD2DEG;
    res.dvDepart = burnFromCircular(muD, rPark, res.vInfDepart);
    res.vInfArrive = vinfA;
    res.periapsisAltKm = enc.periKm - RtKm;
    res.closestApproachJD = enc.timeJD;
    res.impact = res.periapsisAltKm < 0.0;
    res.dvCapture = burnFromCircular(muT, std::max(enc.periKm, RtKm), vinfA);
    {
        Simulation s = launchSnap;
        s.advanceTo(arrival);
        res.arrivalPosition = helioPos(s.system(), tgt);
        res.transferAngleDeg = angleBetween(rDep, res.arrivalPosition) * RAD2DEG;
    }

    // ---- 8. trajectory samples for plots ----------------------------------------------------------
    {
        int sc = 0;
        Simulation s = makeFlight(v0, true, sc);
        s.setIntegrator("dopri5");
        s.integrator().setTolerance(1e-11);
        sampleOrbit(s.system(), dep, res.departOrbit);
        sampleOrbit(s.system(), tgt, res.targetOrbit);
        const std::size_t i = static_cast<std::size_t>(sc);
        const double end = enc.timeJD > 0 ? enc.timeJD : arrival;
        while (s.time() < end) {
            res.trajectory.push_back(helioPos(s.system(), i));
            s.advanceTo(std::min(end, s.time() + 2.0));
        }
        res.trajectory.push_back(helioPos(s.system(), i));
    }
    res.success = std::isfinite(res.periapsisAltKm) && !res.impact &&
                  std::fabs(res.periapsisAltKm - cfg.captureAltitudeKm) < 10.0 * cfg.periapsisToleranceKm + 50.0;
    res.message = res.success ? "converged" : (res.impact ? "trajectory impacts the target" : "targeting did not fully converge");
    return res;
}

} // namespace solsim::analysis
