// =============================================================================
//  SolarSim - OrbitalElements.cpp
// =============================================================================
#include "OrbitalElements.hpp"

#include "Kepler.hpp"

#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"
#include "../core/StringUtil.hpp"

#include <cmath>
#include <limits>

namespace solsim::astro {

using namespace solsim::constants;

double OrbitalElements::a() const {
    if (isParabolic()) return std::numeric_limits<double>::infinity();
    return q / (1.0 - e);
}

double OrbitalElements::Q() const {
    if (!isElliptic()) return std::numeric_limits<double>::infinity();
    return a() * (1.0 + e);
}

double OrbitalElements::meanMotion() const {
    if (isParabolic()) return std::sqrt(mu / (2.0 * q * q * q)); // Barker's pseudo mean motion
    const double aa = std::fabs(a());
    return std::sqrt(mu / (aa * aa * aa));
}

double OrbitalElements::period() const {
    if (!isElliptic()) return std::numeric_limits<double>::infinity();
    return TWO_PI / meanMotion();
}

double OrbitalElements::meanAnomalyAt(double jd) const { return meanMotion() * (jd - tp); }

double OrbitalElements::trueAnomalyAt(double jd) const {
    if (isElliptic()) {
        const double E = solveKeplerElliptic(meanAnomalyAt(jd), e);
        return trueFromEccentric(E, e);
    }
    if (isHyperbolic()) {
        const double H = solveKeplerHyperbolic(meanAnomalyAt(jd), e);
        return trueFromHyperbolic(H, e);
    }
    const double W = 3.0 * std::sqrt(mu / (2.0 * q * q * q)) * (jd - tp);
    return solveBarker(W);
}

double OrbitalElements::meanLongitudeAt(double jd) const {
    return math::wrapTwoPi(node + argp + meanAnomalyAt(jd));
}

double OrbitalElements::specificEnergy() const {
    if (isParabolic()) return 0.0;
    return -mu / (2.0 * a());
}

double OrbitalElements::specificAngularMomentum() const { return std::sqrt(mu * p()); }

OrbitalElements OrbitalElements::fromMeanAnomaly(double a, double e, double i, double node, double argp,
                                                 double M, double epoch, double mu) {
    OrbitalElements el;
    el.e = e;
    el.q = a * (1.0 - e);
    el.i = i;
    el.node = node;
    el.argp = argp;
    el.epoch = epoch;
    el.mu = mu;
    const double n = std::sqrt(mu / std::fabs(a * a * a));
    el.tp = epoch - M / n;
    return el;
}

OrbitalElements OrbitalElements::fromMeanLongitude(double a, double e, double i, double node, double varpi,
                                                   double L, double epoch, double mu) {
    return fromMeanAnomaly(a, e, i, node, varpi - node, math::wrapPi(L - varpi), epoch, mu);
}

std::string OrbitalElements::describe() const {
    return str::format("q=%.6f AU e=%.6f i=%.4f deg node=%.4f deg argp=%.4f deg tp=JD%.4f", q, e,
                       i * RAD2DEG, node * RAD2DEG, argp * RAD2DEG, tp);
}

void perifocalBasis(double i, double node, double argp, Vec3& P, Vec3& Q, Vec3& W) {
    const double cO = std::cos(node), sO = std::sin(node);
    const double cw = std::cos(argp), sw = std::sin(argp);
    const double ci = std::cos(i), si = std::sin(i);
    P = {cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si};
    Q = {-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si};
    W = {sO * si, -cO * si, ci};
}

StateVector stateAtTrueAnomaly(const OrbitalElements& el, double nu) {
    const double p = el.p();
    const double cn = std::cos(nu), sn = std::sin(nu);
    const double r = p / (1.0 + el.e * cn);
    const double sq = std::sqrt(el.mu / p);
    Vec3 P, Q, W;
    perifocalBasis(el.i, el.node, el.argp, P, Q, W);
    StateVector s;
    s.r = P * (r * cn) + Q * (r * sn);
    s.v = P * (-sq * sn) + Q * (sq * (el.e + cn));
    return s;
}

StateVector elementsToState(const OrbitalElements& el, double jd) {
    return stateAtTrueAnomaly(el, el.trueAnomalyAt(jd));
}

OrbitalElements stateToElements(const Vec3& r, const Vec3& v, double mu, double jd) {
    OrbitalElements el;
    el.mu = mu;
    el.epoch = jd;

    const double rmag = r.norm();
    const double vmag2 = v.norm2();
    const Vec3 h = r.cross(v);
    const double hmag = h.norm();
    const Vec3 nVec = Vec3::unitZ().cross(h); // points to ascending node
    const double nmag = nVec.norm();
    const Vec3 eVec = (r * (vmag2 - mu / rmag) - v * r.dot(v)) / mu;
    const double e = eVec.norm();

    el.e = e;
    el.i = std::acos(math::clamp(h.z / hmag, -1.0, 1.0));
    const double p = hmag * hmag / mu;
    el.q = p / (1.0 + e);

    const bool equatorial = nmag < 1e-12 * hmag;
    const bool circular = e < 1e-11;

    // Longitude of ascending node.
    el.node = equatorial ? 0.0 : math::wrapTwoPi(std::atan2(nVec.y, nVec.x));

    // Argument of periapsis and true anomaly, handling singular cases.
    double nu;
    if (circular) {
        el.argp = 0.0;
        if (equatorial) {
            nu = std::atan2(r.y, r.x);                           // true longitude
            if (h.z < 0) nu = -nu;
        } else {
            // Argument of latitude u measured from the node.
            const Vec3 nHat = nVec / nmag;
            const Vec3 inPlaneY = h.cross(nVec).normalized();
            nu = std::atan2(r.dot(inPlaneY), r.dot(nHat));
        }
    } else {
        if (equatorial) {
            el.argp = std::atan2(eVec.y, eVec.x);               // longitude of periapsis
            if (h.z < 0) el.argp = -el.argp;
        } else {
            const Vec3 nHat = nVec / nmag;
            const Vec3 inPlaneY = h.cross(nVec).normalized();
            el.argp = std::atan2(eVec.dot(inPlaneY), eVec.dot(nHat));
        }
        el.argp = math::wrapTwoPi(el.argp);
        const Vec3 eHat = eVec / e;
        const Vec3 qHat = h.cross(eVec).normalized();
        nu = std::atan2(r.dot(qHat), r.dot(eHat));
    }

    // Time of periapsis passage from the anomaly.
    if (el.isElliptic()) {
        const double E = eccentricFromTrue(nu, e);
        const double M = meanFromEccentric(E, e);
        el.tp = jd - M / el.meanMotion();
    } else if (el.isHyperbolic()) {
        const double H = hyperbolicFromTrue(nu, e);
        const double M = meanFromHyperbolic(H, e);
        el.tp = jd - M / el.meanMotion();
    } else {
        // Barker: t - tp = sqrt(2 q^3 / mu) * (D + D^3/3), D = tan(nu/2)
        const double D = std::tan(0.5 * nu);
        el.tp = jd - std::sqrt(2.0 * el.q * el.q * el.q / mu) * (D + D * D * D / 3.0);
    }
    return el;
}

} // namespace solsim::astro
