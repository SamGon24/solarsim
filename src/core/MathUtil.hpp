// =============================================================================
//  SolarSim - MathUtil.hpp
//  Angle helpers, interpolation, generic root finders / optimizers and a tiny
//  least-squares line fit. Root finders are templated on the callable so they
//  can be used with lambdas that run whole simulations (event refinement).
// =============================================================================
#pragma once

#include "Constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace solsim::math {

template <typename T> constexpr T sqr(T v) { return v * v; }
template <typename T> constexpr T cube(T v) { return v * v * v; }
template <typename T> constexpr T clamp(T v, T lo, T hi) { return v < lo ? lo : (v > hi ? hi : v); }
template <typename T> constexpr T sign(T v) { return v < T(0) ? T(-1) : T(1); }
template <typename T> constexpr T lerp(T a, T b, T t) { return a + (b - a) * t; }

inline double smoothstep(double e0, double e1, double x) {
    const double t = clamp((x - e0) / (e1 - e0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

inline double deg(double rad) { return rad * constants::RAD2DEG; }
inline double rad(double deg) { return deg * constants::DEG2RAD; }

// Wrap to [0, 2*pi)
inline double wrapTwoPi(double a) {
    double r = std::fmod(a, constants::TWO_PI);
    if (r < 0.0) r += constants::TWO_PI;
    return r;
}
// Wrap to (-pi, pi]
inline double wrapPi(double a) {
    double r = wrapTwoPi(a);
    if (r > constants::PI) r -= constants::TWO_PI;
    return r;
}
// Wrap to [0, 360)
inline double wrap360(double a) {
    double r = std::fmod(a, 360.0);
    if (r < 0.0) r += 360.0;
    return r;
}
// Wrap to (-180, 180]
inline double wrap180(double a) {
    double r = wrap360(a);
    if (r > 180.0) r -= 360.0;
    return r;
}

// Unwraps a sequence of angles (radians) in place so that consecutive samples
// never jump by more than pi.
inline void unwrapAngles(std::vector<double>& a) {
    for (std::size_t i = 1; i < a.size(); ++i) {
        double d = a[i] - a[i - 1];
        while (d > constants::PI) { a[i] -= constants::TWO_PI; d -= constants::TWO_PI; }
        while (d < -constants::PI) { a[i] += constants::TWO_PI; d += constants::TWO_PI; }
    }
}

// -----------------------------------------------------------------------------
// Root finding
// -----------------------------------------------------------------------------
struct RootResult {
    double x = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
    bool converged = false;
};

// Bisection - robust, linear convergence. Requires f(a)*f(b) <= 0.
template <typename F>
RootResult bisect(F&& f, double a, double b, double tol = 1e-12, int maxIter = 200) {
    RootResult res;
    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0.0) return res;
    for (int i = 0; i < maxIter; ++i) {
        const double m = 0.5 * (a + b);
        const double fm = f(m);
        res.iterations = i + 1;
        if (fm == 0.0 || 0.5 * (b - a) < tol) {
            res.x = m;
            res.converged = true;
            return res;
        }
        if (fa * fm < 0.0) { b = m; fb = fm; } else { a = m; fa = fm; }
    }
    res.x = 0.5 * (a + b);
    return res;
}

// Brent's method (inverse quadratic interpolation + secant + bisection).
// Superlinear convergence with the robustness of bisection.
template <typename F>
RootResult brent(F&& f, double a, double b, double tol = 1e-12, int maxIter = 100) {
    RootResult res;
    double fa = f(a), fb = f(b);
    if (fa * fb > 0.0) return res;
    if (std::fabs(fa) < std::fabs(fb)) { std::swap(a, b); std::swap(fa, fb); }
    double c = a, fc = fa, d = b - a, s = b;
    bool mflag = true;
    for (int i = 0; i < maxIter; ++i) {
        res.iterations = i + 1;
        if (fb == 0.0 || std::fabs(b - a) < tol) {
            res.x = b;
            res.converged = true;
            return res;
        }
        if (fa != fc && fb != fc) {
            s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) +
                c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            s = b - fb * (b - a) / (fb - fa);
        }
        const double q = (3.0 * a + b) / 4.0;
        const bool cond1 = !((s > std::min(q, b) && s < std::max(q, b)));
        const bool cond2 = mflag && std::fabs(s - b) >= std::fabs(b - c) / 2.0;
        const bool cond3 = !mflag && std::fabs(s - b) >= std::fabs(c - d) / 2.0;
        const bool cond4 = mflag && std::fabs(b - c) < tol;
        const bool cond5 = !mflag && std::fabs(c - d) < tol;
        if (cond1 || cond2 || cond3 || cond4 || cond5) {
            s = 0.5 * (a + b);
            mflag = true;
        } else {
            mflag = false;
        }
        const double fs = f(s);
        d = c;
        c = b;
        fc = fb;
        if (fa * fs < 0.0) { b = s; fb = fs; } else { a = s; fa = fs; }
        if (std::fabs(fa) < std::fabs(fb)) { std::swap(a, b); std::swap(fa, fb); }
    }
    res.x = b;
    return res;
}

// Newton-Raphson with derivative callback.
template <typename F, typename DF>
RootResult newton(F&& f, DF&& df, double x0, double tol = 1e-14, int maxIter = 60) {
    RootResult res;
    double x = x0;
    for (int i = 0; i < maxIter; ++i) {
        const double fx = f(x);
        const double dfx = df(x);
        res.iterations = i + 1;
        if (dfx == 0.0) break;
        const double step = fx / dfx;
        x -= step;
        if (std::fabs(step) <= tol * std::max(1.0, std::fabs(x))) {
            res.x = x;
            res.converged = true;
            return res;
        }
    }
    res.x = x;
    return res;
}

// Golden-section search for the MAXIMUM of a unimodal function on [a, b].
template <typename F>
double goldenSectionMax(F&& f, double a, double b, double tol = 1e-9, int maxIter = 200) {
    const double invPhi = (std::sqrt(5.0) - 1.0) / 2.0;
    double c = b - invPhi * (b - a);
    double d = a + invPhi * (b - a);
    double fc = f(c), fd = f(d);
    for (int i = 0; i < maxIter && std::fabs(b - a) > tol; ++i) {
        if (fc > fd) {
            b = d; d = c; fd = fc;
            c = b - invPhi * (b - a);
            fc = f(c);
        } else {
            a = c; c = d; fc = fd;
            d = a + invPhi * (b - a);
            fd = f(d);
        }
    }
    return 0.5 * (a + b);
}

// Vertex of the parabola through three equally spaced samples y(-1), y(0), y(+1).
// Returns the offset in units of the sample spacing (range roughly [-1, 1]).
inline double parabolicVertexOffset(double ym1, double y0, double yp1) {
    const double den = ym1 - 2.0 * y0 + yp1;
    if (den == 0.0) return 0.0;
    return 0.5 * (ym1 - yp1) / den;
}

// -----------------------------------------------------------------------------
// Ordinary least squares fit y = intercept + slope * x
// -----------------------------------------------------------------------------
struct LinearFit {
    double slope = 0.0;
    double intercept = 0.0;
    double slopeStdErr = 0.0;
    double r2 = 0.0;
    std::size_t n = 0;
};

inline LinearFit fitLine(const std::vector<double>& x, const std::vector<double>& y) {
    LinearFit fit;
    const std::size_t n = std::min(x.size(), y.size());
    fit.n = n;
    if (n < 2) return fit;
    double mx = 0.0, my = 0.0;
    for (std::size_t i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
    mx /= double(n);
    my /= double(n);
    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = x[i] - mx, dy = y[i] - my;
        sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
    }
    if (sxx == 0.0) return fit;
    fit.slope = sxy / sxx;
    fit.intercept = my - fit.slope * mx;
    double ssr = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double r = y[i] - (fit.intercept + fit.slope * x[i]);
        ssr += r * r;
    }
    fit.r2 = syy > 0.0 ? 1.0 - ssr / syy : 1.0;
    if (n > 2) fit.slopeStdErr = std::sqrt(ssr / double(n - 2) / sxx);
    return fit;
}

} // namespace solsim::math
