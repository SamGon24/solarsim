// =============================================================================
//  SolarSim - Raster.cpp
// =============================================================================
#include "Raster.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace solsim::render {

namespace {

inline double fpart(double x) { return x - std::floor(x); }
inline double rfpart(double x) { return 1.0 - fpart(x); }

inline void plot(Framebuffer& fb, int x, int y, double z, double coverage, const LineStyle& st) {
    if (!fb.inside(x, y) || coverage <= 0.0) return;
    if (st.depthTest && static_cast<float>(z) >= fb.depth(x, y)) return;
    const float a = static_cast<float>(coverage) * st.alpha;
    if (st.blend == BlendMode::Additive) fb.at(x, y) += st.color * a;
    else fb.blend(x, y, st.color, std::min(1.f, a));
}

// Cohen-Sutherland style clipping of a segment to a padded viewport.
bool clipSegment(double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double w, double h) {
    const double xmin = -2, ymin = -2, xmax = w + 2, ymax = h + 2;
    auto code = [&](double x, double y) {
        int c = 0;
        if (x < xmin) c |= 1; else if (x > xmax) c |= 2;
        if (y < ymin) c |= 4; else if (y > ymax) c |= 8;
        return c;
    };
    int c0 = code(x0, y0), c1 = code(x1, y1);
    for (int iter = 0; iter < 8; ++iter) {
        if (!(c0 | c1)) return true;
        if (c0 & c1) return false;
        const int co = c0 ? c0 : c1;
        double x = 0, y = 0, t = 0;
        if (co & 8) { t = (ymax - y0) / (y1 - y0); y = ymax; x = x0 + (x1 - x0) * t; }
        else if (co & 4) { t = (ymin - y0) / (y1 - y0); y = ymin; x = x0 + (x1 - x0) * t; }
        else if (co & 2) { t = (xmax - x0) / (x1 - x0); x = xmax; y = y0 + (y1 - y0) * t; }
        else { t = (xmin - x0) / (x1 - x0); x = xmin; y = y0 + (y1 - y0) * t; }
        const double z = z0 + (z1 - z0) * t;
        if (co == c0) { x0 = x; y0 = y; z0 = z; c0 = code(x0, y0); }
        else { x1 = x; y1 = y; z1 = z; c1 = code(x1, y1); }
    }
    return false;
}

} // namespace

void drawLine(Framebuffer& fb, double x0, double y0, double z0, double x1, double y1, double z1,
              const LineStyle& st) {
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(x1) || !std::isfinite(y1)) return;
    if (!clipSegment(x0, y0, z0, x1, y1, z1, fb.width(), fb.height())) return;

    const bool steep = std::fabs(y1 - y0) > std::fabs(x1 - x0);
    if (steep) { std::swap(x0, y0); std::swap(x1, y1); }
    if (x0 > x1) { std::swap(x0, x1); std::swap(y0, y1); std::swap(z0, z1); }
    const double dx = x1 - x0, dy = y1 - y0;
    const double gradient = dx == 0.0 ? 1.0 : dy / dx;
    // Brightness compensation so diagonal lines are not dimmer than axial ones.
    const double lenFactor = std::sqrt(1.0 + gradient * gradient);

    auto put = [&](int a, int b, double z, double c) {
        if (steep) plot(fb, b, a, z, c * lenFactor, st);
        else plot(fb, a, b, z, c * lenFactor, st);
    };

    // First endpoint
    double xend = std::round(x0);
    double yend = y0 + gradient * (xend - x0);
    double xgap = rfpart(x0 + 0.5);
    const int xpxl1 = static_cast<int>(xend);
    const int ypxl1 = static_cast<int>(std::floor(yend));
    put(xpxl1, ypxl1, z0, rfpart(yend) * xgap);
    put(xpxl1, ypxl1 + 1, z0, fpart(yend) * xgap);
    double intery = yend + gradient;

    // Second endpoint
    xend = std::round(x1);
    yend = y1 + gradient * (xend - x1);
    xgap = fpart(x1 + 0.5);
    const int xpxl2 = static_cast<int>(xend);
    const int ypxl2 = static_cast<int>(std::floor(yend));
    put(xpxl2, ypxl2, z1, rfpart(yend) * xgap);
    put(xpxl2, ypxl2 + 1, z1, fpart(yend) * xgap);

    const double span = std::max(1, xpxl2 - xpxl1);
    for (int x = xpxl1 + 1; x < xpxl2; ++x) {
        const double t = (x - xpxl1) / span;
        const double z = z0 + (z1 - z0) * t;
        const int yi = static_cast<int>(std::floor(intery));
        put(x, yi, z, rfpart(intery));
        put(x, yi + 1, z, fpart(intery));
        intery += gradient;
    }
}

void drawSplat(Framebuffer& fb, double x, double y, double z, const Color& c, double radius, bool depthTest) {
    const double r = std::max(0.35, radius);
    const int reach = static_cast<int>(std::ceil(r * 2.2));
    const int ix = static_cast<int>(std::floor(x)), iy = static_cast<int>(std::floor(y));
    const double inv2s2 = 1.0 / (2.0 * r * r * 0.5);
    // Normalise so the total energy is independent of the radius.
    double norm = 0.0;
    for (int dy = -reach; dy <= reach; ++dy)
        for (int dx = -reach; dx <= reach; ++dx) {
            const double ddx = ix + dx + 0.5 - x, ddy = iy + dy + 0.5 - y;
            norm += std::exp(-(ddx * ddx + ddy * ddy) * inv2s2);
        }
    if (norm <= 0.0) return;
    for (int dy = -reach; dy <= reach; ++dy)
        for (int dx = -reach; dx <= reach; ++dx) {
            const int px = ix + dx, py = iy + dy;
            if (!fb.inside(px, py)) continue;
            if (depthTest && static_cast<float>(z) >= fb.depth(px, py)) continue;
            const double ddx = px + 0.5 - x, ddy = py + 0.5 - y;
            const double w = std::exp(-(ddx * ddx + ddy * ddy) * inv2s2) / norm;
            fb.at(px, py) += c * static_cast<float>(w);
        }
}

void drawGlow(Framebuffer& fb, double cx, double cy, double r0, double maxRadius, const Color& c, double falloff) {
    if (r0 <= 0.0) return;
    const int x0 = std::max(0, static_cast<int>(cx - maxRadius));
    const int x1 = std::min(fb.width() - 1, static_cast<int>(cx + maxRadius));
    const int y0 = std::max(0, static_cast<int>(cy - maxRadius));
    const int y1 = std::min(fb.height() - 1, static_cast<int>(cy + maxRadius));
    const double inv = 1.0 / (r0 * r0);
    const double max2 = maxRadius * maxRadius;
    for (int y = y0; y <= y1; ++y)
        for (int x = x0; x <= x1; ++x) {
            const double dx = x + 0.5 - cx, dy = y + 0.5 - cy;
            const double d2 = dx * dx + dy * dy;
            if (d2 > max2) continue;
            double v = 1.0 / std::pow(1.0 + d2 * inv, falloff);
            // Smooth fade to zero at the cut-off radius.
            v *= 1.0 - d2 / max2;
            fb.at(x, y) += c * static_cast<float>(v);
        }
}

void drawCircle(Framebuffer& fb, double cx, double cy, double radius, const Color& c, float alpha) {
    const int segments = std::clamp(static_cast<int>(radius * 0.8), 16, 256);
    LineStyle st;
    st.color = c;
    st.alpha = alpha;
    st.depthTest = false;
    double px = cx + radius, py = cy;
    for (int k = 1; k <= segments; ++k) {
        const double a = 6.283185307179586 * k / segments;
        const double nx = cx + radius * std::cos(a), ny = cy + radius * std::sin(a);
        drawLine(fb, px, py, 0, nx, ny, 0, st);
        px = nx;
        py = ny;
    }
}

} // namespace solsim::render
