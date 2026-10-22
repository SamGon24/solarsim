// =============================================================================
//  SolarSim - Raster.hpp
//  2D rasterisation primitives on the HDR framebuffer:
//    * anti-aliased lines (Xiaolin Wu) with per-pixel depth interpolation
//    * gaussian point splats (stars, particles)
//    * radial glows (Sun corona / lens halo)
//    * circles
// =============================================================================
#pragma once

#include "Color.hpp"
#include "Image.hpp"

namespace solsim::render {

enum class BlendMode { Additive, Alpha };

struct LineStyle {
    Color color{1, 1, 1};
    float alpha = 1.f;        // multiplies coverage
    BlendMode blend = BlendMode::Additive;
    bool depthTest = true;
};

// Line between two projected points; z0/z1 are camera-space depths.
void drawLine(Framebuffer& fb, double x0, double y0, double z0, double x1, double y1, double z1,
              const LineStyle& style);

// Soft round point with gaussian falloff. radius ~ 1 gives a 3x3 footprint.
void drawSplat(Framebuffer& fb, double x, double y, double z, const Color& c, double radius, bool depthTest);

// Additive radial glow: I(r) = intensity / (1 + (r/r0)^2)^falloff, cut at maxRadius.
void drawGlow(Framebuffer& fb, double cx, double cy, double r0, double maxRadius, const Color& c,
              double falloff = 1.0);

// Thin anti-aliased circle outline (screen space).
void drawCircle(Framebuffer& fb, double cx, double cy, double radius, const Color& c, float alpha);

} // namespace solsim::render
