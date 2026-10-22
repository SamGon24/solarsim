// =============================================================================
//  SolarSim - Color.hpp
//  Linear-light RGB colour (float, HDR capable) plus conversions:
//    * sRGB <-> linear transfer functions
//    * black-body colour temperature -> RGB
//    * stellar B-V colour index -> effective temperature (Ballesteros 2012)
// =============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace solsim::render {

struct Color {
    float r = 0.f, g = 0.f, b = 0.f;

    constexpr Color() = default;
    constexpr Color(float r_, float g_, float b_) : r(r_), g(g_), b(b_) {}
    static constexpr Color gray(float v) { return {v, v, v}; }

    constexpr Color operator+(const Color& o) const { return {r + o.r, g + o.g, b + o.b}; }
    constexpr Color operator-(const Color& o) const { return {r - o.r, g - o.g, b - o.b}; }
    constexpr Color operator*(float s) const { return {r * s, g * s, b * s}; }
    constexpr Color operator*(const Color& o) const { return {r * o.r, g * o.g, b * o.b}; }
    Color& operator+=(const Color& o) { r += o.r; g += o.g; b += o.b; return *this; }
    Color& operator*=(float s) { r *= s; g *= s; b *= s; return *this; }

    float luminance() const { return 0.2126f * r + 0.7152f * g + 0.0722f * b; }
    float maxComponent() const { return std::max(r, std::max(g, b)); }

    static float srgbToLinear(float c) {
        return c <= 0.04045f ? c / 12.92f : std::pow((c + 0.055f) / 1.055f, 2.4f);
    }
    static float linearToSrgb(float c) {
        c = std::max(0.f, c);
        return c <= 0.0031308f ? 12.92f * c : 1.055f * std::pow(c, 1.f / 2.4f) - 0.055f;
    }

    // 0xRRGGBB (sRGB) -> linear colour
    static Color fromHex(std::uint32_t hex) {
        return {srgbToLinear(float((hex >> 16) & 0xFF) / 255.f), srgbToLinear(float((hex >> 8) & 0xFF) / 255.f),
                srgbToLinear(float(hex & 0xFF) / 255.f)};
    }

    static Color lerp(const Color& a, const Color& b, float t) { return a + (b - a) * t; }
};

// Approximate sRGB colour of a black body (Tanner Helland's fit, 1000-40000 K),
// returned in linear light and normalised so the brightest channel is 1.
inline Color blackbody(double kelvin) {
    const double t = std::clamp(kelvin, 1000.0, 40000.0) / 100.0;
    double r, g, b;
    if (t <= 66.0) {
        r = 255.0;
        g = 99.4708025861 * std::log(t) - 161.1195681661;
    } else {
        r = 329.698727446 * std::pow(t - 60.0, -0.1332047592);
        g = 288.1221695283 * std::pow(t - 60.0, -0.0755148492);
    }
    if (t >= 66.0) b = 255.0;
    else if (t <= 19.0) b = 0.0;
    else b = 138.5177312231 * std::log(t - 10.0) - 305.0447927307;
    auto ch = [](double v) { return Color::srgbToLinear(float(std::clamp(v, 0.0, 255.0) / 255.0)); };
    Color c{ch(r), ch(g), ch(b)};
    const float m = c.maxComponent();
    return m > 0 ? c * (1.f / m) : c;
}

// B-V colour index -> effective temperature (Ballesteros 2012, EPL 97, 34008)
inline double bvToTemperature(double bv) {
    return 4600.0 * (1.0 / (0.92 * bv + 1.7) + 1.0 / (0.92 * bv + 0.62));
}

inline Color starColor(double bv) { return blackbody(bvToTemperature(bv)); }

} // namespace solsim::render
