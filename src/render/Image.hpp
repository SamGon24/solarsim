// =============================================================================
//  SolarSim - Image.hpp
//  HDR framebuffer (linear float RGB + depth buffer) and the final 8-bit sRGB
//  image produced by tone mapping.
// =============================================================================
#pragma once

#include "Color.hpp"

#include <cstdint>
#include <limits>
#include <vector>

namespace solsim::render {

class Framebuffer {
public:
    Framebuffer() = default;
    Framebuffer(int w, int h) { resize(w, h); }

    void resize(int w, int h);
    int width() const { return m_w; }
    int height() const { return m_h; }
    bool inside(int x, int y) const { return x >= 0 && y >= 0 && x < m_w && y < m_h; }

    void clear(const Color& c = Color{});
    Color& at(int x, int y) { return m_color[static_cast<std::size_t>(y) * m_w + x]; }
    const Color& at(int x, int y) const { return m_color[static_cast<std::size_t>(y) * m_w + x]; }
    float& depth(int x, int y) { return m_depth[static_cast<std::size_t>(y) * m_w + x]; }
    float depth(int x, int y) const { return m_depth[static_cast<std::size_t>(y) * m_w + x]; }

    // Additive light (stars, glows, particles, trails).
    void add(int x, int y, const Color& c) {
        if (inside(x, y)) at(x, y) += c;
    }
    // Alpha blend over existing colour.
    void blend(int x, int y, const Color& c, float alpha) {
        if (!inside(x, y)) return;
        Color& d = at(x, y);
        d = d * (1.f - alpha) + c * alpha;
    }
    // Depth-tested variants; `z` is camera-space depth (smaller = closer).
    void addDepth(int x, int y, float z, const Color& c) {
        if (inside(x, y) && z < depth(x, y)) at(x, y) += c;
    }

    std::vector<Color>& pixels() { return m_color; }
    const std::vector<Color>& pixels() const { return m_color; }

    static constexpr float kFar = std::numeric_limits<float>::infinity();

private:
    int m_w = 0, m_h = 0;
    std::vector<Color> m_color;
    std::vector<float> m_depth;
};

// 8-bit sRGB image, row-major RGBRGB...
struct Image8 {
    int width = 0;
    int height = 0;
    std::vector<std::uint8_t> rgb;

    void resize(int w, int h) {
        width = w;
        height = h;
        rgb.assign(static_cast<std::size_t>(w) * h * 3, 0);
    }
    std::uint8_t* pixel(int x, int y) { return &rgb[(static_cast<std::size_t>(y) * width + x) * 3]; }
    const std::uint8_t* pixel(int x, int y) const { return &rgb[(static_cast<std::size_t>(y) * width + x) * 3]; }
};

} // namespace solsim::render
