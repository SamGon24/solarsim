// =============================================================================
//  SolarSim - PostProcess.cpp
// =============================================================================
#include "PostProcess.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace solsim::render {

namespace {

// 1D gaussian kernel (normalised).
std::vector<float> gaussianKernel(float sigma) {
    const int radius = std::max(1, static_cast<int>(std::ceil(sigma * 3.0f)));
    std::vector<float> k(static_cast<std::size_t>(2 * radius + 1));
    float sum = 0.f;
    for (int i = -radius; i <= radius; ++i) {
        const float v = std::exp(-0.5f * float(i * i) / (sigma * sigma));
        k[static_cast<std::size_t>(i + radius)] = v;
        sum += v;
    }
    for (float& v : k) v /= sum;
    return k;
}

void blurSeparable(std::vector<Color>& img, int w, int h, float sigma) {
    const std::vector<float> k = gaussianKernel(sigma);
    const int r = static_cast<int>(k.size() / 2);
    std::vector<Color> tmp(img.size());
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            Color acc;
            for (int i = -r; i <= r; ++i) {
                const int xx = std::clamp(x + i, 0, w - 1);
                acc += img[static_cast<std::size_t>(y) * w + xx] * k[static_cast<std::size_t>(i + r)];
            }
            tmp[static_cast<std::size_t>(y) * w + x] = acc;
        }
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            Color acc;
            for (int i = -r; i <= r; ++i) {
                const int yy = std::clamp(y + i, 0, h - 1);
                acc += tmp[static_cast<std::size_t>(yy) * w + x] * k[static_cast<std::size_t>(i + r)];
            }
            img[static_cast<std::size_t>(y) * w + x] = acc;
        }
}

inline float aces(float x) {
    const float a = 2.51f, b = 0.03f, c = 2.43f, d = 0.59f, e = 0.14f;
    return std::clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0f, 1.0f);
}

} // namespace

void applyBloom(Framebuffer& fb, const PostSettings& s) {
    if (!s.bloom || s.bloomStrength <= 0.f) return;
    const int W = fb.width(), H = fb.height();
    const int factor = (H >= 240) ? 4 : 2;
    const int w = std::max(1, W / factor), h = std::max(1, H / factor);
    std::vector<Color> small(static_cast<std::size_t>(w) * h);
    // Bright pass + box downsample.
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            Color acc;
            for (int dy = 0; dy < factor; ++dy)
                for (int dx = 0; dx < factor; ++dx) {
                    const int X = std::min(W - 1, x * factor + dx), Y = std::min(H - 1, y * factor + dy);
                    const Color& c = fb.at(X, Y);
                    const float l = c.maxComponent();
                    if (l > s.bloomThreshold) acc += c * ((l - s.bloomThreshold) / std::max(l, 1e-6f));
                }
            small[static_cast<std::size_t>(y) * w + x] = acc * (1.f / float(factor * factor));
        }
    const float sigma = std::max(1.0f, s.bloomRadius * float(H) / float(factor));
    blurSeparable(small, w, h, sigma);
    // A second, wider pass gives a softer halo.
    std::vector<Color> wide = small;
    blurSeparable(wide, w, h, sigma * 3.0f);
    // Bilinear upsample and add.
    for (int Y = 0; Y < H; ++Y)
        for (int X = 0; X < W; ++X) {
            const float fx = (float(X) + 0.5f) / float(factor) - 0.5f;
            const float fy = (float(Y) + 0.5f) / float(factor) - 0.5f;
            const int x0 = std::clamp(static_cast<int>(std::floor(fx)), 0, w - 1);
            const int y0 = std::clamp(static_cast<int>(std::floor(fy)), 0, h - 1);
            const int x1 = std::min(w - 1, x0 + 1), y1 = std::min(h - 1, y0 + 1);
            const float tx = std::clamp(fx - float(x0), 0.f, 1.f), ty = std::clamp(fy - float(y0), 0.f, 1.f);
            auto sample = [&](const std::vector<Color>& img) {
                const Color a = Color::lerp(img[static_cast<std::size_t>(y0) * w + x0], img[static_cast<std::size_t>(y0) * w + x1], tx);
                const Color b = Color::lerp(img[static_cast<std::size_t>(y1) * w + x0], img[static_cast<std::size_t>(y1) * w + x1], tx);
                return Color::lerp(a, b, ty);
            };
            fb.at(X, Y) += (sample(small) * 0.6f + sample(wide) * 0.4f) * s.bloomStrength;
        }
}

void toneMap(const Framebuffer& fb, Image8& out, const PostSettings& s) {
    const int W = fb.width(), H = fb.height();
    out.resize(W, H);
    const float cx = 0.5f * float(W), cy = 0.5f * float(H);
    const float invR2 = 1.f / (cx * cx + cy * cy);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
            Color c = fb.at(x, y) * s.exposure;
            if (s.vignette) {
                const float dx = float(x) - cx, dy = float(y) - cy;
                c *= 1.f - 0.25f * (dx * dx + dy * dy) * invR2;
            }
            std::uint8_t* p = out.pixel(x, y);
            p[0] = static_cast<std::uint8_t>(std::lround(Color::linearToSrgb(aces(c.r)) * 255.f));
            p[1] = static_cast<std::uint8_t>(std::lround(Color::linearToSrgb(aces(c.g)) * 255.f));
            p[2] = static_cast<std::uint8_t>(std::lround(Color::linearToSrgb(aces(c.b)) * 255.f));
        }
}

void finishFrame(Framebuffer& fb, Image8& out, const PostSettings& s) {
    applyBloom(fb, s);
    toneMap(fb, out, s);
}

} // namespace solsim::render
