// =============================================================================
//  SolarSim - Font.cpp
// =============================================================================
#include "Font.hpp"

#include "Image.hpp"

namespace solsim::render {

namespace {
// Decodes the text into glyph indices (handles the UTF-8 degree sign).
template <typename F>
void forEachGlyph(const std::string& s, F&& fn) {
    for (std::size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        if (c == 0xC2 && i + 1 < s.size() && static_cast<unsigned char>(s[i + 1]) == 0xB0) {
            fn(static_cast<char>(127));
            ++i;
            continue;
        }
        if (c >= 0x80) {
            // Skip other multi-byte sequences, render as '?'
            while (i + 1 < s.size() && (static_cast<unsigned char>(s[i + 1]) & 0xC0) == 0x80) ++i;
            fn('?');
            continue;
        }
        fn(static_cast<char>(c));
    }
}
} // namespace

const std::uint8_t* BitmapFont::glyph(char c) {
    int code = static_cast<unsigned char>(c);
    if (code < 32 || code > 127) code = '?';
    return kFont5x7[code - 32];
}

bool BitmapFont::pixel(char c, int col, int row) {
    if (col < 0 || col >= kGlyphW || row < 0 || row >= kGlyphH) return false;
    return (glyph(c)[row] >> (kGlyphW - 1 - col)) & 1u;
}

int BitmapFont::textWidth(const std::string& s, int scale) {
    int n = 0;
    forEachGlyph(s, [&](char) { ++n; });
    return n == 0 ? 0 : (n * kAdvance - 1) * scale;
}

namespace {
// Shared glyph walker: plot(x, y, isShadowPass) is called for every lit pixel.
template <typename Plot>
void walkText(int x, int y, const std::string& text, int scale, bool shadow, Plot&& plot) {
    // Pass 0 draws the drop shadow, pass 1 the glyphs (so shadows never cover glyph pixels).
    for (int pass = shadow ? 0 : 1; pass < 2; ++pass) {
        const int off = pass == 0 ? scale : 0;
        int cx = x;
        forEachGlyph(text, [&](char ch) {
            const std::uint8_t* g = BitmapFont::glyph(ch);
            for (int row = 0; row < BitmapFont::kGlyphH; ++row) {
                for (int col = 0; col < BitmapFont::kGlyphW; ++col) {
                    if (!((g[row] >> (BitmapFont::kGlyphW - 1 - col)) & 1u)) continue;
                    for (int sy = 0; sy < scale; ++sy)
                        for (int sx = 0; sx < scale; ++sx)
                            plot(cx + col * scale + sx + off, y + row * scale + sy + off, pass == 0);
                }
            }
            cx += BitmapFont::kAdvance * scale;
        });
    }
}
} // namespace

void BitmapFont::draw(Framebuffer& fb, int x, int y, const std::string& text, const Color& color, int scale,
                      bool shadow, float alpha) {
    walkText(x, y, text, scale, shadow, [&](int px, int py, bool isShadow) {
        if (isShadow) fb.blend(px, py, Color{}, 0.75f * alpha);
        else fb.blend(px, py, color, alpha);
    });
}

void BitmapFont::draw(Image8& img, int x, int y, const std::string& text, std::uint32_t rgb, int scale, bool shadow,
                      float alpha) {
    const float cr = float((rgb >> 16) & 0xFF), cg = float((rgb >> 8) & 0xFF), cb = float(rgb & 0xFF);
    walkText(x, y, text, scale, shadow, [&](int px, int py, bool isShadow) {
        if (px < 0 || py < 0 || px >= img.width || py >= img.height) return;
        std::uint8_t* p = img.pixel(px, py);
        const float a = isShadow ? 0.7f * alpha : alpha;
        const float tr = isShadow ? 0.f : cr, tg = isShadow ? 0.f : cg, tb = isShadow ? 0.f : cb;
        p[0] = static_cast<std::uint8_t>(p[0] * (1.f - a) + tr * a);
        p[1] = static_cast<std::uint8_t>(p[1] * (1.f - a) + tg * a);
        p[2] = static_cast<std::uint8_t>(p[2] * (1.f - a) + tb * a);
    });
}

} // namespace solsim::render
