// =============================================================================
//  SolarSim - TerminalRenderer.cpp
// =============================================================================
#include "TerminalRenderer.hpp"

#include "../core/StringUtil.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>

namespace solsim::render {

bool parseTerminalColorMode(const std::string& s, TerminalColorMode& out) {
    const std::string k = str::toLower(s);
    if (k == "truecolor" || k == "24bit" || k == "rgb") out = TerminalColorMode::TrueColor;
    else if (k == "256" || k == "ansi256" || k == "xterm") out = TerminalColorMode::Ansi256;
    else if (k == "ascii" || k == "none" || k == "mono") out = TerminalColorMode::Ascii;
    else return false;
    return true;
}

const char* terminalColorModeName(TerminalColorMode m) {
    switch (m) {
        case TerminalColorMode::TrueColor: return "truecolor";
        case TerminalColorMode::Ansi256: return "256";
        case TerminalColorMode::Ascii: return "ascii";
    }
    return "?";
}

int TerminalRenderer::xterm256(int r, int g, int b) {
    // Grayscale ramp (232..255) when the colour is nearly neutral.
    if (std::abs(r - g) < 10 && std::abs(g - b) < 10) {
        const int avg = (r + g + b) / 3;
        if (avg < 8) return 16;
        if (avg > 238) return 231;
        return 232 + (avg - 8) * 24 / 231;
    }
    auto q = [](int v) { return v < 48 ? 0 : (v < 115 ? 1 : (v - 35) / 40); };
    return 16 + 36 * q(r) + 6 * q(g) + q(b);
}

namespace {

struct Cell {
    char ch = 0;               // label character (0 = image cell)
    std::uint32_t color = 0;   // label colour
};

void appendColor(std::string& out, bool fg, int r, int g, int b, TerminalColorMode mode) {
    char buf[32];
    if (mode == TerminalColorMode::TrueColor) {
        std::snprintf(buf, sizeof buf, "\x1b[%d;2;%d;%d;%dm", fg ? 38 : 48, r, g, b);
    } else {
        std::snprintf(buf, sizeof buf, "\x1b[%d;5;%dm", fg ? 38 : 48, TerminalRenderer::xterm256(r, g, b));
    }
    out += buf;
}

std::string fitLine(const std::string& s, int columns) {
    std::string out;
    int n = 0;
    for (std::size_t i = 0; i < s.size() && n < columns; ++i) {
        out.push_back(s[i]);
        if ((static_cast<unsigned char>(s[i]) & 0xC0) != 0x80) ++n;
    }
    // Remove a trailing incomplete UTF-8 lead byte, if any.
    while (!out.empty() && (static_cast<unsigned char>(out.back()) & 0xC0) == 0xC0) out.pop_back();
    if (n < columns) out.append(static_cast<std::size_t>(columns - n), ' ');
    return out;
}

} // namespace

std::string TerminalRenderer::compose(const Image8& img, const std::vector<ScreenLabel>& labels,
                                      const std::vector<std::string>& header, const std::vector<std::string>& footer,
                                      int columns, int rows) const {
    const int imgRows = std::max(0, rows - static_cast<int>(header.size()) - static_cast<int>(footer.size()));
    std::vector<Cell> cells(static_cast<std::size_t>(columns) * static_cast<std::size_t>(std::max(1, imgRows)));
    for (const auto& l : labels) {
        const int cx = static_cast<int>(std::lround(l.x));
        const int cy = static_cast<int>(std::floor(l.y / 2.0 + 0.5));
        if (cy < 0 || cy >= imgRows) continue;
        for (std::size_t k = 0; k < l.text.size(); ++k) {
            const int x = cx + static_cast<int>(k);
            if (x < 0 || x >= columns) continue;
            const char ch = l.text[k];
            if (static_cast<unsigned char>(ch) < 32 || static_cast<unsigned char>(ch) > 126) continue;
            Cell& c = cells[static_cast<std::size_t>(cy) * columns + x];
            c.ch = ch;
            c.color = l.color;
        }
    }

    std::string out;
    out.reserve(static_cast<std::size_t>(columns) * rows * (m_mode == TerminalColorMode::Ascii ? 2 : 24));
    out += "\x1b[H";
    const bool color = m_mode != TerminalColorMode::Ascii;

    auto textLine = [&](const std::string& s, bool last) {
        if (color) out += "\x1b[0m\x1b[38;5;252m\x1b[48;5;235m";
        out += fitLine(s, columns);
        if (color) out += "\x1b[0m";
        if (!last) out += "\r\n";
    };

    int lineIndex = 0;
    for (const auto& h : header) textLine(h, ++lineIndex == rows);

    static const char ramp[] = " .,:;-=+*#%@";
    for (int cy = 0; cy < imgRows; ++cy) {
        int lastFg[3] = {-1, -1, -1}, lastBg[3] = {-1, -1, -1};
        for (int cx = 0; cx < columns; ++cx) {
            const int ty = std::min(img.height - 1, cy * 2), by = std::min(img.height - 1, cy * 2 + 1);
            const int ix = std::min(img.width - 1, cx);
            const std::uint8_t* top = img.pixel(ix, ty);
            const std::uint8_t* bot = img.pixel(ix, by);
            const Cell& cell = cells[static_cast<std::size_t>(cy) * columns + cx];
            if (!color) {
                if (cell.ch) {
                    out.push_back(cell.ch);
                } else {
                    const double lum = (0.2126 * (top[0] + bot[0]) + 0.7152 * (top[1] + bot[1]) + 0.0722 * (top[2] + bot[2])) / 510.0;
                    const int idx = std::clamp(static_cast<int>(std::pow(lum, 0.8) * 12.0), 0, 11);
                    out.push_back(ramp[idx]);
                }
                continue;
            }
            auto setFg = [&](int r, int g, int b) {
                if (r == lastFg[0] && g == lastFg[1] && b == lastFg[2]) return;
                appendColor(out, true, r, g, b, m_mode);
                lastFg[0] = r; lastFg[1] = g; lastFg[2] = b;
            };
            auto setBg = [&](int r, int g, int b) {
                if (r == lastBg[0] && g == lastBg[1] && b == lastBg[2]) return;
                appendColor(out, false, r, g, b, m_mode);
                lastBg[0] = r; lastBg[1] = g; lastBg[2] = b;
            };
            if (cell.ch) {
                setBg(bot[0] / 3, bot[1] / 3, bot[2] / 3);
                setFg(static_cast<int>((cell.color >> 16) & 0xFF), static_cast<int>((cell.color >> 8) & 0xFF),
                      static_cast<int>(cell.color & 0xFF));
                out.push_back(cell.ch);
            } else if (top[0] == bot[0] && top[1] == bot[1] && top[2] == bot[2]) {
                setBg(bot[0], bot[1], bot[2]);
                out.push_back(' ');
            } else {
                setFg(top[0], top[1], top[2]);
                setBg(bot[0], bot[1], bot[2]);
                out += "\xE2\x96\x80"; // U+2580 UPPER HALF BLOCK
            }
        }
        if (color) out += "\x1b[0m";
        if (++lineIndex < rows) out += "\r\n";
    }
    for (const auto& f : footer) textLine(f, ++lineIndex >= rows);
    return out;
}

} // namespace solsim::render
