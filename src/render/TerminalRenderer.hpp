// =============================================================================
//  SolarSim - TerminalRenderer.hpp
//  Converts a tone-mapped frame to terminal output. Each character cell shows
//  two vertically stacked pixels using the UPPER HALF BLOCK (U+2580): the
//  foreground colour paints the top pixel, the background the bottom one.
//
//  Colour modes:  truecolor (24-bit SGR), 256 (xterm palette), ascii (no
//  colour, luminance ramp - works on any terminal, including legacy consoles)
// =============================================================================
#pragma once

#include "Image.hpp"
#include "SceneRenderer.hpp"

#include <string>
#include <vector>

namespace solsim::render {

enum class TerminalColorMode { TrueColor, Ansi256, Ascii };
bool parseTerminalColorMode(const std::string& s, TerminalColorMode& out);
const char* terminalColorModeName(TerminalColorMode m);

class TerminalRenderer {
public:
    void setMode(TerminalColorMode m) { m_mode = m; }
    TerminalColorMode mode() const { return m_mode; }

    // Builds one full screen update (cursor home + all rows). `img` must be
    // columns wide and 2*(rows - header - footer) tall.
    std::string compose(const Image8& img, const std::vector<ScreenLabel>& labels,
                        const std::vector<std::string>& header, const std::vector<std::string>& footer,
                        int columns, int rows) const;

    static int xterm256(int r, int g, int b);

private:
    TerminalColorMode m_mode = TerminalColorMode::TrueColor;
};

} // namespace solsim::render
