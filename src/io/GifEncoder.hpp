// =============================================================================
//  SolarSim - GifEncoder.hpp
//  Animated GIF89a writer:
//    * per-frame 256-colour palette by MEDIAN CUT (Heckbert 1982) on a
//      15-bit colour histogram
//    * optional Floyd-Steinberg error diffusion dithering
//    * variable-width LZW compression (GIF flavour: 9..12 bit codes,
//      clear/end codes, 255-byte data sub-blocks)
//    * NETSCAPE2.0 application extension for looping
// =============================================================================
#pragma once

#include "../render/Image.hpp"

#include <array>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace solsim::io {

struct Palette {
    std::vector<std::array<std::uint8_t, 3>> colors; // up to 256 entries
};

class MedianCutQuantizer {
public:
    // Builds a palette of at most `maxColors` from the image.
    Palette build(const render::Image8& img, int maxColors = 256);
    // Maps pixels to palette indices (optionally dithered).
    void map(const render::Image8& img, const Palette& pal, bool dither, std::vector<std::uint8_t>& indices);

private:
    std::vector<std::uint16_t> m_binToIndex; // 32768 entries, 0xFFFF = unknown
    std::uint8_t nearest(const Palette& pal, int r, int g, int b);
};

// LZW compression of 8-bit indices with GIF code-size semantics.
std::vector<std::uint8_t> gifLzwEncode(const std::vector<std::uint8_t>& indices, int minCodeSize);

class GifEncoder {
public:
    GifEncoder() = default;
    ~GifEncoder();
    GifEncoder(const GifEncoder&) = delete;
    GifEncoder& operator=(const GifEncoder&) = delete;

    // loopCount 0 = loop forever.
    bool open(const std::string& path, int width, int height, int loopCount = 0);
    // delay in hundredths of a second.
    bool addFrame(const render::Image8& frame, int delayCs, bool dither = false);
    bool close();

    bool isOpen() const { return m_file != nullptr; }
    int frameCount() const { return m_frames; }
    std::size_t bytesWritten() const { return m_bytes; }
    const std::string& error() const { return m_error; }

private:
    void put(const void* data, std::size_t n);
    void putByte(std::uint8_t b) { put(&b, 1); }
    void putU16(std::uint16_t v) {
        putByte(static_cast<std::uint8_t>(v & 0xFF));
        putByte(static_cast<std::uint8_t>(v >> 8));
    }

    std::FILE* m_file = nullptr;
    int m_w = 0, m_h = 0;
    int m_frames = 0;
    std::size_t m_bytes = 0;
    std::string m_error;
    MedianCutQuantizer m_quant;
    std::vector<std::uint8_t> m_indices;
};

} // namespace solsim::io
