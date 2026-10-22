// =============================================================================
//  SolarSim - ImageWriter.cpp
// =============================================================================
#include "ImageWriter.hpp"

#include "Checksum.hpp"
#include "Deflate.hpp"

#include "../core/StringUtil.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace solsim::io {

ImageFormat formatFromPath(const std::string& path) {
    const std::string p = str::toLower(path);
    if (str::endsWith(p, ".png")) return ImageFormat::PNG;
    if (str::endsWith(p, ".bmp")) return ImageFormat::BMP;
    if (str::endsWith(p, ".ppm") || str::endsWith(p, ".pnm")) return ImageFormat::PPM;
    return ImageFormat::Unknown;
}

namespace {

void putU32BE(std::vector<std::uint8_t>& v, std::uint32_t x) {
    v.push_back(static_cast<std::uint8_t>(x >> 24));
    v.push_back(static_cast<std::uint8_t>(x >> 16));
    v.push_back(static_cast<std::uint8_t>(x >> 8));
    v.push_back(static_cast<std::uint8_t>(x));
}
void putU16LE(std::vector<std::uint8_t>& v, std::uint32_t x) {
    v.push_back(static_cast<std::uint8_t>(x));
    v.push_back(static_cast<std::uint8_t>(x >> 8));
}
void putU32LE(std::vector<std::uint8_t>& v, std::uint32_t x) {
    putU16LE(v, x & 0xFFFF);
    putU16LE(v, x >> 16);
}

void pngChunk(std::vector<std::uint8_t>& out, const char type[4], const std::vector<std::uint8_t>& data) {
    putU32BE(out, static_cast<std::uint32_t>(data.size()));
    const std::size_t start = out.size();
    out.insert(out.end(), type, type + 4);
    out.insert(out.end(), data.begin(), data.end());
    putU32BE(out, Crc32::compute(out.data() + start, out.size() - start));
}

inline int paeth(int a, int b, int c) {
    const int p = a + b - c;
    const int pa = std::abs(p - a), pb = std::abs(p - b), pc = std::abs(p - c);
    if (pa <= pb && pa <= pc) return a;
    if (pb <= pc) return b;
    return c;
}

} // namespace

std::vector<std::uint8_t> encodePNG(const render::Image8& img, int level) {
    const int w = img.width, h = img.height, bpp = 3;
    const std::size_t stride = static_cast<std::size_t>(w) * bpp;

    // Adaptive filtering: per row choose the filter with the smallest sum of
    // absolute (signed) residuals - the heuristic recommended by the PNG spec.
    std::vector<std::uint8_t> raw;
    raw.reserve((stride + 1) * static_cast<std::size_t>(h));
    std::vector<std::uint8_t> cand[5];
    for (auto& c : cand) c.resize(stride);
    const std::vector<std::uint8_t> zeros(stride, 0);
    for (int y = 0; y < h; ++y) {
        const std::uint8_t* cur = img.rgb.data() + static_cast<std::size_t>(y) * stride;
        const std::uint8_t* up = y > 0 ? cur - stride : zeros.data();
        for (std::size_t i = 0; i < stride; ++i) {
            const int a = i >= static_cast<std::size_t>(bpp) ? cur[i - bpp] : 0;
            const int b = up[i];
            const int c = i >= static_cast<std::size_t>(bpp) ? up[i - bpp] : 0;
            cand[0][i] = cur[i];
            cand[1][i] = static_cast<std::uint8_t>(cur[i] - a);
            cand[2][i] = static_cast<std::uint8_t>(cur[i] - b);
            cand[3][i] = static_cast<std::uint8_t>(cur[i] - ((a + b) >> 1));
            cand[4][i] = static_cast<std::uint8_t>(cur[i] - paeth(a, b, c));
        }
        int best = 0;
        long bestScore = -1;
        for (int f = 0; f < 5; ++f) {
            long s = 0;
            for (std::size_t i = 0; i < stride; ++i) s += std::abs(static_cast<int>(static_cast<std::int8_t>(cand[f][i])));
            if (bestScore < 0 || s < bestScore) { bestScore = s; best = f; }
        }
        raw.push_back(static_cast<std::uint8_t>(best));
        raw.insert(raw.end(), cand[best].begin(), cand[best].end());
    }

    std::vector<std::uint8_t> out = {0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n'};
    std::vector<std::uint8_t> ihdr;
    putU32BE(ihdr, static_cast<std::uint32_t>(w));
    putU32BE(ihdr, static_cast<std::uint32_t>(h));
    ihdr.push_back(8);  // bit depth
    ihdr.push_back(2);  // colour type: truecolour RGB
    ihdr.push_back(0);  // compression: deflate
    ihdr.push_back(0);  // filter method
    ihdr.push_back(0);  // no interlace
    pngChunk(out, "IHDR", ihdr);
    const char* sw = "Software\0SolarSim (C++ N-body Solar System simulator)";
    pngChunk(out, "tEXt", std::vector<std::uint8_t>(sw, sw + 9 + std::strlen(sw + 9)));
    DeflateOptions dopt;
    dopt.level = level;
    pngChunk(out, "IDAT", zlibCompress(raw.data(), raw.size(), dopt));
    pngChunk(out, "IEND", {});
    return out;
}

std::vector<std::uint8_t> encodeBMP(const render::Image8& img) {
    const int w = img.width, h = img.height;
    const std::uint32_t rowSize = (static_cast<std::uint32_t>(w) * 3 + 3) & ~3u;
    const std::uint32_t dataSize = rowSize * static_cast<std::uint32_t>(h);
    std::vector<std::uint8_t> out;
    out.reserve(54 + dataSize);
    out.push_back('B');
    out.push_back('M');
    putU32LE(out, 54 + dataSize); // file size
    putU32LE(out, 0);             // reserved
    putU32LE(out, 54);            // pixel data offset
    putU32LE(out, 40);            // BITMAPINFOHEADER size
    putU32LE(out, static_cast<std::uint32_t>(w));
    putU32LE(out, static_cast<std::uint32_t>(h)); // positive = bottom-up
    putU16LE(out, 1);   // planes
    putU16LE(out, 24);  // bits per pixel
    putU32LE(out, 0);   // BI_RGB
    putU32LE(out, dataSize);
    putU32LE(out, 2835); // 72 dpi
    putU32LE(out, 2835);
    putU32LE(out, 0);
    putU32LE(out, 0);
    for (int y = h - 1; y >= 0; --y) {
        for (int x = 0; x < w; ++x) {
            const std::uint8_t* p = img.pixel(x, y);
            out.push_back(p[2]);
            out.push_back(p[1]);
            out.push_back(p[0]);
        }
        for (std::uint32_t pad = static_cast<std::uint32_t>(w) * 3; pad < rowSize; ++pad) out.push_back(0);
    }
    return out;
}

std::vector<std::uint8_t> encodePPM(const render::Image8& img) {
    const std::string header = "P6\n" + std::to_string(img.width) + " " + std::to_string(img.height) + "\n255\n";
    std::vector<std::uint8_t> out(header.begin(), header.end());
    out.insert(out.end(), img.rgb.begin(), img.rgb.end());
    return out;
}

bool writeFile(const std::string& path, const std::vector<std::uint8_t>& bytes, std::string& error) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) {
        error = "cannot open '" + path + "' for writing";
        return false;
    }
    const std::size_t n = bytes.empty() ? 0 : std::fwrite(bytes.data(), 1, bytes.size(), f);
    const bool ok = n == bytes.size() && std::fclose(f) == 0;
    if (!ok) error = "write error on '" + path + "'";
    return ok;
}

bool writeImage(const std::string& path, const render::Image8& img, std::string& error) {
    switch (formatFromPath(path)) {
        case ImageFormat::BMP: return writeFile(path, encodeBMP(img), error);
        case ImageFormat::PPM: return writeFile(path, encodePPM(img), error);
        case ImageFormat::PNG:
        case ImageFormat::Unknown:
        default: return writeFile(path, encodePNG(img), error);
    }
}

} // namespace solsim::io
