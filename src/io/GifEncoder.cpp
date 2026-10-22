// =============================================================================
//  SolarSim - GifEncoder.cpp
// =============================================================================
#include "GifEncoder.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace solsim::io {

namespace {

inline int bin15(int r, int g, int b) { return ((r >> 3) << 10) | ((g >> 3) << 5) | (b >> 3); }

struct HistEntry {
    std::uint16_t bin;
    std::uint32_t count;
    std::uint64_t sr, sg, sb; // exact 8-bit sums for accurate averages
};

struct Box {
    std::size_t begin, end; // range in the entry array
    std::uint64_t total;
    int lo[3], hi[3];
    int longestAxis() const {
        int best = 0;
        for (int a = 1; a < 3; ++a)
            if (hi[a] - lo[a] > hi[best] - lo[best]) best = a;
        return best;
    }
    int range() const {
        const int a = longestAxis();
        return hi[a] - lo[a];
    }
};

inline int channelOf(std::uint16_t bin, int axis) {
    return axis == 0 ? (bin >> 10) & 31 : (axis == 1 ? (bin >> 5) & 31 : bin & 31);
}

void computeBounds(Box& b, const std::vector<HistEntry>& e) {
    for (int a = 0; a < 3; ++a) { b.lo[a] = 31; b.hi[a] = 0; }
    b.total = 0;
    for (std::size_t i = b.begin; i < b.end; ++i) {
        for (int a = 0; a < 3; ++a) {
            const int c = channelOf(e[i].bin, a);
            b.lo[a] = std::min(b.lo[a], c);
            b.hi[a] = std::max(b.hi[a], c);
        }
        b.total += e[i].count;
    }
}

} // namespace

Palette MedianCutQuantizer::build(const render::Image8& img, int maxColors) {
    std::vector<HistEntry> hist(32768);
    for (std::size_t i = 0; i < hist.size(); ++i) hist[i] = {static_cast<std::uint16_t>(i), 0, 0, 0, 0};
    const std::size_t npx = static_cast<std::size_t>(img.width) * img.height;
    for (std::size_t p = 0; p < npx; ++p) {
        const std::uint8_t* c = &img.rgb[p * 3];
        HistEntry& h = hist[static_cast<std::size_t>(bin15(c[0], c[1], c[2]))];
        ++h.count;
        h.sr += c[0];
        h.sg += c[1];
        h.sb += c[2];
    }
    std::vector<HistEntry> entries;
    for (const auto& h : hist)
        if (h.count) entries.push_back(h);

    std::vector<Box> boxes;
    if (!entries.empty()) {
        Box root{0, entries.size(), 0, {0, 0, 0}, {0, 0, 0}};
        computeBounds(root, entries);
        boxes.push_back(root);
    }
    while (static_cast<int>(boxes.size()) < maxColors) {
        // Split the box with the largest (range x population) score.
        int pick = -1;
        double bestScore = 0.0;
        for (std::size_t k = 0; k < boxes.size(); ++k) {
            const Box& b = boxes[k];
            if (b.end - b.begin < 2 || b.range() == 0) continue;
            const double score = double(b.range()) * std::sqrt(double(b.total));
            if (score > bestScore) { bestScore = score; pick = static_cast<int>(k); }
        }
        if (pick < 0) break;
        Box b = boxes[static_cast<std::size_t>(pick)];
        const int axis = b.longestAxis();
        std::sort(entries.begin() + static_cast<std::ptrdiff_t>(b.begin), entries.begin() + static_cast<std::ptrdiff_t>(b.end),
                  [axis](const HistEntry& x, const HistEntry& y) { return channelOf(x.bin, axis) < channelOf(y.bin, axis); });
        // Weighted median split.
        std::uint64_t acc = 0;
        std::size_t split = b.begin + 1;
        for (std::size_t i = b.begin; i < b.end; ++i) {
            acc += entries[i].count;
            if (acc * 2 >= b.total) {
                split = i + 1;
                break;
            }
        }
        split = std::clamp(split, b.begin + 1, b.end - 1);
        Box left{b.begin, split, 0, {0, 0, 0}, {0, 0, 0}};
        Box right{split, b.end, 0, {0, 0, 0}, {0, 0, 0}};
        computeBounds(left, entries);
        computeBounds(right, entries);
        boxes[static_cast<std::size_t>(pick)] = left;
        boxes.push_back(right);
    }

    Palette pal;
    m_binToIndex.assign(32768, 0xFFFF);
    for (std::size_t k = 0; k < boxes.size(); ++k) {
        std::uint64_t r = 0, g = 0, bl = 0, n = 0;
        for (std::size_t i = boxes[k].begin; i < boxes[k].end; ++i) {
            r += entries[i].sr;
            g += entries[i].sg;
            bl += entries[i].sb;
            n += entries[i].count;
            m_binToIndex[entries[i].bin] = static_cast<std::uint16_t>(k);
        }
        if (n == 0) n = 1;
        pal.colors.push_back({static_cast<std::uint8_t>(r / n), static_cast<std::uint8_t>(g / n),
                              static_cast<std::uint8_t>(bl / n)});
    }
    if (pal.colors.empty()) pal.colors.push_back({0, 0, 0});
    return pal;
}

std::uint8_t MedianCutQuantizer::nearest(const Palette& pal, int r, int g, int b) {
    const int key = bin15(r, g, b);
    if (m_binToIndex.size() == 32768 && m_binToIndex[static_cast<std::size_t>(key)] != 0xFFFF)
        return static_cast<std::uint8_t>(m_binToIndex[static_cast<std::size_t>(key)]);
    int best = 0;
    long bestD = -1;
    for (std::size_t k = 0; k < pal.colors.size(); ++k) {
        const long dr = r - pal.colors[k][0], dg = g - pal.colors[k][1], db = b - pal.colors[k][2];
        // Perceptual-ish weights (green most visible).
        const long d = 2 * dr * dr + 4 * dg * dg + 3 * db * db;
        if (bestD < 0 || d < bestD) { bestD = d; best = static_cast<int>(k); }
    }
    if (m_binToIndex.size() == 32768) m_binToIndex[static_cast<std::size_t>(key)] = static_cast<std::uint16_t>(best);
    return static_cast<std::uint8_t>(best);
}

void MedianCutQuantizer::map(const render::Image8& img, const Palette& pal, bool dither,
                             std::vector<std::uint8_t>& indices) {
    const int w = img.width, h = img.height;
    indices.resize(static_cast<std::size_t>(w) * h);
    if (!dither) {
        for (std::size_t p = 0; p < indices.size(); ++p) {
            const std::uint8_t* c = &img.rgb[p * 3];
            indices[p] = nearest(pal, c[0], c[1], c[2]);
        }
        return;
    }
    // Floyd-Steinberg error diffusion (serpentine scan).
    std::vector<int> err0(static_cast<std::size_t>(w + 2) * 3, 0), err1(static_cast<std::size_t>(w + 2) * 3, 0);
    for (int y = 0; y < h; ++y) {
        std::fill(err1.begin(), err1.end(), 0);
        const bool rev = (y & 1) != 0;
        for (int i = 0; i < w; ++i) {
            const int x = rev ? w - 1 - i : i;
            const std::uint8_t* c = img.pixel(x, y);
            int v[3];
            for (int a = 0; a < 3; ++a) v[a] = std::clamp(c[a] + err0[static_cast<std::size_t>(x + 1) * 3 + a] / 16, 0, 255);
            const std::uint8_t idx = nearest(pal, v[0], v[1], v[2]);
            indices[static_cast<std::size_t>(y) * w + x] = idx;
            const int dx = rev ? -1 : 1;
            for (int a = 0; a < 3; ++a) {
                const int e = v[a] - pal.colors[idx][static_cast<std::size_t>(a)];
                const auto at = [&](int xx) { return static_cast<std::size_t>(xx + 1) * 3 + static_cast<std::size_t>(a); };
                if (x + dx >= 0 && x + dx < w) err0[at(x + dx)] += e * 7;
                if (x - dx >= 0 && x - dx < w) err1[at(x - dx)] += e * 3;
                err1[at(x)] += e * 5;
                if (x + dx >= 0 && x + dx < w) err1[at(x + dx)] += e * 1;
            }
        }
        std::swap(err0, err1);
    }
}

// -----------------------------------------------------------------------------
std::vector<std::uint8_t> gifLzwEncode(const std::vector<std::uint8_t>& indices, int minCodeSize) {
    const int clearCode = 1 << minCodeSize;
    const int eoiCode = clearCode + 1;
    int codeSize = minCodeSize + 1;
    int nextCode = eoiCode + 1;
    // Dictionary as a dense trie: child[code * 256 + byte] -> code (0 = none).
    std::vector<std::uint16_t> child(4096 * 256, 0);

    std::vector<std::uint8_t> out;
    std::uint32_t acc = 0;
    int bits = 0;
    auto emit = [&](int code) {
        acc |= static_cast<std::uint32_t>(code) << bits;
        bits += codeSize;
        while (bits >= 8) {
            out.push_back(static_cast<std::uint8_t>(acc & 0xFF));
            acc >>= 8;
            bits -= 8;
        }
    };
    auto resetDictionary = [&]() {
        std::fill(child.begin(), child.end(), 0);
        codeSize = minCodeSize + 1;
        nextCode = eoiCode + 1;
    };

    emit(clearCode);
    if (indices.empty()) {
        emit(eoiCode);
        if (bits > 0) out.push_back(static_cast<std::uint8_t>(acc & 0xFF));
        return out;
    }
    int cur = indices[0];
    for (std::size_t i = 1; i < indices.size(); ++i) {
        const int k = indices[i];
        const std::uint16_t nxt = child[static_cast<std::size_t>(cur) * 256 + static_cast<std::size_t>(k)];
        if (nxt != 0) {
            cur = nxt;
            continue;
        }
        emit(cur);
        if (nextCode < 4096) {
            child[static_cast<std::size_t>(cur) * 256 + static_cast<std::size_t>(k)] = static_cast<std::uint16_t>(nextCode++);
            if (nextCode > (1 << codeSize) && codeSize < 12) ++codeSize;
        } else {
            emit(clearCode);
            resetDictionary();
        }
        cur = k;
    }
    emit(cur);
    emit(eoiCode);
    if (bits > 0) out.push_back(static_cast<std::uint8_t>(acc & 0xFF));
    return out;
}

// -----------------------------------------------------------------------------
GifEncoder::~GifEncoder() {
    if (m_file) close();
}

void GifEncoder::put(const void* data, std::size_t n) {
    if (!m_file) return;
    if (std::fwrite(data, 1, n, m_file) != n) m_error = "write error";
    m_bytes += n;
}

bool GifEncoder::open(const std::string& path, int width, int height, int loopCount) {
    if (width <= 0 || height <= 0 || width > 65535 || height > 65535) {
        m_error = "invalid GIF dimensions";
        return false;
    }
    m_file = std::fopen(path.c_str(), "wb");
    if (!m_file) {
        m_error = "cannot open '" + path + "' for writing";
        return false;
    }
    m_w = width;
    m_h = height;
    m_frames = 0;
    m_bytes = 0;
    put("GIF89a", 6);
    putU16(static_cast<std::uint16_t>(width));
    putU16(static_cast<std::uint16_t>(height));
    putByte(0x70); // no global colour table, 8-bit colour resolution
    putByte(0);    // background colour index
    putByte(0);    // pixel aspect ratio
    // NETSCAPE2.0 looping extension
    const std::uint8_t app[] = {0x21, 0xFF, 0x0B, 'N', 'E', 'T', 'S', 'C', 'A', 'P', 'E', '2', '.', '0', 0x03, 0x01};
    put(app, sizeof app);
    putU16(static_cast<std::uint16_t>(loopCount));
    putByte(0);
    return m_error.empty();
}

bool GifEncoder::addFrame(const render::Image8& frame, int delayCs, bool dither) {
    if (!m_file) return false;
    if (frame.width != m_w || frame.height != m_h) {
        m_error = "frame size mismatch";
        return false;
    }
    const Palette pal = m_quant.build(frame, 256);
    m_quant.map(frame, pal, dither, m_indices);

    // Graphic Control Extension: disposal = 1 (do not dispose), no transparency.
    const std::uint8_t gce[] = {0x21, 0xF9, 0x04, 0x04};
    put(gce, sizeof gce);
    putU16(static_cast<std::uint16_t>(std::clamp(delayCs, 0, 65535)));
    putByte(0);
    putByte(0);

    // Image descriptor with a 256-entry local colour table.
    putByte(0x2C);
    putU16(0);
    putU16(0);
    putU16(static_cast<std::uint16_t>(m_w));
    putU16(static_cast<std::uint16_t>(m_h));
    putByte(0x87); // local table present, size 2^(7+1) = 256
    for (std::size_t k = 0; k < 256; ++k) {
        const auto c = k < pal.colors.size() ? pal.colors[k] : std::array<std::uint8_t, 3>{0, 0, 0};
        put(c.data(), 3);
    }

    const int minCodeSize = 8;
    putByte(static_cast<std::uint8_t>(minCodeSize));
    const std::vector<std::uint8_t> lzw = gifLzwEncode(m_indices, minCodeSize);
    for (std::size_t off = 0; off < lzw.size(); off += 255) {
        const std::size_t n = std::min<std::size_t>(255, lzw.size() - off);
        putByte(static_cast<std::uint8_t>(n));
        put(lzw.data() + off, n);
    }
    putByte(0); // block terminator
    ++m_frames;
    return m_error.empty();
}

bool GifEncoder::close() {
    if (!m_file) return false;
    putByte(0x3B); // trailer
    const bool ok = std::fclose(m_file) == 0 && m_error.empty();
    m_file = nullptr;
    return ok;
}

} // namespace solsim::io
