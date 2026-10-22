// =============================================================================
//  SolarSim - test_io.cpp : checksums, DEFLATE (verified with an independent
//  inflater written here), GIF LZW (verified with an independent decoder),
//  PNG/BMP containers, body-file loader and CSV escaping.
// =============================================================================
#include "TestFramework.hpp"

#include "io/BodyFileLoader.hpp"
#include "io/Checksum.hpp"
#include "io/CsvWriter.hpp"
#include "io/Deflate.hpp"
#include "io/GifEncoder.hpp"
#include "io/ImageWriter.hpp"

#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <vector>

using namespace solsim;

namespace {

// ---- minimal RFC 1951 inflater (stored + fixed Huffman blocks only) ------------------------------
class BitReader {
public:
    BitReader(const std::uint8_t* d, std::size_t n) : m_d(d), m_n(n) {}
    std::uint32_t bits(int count) {
        std::uint32_t v = 0;
        for (int i = 0; i < count; ++i) v |= static_cast<std::uint32_t>(bit()) << i;
        return v;
    }
    int bit() {
        if (m_pos >= m_n) throw std::runtime_error("inflate: out of data");
        const int b = (m_d[m_pos] >> m_bit) & 1;
        if (++m_bit == 8) { m_bit = 0; ++m_pos; }
        return b;
    }
    void align() { if (m_bit) { m_bit = 0; ++m_pos; } }
    std::uint8_t byte() { align(); return m_d[m_pos++]; }

private:
    const std::uint8_t* m_d;
    std::size_t m_n, m_pos = 0;
    int m_bit = 0;
};

int fixedLitLen(BitReader& br) {
    int code = 0;
    for (int len = 1; len <= 9; ++len) {
        code = (code << 1) | br.bit();
        if (len == 7 && code <= 23) return 256 + code;
        if (len == 8 && code >= 0x30 && code <= 0xBF) return code - 0x30;
        if (len == 8 && code >= 0xC0 && code <= 0xC7) return 280 + code - 0xC0;
        if (len == 9 && code >= 0x190 && code <= 0x1FF) return 144 + code - 0x190;
    }
    throw std::runtime_error("inflate: bad literal/length code");
}

std::vector<std::uint8_t> inflateRaw(const std::uint8_t* d, std::size_t n) {
    static const int lbase[] = {3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258};
    static const int lext[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
    static const int dbase[] = {1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577};
    static const int dext[] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};
    BitReader br(d, n);
    std::vector<std::uint8_t> out;
    for (;;) {
        const int final = br.bit();
        const int type = static_cast<int>(br.bits(2));
        if (type == 0) {
            br.align();
            // Separate statements: the evaluation order of operands of '|' is unspecified.
            const int l0 = br.byte();
            const int l1 = br.byte();
            const int n0 = br.byte();
            const int n1 = br.byte();
            const int len = l0 | (l1 << 8);
            const int nlen = n0 | (n1 << 8);
            if ((len ^ 0xFFFF) != nlen) throw std::runtime_error("inflate: stored length mismatch");
            for (int i = 0; i < len; ++i) out.push_back(br.byte());
        } else if (type == 1) {
            for (;;) {
                const int sym = fixedLitLen(br);
                if (sym < 256) { out.push_back(static_cast<std::uint8_t>(sym)); continue; }
                if (sym == 256) break;
                const int li = sym - 257;
                const int len = lbase[li] + static_cast<int>(br.bits(lext[li]));
                int dcode = 0;
                for (int k = 0; k < 5; ++k) dcode = (dcode << 1) | br.bit();
                const int dist = dbase[dcode] + static_cast<int>(br.bits(dext[dcode]));
                if (dist > static_cast<int>(out.size())) throw std::runtime_error("inflate: distance too far");
                for (int k = 0; k < len; ++k) out.push_back(out[out.size() - static_cast<std::size_t>(dist)]);
            }
        } else {
            throw std::runtime_error("inflate: unsupported block type");
        }
        if (final) break;
    }
    return out;
}

// ---- independent GIF LZW decoder --------------------------------------------------------------------
std::vector<std::uint8_t> gifLzwDecode(const std::vector<std::uint8_t>& data, int minCode) {
    const int clear = 1 << minCode, eoi = clear + 1;
    int codeSize = minCode + 1;
    std::vector<std::vector<std::uint8_t>> dict;
    auto reset = [&] {
        dict.clear();
        for (int i = 0; i < clear; ++i) dict.push_back({static_cast<std::uint8_t>(i)});
        dict.push_back({});
        dict.push_back({});
        codeSize = minCode + 1;
    };
    reset();
    std::vector<std::uint8_t> out;
    std::size_t bitPos = 0;
    auto read = [&](int n) {
        int v = 0;
        for (int i = 0; i < n; ++i, ++bitPos) {
            if (bitPos / 8 >= data.size()) throw std::runtime_error("lzw: out of data");
            v |= ((data[bitPos / 8] >> (bitPos % 8)) & 1) << i;
        }
        return v;
    };
    int prev = -1;
    for (;;) {
        const int code = read(codeSize);
        if (code == clear) { reset(); prev = -1; continue; }
        if (code == eoi) break;
        std::vector<std::uint8_t> entry;
        if (prev < 0) {
            entry = dict.at(static_cast<std::size_t>(code));
        } else {
            if (code < static_cast<int>(dict.size())) {
                entry = dict[static_cast<std::size_t>(code)];
            } else if (code == static_cast<int>(dict.size())) {
                entry = dict[static_cast<std::size_t>(prev)];
                entry.push_back(entry.front());
            } else {
                throw std::runtime_error("lzw: invalid code");
            }
            if (dict.size() < 4096) {
                auto ne = dict[static_cast<std::size_t>(prev)];
                ne.push_back(entry.front());
                dict.push_back(ne);
                if (static_cast<int>(dict.size()) == (1 << codeSize) && codeSize < 12) ++codeSize;
            }
        }
        out.insert(out.end(), entry.begin(), entry.end());
        prev = code;
    }
    return out;
}

} // namespace

TEST(IO, Checksums) {
    const char* s = "123456789";
    CHECK_EQ(io::Crc32::compute(reinterpret_cast<const std::uint8_t*>(s), 9), 0xCBF43926u);
    CHECK_EQ(io::Adler32::compute(reinterpret_cast<const std::uint8_t*>("Wikipedia"), 9), 0x11E60398u);
    // Incremental CRC equals one-shot CRC.
    io::Crc32 c;
    c.update(reinterpret_cast<const std::uint8_t*>(s), 4);
    c.update(reinterpret_cast<const std::uint8_t*>(s) + 4, 5);
    CHECK_EQ(c.value(), 0xCBF43926u);
}

TEST(IO, DeflateRoundTrip) {
    std::vector<std::vector<std::uint8_t>> inputs;
    inputs.push_back({});
    inputs.push_back({42});
    std::vector<std::uint8_t> text;
    const std::string phrase = "The Solar System formed 4.6 billion years ago. ";
    for (int i = 0; i < 300; ++i) text.insert(text.end(), phrase.begin(), phrase.end());
    inputs.push_back(text);
    std::vector<std::uint8_t> noise(70000);
    std::uint32_t x = 12345;
    for (auto& b : noise) { x = x * 1103515245u + 12345u; b = static_cast<std::uint8_t>(x >> 24); }
    inputs.push_back(noise);
    std::vector<std::uint8_t> runs(100000, 7);
    inputs.push_back(runs);
    for (const auto& in : inputs) {
        for (int level : {0, 1, 6, 9}) {
            io::DeflateOptions o;
            o.level = level;
            const auto z = io::deflateRaw(in.data(), in.size(), o);
            const auto back = inflateRaw(z.data(), z.size());
            CHECK(back == in);
        }
    }
    const auto z = io::zlibCompress(text.data(), text.size());
    CHECK_EQ(static_cast<int>(z[0]), 0x78);
    CHECK_EQ(((z[0] << 8) | z[1]) % 31, 0);
    CHECK(z.size() < text.size() / 10);
}

TEST(IO, GifLzwRoundTrip) {
    std::vector<std::vector<std::uint8_t>> inputs;
    std::vector<std::uint8_t> a(50000);
    for (std::size_t i = 0; i < a.size(); ++i) a[i] = static_cast<std::uint8_t>((i * i) % 251);
    inputs.push_back(a); // forces dictionary resets
    inputs.push_back(std::vector<std::uint8_t>(9000, 3));
    inputs.push_back({1, 2, 3});
    for (const auto& in : inputs) {
        const auto enc = io::gifLzwEncode(in, 8);
        CHECK(gifLzwDecode(enc, 8) == in);
    }
}

TEST(IO, PngContainer) {
    render::Image8 img;
    img.resize(37, 21);
    for (int y = 0; y < 21; ++y)
        for (int x = 0; x < 37; ++x) {
            auto* p = img.pixel(x, y);
            p[0] = static_cast<std::uint8_t>(x * 7);
            p[1] = static_cast<std::uint8_t>(y * 11);
            p[2] = static_cast<std::uint8_t>((x ^ y) * 3);
        }
    const auto png = io::encodePNG(img);
    CHECK(png.size() > 60);
    CHECK_EQ(static_cast<int>(png[0]), 0x89);
    CHECK_EQ(static_cast<int>(png[1]), 'P');
    // Walk the chunks, verify CRCs, collect IDAT and decode the pixels again.
    std::size_t pos = 8;
    std::vector<std::uint8_t> idat;
    bool sawEnd = false;
    while (pos + 12 <= png.size()) {
        const std::uint32_t len = (png[pos] << 24) | (png[pos + 1] << 16) | (png[pos + 2] << 8) | png[pos + 3];
        const std::string type(reinterpret_cast<const char*>(&png[pos + 4]), 4);
        const std::uint32_t crc = (png[pos + 8 + len] << 24) | (png[pos + 9 + len] << 16) | (png[pos + 10 + len] << 8) | png[pos + 11 + len];
        CHECK_EQ(io::Crc32::compute(&png[pos + 4], len + 4), crc);
        if (type == "IDAT") idat.insert(idat.end(), png.begin() + static_cast<std::ptrdiff_t>(pos + 8), png.begin() + static_cast<std::ptrdiff_t>(pos + 8 + len));
        if (type == "IEND") sawEnd = true;
        pos += 12 + len;
    }
    CHECK(sawEnd);
    const auto raw = inflateRaw(idat.data() + 2, idat.size() - 6);
    CHECK_EQ(raw.size(), static_cast<std::size_t>(21 * (37 * 3 + 1)));
    // Undo the scanline filters.
    std::vector<std::uint8_t> prev(37 * 3, 0), cur(37 * 3);
    bool same = true;
    for (int y = 0; y < 21; ++y) {
        const std::uint8_t f = raw[static_cast<std::size_t>(y) * (37 * 3 + 1)];
        const std::uint8_t* s = &raw[static_cast<std::size_t>(y) * (37 * 3 + 1) + 1];
        for (int i = 0; i < 37 * 3; ++i) {
            const int a = i >= 3 ? cur[static_cast<std::size_t>(i - 3)] : 0, b = prev[static_cast<std::size_t>(i)],
                      c = i >= 3 ? prev[static_cast<std::size_t>(i - 3)] : 0;
            int pred = 0;
            if (f == 1) pred = a;
            else if (f == 2) pred = b;
            else if (f == 3) pred = (a + b) / 2;
            else if (f == 4) {
                const int p = a + b - c, pa = std::abs(p - a), pb = std::abs(p - b), pc = std::abs(p - c);
                pred = (pa <= pb && pa <= pc) ? a : (pb <= pc ? b : c);
            }
            cur[static_cast<std::size_t>(i)] = static_cast<std::uint8_t>(s[i] + pred);
        }
        for (int i = 0; i < 37 * 3; ++i)
            if (cur[static_cast<std::size_t>(i)] != img.rgb[static_cast<std::size_t>(y) * 37 * 3 + static_cast<std::size_t>(i)]) same = false;
        prev = cur;
    }
    CHECK(same);
}

TEST(IO, BmpLayout) {
    render::Image8 img;
    img.resize(3, 2); // row size 9 -> padded to 12
    const auto bmp = io::encodeBMP(img);
    CHECK_EQ(bmp.size(), static_cast<std::size_t>(54 + 12 * 2));
    CHECK_EQ(static_cast<int>(bmp[0]), 'B');
    CHECK_EQ(static_cast<int>(bmp[28]), 24);
}

TEST(IO, MedianCutPalette) {
    render::Image8 img;
    img.resize(64, 64);
    for (int y = 0; y < 64; ++y)
        for (int x = 0; x < 64; ++x) {
            auto* p = img.pixel(x, y);
            p[0] = static_cast<std::uint8_t>(x * 4);
            p[1] = static_cast<std::uint8_t>(y * 4);
            p[2] = 128;
        }
    io::MedianCutQuantizer q;
    const auto pal = q.build(img, 256);
    CHECK(pal.colors.size() <= 256);
    CHECK(pal.colors.size() >= 200);
    std::vector<std::uint8_t> idx;
    q.map(img, pal, false, idx);
    double err = 0;
    for (std::size_t i = 0; i < idx.size(); ++i)
        for (int c = 0; c < 3; ++c) err += std::fabs(double(pal.colors[idx[i]][static_cast<std::size_t>(c)]) - img.rgb[i * 3 + static_cast<std::size_t>(c)]);
    CHECK(err / (idx.size() * 3.0) < 6.0);
}

TEST(IO, BodyFileLoader) {
    const auto r = io::loadBodyString(
        "[body]\nname = TestComet\nkind = comet\nq = 0.5\ne = 0.9\ni = 10\nnode = 20\nargp = 30\ntp = 2030-01-01\n"
        "[body]\nname = Rock\nkind = asteroid\na = 2.5\ne = 0.1\nM = 45\nepoch = J2000\nmass_kg = 1e18\ncolor = #FF8800\n");
    CHECK(r.ok());
    CHECK_EQ(r.bodies.size(), static_cast<std::size_t>(2));
    CHECK_NEAR(r.bodies[0].q, 0.5, 0);
    CHECK_NEAR(r.bodies[1].q, 2.25, 1e-12);
    CHECK_EQ(r.bodies[1].color, 0xFF8800u);
    CHECK(r.bodies[1].gmKm3s2 > 0.0);
    const auto bad = io::loadBodyString("[body]\nkind = comet\n");
    CHECK(!bad.ok());
}

TEST(IO, CsvEscaping) {
    CHECK_EQ(io::CsvWriter::escape("plain"), std::string("plain"));
    CHECK_EQ(io::CsvWriter::escape("a,b"), std::string("\"a,b\""));
    CHECK_EQ(io::CsvWriter::escape("say \"hi\""), std::string("\"say \"\"hi\"\"\""));
}
