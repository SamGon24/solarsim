// =============================================================================
//  SolarSim - Deflate.cpp
// =============================================================================
#include "Deflate.hpp"

#include "Checksum.hpp"

#include <algorithm>
#include <cstring>

namespace solsim::io {

// -----------------------------------------------------------------------------
void BitWriter::writeBits(std::uint32_t value, int count) {
    for (int i = 0; i < count; ++i) {
        m_acc |= ((value >> i) & 1u) << m_bits;
        if (++m_bits == 8) {
            m_out.push_back(static_cast<std::uint8_t>(m_acc));
            m_acc = 0;
            m_bits = 0;
        }
    }
}

void BitWriter::writeHuffman(std::uint32_t code, int length) {
    // Huffman codes are packed starting with the most significant bit.
    for (int i = length - 1; i >= 0; --i) writeBits((code >> i) & 1u, 1);
}

void BitWriter::alignToByte() {
    if (m_bits > 0) {
        m_out.push_back(static_cast<std::uint8_t>(m_acc));
        m_acc = 0;
        m_bits = 0;
    }
}

// -----------------------------------------------------------------------------
namespace {

constexpr int kWindow = 32768;
constexpr int kMinMatch = 3;
constexpr int kMaxMatch = 258;
constexpr int kHashBits = 15;
constexpr int kHashSize = 1 << kHashBits;

const std::uint16_t kLengthBase[29] = {3,  4,  5,  6,  7,  8,  9,  10, 11,  13,  15,  17,  19,  23, 27,
                                       31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258};
const std::uint8_t kLengthExtra[29] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
const std::uint16_t kDistBase[30] = {1,   2,   3,   4,    5,    7,    9,    13,   17,   25,    33,    49,    65,    97,    129,
                                     193, 257, 385, 513,  769,  1025, 1537, 2049, 3073, 4097,  6145,  8193,  12289, 16385, 24577};
const std::uint8_t kDistExtra[30] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};

// Fixed Huffman literal/length codes (RFC 1951 section 3.2.6).
void writeFixedLiteral(BitWriter& bw, int sym) {
    if (sym < 144) bw.writeHuffman(0x30u + static_cast<std::uint32_t>(sym), 8);
    else if (sym < 256) bw.writeHuffman(0x190u + static_cast<std::uint32_t>(sym - 144), 9);
    else if (sym < 280) bw.writeHuffman(static_cast<std::uint32_t>(sym - 256), 7);
    else bw.writeHuffman(0xC0u + static_cast<std::uint32_t>(sym - 280), 8);
}

void writeMatch(BitWriter& bw, int length, int distance) {
    int li = 28;
    while (li > 0 && kLengthBase[li] > length) --li;
    writeFixedLiteral(bw, 257 + li);
    bw.writeBits(static_cast<std::uint32_t>(length - kLengthBase[li]), kLengthExtra[li]);
    int di = 29;
    while (di > 0 && kDistBase[di] > distance) --di;
    bw.writeHuffman(static_cast<std::uint32_t>(di), 5);
    bw.writeBits(static_cast<std::uint32_t>(distance - kDistBase[di]), kDistExtra[di]);
}

inline std::uint32_t hash3(const std::uint8_t* p) {
    return ((static_cast<std::uint32_t>(p[0]) << 10) ^ (static_cast<std::uint32_t>(p[1]) << 5) ^ p[2]) & (kHashSize - 1);
}

std::vector<std::uint8_t> storedBlocks(const std::uint8_t* data, std::size_t n) {
    std::vector<std::uint8_t> out;
    std::size_t pos = 0;
    do {
        const std::size_t len = std::min<std::size_t>(65535, n - pos);
        const bool final = pos + len >= n;
        out.push_back(final ? 1 : 0); // BFINAL + BTYPE=00, already byte aligned
        out.push_back(static_cast<std::uint8_t>(len & 0xFF));
        out.push_back(static_cast<std::uint8_t>(len >> 8));
        out.push_back(static_cast<std::uint8_t>(~len & 0xFF));
        out.push_back(static_cast<std::uint8_t>((~len >> 8) & 0xFF));
        out.insert(out.end(), data + pos, data + pos + len);
        pos += len;
    } while (pos < n);
    return out;
}

} // namespace

std::vector<std::uint8_t> deflateRaw(const std::uint8_t* data, std::size_t n, const DeflateOptions& opt) {
    if (opt.level <= 0) return storedBlocks(data, n);
    const int maxChain = std::clamp(opt.level, 1, 9) * 16;

    BitWriter bw;
    bw.writeBits(1, 1); // BFINAL
    bw.writeBits(1, 2); // BTYPE = 01 (fixed Huffman)

    std::vector<int> head(kHashSize, -1);
    std::vector<int> prev(kWindow, -1);
    std::size_t i = 0;
    auto insert = [&](std::size_t pos) {
        if (pos + kMinMatch > n) return;
        const std::uint32_t h = hash3(data + pos);
        prev[pos & (kWindow - 1)] = head[h];
        head[h] = static_cast<int>(pos);
    };

    while (i < n) {
        int bestLen = 0, bestDist = 0;
        if (i + kMinMatch <= n) {
            const std::uint32_t h = hash3(data + i);
            int cand = head[h];
            int chain = maxChain;
            const std::size_t maxLen = std::min<std::size_t>(kMaxMatch, n - i);
            while (cand >= 0 && chain-- > 0) {
                const std::size_t dist = i - static_cast<std::size_t>(cand);
                if (dist == 0 || dist > static_cast<std::size_t>(kWindow - 1)) break;
                const std::uint8_t* a = data + cand;
                const std::uint8_t* b = data + i;
                if (a[bestLen] == b[bestLen]) {
                    std::size_t len = 0;
                    while (len < maxLen && a[len] == b[len]) ++len;
                    if (static_cast<int>(len) > bestLen) {
                        bestLen = static_cast<int>(len);
                        bestDist = static_cast<int>(dist);
                        if (len == maxLen) break;
                    }
                }
                cand = prev[static_cast<std::size_t>(cand) & (kWindow - 1)];
            }
        }
        if (bestLen >= kMinMatch) {
            writeMatch(bw, bestLen, bestDist);
            for (int k = 0; k < bestLen; ++k) insert(i + static_cast<std::size_t>(k));
            i += static_cast<std::size_t>(bestLen);
        } else {
            writeFixedLiteral(bw, data[i]);
            insert(i);
            ++i;
        }
    }
    writeFixedLiteral(bw, 256); // end of block
    bw.alignToByte();
    return std::move(bw.bytes());
}

std::vector<std::uint8_t> zlibCompress(const std::uint8_t* data, std::size_t n, const DeflateOptions& opt) {
    std::vector<std::uint8_t> out;
    out.push_back(0x78); // CMF: deflate, 32K window
    out.push_back(0x9C); // FLG: default level, check bits make (CMF*256+FLG) % 31 == 0
    const std::vector<std::uint8_t> body = deflateRaw(data, n, opt);
    out.insert(out.end(), body.begin(), body.end());
    const std::uint32_t adler = Adler32::compute(data, n);
    out.push_back(static_cast<std::uint8_t>(adler >> 24));
    out.push_back(static_cast<std::uint8_t>(adler >> 16));
    out.push_back(static_cast<std::uint8_t>(adler >> 8));
    out.push_back(static_cast<std::uint8_t>(adler));
    return out;
}

} // namespace solsim::io
