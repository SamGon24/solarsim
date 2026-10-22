// =============================================================================
//  SolarSim - Deflate.hpp
//  A from-scratch DEFLATE (RFC 1951) compressor with a zlib (RFC 1950)
//  wrapper, used by the PNG writer:
//    * LZ77 with hash chains over a 32 KiB window (3..258 byte matches)
//    * fixed Huffman code blocks (BTYPE = 01)
//    * optional stored blocks (level 0)
// =============================================================================
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace solsim::io {

class BitWriter {
public:
    void writeBits(std::uint32_t value, int count); // LSB first
    void writeHuffman(std::uint32_t code, int length); // code given MSB first
    void alignToByte();
    std::vector<std::uint8_t>& bytes() { return m_out; }

private:
    std::vector<std::uint8_t> m_out;
    std::uint32_t m_acc = 0;
    int m_bits = 0;
};

struct DeflateOptions {
    int level = 6;           // 0 = stored, 1..9 = chain length / effort
};

// Raw DEFLATE stream.
std::vector<std::uint8_t> deflateRaw(const std::uint8_t* data, std::size_t n, const DeflateOptions& opt = {});
// zlib stream: 2-byte header + DEFLATE + Adler-32.
std::vector<std::uint8_t> zlibCompress(const std::uint8_t* data, std::size_t n, const DeflateOptions& opt = {});

} // namespace solsim::io
