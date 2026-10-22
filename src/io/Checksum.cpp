// =============================================================================
//  SolarSim - Checksum.cpp
// =============================================================================
#include "Checksum.hpp"

#include <array>

namespace solsim::io {

namespace {
std::array<std::uint32_t, 256> makeCrcTable() {
    std::array<std::uint32_t, 256> t{};
    for (std::uint32_t n = 0; n < 256; ++n) {
        std::uint32_t c = n;
        for (int k = 0; k < 8; ++k) c = (c & 1u) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
        t[n] = c;
    }
    return t;
}
const std::array<std::uint32_t, 256>& crcTable() {
    static const std::array<std::uint32_t, 256> table = makeCrcTable();
    return table;
}
} // namespace

void Crc32::update(const std::uint8_t* data, std::size_t n) {
    const auto& t = crcTable();
    std::uint32_t c = m_crc;
    for (std::size_t i = 0; i < n; ++i) c = t[(c ^ data[i]) & 0xFFu] ^ (c >> 8);
    m_crc = c;
}

std::uint32_t Crc32::compute(const std::uint8_t* data, std::size_t n) {
    Crc32 c;
    c.update(data, n);
    return c.value();
}

void Adler32::update(const std::uint8_t* data, std::size_t n) {
    constexpr std::uint32_t MOD = 65521u;
    // Process in blocks of 5552 bytes - the largest n for which b cannot overflow.
    while (n > 0) {
        const std::size_t block = n < 5552 ? n : 5552;
        for (std::size_t i = 0; i < block; ++i) {
            m_a += data[i];
            m_b += m_a;
        }
        m_a %= MOD;
        m_b %= MOD;
        data += block;
        n -= block;
    }
}

std::uint32_t Adler32::compute(const std::uint8_t* data, std::size_t n) {
    Adler32 a;
    a.update(data, n);
    return a.value();
}

} // namespace solsim::io
