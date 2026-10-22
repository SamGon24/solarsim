// =============================================================================
//  SolarSim - Checksum.hpp
//  CRC-32 (IEEE 802.3, reflected polynomial 0xEDB88320) used by PNG chunks and
//  Adler-32 used by the zlib container.
// =============================================================================
#pragma once

#include <cstddef>
#include <cstdint>

namespace solsim::io {

class Crc32 {
public:
    Crc32() = default;
    void update(const std::uint8_t* data, std::size_t n);
    std::uint32_t value() const { return m_crc ^ 0xFFFFFFFFu; }
    static std::uint32_t compute(const std::uint8_t* data, std::size_t n);

private:
    std::uint32_t m_crc = 0xFFFFFFFFu;
};

class Adler32 {
public:
    void update(const std::uint8_t* data, std::size_t n);
    std::uint32_t value() const { return (m_b << 16) | m_a; }
    static std::uint32_t compute(const std::uint8_t* data, std::size_t n);

private:
    std::uint32_t m_a = 1, m_b = 0;
};

} // namespace solsim::io
