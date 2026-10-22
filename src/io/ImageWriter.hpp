// =============================================================================
//  SolarSim - ImageWriter.hpp
//  Still image encoders: PNG (own DEFLATE, adaptive scanline filters),
//  BMP (24-bit, bottom-up) and binary PPM (P6).
// =============================================================================
#pragma once

#include "../render/Image.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace solsim::io {

enum class ImageFormat { PNG, BMP, PPM, Unknown };

ImageFormat formatFromPath(const std::string& path);

std::vector<std::uint8_t> encodePNG(const render::Image8& img, int compressionLevel = 6);
std::vector<std::uint8_t> encodeBMP(const render::Image8& img);
std::vector<std::uint8_t> encodePPM(const render::Image8& img);

// Writes using the format implied by the extension (PNG if unknown).
// Returns false and fills `error` on failure.
bool writeImage(const std::string& path, const render::Image8& img, std::string& error);

bool writeFile(const std::string& path, const std::vector<std::uint8_t>& bytes, std::string& error);

} // namespace solsim::io
