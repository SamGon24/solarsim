// =============================================================================
//  SolarSim - PostProcess.hpp
//  HDR post-processing: bloom (bright-pass, downsample, separable gaussian
//  blur, additive upsample), exposure, filmic tone mapping (ACES fit by
//  K. Narkowicz) and conversion to 8-bit sRGB.
// =============================================================================
#pragma once

#include "Image.hpp"

namespace solsim::render {

struct PostSettings {
    bool bloom = true;
    float bloomThreshold = 0.9f;
    float bloomStrength = 0.55f;
    float bloomRadius = 0.012f;   // gaussian sigma as a fraction of image height
    float exposure = 1.0f;
    bool vignette = true;
};

void applyBloom(Framebuffer& fb, const PostSettings& s);
void toneMap(const Framebuffer& fb, Image8& out, const PostSettings& s);

// Convenience: bloom (in place) + tone map.
void finishFrame(Framebuffer& fb, Image8& out, const PostSettings& s);

} // namespace solsim::render
