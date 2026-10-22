// =============================================================================
//  SolarSim - Image.cpp
// =============================================================================
#include "Image.hpp"

#include <algorithm>

namespace solsim::render {

void Framebuffer::resize(int w, int h) {
    m_w = std::max(1, w);
    m_h = std::max(1, h);
    m_color.assign(static_cast<std::size_t>(m_w) * m_h, Color{});
    m_depth.assign(static_cast<std::size_t>(m_w) * m_h, kFar);
}

void Framebuffer::clear(const Color& c) {
    std::fill(m_color.begin(), m_color.end(), c);
    std::fill(m_depth.begin(), m_depth.end(), kFar);
}

} // namespace solsim::render
