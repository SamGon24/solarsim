// =============================================================================
//  SolarSim - Surface.hpp
//  Procedural surface appearance of bodies. Given a unit normal in the body's
//  rotating frame (z = spin axis) it returns albedo colour and atmospheric
//  properties used by the shader in SceneRenderer.
// =============================================================================
#pragma once

#include "Color.hpp"
#include "Noise.hpp"

#include "../astro/BodyCatalog.hpp"
#include "../core/Vec3.hpp"

#include <cstdint>

namespace solsim::render {

struct SurfaceSample {
    Color albedo;        // diffuse reflectance (linear)
    float specular = 0;  // ocean glint strength
    float emissive = 0;  // self illumination (stars, lava)
};

struct AtmosphereParams {
    bool present = false;
    Color color{0, 0, 0};
    float strength = 0.f;
};

class SurfaceShader {
public:
    SurfaceShader(astro::SurfaceStyle style, std::uint32_t c1, std::uint32_t c2, std::uint32_t seed,
                  const std::string& name);

    // nBody: unit normal in the body frame (z = north pole)
    SurfaceSample sample(const Vec3& nBody) const;
    AtmosphereParams atmosphere() const { return m_atmo; }
    astro::SurfaceStyle style() const { return m_style; }

private:
    astro::SurfaceStyle m_style;
    Color m_c1, m_c2;
    ValueNoise m_noise;
    AtmosphereParams m_atmo;
    bool m_greatSpot = false; // Jupiter's Great Red Spot
};

} // namespace solsim::render
