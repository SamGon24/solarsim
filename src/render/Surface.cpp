// =============================================================================
//  SolarSim - Surface.cpp
// =============================================================================
#include "Surface.hpp"

#include "../core/MathUtil.hpp"
#include "../core/StringUtil.hpp"

#include <cmath>

namespace solsim::render {

using astro::SurfaceStyle;

SurfaceShader::SurfaceShader(SurfaceStyle style, std::uint32_t c1, std::uint32_t c2, std::uint32_t seed,
                             const std::string& name)
    : m_style(style), m_c1(Color::fromHex(c1)), m_c2(Color::fromHex(c2)), m_noise(seed) {
    switch (style) {
        case SurfaceStyle::Terran: m_atmo = {true, Color{0.30f, 0.55f, 1.0f}, 0.9f}; break;
        case SurfaceStyle::Cloudy: m_atmo = {true, m_c1 * 1.1f, 0.6f}; break;
        case SurfaceStyle::GasBanded: m_atmo = {true, m_c1, 0.35f}; break;
        case SurfaceStyle::IceGiant: m_atmo = {true, m_c1 * 1.2f, 0.55f}; break;
        case SurfaceStyle::Desert: m_atmo = {true, Color{0.9f, 0.55f, 0.35f}, 0.25f}; break;
        default: m_atmo = {}; break;
    }
    m_greatSpot = str::iequals(name, "Jupiter");
}

SurfaceSample SurfaceShader::sample(const Vec3& n) const {
    SurfaceSample s;
    const double lat = std::asin(math::clamp(n.z, -1.0, 1.0));
    const double absLat = std::fabs(lat);
    switch (m_style) {
        case SurfaceStyle::Star: {
            // Granulation + faint sunspot belts (+-10..30 deg latitude)
            const double g = m_noise.fbm(n * 24.0, 3);
            double spots = 0.0;
            if (absLat > 0.15 && absLat < 0.55) spots = math::smoothstep(0.80, 0.86, m_noise.fbm(n * 6.0, 2));
            s.albedo = m_c1 * static_cast<float>(0.92 + 0.16 * g) * static_cast<float>(1.0 - 0.6 * spots);
            s.emissive = 1.f;
            break;
        }
        case SurfaceStyle::Terran: {
            const double h = m_noise.fbm(n * 2.2, 5);
            const double clouds = m_noise.fbm(n * 3.5 + Vec3{5, 5, 5}, 4);
            const bool land = h > 0.54;
            Color base;
            if (land) {
                const double dry = math::smoothstep(0.0, 0.45, 1.0 - std::fabs(absLat - 0.4) * 2.5) *
                                   m_noise.noise(n * 6.0);
                base = Color::lerp(m_c2, Color{0.55f, 0.45f, 0.28f}, static_cast<float>(dry));
            } else {
                base = m_c1 * static_cast<float>(0.55 + 0.6 * h);
                s.specular = 0.6f;
            }
            if (absLat > 1.18 + 0.08 * (h - 0.5)) { base = Color::gray(0.85f); s.specular = 0.f; }
            const double c = math::smoothstep(0.55, 0.75, clouds);
            s.albedo = Color::lerp(base, Color::gray(0.95f), static_cast<float>(c * 0.85));
            s.specular *= static_cast<float>(1.0 - c);
            break;
        }
        case SurfaceStyle::Desert: {
            const double h = m_noise.fbm(n * 3.0, 5);
            Color base = Color::lerp(m_c2, m_c1, static_cast<float>(math::smoothstep(0.35, 0.65, h)));
            if (absLat > 1.35) base = Color::gray(0.9f);
            s.albedo = base;
            break;
        }
        case SurfaceStyle::GasBanded: {
            // Noise stretched along latitude circles -> zonal turbulence.
            const double turb = m_noise.fbm(Vec3{n.x * 3.0, n.y * 3.0, n.z * 10.0}, 4);
            const double band = std::sin(lat * 15.0 + (turb - 0.5) * 2.4);
            const double fine = std::sin(lat * 43.0 + (turb - 0.5) * 3.0);
            double t = 0.5 + 0.38 * band + 0.12 * fine;
            Color base = Color::lerp(m_c2, m_c1, static_cast<float>(math::clamp(t, 0.0, 1.0)));
            if (m_greatSpot) {
                // Great Red Spot near -22 deg latitude, fixed System II longitude.
                const double lon = std::atan2(n.y, n.x);
                const double dl = math::wrapPi(lon - 1.0) / 0.20;
                const double db = (lat + 0.384) / 0.09;
                const double r2 = dl * dl + db * db;
                if (r2 < 1.0) base = Color::lerp(base, Color{0.70f, 0.28f, 0.14f}, static_cast<float>(1.0 - r2));
            }
            s.albedo = base;
            break;
        }
        case SurfaceStyle::IceGiant: {
            const double t = 0.5 + 0.1 * std::sin(lat * 9.0) + 0.1 * (m_noise.fbm(n * 3.0, 3) - 0.5);
            s.albedo = Color::lerp(m_c2, m_c1, static_cast<float>(t));
            break;
        }
        case SurfaceStyle::Cloudy: {
            const double swirl = m_noise.fbm(n * 2.5 + Vec3{0, 0, lat * 3.0}, 4);
            s.albedo = Color::lerp(m_c2, m_c1, static_cast<float>(0.55 + 0.45 * swirl));
            break;
        }
        case SurfaceStyle::Icy: {
            const double r = m_noise.ridged(n * 3.0, 4);
            s.albedo = Color::lerp(m_c1, m_c2, static_cast<float>(math::smoothstep(0.55, 0.9, r)));
            break;
        }
        case SurfaceStyle::Volcanic: {
            const double h = m_noise.fbm(n * 4.0, 4);
            const double spot = math::smoothstep(0.72, 0.78, m_noise.noise(n * 9.0));
            s.albedo = Color::lerp(m_c1, m_c2, static_cast<float>(h * 0.6 + spot * 0.4));
            s.emissive = static_cast<float>(spot * 0.15);
            break;
        }
        case SurfaceStyle::Comet: {
            s.albedo = m_c2 * static_cast<float>(0.4 + 0.3 * m_noise.fbm(n * 4.0, 3));
            break;
        }
        case SurfaceStyle::Cratered:
        case SurfaceStyle::Rocky:
        case SurfaceStyle::Metallic:
        default: {
            const double h = m_noise.fbm(n * 3.0, 5);
            const double craters = m_noise.ridged(n * 7.0, 3);
            double t = math::smoothstep(0.3, 0.7, h);
            if (m_style == SurfaceStyle::Cratered) t = t * 0.7 + 0.3 * math::smoothstep(0.6, 0.95, craters);
            s.albedo = Color::lerp(m_c2, m_c1, static_cast<float>(t));
            if (m_style == SurfaceStyle::Metallic) s.specular = 0.4f;
            break;
        }
    }
    return s;
}

} // namespace solsim::render
