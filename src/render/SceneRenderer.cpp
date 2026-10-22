// =============================================================================
//  SolarSim - SceneRenderer.cpp
// =============================================================================
#include "SceneRenderer.hpp"

#include "Font.hpp"
#include "Raster.hpp"

#include "../astro/OrbitalElements.hpp"
#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"
#include "../core/StringUtil.hpp"

#include <algorithm>
#include <cmath>

namespace solsim::render {

using namespace solsim::constants;
using astro::BodyKind;

namespace {

constexpr float kSunIntensity = 1.55f;
constexpr float kAmbient = 0.018f;

std::uint32_t lighten(std::uint32_t rgb, float f) {
    auto ch = [&](int shift) {
        const float v = float((rgb >> shift) & 0xFF);
        return static_cast<std::uint32_t>(std::min(255.f, v + (255.f - v) * f));
    };
    return (ch(16) << 16) | (ch(8) << 8) | ch(0);
}

int kindPriority(BodyKind k) {
    switch (k) {
        case BodyKind::Star: return 0;
        case BodyKind::Planet: return 1;
        case BodyKind::DwarfPlanet: return 2;
        case BodyKind::Spacecraft: return 2;
        case BodyKind::Comet: return 3;
        case BodyKind::Moon: return 3;
        default: return 4;
    }
}

} // namespace

SceneRenderer::SceneRenderer() : m_stars(9000, 7) {}

void SceneRenderer::prepare(const physics::NBodySystem& sys) {
    m_map.prepare(sys);
    if (m_revision == sys.revision() && m_shaders.size() == sys.namedCount()) return;
    m_revision = sys.revision();
    m_shaders.clear();
    for (std::size_t i = 0; i < sys.namedCount(); ++i) {
        const auto& inf = sys.info(i);
        m_shaders.emplace_back(inf.style, inf.color, inf.color2, static_cast<std::uint32_t>(1234 + 977 * i), inf.name);
    }
}

Vec3 SceneRenderer::displayPosition(const physics::NBodySystem& sys, std::size_t index) {
    m_map.prepare(sys);
    return m_map.position(sys, index);
}

double SceneRenderer::displayRadius(const physics::NBodySystem& sys, std::size_t index) {
    m_map.prepare(sys);
    return m_map.displayRadius(sys, index);
}

void SceneRenderer::record(const physics::Simulation& sim) { m_trails.record(sim.system(), sim.time()); }

// -----------------------------------------------------------------------------
void SceneRenderer::render(const physics::Simulation& sim, const Camera& cam, Framebuffer& fb) {
    const auto& sys = sim.system();
    prepare(sys);
    m_labels.clear();
    fb.clear(Color{0.0015f, 0.0018f, 0.003f});
    m_sunDisplay = sys.size() ? m_map.position(sys, 0) : Vec3{};

    if (m_opt.stars) drawStars(cam, fb);
    computeBodyDraws(sys, cam);

    // Opaque bodies first so lines and particles can be depth tested against them.
    for (const auto& bd : m_draws)
        if (bd.visible) drawSphere(sim, bd, cam, fb);
    if (m_opt.grid) drawGrid(cam, fb);
    if (m_opt.orbits) drawOrbits(sys, cam, fb);
    if (m_opt.trails) drawTrails(sys, cam, fb);
    if (m_opt.particles) drawParticles(sys, cam, fb);
    if (m_opt.rings)
        for (const auto& bd : m_draws)
            if (bd.visible && sys.info(bd.index).rings.present) drawRings(sys, bd, cam, fb);
    if (m_opt.cometTails)
        for (const auto& bd : m_draws)
            if (bd.visible && sys.info(bd.index).kind == BodyKind::Comet) drawCometTail(sys, bd, cam, fb);

    // Stellar glow (corona + optical halo).
    if (m_opt.glow) {
        for (const auto& bd : m_draws) {
            if (!bd.visible || sys.info(bd.index).kind != BodyKind::Star) continue;
            const Color c = Color::fromHex(sys.info(bd.index).color);
            const double r0 = std::max(1.5, bd.rpx * 1.05);
            drawGlow(fb, bd.proj.x, bd.proj.y, r0, std::min(r0 * 12.0, 0.45 * fb.height()), c * 0.6f, 1.45);
        }
    }

    // Focus marker.
    if (m_opt.focus >= 0) {
        for (const auto& bd : m_draws) {
            if (static_cast<int>(bd.index) != m_opt.focus || !bd.visible) continue;
            drawCircle(fb, bd.proj.x, bd.proj.y, bd.rpx + 4.0 + 0.004 * fb.height(), Color{0.4f, 0.8f, 1.0f}, 0.6f);
        }
    }
    collectLabels(sys);
}

// -----------------------------------------------------------------------------
void SceneRenderer::drawStars(const Camera& cam, Framebuffer& fb) {
    const float k = 0.55f * m_opt.starBrightness;
    const bool lowRes = fb.height() < 200;
    for (const Star& s : m_stars.stars()) {
        if (lowRes && s.magnitude > 5.5f) continue;
        Projected p;
        if (!cam.projectDirection(s.direction, p)) continue;
        if (p.x < -3 || p.y < -3 || p.x > fb.width() + 3 || p.y > fb.height() + 3) continue;
        // Perceptually compressed flux: 10^(-0.4 * 0.6 * (m - 1))
        const float flux = k * std::pow(10.f, -0.24f * (s.magnitude - 1.0f));
        const double radius = 0.55 + 0.22 * std::max(0.0, 2.5 - double(s.magnitude));
        drawSplat(fb, p.x, p.y, 1e30, s.color * flux, radius, false);
    }
}

void SceneRenderer::drawGrid(const Camera& cam, Framebuffer& fb) {
    static const double radii[] = {0.5, 1, 2, 5, 10, 20, 30, 50, 100};
    const Color c{0.25f, 0.35f, 0.55f};
    for (double r : radii) {
        const double R = m_map.radial(r);
        std::vector<Vec3> pts;
        for (int k = 0; k < 160; ++k) {
            const double a = TWO_PI * k / 160.0;
            pts.emplace_back(R * std::cos(a), R * std::sin(a), 0.0);
        }
        polyline3D(pts, cam, fb, c, 0.22f, false, true);
        Projected p;
        if (cam.project(Vec3{0.0, -R, 0.0}, p)) {
            ScreenLabel l;
            l.x = p.x + 2;
            l.y = p.y + 2;
            l.text = str::format("%g AU", r);
            l.color = 0x5A7AA8;
            l.priority = 6;
            m_labels.push_back(l);
        }
    }
    // Direction of the vernal equinox (+x axis of the ecliptic frame).
    std::vector<Vec3> axis{{0, 0, 0}, {m_map.radial(60.0), 0, 0}};
    polyline3D(axis, cam, fb, Color{0.5f, 0.4f, 0.2f}, 0.25f, false, false);
}

void SceneRenderer::polyline3D(const std::vector<Vec3>& pts, const Camera& cam, Framebuffer& fb, const Color& c,
                               float alpha, bool fade, bool closed) {
    if (pts.size() < 2) return;
    LineStyle st;
    st.color = c;
    st.blend = BlendMode::Additive;
    st.depthTest = true;
    const std::size_t n = pts.size();
    const std::size_t segs = closed ? n : n - 1;
    Projected a, b;
    bool aOk = cam.project(pts[0], a);
    for (std::size_t k = 0; k < segs; ++k) {
        const Vec3& pb = pts[(k + 1) % n];
        const bool bOk = cam.project(pb, b);
        if (aOk && bOk) {
            st.alpha = fade ? alpha * float(k + 1) / float(segs) : alpha;
            drawLine(fb, a.x, a.y, a.depth, b.x, b.y, b.depth, st);
        }
        a = b;
        aOk = bOk;
    }
}

void SceneRenderer::drawOrbits(const physics::NBodySystem& sys, const Camera& cam, Framebuffer& fb) {
    for (const auto& bd : m_draws) {
        const std::size_t i = bd.index;
        const auto& inf = sys.info(i);
        if (inf.parent < 0 || inf.kind == BodyKind::Star) continue;
        // Satellites: skip orbits while the system is too small on screen.
        if (inf.parent > 0) {
            const BodyDraw& pd = m_draws[static_cast<std::size_t>(inf.parent)];
            if (!pd.visible || !bd.visible) continue;
        }
        const double mu = sys.twoBodyMu(i);
        const Vec3 r = sys.relativePosition(i, inf.parent);
        const Vec3 v = sys.relativeVelocity(i, inf.parent);
        const astro::OrbitalElements el = astro::stateToElements(r, v, mu, 0.0);
        Vec3 P, Q, W;
        astro::perifocalBasis(el.i, el.node, el.argp, P, Q, W);
        std::vector<Vec3> pts;
        const int N = 240;
        bool closed = false;
        if (el.isElliptic()) {
            const double a = el.a(), b = a * std::sqrt(1.0 - el.e * el.e);
            for (int k = 0; k < N; ++k) {
                const double E = TWO_PI * k / N;
                const Vec3 rel = P * (a * (std::cos(E) - el.e)) + Q * (b * std::sin(E));
                pts.push_back(m_map.relativeToParent(sys, inf.parent, rel));
            }
            closed = true;
        } else {
            // Hyperbolic / parabolic arc out to a few times the current distance.
            const double rmax = std::max(3.0 * r.norm(), 4.0 * el.q);
            const double p = el.p();
            const double nuInf = el.isHyperbolic() ? std::acos(-1.0 / el.e) : PI;
            double nuMax = nuInf * 0.999;
            // limit by r <= rmax: r = p / (1 + e cos nu)
            const double cosLim = (p / rmax - 1.0) / el.e;
            if (cosLim > -1.0 && cosLim < 1.0) nuMax = std::min(nuMax, std::acos(cosLim));
            for (int k = 0; k <= N; ++k) {
                const double nu = -nuMax + 2.0 * nuMax * k / N;
                const double rr = p / (1.0 + el.e * std::cos(nu));
                const Vec3 rel = P * (rr * std::cos(nu)) + Q * (rr * std::sin(nu));
                pts.push_back(m_map.relativeToParent(sys, inf.parent, rel));
            }
        }
        Color c = Color::fromHex(inf.color);
        c = c * (0.7f / std::max(0.05f, c.maxComponent()));
        const float alpha = (inf.kind == BodyKind::Moon ? 0.22f : 0.30f) * m_opt.orbitBrightness;
        polyline3D(pts, cam, fb, c, alpha, false, closed);
    }
}

void SceneRenderer::drawTrails(const physics::NBodySystem& sys, const Camera& cam, Framebuffer& fb) {
    std::vector<Vec3> pts;
    for (const auto& bd : m_draws) {
        const std::size_t i = bd.index;
        const std::size_t n = m_trails.count(i);
        if (n < 2) continue;
        const int parent = m_trails.parent(i);
        if (parent > 0) {
            const BodyDraw& pd = m_draws[static_cast<std::size_t>(parent)];
            if (!pd.visible || !bd.visible) continue;
        }
        pts.clear();
        for (std::size_t k = 0; k < n; ++k) pts.push_back(m_map.relativeToParent(sys, parent, m_trails.point(i, k)));
        // Connect the trail to the current position.
        pts.push_back(m_map.position(sys, i));
        Color c = Color::fromHex(sys.info(i).color);
        c = c * (1.0f / std::max(0.05f, c.maxComponent()));
        polyline3D(pts, cam, fb, c, 0.85f * m_opt.orbitBrightness, true, false);
    }
}

void SceneRenderer::drawParticles(const physics::NBodySystem& sys, const Camera& cam, Framebuffer& fb) {
    for (const auto& g : sys.groups()) {
        const std::size_t count = g.end - g.begin;
        if (count == 0) continue;
        // Keep the integrated brightness of a population roughly constant.
        const float per = m_opt.particleBrightness * std::clamp(22.0f / std::sqrt(float(count)), 0.08f, 0.6f);
        const Color c = Color::fromHex(g.color) * per;
        for (std::size_t i = g.begin; i < g.end; ++i) {
            Projected p;
            if (!cam.project(m_map.position(sys, i), p)) continue;
            if (p.x < -2 || p.y < -2 || p.x > fb.width() + 2 || p.y > fb.height() + 2) continue;
            drawSplat(fb, p.x, p.y, p.depth, c, 0.55, true);
        }
    }
}

// -----------------------------------------------------------------------------
void SceneRenderer::computeBodyDraws(const physics::NBodySystem& sys, const Camera& cam) {
    m_draws.assign(sys.namedCount(), BodyDraw{});
    for (std::size_t i = 0; i < sys.namedCount(); ++i) {
        BodyDraw& bd = m_draws[i];
        bd.index = i;
        bd.center = m_map.position(sys, i);
        bd.baseRadius = m_map.displayRadius(sys, i);
        if (!cam.project(bd.center, bd.proj)) continue;
        const auto kind = sys.info(i).kind;
        double minPx = m_opt.minBodyPixels;
        if (kind == BodyKind::Moon || kind == BodyKind::Asteroid || kind == BodyKind::Comet ||
            kind == BodyKind::Spacecraft)
            minPx *= 0.6;
        if (kind == BodyKind::Star) minPx *= 1.6;
        double rpx = cam.projectedRadius(bd.baseRadius, bd.proj.depth);
        bd.radius = bd.baseRadius;
        if (rpx < minPx) {
            bd.radius = minPx * bd.proj.depth / cam.focalPixels();
            rpx = minPx;
        }
        bd.rpx = rpx;
        bd.visible = bd.proj.x > -rpx - 50 && bd.proj.y > -rpx - 50 && bd.proj.x < cam.width() + rpx + 50 &&
                     bd.proj.y < cam.height() + rpx + 50;
    }
    // Hide satellites that would sit on top of their (inflated) parent disc.
    for (std::size_t i = 0; i < m_draws.size(); ++i) {
        const int p = sys.info(i).parent;
        if (p <= 0 || !m_draws[i].visible) continue;
        const BodyDraw& pd = m_draws[static_cast<std::size_t>(p)];
        const double dx = m_draws[i].proj.x - pd.proj.x, dy = m_draws[i].proj.y - pd.proj.y;
        if (std::sqrt(dx * dx + dy * dy) < pd.rpx + m_draws[i].rpx + 2.0) m_draws[i].visible = false;
    }
}

void SceneRenderer::drawSphere(const physics::Simulation& sim, const BodyDraw& bd, const Camera& cam,
                               Framebuffer& fb) {
    const auto& sys = sim.system();
    const auto& inf = sys.info(bd.index);
    const SurfaceShader& shader = m_shaders[bd.index];
    const AtmosphereParams atmo = shader.atmosphere();
    const bool isStar = inf.kind == BodyKind::Star;
    const Vec3 C = bd.center;
    const double R = bd.radius;
    const Vec3 eye = cam.eye();
    const Vec3 L = (m_sunDisplay - C).normalized();
    const double pixelWorld = bd.proj.depth / cam.focalPixels(); // world size of a pixel at the body

    // Body-fixed frame for texturing (spin axis = pole).
    const Vec3 pole = inf.pole.normalized();
    Vec3 e1 = Vec3::unitZ().cross(pole);
    if (e1.norm() < 1e-9) e1 = Vec3::unitX();
    e1 = e1.normalized();
    const Vec3 e2 = pole.cross(e1);
    double theta = inf.rotationPhase0;
    if (std::fabs(inf.rotationPeriod) > 1e-9) theta += TWO_PI * (sim.time() - J2000_JD) / inf.rotationPeriod;
    const double cT = std::cos(theta), sT = std::sin(theta);
    const bool hasRings = inf.rings.present && m_opt.rings;

    const double reach = bd.rpx * 1.08 + 1.5;
    const int x0 = std::max(0, static_cast<int>(std::floor(bd.proj.x - reach)));
    const int x1 = std::min(fb.width() - 1, static_cast<int>(std::ceil(bd.proj.x + reach)));
    const int y0 = std::max(0, static_cast<int>(std::floor(bd.proj.y - reach)));
    const int y1 = std::min(fb.height() - 1, static_cast<int>(std::ceil(bd.proj.y + reach)));
    const Vec3 oc = eye - C;

    for (int py = y0; py <= y1; ++py) {
        for (int px = x0; px <= x1; ++px) {
            const Vec3 dir = cam.rayDirection(px + 0.5, py + 0.5);
            const double b = oc.dot(dir);
            if (b > 0.0) continue; // sphere behind the ray origin
            const Vec3 closest = oc - dir * b; // perpendicular offset from the centre
            const double dperp = closest.norm();
            const double coverage = math::clamp((R - dperp) / pixelWorld + 0.5, 0.0, 1.0);
            if (coverage <= 0.0) continue;
            const double h = std::sqrt(std::max(0.0, R * R - dperp * dperp));
            const double t = -b - h;
            const Vec3 P = eye + dir * t;
            Vec3 N = (P - C) / R;
            if (dperp >= R) N = closest.normalized();
            N = N.normalized();
            const float depth = static_cast<float>(t * dir.dot(cam.forward()));
            if (depth >= fb.depth(px, py)) continue;

            // Surface coordinates in the rotating body frame.
            const double bx = N.dot(e1), by = N.dot(e2), bz = N.dot(pole);
            const Vec3 nb{bx * cT + by * sT, -bx * sT + by * cT, bz};
            const SurfaceSample s = shader.sample(nb);
            const Vec3 V = -dir;
            const double mu = std::max(0.0, N.dot(V));

            Color col;
            if (isStar) {
                // Eddington limb darkening I(mu) = 1 - u (1 - mu), u ~ 0.6 for the Sun in V.
                const float limb = static_cast<float>(1.0 - 0.6 * (1.0 - mu));
                col = s.albedo * (3.2f * limb);
            } else {
                const double ndl = N.dot(L);
                double diffuse = std::max(0.0, ndl);
                if (atmo.present) diffuse = math::clamp((ndl + 0.06) / 1.06, 0.0, 1.0);
                // Ring shadow on the planet.
                if (hasRings && ndl > 0.0) {
                    const double denom = L.dot(pole);
                    if (std::fabs(denom) > 1e-9) {
                        const double ts = (C - P).dot(pole) / denom;
                        if (ts > 0.0) {
                            const double rr = (P + L * ts - C).norm() / R * (inf.radius * AU_KM);
                            diffuse *= 1.0 - 0.85 * ringOpacity(inf.rings, inf.name, rr);
                        }
                    }
                }
                col = s.albedo * static_cast<float>(diffuse * kSunIntensity + kAmbient);
                if (s.specular > 0.f && ndl > 0.0) {
                    const Vec3 H = (L + V).normalized();
                    const double spec = std::pow(std::max(0.0, N.dot(H)), 60.0);
                    col += Color{1.0f, 0.95f, 0.85f} * static_cast<float>(s.specular * spec * 1.4);
                }
                if (atmo.present) {
                    const double rim = std::pow(1.0 - mu, 3.0) * math::clamp(ndl + 0.35, 0.0, 1.0);
                    col += atmo.color * static_cast<float>(atmo.strength * rim * 0.9);
                }
                if (s.emissive > 0.f) col += s.albedo * s.emissive;
            }
            if (coverage >= 0.999) {
                fb.at(px, py) = col;
                fb.depth(px, py) = depth;
            } else {
                fb.blend(px, py, col, static_cast<float>(coverage));
                if (coverage > 0.5) fb.depth(px, py) = depth;
            }
        }
    }
}

double SceneRenderer::ringOpacity(const astro::RingSystem& rings, const std::string& planet, double r) {
    if (!rings.present || r < rings.innerKm || r > rings.outerKm) return 0.0;
    double o;
    if (str::iequals(planet, "Saturn")) {
        // Radial structure of Saturn's main rings (km from the planet centre).
        if (r < 92000.0) o = 0.12 + 0.08 * math::smoothstep(74658.0, 92000.0, r);        // C ring
        else if (r < 117580.0) o = 0.55 + 0.4 * math::smoothstep(92000.0, 105000.0, r); // B ring
        else if (r < 122170.0) o = 0.06;                                                  // Cassini division
        else if (r > 133410.0 && r < 133740.0) o = 0.04;                                  // Encke gap
        else o = 0.55 - 0.1 * math::smoothstep(122170.0, 136775.0, r);                    // A ring
        // gentle ringlet modulation (low frequency to avoid aliasing)
        o *= 0.9 + 0.1 * std::sin(r * 0.00045);
    } else {
        // Narrow ringlets for the other giants.
        const double t = (r - rings.innerKm) / (rings.outerKm - rings.innerKm);
        o = 0.35 * std::pow(std::fabs(std::sin(t * 31.0)), 12.0) + 0.05;
    }
    return math::clamp(o * rings.opacity, 0.0, 1.0);
}

void SceneRenderer::drawRings(const physics::NBodySystem& sys, const BodyDraw& bd, const Camera& cam,
                              Framebuffer& fb) {
    const auto& inf = sys.info(bd.index);
    const Vec3 pole = inf.pole.normalized();
    const Vec3 C = bd.center;
    const double R = bd.radius;
    const double planetKm = inf.radius * AU_KM;
    const double outer = R * inf.rings.outerKm / planetKm;
    const double outerPx = bd.rpx * inf.rings.outerKm / planetKm;
    if (outerPx < 1.5) return;
    const Vec3 eye = cam.eye();
    const Vec3 L = (m_sunDisplay - C).normalized();
    const Color base = str::iequals(inf.name, "Saturn") ? Color::fromHex(0xD9C9A3) : Color::fromHex(0x8899AA);
    const double reach = outerPx * 1.05 + 2.0;
    const int x0 = std::max(0, static_cast<int>(bd.proj.x - reach)), x1 = std::min(fb.width() - 1, static_cast<int>(bd.proj.x + reach));
    const int y0 = std::max(0, static_cast<int>(bd.proj.y - reach)), y1 = std::min(fb.height() - 1, static_cast<int>(bd.proj.y + reach));
    const double lit = 0.30 + 0.70 * std::pow(std::fabs(L.dot(pole)), 0.5);

    for (int py = y0; py <= y1; ++py)
        for (int px = x0; px <= x1; ++px) {
            const Vec3 dir = cam.rayDirection(px + 0.5, py + 0.5);
            const double denom = dir.dot(pole);
            if (std::fabs(denom) < 1e-12) continue;
            const double t = (C - eye).dot(pole) / denom;
            if (t <= 0.0) continue;
            const Vec3 Q = eye + dir * t;
            const double rr = (Q - C).norm();
            if (rr > outer * 1.001) continue;
            const double op = ringOpacity(inf.rings, inf.name, rr / R * planetKm);
            if (op <= 0.0) continue;
            const float depth = static_cast<float>(t * dir.dot(cam.forward()));
            if (depth >= fb.depth(px, py)) continue;
            // Planet shadow on the ring: does the ray Q -> Sun hit the planet sphere?
            double shade = 1.0;
            const Vec3 oq = Q - C;
            const double bq = oq.dot(L);
            if (bq < 0.0 && (oq - L * bq).norm() < R) shade = 0.12;
            const Color col = base * static_cast<float>(lit * shade * 0.95);
            fb.blend(px, py, col, static_cast<float>(op));
        }
}

void SceneRenderer::drawCometTail(const physics::NBodySystem& sys, const BodyDraw& bd, const Camera& cam,
                                  Framebuffer& fb) {
    const std::size_t i = bd.index;
    const double r = (sys.pos(i) - sys.pos(0)).norm();
    if (r > 5.0) return;
    // Activity rises steeply inside ~3 AU (water ice sublimation).
    const double activity = math::clamp(2.0 / (r * r), 0.0, 6.0);
    const Vec3 C = bd.center;
    const Vec3 anti = (C - m_sunDisplay).normalized();
    const Vec3 vel = (sys.vel(i) - sys.vel(0)).normalized();
    const double len = 0.12 * std::sqrt(activity) * std::max(0.4, m_map.radial(1.0));
    const Vec3 dustDir = (anti * 0.85 - vel * 0.35).normalized();

    LineStyle ion;
    ion.color = Color{0.35f, 0.55f, 1.0f};
    ion.blend = BlendMode::Additive;
    LineStyle dust = ion;
    dust.color = Color{1.0f, 0.88f, 0.65f};
    const int N = 40;
    Projected a, b, c2, d;
    for (int k = 0; k < N; ++k) {
        const double s0 = double(k) / N, s1 = double(k + 1) / N;
        const float fade = static_cast<float>(std::pow(1.0 - s0, 1.6) * std::min(1.0, activity));
        // Ion tail: straight, anti-solar.
        if (cam.project(C + anti * (len * 1.3 * s0), a) && cam.project(C + anti * (len * 1.3 * s1), b)) {
            ion.alpha = 0.9f * fade;
            drawLine(fb, a.x, a.y, a.depth, b.x, b.y, b.depth, ion);
        }
        // Dust tail: curved (quadratic Bezier from anti-solar to the lagging direction).
        auto dustPoint = [&](double s) {
            const Vec3 p0 = C, p1 = C + anti * (0.5 * len), p2 = C + dustDir * len;
            return p0 * ((1 - s) * (1 - s)) + p1 * (2 * s * (1 - s)) + p2 * (s * s);
        };
        if (cam.project(dustPoint(s0), c2) && cam.project(dustPoint(s1), d)) {
            dust.alpha = 0.7f * fade;
            drawLine(fb, c2.x, c2.y, c2.depth, d.x, d.y, d.depth, dust);
        }
    }
    // Coma
    drawGlow(fb, bd.proj.x, bd.proj.y, std::max(1.2, bd.rpx * 1.5), std::max(6.0, bd.rpx * 10.0 * std::sqrt(activity)),
             Color{0.55f, 0.75f, 1.0f} * static_cast<float>(0.4 * std::min(1.0, activity)), 1.2);
}

// -----------------------------------------------------------------------------
void SceneRenderer::collectLabels(const physics::NBodySystem& sys) {
    struct Candidate {
        ScreenLabel label;
        double cx = 0, cy = 0, r = -1; // body centre and radius (r < 0: fixed position)
        double size = 0;               // physical radius, larger bodies win ties
    };
    std::vector<Candidate> candidates;
    for (auto& l : m_labels) candidates.push_back({l, 0, 0, -1, 0}); // grid labels added earlier this frame
    m_labels.clear();
    if (!m_opt.labels) return;
    for (const auto& bd : m_draws) {
        if (!bd.visible) continue;
        const auto& inf = sys.info(bd.index);
        Candidate c;
        c.label.text = inf.name;
        c.label.color = lighten(inf.color, 0.45f);
        c.label.priority = kindPriority(inf.kind);
        if (static_cast<int>(bd.index) == m_opt.focus) c.label.priority = -1;
        c.cx = bd.proj.x;
        c.cy = bd.proj.y;
        c.r = bd.rpx;
        c.size = inf.radius;
        candidates.push_back(c);
    }
    std::stable_sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
        if (a.label.priority != b.label.priority) return a.label.priority < b.label.priority;
        return a.size > b.size;
    });
    struct Box { double x0, y0, x1, y1; };
    std::vector<Box> taken;
    const double cw = m_opt.labelCharWidth, ch = m_opt.labelCharHeight;
    auto overlaps = [&](const Box& bx) {
        for (const Box& t : taken)
            if (!(bx.x1 < t.x0 || bx.x0 > t.x1 || bx.y1 < t.y0 || bx.y0 > t.y1)) return true;
        return false;
    };
    for (auto& c : candidates) {
        const double w = cw * double(str::utf8Length(c.label.text));
        std::vector<std::pair<double, double>> spots;
        if (c.r < 0) {
            spots.emplace_back(c.label.x, c.label.y);
        } else {
            const double gap = 0.6 * cw;
            spots.emplace_back(c.cx + c.r + gap, c.cy - 0.5 * ch);      // right
            spots.emplace_back(c.cx - c.r - gap - w, c.cy - 0.5 * ch);  // left
            spots.emplace_back(c.cx - 0.5 * w, c.cy - c.r - ch - 1.0);  // above
            spots.emplace_back(c.cx - 0.5 * w, c.cy + c.r + 1.0);       // below
        }
        for (const auto& sp : spots) {
            const Box bx{sp.first - 1, sp.second - 1, sp.first + w + 1, sp.second + ch + 1};
            if (overlaps(bx)) continue;
            taken.push_back(bx);
            c.label.x = sp.first;
            c.label.y = sp.second;
            m_labels.push_back(c.label);
            break;
        }
    }
}

void SceneRenderer::drawOverlay(Image8& img, const std::vector<std::string>& hud) const {
    const int s = std::max(1, m_opt.textScale);
    for (const auto& l : m_labels)
        BitmapFont::draw(img, static_cast<int>(std::lround(l.x)), static_cast<int>(std::lround(l.y)), l.text, l.color, s);
    int y = 6 * s;
    for (const auto& line : hud) {
        BitmapFont::draw(img, 6 * s, y, line, 0xDDE6F0, s);
        y += 10 * s;
    }
}

} // namespace solsim::render
