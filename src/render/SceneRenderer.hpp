// =============================================================================
//  SolarSim - SceneRenderer.hpp
//  Software renderer for the simulated system:
//    1. star field (real bright stars + Milky Way population)
//    2. ecliptic distance grid
//    3. osculating orbits (conics from the instantaneous state vectors)
//    4. trails (recent history, relative to the parent body)
//    5. test particles (asteroid / Kuiper belts)
//    6. bodies: per-pixel ray casting of spheres with procedural surfaces,
//       Lambert lighting from the Sun, limb darkening, atmospheres, ocean glint
//    7. planetary rings with ring shadows on the planet and planet shadows on
//       the rings
//    8. comet comae and ion/dust tails (pointing away from the Sun)
//    9. labels (collected for overlay drawing after tone mapping)
// =============================================================================
#pragma once

#include "Camera.hpp"
#include "DisplayMapping.hpp"
#include "Image.hpp"
#include "StarField.hpp"
#include "Surface.hpp"
#include "Trails.hpp"

#include "../physics/Simulation.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace solsim::render {

struct SceneOptions {
    bool stars = true;
    bool grid = false;
    bool orbits = true;
    bool trails = true;
    bool particles = true;
    bool labels = true;
    bool rings = true;
    bool glow = true;
    bool cometTails = true;
    double minBodyPixels = 2.5;   // planets never shrink below this radius
    float starBrightness = 1.0f;
    float orbitBrightness = 1.0f;
    float particleBrightness = 1.0f;
    int focus = -1;               // index of the focused body (marker), -1 = none
    // Label layout metrics (pixels per character): 6x7*scale for images,
    // 1x2 for terminal half-block rendering.
    double labelCharWidth = 6.0;
    double labelCharHeight = 8.0;
    int textScale = 1;
};

struct ScreenLabel {
    double x = 0, y = 0;   // top-left anchor in framebuffer pixels
    std::string text;
    std::uint32_t color = 0xFFFFFF;
    int priority = 0;      // lower = more important
};

class SceneRenderer {
public:
    SceneRenderer();

    SceneOptions& options() { return m_opt; }
    const SceneOptions& options() const { return m_opt; }
    DisplayMapping& mapping() { return m_map; }
    const DisplayMapping& mapping() const { return m_map; }
    TrailRecorder& trails() { return m_trails; }

    // Record trail samples (call after advancing the simulation).
    void record(const physics::Simulation& sim);
    void render(const physics::Simulation& sim, const Camera& cam, Framebuffer& fb);

    // Labels from the last render() call (already de-overlapped).
    const std::vector<ScreenLabel>& labels() const { return m_labels; }
    // Draws labels + HUD text on the tone-mapped image.
    void drawOverlay(Image8& img, const std::vector<std::string>& hud) const;

    Vec3 displayPosition(const physics::NBodySystem& sys, std::size_t index);
    double displayRadius(const physics::NBodySystem& sys, std::size_t index);

private:
    struct BodyDraw {
        std::size_t index = 0;
        Vec3 center;
        double radius = 0;   // display radius actually drawn (after min-pixel inflation)
        double baseRadius = 0;
        Projected proj;
        double rpx = 0;      // projected radius in pixels
        bool visible = false;
    };

    void prepare(const physics::NBodySystem& sys);
    void drawStars(const Camera& cam, Framebuffer& fb);
    void drawGrid(const Camera& cam, Framebuffer& fb);
    void drawOrbits(const physics::NBodySystem& sys, const Camera& cam, Framebuffer& fb);
    void drawTrails(const physics::NBodySystem& sys, const Camera& cam, Framebuffer& fb);
    void drawParticles(const physics::NBodySystem& sys, const Camera& cam, Framebuffer& fb);
    void computeBodyDraws(const physics::NBodySystem& sys, const Camera& cam);
    void drawSphere(const physics::Simulation& sim, const BodyDraw& bd, const Camera& cam, Framebuffer& fb);
    void drawRings(const physics::NBodySystem& sys, const BodyDraw& bd, const Camera& cam, Framebuffer& fb);
    void drawCometTail(const physics::NBodySystem& sys, const BodyDraw& bd, const Camera& cam, Framebuffer& fb);
    void collectLabels(const physics::NBodySystem& sys);
    void polyline3D(const std::vector<Vec3>& pts, const Camera& cam, Framebuffer& fb, const Color& c, float alpha,
                    bool fade, bool closed);

    static double ringOpacity(const astro::RingSystem& rings, const std::string& planet, double rKm);

    SceneOptions m_opt;
    DisplayMapping m_map;
    TrailRecorder m_trails;
    StarField m_stars;
    std::vector<SurfaceShader> m_shaders;
    std::vector<BodyDraw> m_draws;
    std::vector<ScreenLabel> m_labels;
    std::uint64_t m_revision = ~0ULL;
    Vec3 m_sunDisplay;
};

} // namespace solsim::render
