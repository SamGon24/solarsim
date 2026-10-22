// =============================================================================
//  SolarSim - CmdView.cpp
//  `solarsim view` - real-time interactive 3D view rendered into the terminal
//  with 24-bit colour half-block characters (or 256 colours / plain ASCII).
// =============================================================================
#include "AppContext.hpp"

#include "../astro/JulianDate.hpp"
#include "../core/Constants.hpp"
#include "../core/Log.hpp"
#include "../core/StringUtil.hpp"
#include "../core/Timer.hpp"
#include "../io/ImageWriter.hpp"
#include "../platform/Terminal.hpp"
#include "../render/TerminalRenderer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <memory>

namespace solsim::app {

using namespace solsim::constants;
using platform::Key;

namespace {

struct ViewState {
    bool paused = false;
    double warp = 10.0;       // days per second
    double direction = 1.0;   // +1 forward, -1 backward in time
    int focus = 0;            // body index the camera follows
    bool showHelp = false;
    std::string message;
    double messageUntil = 0.0;
    render::TerminalColorMode colorMode = render::TerminalColorMode::TrueColor;
};

const char* kIntegratorCycle[] = {"yoshida4", "leapfrog", "rk4", "dopri5", "wh", "symplectic-euler", "euler"};

double frameDistance(render::SceneRenderer& r, const physics::NBodySystem& sys, int focus, double outerFallback) {
    if (focus <= 0) return outerFallback;
    const Vec3 c = r.displayPosition(sys, static_cast<std::size_t>(focus));
    double span = r.displayRadius(sys, static_cast<std::size_t>(focus)) * 14.0;
    for (std::size_t i = 0; i < sys.namedCount(); ++i)
        if (sys.info(i).parent == focus) span = std::max(span, 2.2 * (r.displayPosition(sys, i) - c).norm());
    return span;
}

std::string screenshot(const physics::Simulation& sim, render::SceneRenderer& r, const render::Camera& cam,
                       const render::PostSettings& post) {
    const int W = 1920, H = 1080;
    render::SceneRenderer hi;
    hi.options() = r.options();
    hi.options().textScale = 2;
    hi.options().labelCharWidth = 12;
    hi.options().labelCharHeight = 16;
    hi.options().minBodyPixels = 5.0;
    hi.mapping() = r.mapping();
    hi.trails() = r.trails();
    render::Camera c = cam;
    c.setViewport(W, H);
    render::Framebuffer fb(W, H);
    render::Image8 img;
    hi.render(sim, c, fb);
    render::finishFrame(fb, img, post);
    hi.drawOverlay(img, {"SolarSim  " + astro::formatUTC(sim.time(), false)});
    const std::time_t now = std::time(nullptr);
    char stamp[32] = "screenshot";
    if (const std::tm* tm = std::localtime(&now)) std::strftime(stamp, sizeof stamp, "%Y%m%d_%H%M%S", tm);
    const std::string path = std::string("solarsim_") + stamp + ".png";
    std::string err;
    if (!io::writeImage(path, img, err)) return "screenshot failed: " + err;
    return "saved " + path;
}

} // namespace

int cmdView(AppContext& ctx) {
    const bool tty = platform::stdoutIsTerminal();
    const long long maxFrames = ctx.args().getInt("max-frames", -1);
    if (!tty && maxFrames < 0) {
        std::fprintf(stderr, "error: 'view' needs an interactive terminal.\n"
                             "Use 'solarsim render --out solar.gif' for files, or add --max-frames N to force output.\n");
        return 2;
    }

    ViewState st;
    st.warp = ctx.args().getDouble("warp", 10.0);
    if (!(st.warp > 0.0)) throw UsageError("--warp must be positive");
    const std::string colorArg = str::toLower(ctx.args().getString("color", "auto"));
    if (colorArg == "auto") {
        const bool vt = platform::initConsole();
        st.colorMode = !vt ? render::TerminalColorMode::Ascii
                           : (platform::supportsTrueColor() ? render::TerminalColorMode::TrueColor : render::TerminalColorMode::Ansi256);
    } else if (!render::parseTerminalColorMode(colorArg, st.colorMode)) {
        throw UsageError("--color must be truecolor, 256 or ascii");
    }
    const double targetFps = std::clamp(ctx.args().getDouble("target-fps", 30.0), 1.0, 120.0);

    auto sim = ctx.makeSimulation(ctx.epoch());
    render::SceneRenderer renderer;
    ctx.configureRenderer(renderer, 80);
    auto& opt = renderer.options();
    opt.labelCharWidth = 1.0;   // terminal cells: 1 column = 1 px, 1 row = 2 px
    opt.labelCharHeight = 2.0;
    opt.minBodyPixels = 1.0;
    opt.starBrightness = 1.0f;
    opt.particleBrightness = 1.8f;
    opt.orbitBrightness = 0.9f;
    render::PostSettings post = ctx.postSettings();
    post.vignette = false;
    post.bloomStrength = 0.35f;

    render::Camera cam;
    platform::TerminalSize ts = platform::terminalSize();
    const int headerRows = 1;
    auto footerRows = [&]() { return st.showHelp ? 4 : 1; };
    int fbW = std::max(20, ts.columns), fbH = std::max(10, (ts.rows - headerRows - footerRows()) * 2);
    ctx.configureCamera(cam, renderer, sim->system(), fbW, fbH);
    const double homeDistance = cam.distance();
    const double homeYaw = cam.yaw(), homePitch = cam.pitch();
    st.focus = std::max(0, ctx.focusIndex(sim->system()));
    render::TerminalRenderer term;
    term.setMode(st.colorMode);
    render::Framebuffer fb(fbW, fbH);
    render::Image8 img;
    renderer.record(*sim);

    std::unique_ptr<platform::RawTerminal> raw;
    if (tty && platform::stdinIsTerminal()) raw = std::make_unique<platform::RawTerminal>();

    Stopwatch clock;
    RateMeter fpsMeter(0.1);
    double last = clock.seconds();
    long long frame = 0;
    bool quit = false;
    auto say = [&](const std::string& m) {
        st.message = m;
        st.messageUntil = clock.seconds() + 3.0;
    };

    while (!quit && !platform::interruptRequested()) {
        // ---- input -----------------------------------------------------------------------
        if (raw) {
            for (int guard = 0; guard < 32; ++guard) {
                const platform::KeyEvent ev = raw->poll();
                if (ev.key == Key::None) break;
                const char c = ev.key == Key::Char ? ev.ch : 0;
                if (ev.key == Key::Escape || c == 'q' || c == 'Q') quit = true;
                else if (c == ' ') st.paused = !st.paused;
                else if (ev.key == Key::Left || c == 'a') cam.orbit(-5.0 * DEG2RAD, 0.0);
                else if (ev.key == Key::Right || c == 'd') cam.orbit(5.0 * DEG2RAD, 0.0);
                else if (ev.key == Key::Up || c == 'w') cam.orbit(0.0, 4.0 * DEG2RAD);
                else if (ev.key == Key::Down || c == 's') cam.orbit(0.0, -4.0 * DEG2RAD);
                else if (c == '+' || c == '=' || ev.key == Key::PageUp) cam.zoom(0.85);
                else if (c == '-' || c == '_' || ev.key == Key::PageDown) cam.zoom(1.0 / 0.85);
                else if (c == ']') { st.warp = std::min(st.warp * 2.0, 1e6); say(str::format("warp %.3g days/s", st.warp)); }
                else if (c == '[') { st.warp = std::max(st.warp / 2.0, 1e-4); say(str::format("warp %.3g days/s", st.warp)); }
                else if (c == 'r') { st.direction = -st.direction; say(st.direction > 0 ? "time runs forward" : "time runs BACKWARD"); }
                else if (ev.key == Key::Tab || (c >= '0' && c <= '9')) {
                    const auto& sys = sim->system();
                    if (ev.key == Key::Tab) {
                        st.focus = (st.focus + 1) % static_cast<int>(sys.namedCount());
                    } else {
                        static const char* byDigit[] = {"Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};
                        const int idx = sys.indexOf(byDigit[c - '0']);
                        if (idx >= 0) st.focus = idx;
                    }
                    cam.setDistance(frameDistance(renderer, sys, st.focus, homeDistance));
                    say("focus: " + sys.info(static_cast<std::size_t>(st.focus)).name);
                }
                else if (c == 'z' || ev.key == Key::Home) { st.focus = 0; cam.setDistance(homeDistance); cam.setAngles(homeYaw, homePitch); say("camera reset"); }
                else if (c == 't') { opt.trails = !opt.trails; say(opt.trails ? "trails on" : "trails off"); }
                else if (c == 'o') { opt.orbits = !opt.orbits; say(opt.orbits ? "orbits on" : "orbits off"); }
                else if (c == 'l') { opt.labels = !opt.labels; say(opt.labels ? "labels on" : "labels off"); }
                else if (c == 'g') { opt.grid = !opt.grid; say(opt.grid ? "grid on" : "grid off"); }
                else if (c == 'p') { opt.particles = !opt.particles; say(opt.particles ? "particles on" : "particles off"); }
                else if (c == 'b') { post.bloom = !post.bloom; say(post.bloom ? "bloom on" : "bloom off"); }
                else if (c == 'm') {
                    auto& m = renderer.mapping();
                    m.setMode(m.mode() == render::ScaleMode::Sqrt ? render::ScaleMode::Log
                                                                  : (m.mode() == render::ScaleMode::Log ? render::ScaleMode::Linear : render::ScaleMode::Sqrt));
                    say(std::string("radial scale: ") + render::scaleModeName(m.mode()));
                }
                else if (c == 'i') {
                    std::size_t k = 0;
                    for (; k < 7; ++k)
                        if (sim->integrator().name() == kIntegratorCycle[k]) break;
                    sim->setIntegrator(kIntegratorCycle[(k + 1) % 7]);
                    sim->rebaseInvariants();
                    say("integrator: " + sim->integrator().displayName() + " (energy error reset)");
                }
                else if (c == 'e') { sim->setRelativity(!sim->relativity()); sim->rebaseInvariants(); say(sim->relativity() ? "general relativity ON" : "general relativity off"); }
                else if (c == 'c') {
                    st.colorMode = st.colorMode == render::TerminalColorMode::TrueColor ? render::TerminalColorMode::Ansi256
                                   : st.colorMode == render::TerminalColorMode::Ansi256 ? render::TerminalColorMode::Ascii
                                                                                         : render::TerminalColorMode::TrueColor;
                    term.setMode(st.colorMode);
                    say(std::string("colour mode: ") + render::terminalColorModeName(st.colorMode));
                }
                else if (c == 'x') say(screenshot(*sim, renderer, cam, post));
                else if (c == 'h' || c == '?') st.showHelp = !st.showHelp;
            }
        }

        // ---- layout ----------------------------------------------------------------------
        ts = platform::terminalSize(); // falls back to $COLUMNS/$LINES when not a terminal
        const int rows = std::max(12, ts.rows), cols = std::max(30, ts.columns);
        const int w = cols, h = std::max(4, (rows - headerRows - footerRows()) * 2);
        if (w != fb.width() || h != fb.height()) {
            fb.resize(w, h);
            cam.setViewport(w, h);
        }

        // ---- simulation ---------------------------------------------------------------
        const double now = clock.seconds();
        const double realDt = std::min(0.1, now - last);
        last = now;
        if (!st.paused) {
            const double days = st.warp * st.direction * (maxFrames >= 0 && !tty ? 1.0 / targetFps : realDt);
            const double chunk = std::max(0.25, 4.0 * sim->timeStep());
            const int sub = std::clamp(static_cast<int>(std::ceil(std::fabs(days) / chunk)), 1, 64);
            for (int k = 0; k < sub; ++k) {
                sim->advance(days / sub);
                renderer.record(*sim);
            }
        }
        const auto& sys = sim->system();
        if (st.focus >= static_cast<int>(sys.namedCount())) st.focus = 0;
        cam.setTarget(st.focus > 0 ? renderer.displayPosition(sys, static_cast<std::size_t>(st.focus)) : Vec3{});
        opt.focus = st.focus > 0 ? st.focus : -1;

        // ---- render -------------------------------------------------------------------
        renderer.render(*sim, cam, fb);
        render::finishFrame(fb, img, post);
        const auto inv = sim->invariantErrors();
        std::vector<std::string> header = {str::format(
            " SolarSim | %s | %s%.3g d/s%s | %s dt=%.3gd | dE/E %+.1e | GR %s | focus %s | %.0f fps",
            astro::formatUTC(sim->time(), false).c_str(), st.direction < 0 ? "-" : "", st.warp, st.paused ? " (paused)" : "",
            sim->integrator().name().c_str(), sim->timeStep(), inv.energyRel, sim->relativity() ? "on" : "off",
            sys.info(static_cast<std::size_t>(st.focus)).name.c_str(), fpsMeter.value())};
        std::vector<std::string> footer;
        if (st.showHelp) {
            footer.push_back(" q/Esc quit | space pause | arrows/WASD rotate | +/- zoom | [ ] slower/faster | r reverse time | z reset");
            footer.push_back(" Tab next body | 0 Sun 1 Mercury 2 Venus 3 Earth 4 Mars 5 Jupiter 6 Saturn 7 Uranus 8 Neptune 9 Pluto");
            footer.push_back(" t trails | o orbits | l labels | g grid | p particles | b bloom | m scale | c colours | i integrator | e relativity");
            footer.push_back(" x screenshot (1920x1080 PNG) | h hide help");
        } else {
            footer.push_back(clock.seconds() < st.messageUntil ? " " + st.message
                                                              : " q quit | space pause | arrows rotate | +/- zoom | [ ] warp | Tab/0-9 focus | h help");
        }
        platform::writeStdout(term.compose(img, renderer.labels(), header, footer, cols, rows));

        // ---- pacing ---------------------------------------------------------------------
        ++frame;
        if (maxFrames >= 0 && frame >= maxFrames) break;
        const double frameTime = clock.seconds() - now;
        fpsMeter.sample(1.0 / std::max(1e-3, frameTime));
        const double sleep = 1.0 / targetFps - frameTime;
        if (sleep > 0.0 && tty) platform::sleepMs(static_cast<int>(sleep * 1000.0));
    }
    raw.reset(); // restore the terminal before printing
    if (!tty) std::printf("\n");
    const auto inv = sim->invariantErrors();
    std::printf("SolarSim: %lld frames, simulated %s .. %s (%.1f days), %llu steps, final dE/E = %+.2e\n", frame,
                astro::formatDateOnly(sim->epoch()).c_str(), astro::formatDateOnly(sim->time()).c_str(), sim->elapsed(),
                static_cast<unsigned long long>(sim->totalSteps()), inv.energyRel);
    return 0;
}

} // namespace solsim::app
