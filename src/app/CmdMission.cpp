// =============================================================================
//  SolarSim - CmdMission.cpp
//  `solarsim mission` - plans an interplanetary transfer and renders the
//  porkchop plot (PNG) and a trajectory diagram (SVG).
// =============================================================================
#include "AppContext.hpp"

#include "../analysis/Transfer.hpp"
#include "../astro/JulianDate.hpp"
#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"
#include "../core/StringUtil.hpp"
#include "../core/Timer.hpp"
#include "../io/CsvWriter.hpp"
#include "../io/ImageWriter.hpp"
#include "../io/SvgWriter.hpp"
#include "../render/Font.hpp"
#include "../render/Image.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace solsim::app {

using namespace solsim::constants;

namespace {

// Perceptually ordered colour map (viridis-like anchors).
void colormap(double t, std::uint8_t out[3]) {
    static const double anchors[][3] = {{68, 1, 84}, {59, 82, 139}, {33, 145, 140}, {94, 201, 98}, {253, 231, 37}};
    t = math::clamp(t, 0.0, 1.0) * 4.0;
    const int i = std::min(3, static_cast<int>(t));
    const double f = t - i;
    for (int c = 0; c < 3; ++c)
        out[c] = static_cast<std::uint8_t>(anchors[i][c] + (anchors[i + 1][c] - anchors[i][c]) * f);
}

void fillRect(render::Image8& img, int x0, int y0, int x1, int y1, std::uint8_t r, std::uint8_t g, std::uint8_t b) {
    for (int y = std::max(0, y0); y <= std::min(img.height - 1, y1); ++y)
        for (int x = std::max(0, x0); x <= std::min(img.width - 1, x1); ++x) {
            std::uint8_t* p = img.pixel(x, y);
            p[0] = r; p[1] = g; p[2] = b;
        }
}

bool writePorkchop(const analysis::MissionResult& res, const std::string& from, const std::string& to,
                   const std::string& path, std::string& err) {
    const auto& pc = res.porkchop;
    if (!pc.valid()) {
        err = "no porkchop data";
        return false;
    }
    const int W = 1000, H = 640, L = 90, R = 150, T = 60, B = 80;
    const int pw = W - L - R, ph = H - T - B;
    render::Image8 img;
    img.resize(W, H);
    fillRect(img, 0, 0, W - 1, H - 1, 12, 15, 24);
    const std::size_t nL = pc.launchJD.size(), nT = pc.tofDays.size();
    const double best = pc.bestC3;
    const double cmax = std::max(best * 4.0, best + 25.0);
    auto value = [&](double fx, double fy) {
        // bilinear interpolation in the grid (fx: launch index, fy: tof index)
        const std::size_t x0 = static_cast<std::size_t>(math::clamp(std::floor(fx), 0.0, double(nL - 2)));
        const std::size_t y0 = static_cast<std::size_t>(math::clamp(std::floor(fy), 0.0, double(nT - 2)));
        const double tx = math::clamp(fx - double(x0), 0.0, 1.0), ty = math::clamp(fy - double(y0), 0.0, 1.0);
        auto at = [&](std::size_t x, std::size_t y) { return pc.c3[y * nL + x]; };
        const double a = at(x0, y0), b = at(x0 + 1, y0), c = at(x0, y0 + 1), d = at(x0 + 1, y0 + 1);
        if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(d))
            return std::numeric_limits<double>::quiet_NaN();
        return (a * (1 - tx) + b * tx) * (1 - ty) + (c * (1 - tx) + d * tx) * ty;
    };
    // Contour levels (km^2/s^2).
    std::vector<double> levels;
    for (double lv : {1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 16.0, 20.0, 30.0, 40.0, 60.0, 80.0, 120.0})
        levels.push_back(best + lv);
    auto band = [&](double v) {
        int k = 0;
        while (k < static_cast<int>(levels.size()) && v > levels[static_cast<std::size_t>(k)]) ++k;
        return k;
    };
    std::vector<int> bands(static_cast<std::size_t>(pw * ph), -1);
    for (int py = 0; py < ph; ++py) {
        for (int px = 0; px < pw; ++px) {
            const double fx = double(px) / (pw - 1) * double(nL - 1);
            const double fy = double(ph - 1 - py) / (ph - 1) * double(nT - 1);
            const double v = value(fx, fy);
            std::uint8_t c[3] = {30, 30, 40};
            if (std::isfinite(v)) {
                if (v <= cmax) colormap(1.0 - std::log(v / best) / std::log(cmax / best), c);
                else { c[0] = 40; c[1] = 20; c[2] = 50; }
                bands[static_cast<std::size_t>(py * pw + px)] = band(v);
            }
            std::uint8_t* p = img.pixel(L + px, T + py);
            p[0] = c[0]; p[1] = c[1]; p[2] = c[2];
        }
    }
    // Contour lines where the band index changes.
    for (int py = 1; py < ph; ++py)
        for (int px = 1; px < pw; ++px) {
            const int b0 = bands[static_cast<std::size_t>(py * pw + px)];
            if (b0 < 0) continue;
            const int bx = bands[static_cast<std::size_t>(py * pw + px - 1)], by = bands[static_cast<std::size_t>((py - 1) * pw + px)];
            if ((bx >= 0 && bx != b0) || (by >= 0 && by != b0)) {
                std::uint8_t* p = img.pixel(L + px, T + py);
                p[0] = static_cast<std::uint8_t>(p[0] / 3 + 150);
                p[1] = static_cast<std::uint8_t>(p[1] / 3 + 150);
                p[2] = static_cast<std::uint8_t>(p[2] / 3 + 150);
            }
        }
    // Best point marker.
    const double bx = L + (res.launchJD - pc.launchJD.front()) / (pc.launchJD.back() - pc.launchJD.front()) * (pw - 1);
    const double by = T + (ph - 1) - (res.timeOfFlight - pc.tofDays.front()) / (pc.tofDays.back() - pc.tofDays.front()) * (ph - 1);
    for (int k = -7; k <= 7; ++k) {
        fillRect(img, static_cast<int>(bx) + k, static_cast<int>(by) - 1, static_cast<int>(bx) + k, static_cast<int>(by) + 1, 255, 60, 60);
        fillRect(img, static_cast<int>(bx) - 1, static_cast<int>(by) + k, static_cast<int>(bx) + 1, static_cast<int>(by) + k, 255, 60, 60);
    }
    // Axes and labels.
    using render::BitmapFont;
    const int s = 2;
    BitmapFont::draw(img, L, 16, from + " -> " + to + "  departure C3 (km^2/s^2), Lambert arcs", 0xE8EEF6, s);
    for (int k = 0; k <= 6; ++k) {
        const double jd = pc.launchJD.front() + (pc.launchJD.back() - pc.launchJD.front()) * k / 6.0;
        const int x = L + static_cast<int>((pw - 1) * k / 6.0);
        fillRect(img, x, T + ph, x, T + ph + 6, 200, 200, 210);
        const std::string lbl = astro::formatDateOnly(jd).substr(2);
        BitmapFont::draw(img, x - BitmapFont::textWidth(lbl, s) / 2, T + ph + 12, lbl, 0xC8D2E0, s);
    }
    BitmapFont::draw(img, L + pw / 2 - 100, H - 26, "launch date (UTC)", 0x9FB0C8, s);
    for (int k = 0; k <= 5; ++k) {
        const double tof = pc.tofDays.front() + (pc.tofDays.back() - pc.tofDays.front()) * k / 5.0;
        const int y = T + ph - 1 - static_cast<int>((ph - 1) * k / 5.0);
        fillRect(img, L - 6, y, L, y, 200, 200, 210);
        const std::string lbl = str::format("%.0f", tof);
        BitmapFont::draw(img, L - 10 - BitmapFont::textWidth(lbl, s), y - 7, lbl, 0xC8D2E0, s);
    }
    BitmapFont::draw(img, 8, T - 22, "TOF (days)", 0x9FB0C8, s);
    // Colour bar.
    const int cbx = L + pw + 30;
    for (int y = 0; y < ph; ++y) {
        std::uint8_t c[3];
        colormap(1.0 - double(y) / (ph - 1), c);
        fillRect(img, cbx, T + y, cbx + 22, T + y, c[0], c[1], c[2]);
    }
    BitmapFont::draw(img, cbx + 28, T - 4, str::format("%.1f", best), 0xC8D2E0, s);
    BitmapFont::draw(img, cbx + 28, T + ph - 12, str::format("%.0f", cmax), 0xC8D2E0, s);
    BitmapFont::draw(img, L + 8, T + 8,
                     str::format("best: %s  TOF %.0f d  C3 %.2f", astro::formatDateOnly(res.launchJD).c_str(),
                                 res.timeOfFlight, pc.bestC3),
                     0xFFFFFF, s);
    return io::writeImage(path, img, err);
}

bool writeMissionSvg(const analysis::MissionResult& res, const std::string& from, const std::string& to,
                     const std::string& path, std::string& err) {
    const double W = 900, H = 1000, C = 450;
    double extent = 0.1;
    for (const auto& p : res.departOrbit) extent = std::max(extent, std::hypot(p.x, p.y));
    for (const auto& p : res.targetOrbit) extent = std::max(extent, std::hypot(p.x, p.y));
    for (const auto& p : res.trajectory) extent = std::max(extent, std::hypot(p.x, p.y));
    const double scale = 380.0 / extent;
    auto P = [&](const Vec3& v) { return std::make_pair(C + v.x * scale, C - v.y * scale); };
    io::SvgWriter svg(W, H);
    svg.text(C, 34, from + " → " + to + " transfer (ecliptic plane, view from north)", "#e8eef6", 18, "middle", "bold");
    for (double r = 0.5; r < extent * 1.05; r += 0.5)
        svg.circle(C, C, r * scale, "none", "#1f2a40", 1.0);
    std::vector<std::pair<double, double>> pts;
    for (const auto& p : res.departOrbit) pts.push_back(P(p));
    svg.polyline(pts, "#4dabf7", 1.5, 0.8, true);
    pts.clear();
    for (const auto& p : res.targetOrbit) pts.push_back(P(p));
    svg.polyline(pts, "#ff8a5c", 1.5, 0.8, true);
    pts.clear();
    for (const auto& p : res.trajectory) pts.push_back(P(p));
    svg.polyline(pts, "#7cf2c2", 2.5, 1.0, false);
    svg.circle(C, C, 9, "#ffd43b");
    svg.text(C + 12, C - 10, "Sun", "#ffd43b", 13);
    const auto lp = P(res.launchPosition), ap = P(res.arrivalPosition);
    svg.circle(lp.first, lp.second, 6, "#4dabf7");
    svg.text(lp.first + 10, lp.second + 4, from + " at launch " + astro::formatDateOnly(res.launchJD), "#4dabf7", 13);
    svg.circle(ap.first, ap.second, 6, "#ff8a5c");
    svg.text(ap.first + 10, ap.second + 4, to + " at arrival " + astro::formatDateOnly(res.closestApproachJD), "#ff8a5c", 13);
    double y = H - 120;
    for (const std::string& line :
         {str::format("time of flight %.1f days, transfer angle %.1f deg", res.timeOfFlight, res.transferAngleDeg),
          str::format("departure: C3 %.2f km²/s², v∞ %.3f km/s, DLA %+.1f deg, Δv %.3f km/s", res.c3,
                      res.vInfDepart, res.declinationDeg, res.dvDepart),
          str::format("arrival: v∞ %.3f km/s, periapsis %.0f km altitude, capture Δv %.3f km/s", res.vInfArrive,
                      res.periapsisAltKm, res.dvCapture),
          str::format("total Δv %.3f km/s (N-body verified)", res.dvDepart + res.dvCapture)}) {
        svg.text(40, y, line, "#c8d2e0", 14);
        y += 22;
    }
    svg.text(W - 40, H - 20, str::format("grid: 0.5 AU"), "#5a6a88", 12, "end");
    return svg.save(path, err);
}

} // namespace

int cmdMission(AppContext& ctx) {
    const double start = ctx.args().has("from") ? ctx.dateOption("from", 0) : ctx.epoch();
    analysis::MissionConfig cfg;
    cfg.from = ctx.args().getString("origin", "Earth");
    cfg.to = ctx.args().getString("target", "Mars");
    cfg.searchStartJD = start;
    cfg.parkingAltitudeKm = ctx.args().getDouble("park-alt", 300.0);
    cfg.captureAltitudeKm = ctx.args().getDouble("capture-alt", 400.0);
    if (cfg.parkingAltitudeKm < 100.0 || cfg.captureAltitudeKm < 50.0)
        throw UsageError("parking/capture altitudes are too low (min 100 km / 50 km)");

    physics::BuildOptions bo = ctx.buildOptions(start);
    bo.populations = {};
    physics::SimulationConfig sc = ctx.simulationConfig();
    physics::Simulation sim(physics::buildSystem(bo), start, sc);
    if (sim.system().indexOf(cfg.from) < 0 || sim.system().indexOf(cfg.to) < 0)
        throw UsageError("unknown origin/target body");

    if (!ctx.quiet())
        std::printf("Mission design %s -> %s, searching from %s\n\n", cfg.from.c_str(), cfg.to.c_str(),
                    astro::formatDateOnly(start).c_str());
    Stopwatch sw;
    const analysis::MissionResult res = analysis::planMission(sim, cfg, [&](const std::string& m) {
        if (!ctx.quiet()) std::printf("  %s\n", m.c_str());
    });
    if (res.launchJD == 0.0) {
        std::fprintf(stderr, "mission planning failed: %s\n", res.message.c_str());
        return 1;
    }

    const double t = res.timeOfFlight;
    std::printf("\nResult (%s, %s):\n", res.message.c_str(), str::formatDuration(sw.seconds()).c_str());
    std::printf("  launch                %s\n", astro::formatUTC(res.launchJD, false).c_str());
    std::printf("  arrival (periapsis)   %s   time of flight %.1f days (%.2f months)\n",
                astro::formatUTC(res.closestApproachJD, false).c_str(), t, t / 30.44);
    std::printf("  transfer angle        %.1f deg  (Hohmann would be 180)\n", res.transferAngleDeg);
    std::printf("  departure C3          %.3f km^2/s^2   (v_inf %.3f km/s, DLA %+.2f deg)\n", res.c3, res.vInfDepart,
                res.declinationDeg);
    std::printf("  departure burn        %.3f km/s from a %.0f km circular parking orbit\n", res.dvDepart, cfg.parkingAltitudeKm);
    std::printf("  arrival v_inf         %.3f km/s\n", res.vInfArrive);
    std::printf("  periapsis altitude    %.1f km (requested %.0f km) via %.0f km B-plane offset%s\n", res.periapsisAltKm,
                cfg.captureAltitudeKm, res.aimOffsetKm, res.impact ? "  ** IMPACT **" : "");
    std::printf("  capture burn          %.3f km/s into a circular orbit\n", res.dvCapture);
    std::printf("  total delta-v         %.3f km/s\n", res.dvDepart + res.dvCapture);
    std::printf("  analytic Hohmann      %.3f + %.3f km/s heliocentric, TOF %.1f d\n", res.analytic.dv1, res.analytic.dv2,
                res.analytic.timeOfFlight);

    std::string err;
    if (ctx.args().has("porkchop")) {
        if (!writePorkchop(res, cfg.from, cfg.to, ctx.args().getString("porkchop"), err)) throw std::runtime_error(err);
        std::printf("wrote porkchop plot %s\n", ctx.args().getString("porkchop").c_str());
    }
    if (ctx.args().has("plot")) {
        if (!writeMissionSvg(res, cfg.from, cfg.to, ctx.args().getString("plot"), err)) throw std::runtime_error(err);
        std::printf("wrote trajectory diagram %s\n", ctx.args().getString("plot").c_str());
    }
    if (ctx.args().has("csv")) {
        io::CsvWriter csv(ctx.args().getString("csv"));
        csv.header({"sample", "x_au", "y_au", "z_au", "r_au"});
        for (std::size_t i = 0; i < res.trajectory.size(); ++i) {
            const Vec3& p = res.trajectory[i];
            csv.row({std::to_string(i), io::CsvWriter::num(p.x), io::CsvWriter::num(p.y), io::CsvWriter::num(p.z),
                     io::CsvWriter::num(p.norm())});
        }
        std::printf("wrote trajectory %s\n", ctx.args().getString("csv").c_str());
    }
    return res.success ? 0 : 1;
}

} // namespace solsim::app
