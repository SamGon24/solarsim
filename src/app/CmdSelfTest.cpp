// =============================================================================
//  SolarSim - CmdSelfTest.cpp
//  `solarsim selftest` - fast built-in checks of the numerical building blocks
//  (the full unit-test suite lives in tests/ and is built as solarsim_tests).
// =============================================================================
#include "AppContext.hpp"

#include "../astro/JulianDate.hpp"
#include "../astro/Kepler.hpp"
#include "../astro/Lambert.hpp"
#include "../astro/LunarTheory.hpp"
#include "../astro/OrbitalElements.hpp"
#include "../core/Constants.hpp"
#include "../core/StringUtil.hpp"
#include "../io/Checksum.hpp"
#include "../io/Deflate.hpp"
#include "../io/GifEncoder.hpp"

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace solsim::app {

using namespace solsim::constants;

int cmdSelfTest(AppContext& ctx) {
    struct T {
        const char* name;
        std::function<bool(std::string&)> fn;
    };
    const std::vector<T> tests = {
        {"Julian date of J2000 epoch", [](std::string& d) {
             const double jd = astro::calendarToJD(2000, 1, 1.5);
             d = str::format("%.4f", jd);
             return std::fabs(jd - 2451545.0) < 1e-9;
         }},
        {"Julian date of Sputnik launch (Meeus 7.a)", [](std::string& d) {
             const double jd = astro::calendarToJD(1957, 10, 4.81);
             d = str::format("%.2f", jd);
             return std::fabs(jd - 2436116.31) < 1e-6;
         }},
        {"calendar round trip", [](std::string& d) {
             const auto c = astro::jdToCalendar(2460604.75);
             d = str::format("%04d-%02d-%02d %02d:%02d", c.year, c.month, c.day, c.hour, c.minute);
             return c.year == 2024 && c.month == 10 && c.day == 21 && c.hour == 6;
         }},
        {"Kepler equation e = 0.9", [](std::string& d) {
             const double E = astro::solveKeplerElliptic(0.3, 0.9);
             d = str::format("residual %.1e", E - 0.9 * std::sin(E) - 0.3);
             return std::fabs(E - 0.9 * std::sin(E) - 0.3) < 1e-13;
         }},
        {"hyperbolic Kepler e = 2.5", [](std::string& d) {
             const double H = astro::solveKeplerHyperbolic(10.0, 2.5);
             d = str::format("residual %.1e", 2.5 * std::sinh(H) - H - 10.0);
             return std::fabs(2.5 * std::sinh(H) - H - 10.0) < 1e-10;
         }},
        {"elements <-> state round trip", [](std::string& d) {
             const auto el = astro::OrbitalElements::fromMeanAnomaly(2.5, 0.3, 0.4, 1.0, 2.0, 0.7, J2000_JD, GM_SUN);
             const auto s = astro::elementsToState(el, J2000_JD + 123.0);
             const auto el2 = astro::stateToElements(s.r, s.v, GM_SUN, J2000_JD + 123.0);
             d = str::format("da = %.1e AU", std::fabs(el2.a() - 2.5));
             return std::fabs(el2.a() - 2.5) < 1e-10 && std::fabs(el2.e - 0.3) < 1e-12 && std::fabs(el2.tp - el.tp) < 1e-6;
         }},
        {"universal Kepler propagation (1 orbit)", [](std::string& d) {
             const Vec3 r{1, 0, 0}, v{0, std::sqrt(GM_SUN) * 1.1, 0};
             const auto el = astro::stateToElements(r, v, GM_SUN, 0.0);
             const auto res = astro::propagateKepler(r, v, GM_SUN, el.period());
             d = str::format("error %.2e AU", (res.r - r).norm());
             return (res.r - r).norm() < 1e-11;
         }},
        {"Lambert vs Kepler consistency", [](std::string& d) {
             const Vec3 r1{1, 0, 0}, v1{0.002, 0.0185, 0.0005};
             const auto p = astro::propagateKepler(r1, v1, GM_SUN, 150.0);
             const auto l = astro::solveLambert(r1, p.r, 150.0, GM_SUN, true);
             d = str::format("dv = %.1e AU/d", (l.v1 - v1).norm());
             return l.ok && (l.v1 - v1).norm() < 1e-9;
         }},
        {"ELP lunar distance in range", [](std::string& d) {
             double mn = 1e9, mx = 0;
             for (int k = 0; k < 400; ++k) {
                 const double r = astro::lunarPositionOfDate(J2000_JD + k * 0.5).distanceKm;
                 mn = std::min(mn, r);
                 mx = std::max(mx, r);
             }
             d = str::format("%.0f..%.0f km", mn, mx);
             return mn > 356000 && mx < 407000;
         }},
        {"CRC-32 check value", [](std::string& d) {
             const char* s = "123456789";
             const auto c = io::Crc32::compute(reinterpret_cast<const std::uint8_t*>(s), 9);
             d = str::format("0x%08X", c);
             return c == 0xCBF43926u;
         }},
        {"Adler-32 check value", [](std::string& d) {
             const char* s = "Wikipedia";
             const auto a = io::Adler32::compute(reinterpret_cast<const std::uint8_t*>(s), 9);
             d = str::format("0x%08X", a);
             return a == 0x11E60398u;
         }},
        {"DEFLATE compresses repetitive data", [](std::string& d) {
             std::vector<std::uint8_t> data(20000);
             for (std::size_t i = 0; i < data.size(); ++i) data[i] = static_cast<std::uint8_t>((i % 37) * 7);
             const auto z = io::zlibCompress(data.data(), data.size());
             d = str::format("%zu -> %zu bytes", data.size(), z.size());
             return z.size() < data.size() / 20;
         }},
        {"GIF LZW stream sanity", [](std::string& d) {
             std::vector<std::uint8_t> idx(5000);
             for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = static_cast<std::uint8_t>((i / 7) % 256);
             const auto out = io::gifLzwEncode(idx, 8);
             d = str::format("%zu bytes", out.size());
             return !out.empty() && out.size() < idx.size();
         }},
        {"N-body energy conservation (1 yr)", [&ctx](std::string& d) {
             const auto sys = ctx.buildSystem(J2000_JD);
             physics::Simulation sim(sys, J2000_JD);
             sim.advance(365.25);
             const double e = sim.invariantErrors().energyRel;
             d = str::format("dE/E = %.1e", e);
             return std::fabs(e) < 1e-9;
         }},
    };
    int failed = 0;
    for (const auto& t : tests) {
        std::string detail;
        bool ok = false;
        try {
            ok = t.fn(detail);
        } catch (const std::exception& e) {
            detail = std::string("exception: ") + e.what();
        }
        failed += !ok;
        std::printf("  [%s] %-42s %s\n", ok ? "PASS" : "FAIL", t.name, detail.c_str());
    }
    std::printf("%zu tests, %d failed\n", tests.size(), failed);
    return failed ? 1 : 0;
}

} // namespace solsim::app
