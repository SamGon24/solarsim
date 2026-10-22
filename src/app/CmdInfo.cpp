// =============================================================================
//  SolarSim - CmdInfo.cpp
//  `solarsim info` - catalogue of bodies with physical and derived quantities.
// =============================================================================
#include "AppContext.hpp"

#include "../analysis/Lagrange.hpp"
#include "../astro/Coordinates.hpp"
#include "../astro/JulianDate.hpp"
#include "../astro/OrbitalElements.hpp"
#include "../core/Constants.hpp"
#include "../core/StringUtil.hpp"
#include "../core/Table.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>

namespace solsim::app {

using namespace solsim::constants;

namespace {

struct Derived {
    double massKg = 0, densityGcm3 = 0, surfaceG = 0, escapeKms = 0;
    astro::OrbitalElements el;
    double periodDays = 0;
    double hillAU = 0, soiAU = 0;
    double speedKms = 0;
    double teqK = 0;
};

Derived derive(const physics::NBodySystem& sys, std::size_t i, double jd) {
    Derived d;
    const auto& inf = sys.info(i);
    const double gmKm = sys.gm()[i] * AU3D2_TO_KM3S2;
    const double rKm = inf.radius * AU_KM;
    d.massKg = gmKm * 1e9 / G_SI;
    if (rKm > 0 && gmKm > 0) {
        const double vol = 4.0 / 3.0 * PI * std::pow(rKm * 1e5, 3.0); // cm^3
        d.densityGcm3 = d.massKg * 1e3 / vol;
        d.surfaceG = gmKm / (rKm * rKm) * 1e3;        // m/s^2
        d.escapeKms = std::sqrt(2.0 * gmKm / rKm);
    }
    const int p = inf.parent;
    if (p >= 0) {
        // A planet with massive moons orbits through the barycentre of its subsystem
        // (e.g. the Earth-Moon barycentre): use it for meaningful osculating elements.
        Vec3 r = sys.relativePosition(i, p), v = sys.relativeVelocity(i, p);
        double muSys = sys.twoBodyMu(i);
        if (p == 0) {
            double m = sys.gm()[i];
            Vec3 rw = r * m, vw = v * m;
            for (std::size_t c = 0; c < sys.namedCount(); ++c) {
                if (sys.info(c).parent != static_cast<int>(i) || sys.gm()[c] <= 0.0) continue;
                rw += (sys.pos(c) - sys.pos(0)) * sys.gm()[c];
                vw += (sys.vel(c) - sys.vel(0)) * sys.gm()[c];
                m += sys.gm()[c];
                muSys += sys.gm()[c];
            }
            if (m > 0.0) {
                r = rw / m;
                v = vw / m;
            }
        }
        d.el = astro::stateToElements(r, v, muSys, jd);
        d.periodDays = d.el.period();
        d.speedKms = v.norm() * AU_KM / DAY_S;
        if (sys.gm()[i] > 0 && sys.gm()[static_cast<std::size_t>(p)] > 0) {
            const double ratio = sys.gm()[i] / sys.gm()[static_cast<std::size_t>(p)];
            const double a = d.el.isElliptic() ? d.el.a() : r.norm();
            d.hillAU = a * (1.0 - std::min(d.el.e, 0.99)) * std::cbrt(ratio / 3.0);
            d.soiAU = a * std::pow(ratio, 0.4);
        }
        // Equilibrium temperature (fast rotator): T = 278.6 K (1 - A)^0.25 / sqrt(d_AU)
        const double dSun = (sys.pos(i) - sys.pos(0)).norm();
        if (dSun > 0) d.teqK = 278.6 * std::pow(1.0 - std::min(0.95, inf.albedo), 0.25) / std::sqrt(dSun);
    }
    return d;
}

std::string fmtPeriod(double days) {
    if (!std::isfinite(days)) return "unbound";
    if (days < 2.0) return str::format("%.2f h", days * 24.0);
    if (days < 1000.0) return str::format("%.2f d", days);
    return str::format("%.2f yr", days / JULIAN_YEAR_D);
}

std::string fmtSci(double v) {
    if (v == 0.0) return "-";
    return str::format("%.3e", v);
}

void bodyCard(const physics::NBodySystem& sys, std::size_t i, double jd) {
    const auto& inf = sys.info(i);
    const Derived d = derive(sys, i, jd);
    std::printf("\n%s (%s)\n", inf.name.c_str(), astro::kindName(inf.kind));
    if (!inf.notes.empty()) std::printf("  %s\n", inf.notes.c_str());
    std::printf("\n  Physical\n");
    std::printf("    GM                 %.6e km^3/s^2\n", sys.gm()[i] * AU3D2_TO_KM3S2);
    std::printf("    mass               %.4e kg  (%.6f Earth masses)\n", d.massKg, d.massKg / 5.9722e24);
    std::printf("    mean radius        %.1f km\n", inf.radius * AU_KM);
    if (d.densityGcm3 > 0) std::printf("    mean density       %.3f g/cm^3\n", d.densityGcm3);
    if (d.surfaceG > 0) std::printf("    surface gravity    %.3f m/s^2 (%.3f g)\n", d.surfaceG, d.surfaceG / 9.80665);
    if (d.escapeKms > 0) std::printf("    escape velocity    %.3f km/s\n", d.escapeKms);
    if (inf.rotationPeriod != 0)
        std::printf("    rotation period    %s%s\n", fmtPeriod(std::fabs(inf.rotationPeriod)).c_str(),
                    inf.rotationPeriod < 0 ? " (retrograde)" : "");
    const astro::Spherical pole = astro::toSpherical(astro::eclipticToEquatorial(inf.pole));
    std::printf("    north pole         RA %.2f deg, Dec %+.2f deg\n", pole.lon * RAD2DEG, pole.lat * RAD2DEG);
    if (inf.parent >= 0) {
        const Vec3 orbitNormal = sys.relativePosition(i, inf.parent).cross(sys.relativeVelocity(i, inf.parent)).normalized();
        std::printf("    obliquity          %.2f deg (spin axis vs orbit normal)\n", angleBetween(inf.pole, orbitNormal) * RAD2DEG);
    }
    if (d.teqK > 0) std::printf("    equilibrium temp.  %.0f K (albedo %.2f, no greenhouse)\n", d.teqK, inf.albedo);
    if (inf.rings.present)
        std::printf("    rings              %.0f - %.0f km\n", inf.rings.innerKm, inf.rings.outerKm);
    if (inf.parent >= 0) {
        std::printf("\n  Osculating orbit around %s at %s\n", sys.info(static_cast<std::size_t>(inf.parent)).name.c_str(),
                    astro::formatUTC(jd).c_str());
        const auto& el = d.el;
        if (el.isElliptic()) std::printf("    semi-major axis a  %.8f AU (%.0f km)\n", el.a(), el.a() * AU_KM);
        std::printf("    eccentricity e     %.8f\n", el.e);
        std::printf("    perihelion q       %.8f AU\n", el.q);
        if (el.isElliptic()) std::printf("    aphelion Q         %.8f AU\n", el.Q());
        std::printf("    inclination i      %.6f deg (to the ecliptic J2000)\n", el.i * RAD2DEG);
        std::printf("    asc. node          %.6f deg\n", el.node * RAD2DEG);
        std::printf("    arg. periapsis     %.6f deg\n", el.argp * RAD2DEG);
        std::printf("    true anomaly       %.6f deg\n", el.trueAnomalyAt(jd) * RAD2DEG);
        std::printf("    periapsis passage  %s\n", astro::formatUTC(el.tp, false).c_str());
        std::printf("    period             %s\n", fmtPeriod(d.periodDays).c_str());
        std::printf("    orbital speed      %.4f km/s\n", d.speedKms);
        if (d.hillAU > 0) std::printf("    Hill sphere        %.0f km\n", d.hillAU * AU_KM);
        if (d.soiAU > 0) std::printf("    sphere of infl.    %.0f km\n", d.soiAU * AU_KM);
    }
    if (inf.approximate) std::printf("\n  note: initial conditions are approximate for this body\n");
}

} // namespace

int cmdInfo(AppContext& ctx) {
    const double jd = ctx.epoch();
    const physics::NBodySystem sys = ctx.buildSystem(jd);

    if (ctx.args().has("body")) {
        const int i = sys.indexOf(ctx.args().getString("body"));
        if (i < 0) throw UsageError("unknown body '" + ctx.args().getString("body") + "' (try --moons / --minor)");
        bodyCard(sys, static_cast<std::size_t>(i), jd);
        return 0;
    }

    if (!ctx.quiet()) {
        printBanner();
        std::printf("Epoch %s (JD %.5f TT), %zu bodies", astro::formatUTC(jd).c_str(), jd, sys.namedCount());
        if (sys.particleCount()) std::printf(" + %zu test particles", sys.particleCount());
        std::printf("\n\n");
    }

    TextTable phys({"Body", "Kind", "Mass (kg)", "Radius (km)", "Density", "g (m/s2)", "v_esc (km/s)", "Rotation"});
    for (std::size_t i = 0; i < sys.namedCount(); ++i) {
        const auto& inf = sys.info(i);
        const Derived d = derive(sys, i, jd);
        phys.addRow({inf.name, astro::kindName(inf.kind), d.massKg > 0 ? fmtSci(d.massKg) : "~0",
                     str::format("%.1f", inf.radius * AU_KM), d.densityGcm3 > 0 ? str::format("%.3f", d.densityGcm3) : "-",
                     d.surfaceG > 0 ? str::format("%.2f", d.surfaceG) : "-",
                     d.escapeKms > 0 ? str::format("%.2f", d.escapeKms) : "-",
                     inf.rotationPeriod != 0 ? fmtPeriod(std::fabs(inf.rotationPeriod)) + (inf.rotationPeriod < 0 ? " R" : "") : "-"});
    }
    phys.setAlign(1, TextTable::Align::Left);
    phys.setTitle("Physical properties (R = retrograde rotation)");
    phys.print(std::cout, TextTable::defaultStyle());

    TextTable orb({"Body", "Around", "a (AU)", "e", "i (deg)", "Period", "Speed (km/s)", "Hill (km)", "SOI (km)", "T_eq (K)"});
    for (std::size_t i = 1; i < sys.namedCount(); ++i) {
        const auto& inf = sys.info(i);
        if (inf.parent < 0) continue;
        const Derived d = derive(sys, i, jd);
        const auto& el = d.el;
        orb.addRow({inf.name, sys.info(static_cast<std::size_t>(inf.parent)).name,
                    el.isElliptic() ? str::format("%.6f", el.a()) : str::format("q=%.3f", el.q), str::format("%.5f", el.e),
                    str::format("%.3f", el.i * RAD2DEG), fmtPeriod(d.periodDays), str::format("%.3f", d.speedKms),
                    d.hillAU > 0 ? str::format("%.0f", d.hillAU * AU_KM) : "-",
                    d.soiAU > 0 ? str::format("%.0f", d.soiAU * AU_KM) : "-",
                    d.teqK > 0 ? str::format("%.0f", d.teqK) : "-"});
    }
    orb.setAlign(1, TextTable::Align::Left);
    orb.setTitle("\nOsculating orbits at the epoch (relative to the parent; planets with moons: subsystem barycentre)");
    orb.print(std::cout, TextTable::defaultStyle());

    if (!sys.groups().empty()) {
        std::printf("\nTest-particle populations:\n");
        for (const auto& g : sys.groups()) std::printf("  %-18s %zu particles\n", g.name.c_str(), g.end - g.begin);
    }
    if (!ctx.quiet()) std::printf("\nTip: 'solarsim info --body Jupiter' for details, 'solarsim help' for all commands.\n");
    return 0;
}

} // namespace solsim::app
