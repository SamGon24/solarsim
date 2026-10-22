// =============================================================================
//  SolarSim - CmdEphemeris.cpp
//  `solarsim ephem` - the sky as seen from Earth at a given date.
// =============================================================================
#include "AppContext.hpp"

#include "../analysis/Observation.hpp"
#include "../astro/Coordinates.hpp"
#include "../astro/JulianDate.hpp"
#include "../core/Constants.hpp"
#include "../core/StringUtil.hpp"
#include "../core/Table.hpp"
#include "../io/CsvWriter.hpp"

#include <cstdio>
#include <iostream>

namespace solsim::app {

using namespace solsim::constants;

int cmdEphemeris(AppContext& ctx) {
    const double date = ctx.dateOption("date", astro::utcToTT(astro::nowJDUtc()));
    // Integrate from --epoch (default: the requested date itself).
    const double epochJD = ctx.args().has("epoch") ? ctx.epoch() : date;
    auto sim = ctx.makeSimulation(epochJD);
    if (std::fabs(date - epochJD) > 1e-9) sim->advanceTo(date);
    const auto& sys = sim->system();
    const int earth = sys.indexOf("Earth");
    if (earth < 0) throw UsageError("ephem needs the Earth in the simulation");

    analysis::ObserverSite site;
    if (ctx.args().has("lat") || ctx.args().has("lon")) {
        site.enabled = true;
        site.latitudeDeg = ctx.args().getDouble("lat", 0.0);
        site.longitudeDeg = ctx.args().getDouble("lon", 0.0);
        site.altitudeM = ctx.args().getDouble("alt-m", 0.0);
        if (std::fabs(site.latitudeDeg) > 90.0) throw UsageError("--lat must be within [-90, 90]");
    }

    if (!ctx.quiet()) {
        std::printf("Apparent positions for %s (JD %.5f TT)\n", astro::formatUTC(date).c_str(), date);
        if (site.enabled) {
            const double lst = astro::gmstDegrees(astro::ttToUTC(date)) + site.longitudeDeg;
            std::printf("Observer: lat %+.4f, lon %+.4f, alt %.0f m - local sidereal time %s\n", site.latitudeDeg,
                        site.longitudeDeg, site.altitudeM, str::formatHoursHMS(lst / 15.0).c_str());
        } else {
            std::printf("Geocentric (add --lat/--lon for topocentric altitude/azimuth)\n");
        }
        if (std::fabs(date - epochJD) > 1e-9)
            std::printf("Integrated %.1f days from epoch %s\n", date - epochJD, astro::formatUTC(epochJD, false).c_str());
        std::printf("RA/Dec: mean equator & equinox of date; positions corrected for light time.\n\n");
    }

    std::vector<std::string> headers = {"Body", "RA", "Dec", "Dist (AU)", "Light time", "r (AU)", "Elong",
                                        "Phase", "Illum", "Mag", "Diam"};
    if (site.enabled) {
        headers.push_back("Alt");
        headers.push_back("Az");
    }
    TextTable t(headers);
    io::CsvWriter csv;
    if (ctx.args().has("csv")) {
        if (!csv.open(ctx.args().getString("csv"))) throw std::runtime_error("cannot write " + ctx.args().getString("csv"));
        csv.header({"body", "ra_deg", "dec_deg", "ra_j2000_deg", "dec_j2000_deg", "distance_au", "r_au", "elongation_deg",
                    "phase_deg", "illuminated", "magnitude", "diameter_arcsec", "alt_deg", "az_deg"});
    }

    for (std::size_t i = 0; i < sys.namedCount(); ++i) {
        if (static_cast<int>(i) == earth) continue;
        const auto& inf = sys.info(i);
        if (inf.kind == astro::BodyKind::Moon && inf.name != "Moon") continue; // outer-planet moons: skip
        const analysis::Observation ob = analysis::observe(sys, i, static_cast<std::size_t>(earth), date, site);
        std::string diam;
        if (ob.angularDiameter >= 60.0) diam = str::format("%.1f'", ob.angularDiameter / 60.0);
        else diam = str::format("%.2f\"", ob.angularDiameter);
        std::vector<std::string> row = {
            inf.name,
            str::formatHoursHMS(ob.raDate / 15.0),
            str::formatAngleDMS(ob.decDate),
            ob.distance < 0.1 ? str::format("%.6f", ob.distance) : str::format("%.4f", ob.distance),
            ob.lightTimeMin < 60 ? str::format("%.2f min", ob.lightTimeMin) : str::format("%.2f h", ob.lightTimeMin / 60),
            i == 0 ? "-" : str::format("%.4f", ob.helioDistance),
            i == 0 ? "-" : str::format("%.1f%s", ob.elongation, ob.eveningSky ? "E" : "W"),
            i == 0 ? "-" : str::format("%.1f", ob.phaseAngle),
            i == 0 ? "-" : str::format("%.0f%%", ob.illuminated * 100.0),
            ob.hasMagnitude ? str::format("%+.1f", ob.magnitude) : "-",
            diam};
        if (site.enabled) {
            row.push_back(str::format("%+.1f", ob.altitude));
            row.push_back(str::format("%.1f %s", ob.azimuth, astro::compassPoint(ob.azimuth * DEG2RAD)));
        }
        t.addRow(row);
        if (csv.isOpen())
            csv.row({inf.name, io::CsvWriter::num(ob.raDate), io::CsvWriter::num(ob.decDate), io::CsvWriter::num(ob.raJ2000),
                     io::CsvWriter::num(ob.decJ2000), io::CsvWriter::num(ob.distance), io::CsvWriter::num(ob.helioDistance),
                     io::CsvWriter::num(ob.elongation), io::CsvWriter::num(ob.phaseAngle), io::CsvWriter::num(ob.illuminated),
                     ob.hasMagnitude ? io::CsvWriter::num(ob.magnitude, 4) : "", io::CsvWriter::num(ob.angularDiameter),
                     site.enabled ? io::CsvWriter::num(ob.altitude) : "", site.enabled ? io::CsvWriter::num(ob.azimuth) : ""});
    }
    t.print(std::cout, TextTable::defaultStyle());
    if (!ctx.quiet()) {
        std::printf("Elong: angle from the Sun (E = evening sky, W = morning sky). Mag: visual magnitude.\n");
        if (site.enabled) std::printf("Alt includes standard atmospheric refraction; negative = below the horizon.\n");
    }
    if (csv.isOpen()) std::printf("wrote %s\n", ctx.args().getString("csv").c_str());
    return 0;
}

} // namespace solsim::app
