// =============================================================================
//  SolarSim - CmdEvents.cpp
//  `solarsim events` - astronomical calendar computed by the N-body model.
// =============================================================================
#include "AppContext.hpp"

#include "../analysis/Events.hpp"
#include "../astro/JulianDate.hpp"
#include "../core/Constants.hpp"
#include "../core/StringUtil.hpp"
#include "../core/Table.hpp"
#include "../core/Timer.hpp"
#include "../io/CsvWriter.hpp"

#include <cstdio>
#include <iostream>

namespace solsim::app {

using namespace solsim::constants;

int cmdEvents(AppContext& ctx) {
    const double from = ctx.dateOption("from", astro::utcToTT(astro::nowJDUtc()));
    const double to = ctx.dateOption("to", from + 365.25);
    if (!(to > from)) throw UsageError("--to must be after --from");
    if (to - from > 200.0 * JULIAN_YEAR_D) throw UsageError("search span limited to 200 years");
    const double epochJD = ctx.args().has("epoch") ? ctx.epoch() : from;
    if (epochJD > from) throw UsageError("--epoch must not be after --from");

    analysis::EventSearchOptions opt;
    opt.startJD = from;
    opt.endJD = to;
    if (ctx.args().has("types")) {
        opt.seasons = opt.planets = opt.moonPhases = opt.apsides = false;
        for (const auto& t : ctx.args().getList("types")) {
            const std::string k = str::toLower(t);
            if (k == "seasons" || k == "season") opt.seasons = true;
            else if (k == "moon" || k == "phases") opt.moonPhases = true;
            else if (k == "planets" || k == "planet") opt.planets = true;
            else if (k == "apsides" || k == "apsis") opt.apsides = true;
            else if (k == "all") opt.seasons = opt.planets = opt.moonPhases = opt.apsides = true;
            else throw UsageError("unknown event type '" + t + "' (seasons, moon, planets, apsides)");
        }
    }
    if (ctx.args().has("body")) opt.bodies = {ctx.args().getString("body")};

    auto sim = ctx.makeSimulation(epochJD);
    Stopwatch sw;
    const auto events = analysis::findEvents(*sim, opt, [&](double f) {
        if (!ctx.quiet()) printProgress("searching", f);
    });
    if (!ctx.quiet()) {
        printProgress("searching", 1.0);
        finishProgress();
        std::printf("%zu events between %s and %s (%.2f s, %s)\n\n", events.size(), astro::formatDateOnly(from).c_str(),
                    astro::formatDateOnly(to).c_str(), sw.seconds(), sim->integrator().displayName().c_str());
    }

    TextTable t({"Date (UTC)", "Day", "Type", "Event"});
    t.setAlign(1, TextTable::Align::Left);
    t.setAlign(2, TextTable::Align::Left);
    t.setAlign(3, TextTable::Align::Left);
    for (const auto& e : events) {
        const std::string day = std::string(astro::dayOfWeekName(astro::dayOfWeek(astro::ttToUTC(e.jd)))).substr(0, 3);
        t.addRow({astro::formatUTC(e.jd, false).substr(0, 16), day, e.category, e.description});
    }
    t.print(std::cout, TextTable::defaultStyle());

    if (ctx.args().has("csv")) {
        io::CsvWriter csv;
        const std::string path = ctx.args().getString("csv");
        if (!csv.open(path)) throw std::runtime_error("cannot write " + path);
        csv.header({"utc", "jd_tt", "category", "body", "description", "value"});
        for (const auto& e : events)
            csv.row({astro::formatISO(e.jd), io::CsvWriter::num(e.jd, 12), e.category, e.body, e.description,
                     io::CsvWriter::num(e.value)});
        std::printf("wrote %zu rows to %s\n", csv.rows() - 1, path.c_str());
    }
    if (!ctx.quiet())
        std::printf("Times are computed by the N-body integration (typ. accuracy: lunar phases ~1-5 min, seasons ~10 min,\n"
                    "planetary events ~1 h). Run 'solarsim validate' to compare against published dates.\n");
    return 0;
}

} // namespace solsim::app
