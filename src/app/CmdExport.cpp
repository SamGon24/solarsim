// =============================================================================
//  SolarSim - CmdExport.cpp
//  `solarsim export` - time series of state vectors (and optional conserved-
//  quantity diagnostics) in CSV form for external analysis / plotting.
// =============================================================================
#include "AppContext.hpp"

#include "../astro/JulianDate.hpp"
#include "../core/Constants.hpp"
#include "../core/StringUtil.hpp"
#include "../core/Timer.hpp"
#include "../io/CsvWriter.hpp"

#include <cmath>
#include <cstdio>

namespace solsim::app {

using namespace solsim::constants;

int cmdExport(AppContext& ctx) {
    const std::string path = ctx.args().getString("out", "states.csv");
    const double duration = ctx.args().getDouble("duration", 365.25);
    const double every = ctx.args().getDouble("every", 1.0);
    if (!(duration > 0.0) || !(every > 0.0)) throw UsageError("--duration and --every must be positive");
    if (duration / every > 5e6) throw UsageError("too many output samples (reduce --duration or increase --every)");
    const std::string frame = str::toLower(ctx.args().getString("frame", "bary"));
    if (frame != "bary" && frame != "helio") throw UsageError("--frame must be 'bary' or 'helio'");

    auto sim = ctx.makeSimulation(ctx.epoch());
    io::CsvWriter csv;
    if (!csv.open(path)) throw std::runtime_error("cannot write " + path);
    csv.header({"jd_tt", "utc", "body", "x_au", "y_au", "z_au", "vx_au_per_day", "vy_au_per_day", "vz_au_per_day",
                "distance_au", "speed_km_s"});
    io::CsvWriter diag;
    if (ctx.args().has("csv")) {
        if (!diag.open(ctx.args().getString("csv"))) throw std::runtime_error("cannot write " + ctx.args().getString("csv"));
        diag.header({"jd_tt", "days", "energy_rel_error", "angular_momentum_rel_error", "virial_ratio", "steps", "force_evaluations"});
    }

    const std::size_t n = sim->system().namedCount();
    const std::size_t samples = static_cast<std::size_t>(std::floor(duration / every + 1e-9)) + 1;
    Stopwatch sw;
    for (std::size_t k = 0; k < samples; ++k) {
        if (k > 0) sim->advance(every);
        const auto& s = sim->system();
        const std::string utc = astro::formatISO(sim->time());
        for (std::size_t i = 0; i < n; ++i) {
            Vec3 r = s.pos(i), v = s.vel(i);
            if (frame == "helio") {
                r -= s.pos(0);
                v -= s.vel(0);
            }
            csv.row({io::CsvWriter::num(sim->time(), 13), utc, s.info(i).name, io::CsvWriter::num(r.x, 15),
                     io::CsvWriter::num(r.y, 15), io::CsvWriter::num(r.z, 15), io::CsvWriter::num(v.x, 15),
                     io::CsvWriter::num(v.y, 15), io::CsvWriter::num(v.z, 15), io::CsvWriter::num(r.norm(), 12),
                     io::CsvWriter::num(v.norm() * AU_KM / DAY_S, 10)});
        }
        if (diag.isOpen()) {
            const auto e = sim->invariantErrors();
            diag.row({io::CsvWriter::num(sim->time(), 13), io::CsvWriter::num(sim->elapsed()), io::CsvWriter::num(e.energyRel),
                      io::CsvWriter::num(e.angularMomentumRel), io::CsvWriter::num(e.virialRatio),
                      std::to_string(sim->totalSteps()), std::to_string(sim->forceEvaluations())});
        }
        if (!ctx.quiet() && (k % 64 == 0 || k + 1 == samples)) printProgress("exporting", double(k + 1) / samples);
    }
    if (!ctx.quiet()) finishProgress();
    std::printf("wrote %zu rows (%zu bodies x %zu samples, %s frame) to %s in %s\n", csv.rows() - 1, n, samples,
                frame == "bary" ? "barycentric" : "heliocentric", path.c_str(), str::formatDuration(sw.seconds()).c_str());
    if (diag.isOpen()) std::printf("wrote diagnostics to %s\n", ctx.args().getString("csv").c_str());
    return 0;
}

} // namespace solsim::app
