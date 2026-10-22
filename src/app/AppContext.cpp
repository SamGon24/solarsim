// =============================================================================
//  SolarSim - AppContext.cpp
// =============================================================================
#include "AppContext.hpp"

#include "../astro/JulianDate.hpp"
#include "../core/Constants.hpp"
#include "../core/Log.hpp"
#include "../core/StringUtil.hpp"
#include "../io/BodyFileLoader.hpp"
#include "../platform/Terminal.hpp"

#include <cmath>
#include <cstdio>

namespace solsim::app {

using namespace solsim::constants;

void registerOptions(ArgParser& p) {
    // ---- commands -------------------------------------------------------------------
    p.addCommand("view", "Interactive real-time 3D view in the terminal (default)",
                 "Keys: q quit | space pause | arrows/WASD rotate | +/- zoom | [ ] time warp | r reverse time\n"
                 "      tab/0-9 focus body | t trails | o orbits | l labels | g grid | p particles | m scale mode\n"
                 "      i cycle integrator | e toggle relativity | c colour mode | x screenshot (PNG) | h help\n\n"
                 "Examples:\n  solarsim view\n  solarsim view --moons --asteroids 3000 --warp 30\n"
                 "  solarsim view --focus Jupiter --moons --warp 0.5\n  solarsim view --color ascii");
    p.addCommand("render", "Render an animation (GIF), image sequence, still image or SVG diagram",
                 "The output format follows the extension of --out: .gif .png .bmp .ppm .svg\n"
                 "A path ending in '/' writes a numbered PNG sequence into that directory.\n\n"
                 "Examples:\n  solarsim render --out solar.gif --frames 240 --days-per-frame 3\n"
                 "  solarsim render --out inner.gif --view inner --width 800 --height 600\n"
                 "  solarsim render --out saturn.png --focus Saturn --frames 1 --width 1920 --height 1080\n"
                 "  solarsim render --out orbits.svg --duration 4000\n"
                 "  solarsim render --out flyby.gif --intruder --kuiper 3000 --view outer --days-per-frame 400");
    p.addCommand("info", "Physical and orbital data of every body (or --body NAME for details)",
                 "Examples:\n  solarsim info\n  solarsim info --body Saturn\n  solarsim info --moons --minor");
    p.addCommand("ephem", "Apparent positions as seen from Earth (RA/Dec, distance, magnitude, alt/az)",
                 "Examples:\n  solarsim ephem\n  solarsim ephem --date 2027-02-19 --lat 51.48 --lon 0.0\n"
                 "  solarsim ephem --date \"2026-08-12 17:00\" --lat 40.4 --lon -3.7");
    p.addCommand("events", "Find equinoxes, solstices, lunar phases, oppositions, conjunctions, elongations",
                 "Examples:\n  solarsim events\n  solarsim events --from 2026-01-01 --to 2028-01-01 --types planets\n"
                 "  solarsim events --types moon --csv phases.csv");
    p.addCommand("mission", "Plan an interplanetary transfer (porkchop + N-body targeting)",
                 "Examples:\n  solarsim mission\n  solarsim mission --origin Earth --target Venus\n"
                 "  solarsim mission --porkchop porkchop.png --plot transfer.svg --csv trajectory.csv");
    p.addCommand("precession", "Measure perihelion precession (Mercury: Newtonian vs general relativity)",
                 "Examples:\n  solarsim precession\n  solarsim precession --body Venus --years 200");
    p.addCommand("lagrange", "Lagrange points, Hill spheres and a stability experiment",
                 "Examples:\n  solarsim lagrange\n  solarsim lagrange --pair Sun,Jupiter --test-years 200");
    p.addCommand("compare", "Benchmark all integrators: energy error, accuracy and speed",
                 "Examples:\n  solarsim compare\n  solarsim compare --years 100 --dt 1.0");
    p.addCommand("validate", "Compare simulated events against real-world dates and reference models",
                 "Example:\n  solarsim validate");
    p.addCommand("export", "Export state vectors / diagnostics as CSV",
                 "Examples:\n  solarsim export --out states.csv --duration 3650 --every 10\n"
                 "  solarsim export --out earth_moon.csv --only Earth,Moon --frame helio --every 0.5");
    p.addCommand("selftest", "Run built-in sanity checks", "Example:\n  solarsim selftest");
    p.addCommand("help", "Show help (solarsim help <command>)");
    p.addCommand("version", "Print version information");

    // ---- simulation -----------------------------------------------------------------
    const std::string S = "Simulation options";
    p.addOption("epoch", "DATE", "start date (UTC): YYYY-MM-DD[THH:MM], JD<n>, J2000, now", S, "now");
    p.addOption("integrator", "NAME", "euler | symplectic-euler | leapfrog | rk4 | yoshida4 | dopri5 | wh", S, "yoshida4");
    p.addOption("dt", "DAYS", "integration time step", S, "auto");
    p.addOption("tolerance", "RTOL", "relative tolerance of adaptive integrators", S, "1e-10");
    p.addFlag("gr", "add the general-relativity (1PN) correction from the Sun", S);
    p.addFlag("moons", "add Io, Europa, Ganymede, Callisto, Titan and Triton", S);
    p.addFlag("minor", "add Ceres, Vesta and comets Halley & Encke", S);
    p.addFlag("no-moon", "omit Earth's Moon (Earth then represents the Earth-Moon barycentre)", S);
    p.addOption("only", "LIST", "simulate only these bodies (comma separated; parents added)", S);
    p.addOption("exclude", "LIST", "remove these bodies", S);
    p.addOption("bodies", "FILES", "load extra bodies from INI file(s) (see data/)", S);
    p.addOption("asteroids", "N", "main-belt test particles", S, "0");
    p.addOption("trojans", "N", "Jupiter Trojan test particles", S, "0");
    p.addOption("kuiper", "N", "Kuiper-belt test particles", S, "0");
    p.addOption("seed", "N", "random seed for particle populations", S, "20261021");
    p.addFlag("intruder", "add a passing star on a hyperbolic flyby", S);
    p.addOption("intruder-mass", "MSUN", "intruder mass in solar masses", S, "0.5");
    p.addOption("intruder-q", "AU", "intruder closest approach to the Sun", S, "80");
    p.addOption("intruder-vinf", "KM/S", "intruder speed at infinity", S, "3");
    p.addOption("intruder-days", "DAYS", "days from the epoch to closest approach", S, "36525");
    p.addOption("calibrate", "MODE", "initial-condition calibration: none | moon | all", S, "all");
    p.addOption("softening", "AU", "Plummer softening length", S, "0");
    p.addOption("threads", "N", "worker threads for particle forces (0 = all cores)", S, "0");

    // ---- rendering ------------------------------------------------------------------
    const std::string R = "Rendering options";
    p.addOption("out", "FILE", "output file (render/export/screenshots)", R);
    p.addOption("width", "PX", "image width", R, "960");
    p.addOption("height", "PX", "image height", R, "540");
    p.addOption("frames", "N", "number of frames to render", R, "240");
    p.addOption("days-per-frame", "DAYS", "simulated days between frames", R, "3");
    p.addOption("duration", "DAYS", "total simulated days (overrides --days-per-frame)", R);
    p.addOption("fps", "N", "playback rate of animations", R, "30");
    p.addOption("view", "PRESET", "overview | inner | outer | top | side", R, "overview");
    p.addOption("focus", "BODY", "follow a body with the camera", R);
    p.addOption("yaw", "DEG", "camera azimuth", R);
    p.addOption("pitch", "DEG", "camera elevation above the ecliptic", R);
    p.addOption("distance", "UNITS", "camera distance (display units, 1 = 1 AU mapped)", R);
    p.addOption("fov", "DEG", "vertical field of view", R, "45");
    p.addOption("orbit-camera", "DEG", "rotate the camera by this many degrees per frame", R, "0");
    p.addOption("scale", "MODE", "radial distance mapping: linear | sqrt | log", R, "sqrt");
    p.addFlag("true-scale", "physical body sizes (use with --scale linear)", R);
    p.addOption("body-scale", "X", "multiply displayed body sizes", R, "1");
    p.addOption("exposure", "X", "tone-mapping exposure", R, "1");
    p.addFlag("grid", "draw ecliptic distance rings", R);
    p.addFlag("no-orbits", "hide osculating orbits", R);
    p.addFlag("no-trails", "hide trails", R);
    p.addFlag("no-labels", "hide labels", R);
    p.addFlag("no-stars", "hide the star field", R);
    p.addFlag("no-bloom", "disable bloom", R);
    p.addFlag("no-hud", "hide the text overlay", R);
    p.addFlag("dither", "Floyd-Steinberg dithering for GIF output", R);
    p.addOption("warmup", "DAYS", "integrate this long before the first frame (fills trails)", R, "0");

    // ---- terminal view ---------------------------------------------------------------
    const std::string V = "Interactive view options";
    p.addOption("warp", "DAYS/S", "simulated days per real second", V, "10");
    p.addOption("color", "MODE", "truecolor | 256 | ascii", V, "auto");
    p.addOption("target-fps", "N", "frame rate limit", V, "30");
    p.addOption("max-frames", "N", "exit after N frames (also allows non-interactive use)", V);

    // ---- analysis ---------------------------------------------------------------------
    const std::string A = "Analysis options";
    p.addOption("date", "DATE", "date for 'ephem' (UTC)", A, "now");
    p.addOption("lat", "DEG", "observer latitude (enables alt/az)", A);
    p.addOption("lon", "DEG", "observer longitude, east positive", A);
    p.addOption("alt-m", "M", "observer altitude in metres", A, "0");
    p.addOption("from", "DATE", "start date for 'events'", A, "now");
    p.addOption("to", "DATE", "end date for 'events'", A, "from + 1 year");
    p.addOption("types", "LIST", "events: seasons,moon,planets,apsides", A, "all");
    p.addOption("body", "NAME", "body for 'info', 'precession' and 'events'", A);
    p.addOption("years", "N", "length of 'precession' / 'compare' runs", A);
    p.addOption("origin", "BODY", "mission departure body", A, "Earth");
    p.addOption("target", "BODY", "mission destination", A, "Mars");
    p.addOption("park-alt", "KM", "departure parking-orbit altitude", A, "300");
    p.addOption("capture-alt", "KM", "arrival periapsis altitude", A, "400");
    p.addOption("porkchop", "FILE", "write the porkchop plot (PNG)", A);
    p.addOption("plot", "FILE", "write an SVG plot (mission, precession)", A);
    p.addOption("csv", "FILE", "write results as CSV", A);
    p.addOption("pair", "A,B", "primary,secondary for 'lagrange'", A);
    p.addOption("test-years", "N", "integrate test particles at L1..L5 for N years", A, "0");
    p.addOption("every", "DAYS", "output interval for 'export'", A, "1");
    p.addOption("frame", "NAME", "export frame: bary | helio", A, "bary");
    p.addOption("integrators", "LIST", "integrators for 'compare'", A, "all");

    // ---- general -----------------------------------------------------------------------
    const std::string G = "General options";
    p.addFlag("help", "show help", G, 'h');
    p.addFlag("verbose", "more output", G, 'v');
    p.addFlag("quiet", "less output", G, 'q');
    p.addFlag("ascii", "plain ASCII tables (no box drawing)", G);
    p.addFlag("markdown", "Markdown tables", G);
}

// -----------------------------------------------------------------------------
bool AppContext::verbose() const { return m_args.flag("verbose"); }
bool AppContext::quiet() const { return m_args.flag("quiet"); }

double AppContext::dateOption(const std::string& name, double fallback) const {
    if (!m_args.has(name)) return fallback;
    double jd = 0.0;
    std::string err;
    if (!astro::parseDate(m_args.getString(name), jd, err)) throw UsageError("--" + name + ": " + err);
    return jd;
}

double AppContext::epoch() const { return dateOption("epoch", astro::utcToTT(astro::nowJDUtc())); }

physics::BuildOptions AppContext::buildOptions(double epochJD) const {
    physics::BuildOptions o;
    o.epochJD = epochJD;
    o.moons = m_args.flag("moons");
    o.minor = m_args.flag("minor");
    o.includeMoon = !m_args.flag("no-moon");
    o.only = m_args.getList("only");
    o.exclude = m_args.getList("exclude");
    for (const auto& file : m_args.getList("bodies")) {
        const io::LoadResult lr = io::loadBodyFile(file);
        if (!lr.ok()) throw UsageError("--bodies: " + lr.error);
        for (const auto& w : lr.warnings) log::warn(w);
        o.extra.insert(o.extra.end(), lr.bodies.begin(), lr.bodies.end());
        if (!quiet()) log::info("loaded " + std::to_string(lr.bodies.size()) + " body definition(s) from " + file);
    }
    o.populations.asteroids = static_cast<int>(m_args.getInt("asteroids", 0));
    o.populations.trojans = static_cast<int>(m_args.getInt("trojans", 0));
    o.populations.kuiper = static_cast<int>(m_args.getInt("kuiper", 0));
    o.populations.seed = static_cast<std::uint64_t>(m_args.getInt("seed", 20261021));
    if (o.populations.asteroids < 0 || o.populations.trojans < 0 || o.populations.kuiper < 0)
        throw UsageError("particle counts must be >= 0");
    if (o.populations.asteroids + o.populations.trojans + o.populations.kuiper > 2000000)
        throw UsageError("too many particles (limit 2,000,000)");
    o.intruder.enabled = m_args.flag("intruder");
    o.intruder.massSolar = m_args.getDouble("intruder-mass", 0.5);
    o.intruder.perihelionAU = m_args.getDouble("intruder-q", 80.0);
    o.intruder.vInfKms = m_args.getDouble("intruder-vinf", 3.0);
    o.intruder.daysToPerihelion = m_args.getDouble("intruder-days", 36525.0);
    if (o.intruder.enabled && (o.intruder.massSolar <= 0 || o.intruder.perihelionAU <= 0 || o.intruder.vInfKms <= 0))
        throw UsageError("intruder mass, perihelion and speed must be positive");
    if (m_args.has("calibrate") && !physics::parseCalibrationMode(m_args.getString("calibrate"), o.calibration))
        throw UsageError("--calibrate must be none, moon or all");
    return o;
}

physics::SimulationConfig AppContext::simulationConfig() const {
    physics::SimulationConfig c;
    c.integrator = m_args.getString("integrator", "yoshida4");
    if (!physics::createIntegrator(c.integrator))
        throw UsageError("unknown integrator '" + c.integrator + "' (choose: " + str::join(physics::integratorNames(), ", ") + ")");
    c.dt = m_args.getDouble("dt", 0.0);
    if (m_args.has("dt") && !(c.dt > 0.0)) throw UsageError("--dt must be positive");
    c.tolerance = m_args.getDouble("tolerance", 1e-10);
    if (!(c.tolerance > 0.0 && c.tolerance < 1.0)) throw UsageError("--tolerance must be in (0, 1)");
    c.relativity = m_args.flag("gr");
    c.softening = m_args.getDouble("softening", 0.0);
    const long long th = m_args.getInt("threads", 0);
    if (th < 0 || th > 256) throw UsageError("--threads must be between 0 and 256");
    c.threads = static_cast<unsigned>(th);
    return c;
}

physics::NBodySystem AppContext::buildSystem(double epochJD, physics::BuildReport* reportOut) const {
    physics::BuildReport report;
    physics::NBodySystem sys = physics::buildSystem(buildOptions(epochJD), &report);
    for (const auto& w : report.warnings) log::warn(w);
    if (verbose())
        for (const auto& n : report.notes) log::info(n);
    if (reportOut) *reportOut = report;
    return sys;
}

std::unique_ptr<physics::Simulation> AppContext::makeSimulation(double epochJD) const {
    auto sim = std::make_unique<physics::Simulation>(buildSystem(epochJD), epochJD, simulationConfig());
    if (verbose())
        log::info(str::format("simulation: %zu bodies (%zu massive, %zu particles), %s, dt = %.4f d", sim->system().size(),
                              sim->system().massiveCount(), sim->system().particleCount(),
                              sim->integrator().displayName().c_str(), sim->timeStep()));
    return sim;
}

int AppContext::focusIndex(const physics::NBodySystem& sys) const {
    if (!m_args.has("focus")) return -1;
    const int i = sys.indexOf(m_args.getString("focus"));
    if (i < 0) throw UsageError("--focus: unknown body '" + m_args.getString("focus") + "'");
    return i;
}

void AppContext::configureRenderer(render::SceneRenderer& r, int imageHeight) const {
    auto& o = r.options();
    o.grid = m_args.flag("grid");
    o.orbits = !m_args.flag("no-orbits");
    o.trails = !m_args.flag("no-trails");
    o.labels = !m_args.flag("no-labels");
    o.stars = !m_args.flag("no-stars");
    o.textScale = std::max(1, imageHeight / 400);
    o.labelCharWidth = 6.0 * o.textScale;
    o.labelCharHeight = 8.0 * o.textScale;
    o.minBodyPixels = std::max(2.0, imageHeight / 220.0);
    render::ScaleMode mode = render::ScaleMode::Sqrt;
    if (m_args.has("scale") && !render::parseScaleMode(m_args.getString("scale"), mode))
        throw UsageError("--scale must be linear, sqrt or log");
    r.mapping().setMode(mode);
    r.mapping().setTrueScale(m_args.flag("true-scale"));
    const double bs = m_args.getDouble("body-scale", 1.0);
    if (!(bs > 0.0)) throw UsageError("--body-scale must be positive");
    r.mapping().setBodyScale(bs);
}

render::PostSettings AppContext::postSettings() const {
    render::PostSettings p;
    p.bloom = !m_args.flag("no-bloom");
    p.exposure = static_cast<float>(m_args.getDouble("exposure", 1.0));
    return p;
}

void AppContext::configureCamera(render::Camera& cam, render::SceneRenderer& r, const physics::NBodySystem& sys,
                                 int width, int height) const {
    cam.setViewport(width, height);
    cam.setFieldOfView(m_args.getDouble("fov", 45.0) * DEG2RAD);
    r.mapping().prepare(sys);

    // Size of the system in display units: farthest planet from the Sun.
    double outer = 1.0, inner = 1.0;
    for (std::size_t i = 1; i < sys.namedCount(); ++i) {
        if (sys.info(i).parent != 0 || sys.info(i).kind == astro::BodyKind::Star) continue;
        const double d = r.mapping().radial((sys.pos(i) - sys.pos(0)).norm());
        if (sys.info(i).kind == astro::BodyKind::Planet || sys.info(i).kind == astro::BodyKind::DwarfPlanet) {
            outer = std::max(outer, d);
            if ((sys.pos(i) - sys.pos(0)).norm() < 2.0) inner = std::max(inner, d);
        }
    }
    const std::string view = str::toLower(m_args.getString("view", "overview"));
    double yaw = -100.0, pitch = 28.0, dist = outer * 2.4;
    if (view == "inner") { dist = inner * 3.0; pitch = 35.0; }
    else if (view == "outer" || view == "overview") { dist = outer * 2.4; }
    else if (view == "top") { pitch = 89.0; dist = outer * 2.6; }
    else if (view == "side" || view == "edge") { pitch = 2.0; dist = outer * 2.6; }
    else throw UsageError("--view must be overview, inner, outer, top or side");

    Vec3 target{};
    const int focus = focusIndex(sys);
    if (focus >= 0) {
        target = r.displayPosition(sys, static_cast<std::size_t>(focus));
        double rad = r.displayRadius(sys, static_cast<std::size_t>(focus));
        // Frame the satellite system if there is one.
        double span = rad * 14.0;
        for (std::size_t i = 0; i < sys.namedCount(); ++i)
            if (sys.info(i).parent == focus)
                span = std::max(span, 2.2 * (r.displayPosition(sys, i) - target).norm());
        if (focus == 0) span = std::max(span, outer * 0.5);
        dist = span;
        pitch = 22.0;
    }
    cam.setTarget(target);
    cam.setDistance(m_args.getDouble("distance", dist));
    cam.setAngles(m_args.getDouble("yaw", yaw) * DEG2RAD, m_args.getDouble("pitch", pitch) * DEG2RAD);
}

// -----------------------------------------------------------------------------
void printBanner() {
    std::printf("SolarSim %s - N-body Solar System simulator (C++17)\n", kVersion);
}

std::string progressBar(double f, int width) {
    f = std::max(0.0, std::min(1.0, f));
    const int full = static_cast<int>(f * width);
    std::string s = "[";
    for (int i = 0; i < width; ++i) s += i < full ? '#' : (i == full ? '>' : '.');
    s += "]";
    return s + str::format(" %5.1f%%", f * 100.0);
}

namespace {
// Progress bars use carriage returns: only draw them on an interactive stderr.
bool progressEnabled() {
    static const bool tty = platform::stderrIsTerminal();
    return tty && log::level() <= log::Level::Info;
}
} // namespace

void printProgress(const std::string& label, double fraction) {
    if (!progressEnabled()) return;
    std::fprintf(stderr, "\r%s %s", label.c_str(), progressBar(fraction).c_str());
    std::fflush(stderr);
}

void finishProgress() {
    if (!progressEnabled()) return;
    std::fprintf(stderr, "\n");
}

} // namespace solsim::app
