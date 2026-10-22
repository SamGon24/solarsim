// =============================================================================
//  SolarSim - AppContext.hpp
//  Shared state and helpers for all sub-commands: the parsed command line and
//  the translation of options into BuildOptions / SimulationConfig / camera
//  and renderer settings.
// =============================================================================
#pragma once

#include "../core/ArgParser.hpp"
#include "../physics/Simulation.hpp"
#include "../physics/SystemBuilder.hpp"
#include "../render/Camera.hpp"
#include "../render/PostProcess.hpp"
#include "../render/SceneRenderer.hpp"

#include <memory>
#include <string>
#include <vector>

namespace solsim::app {

constexpr const char* kVersion = "1.0.0";

class AppContext {
public:
    explicit AppContext(ArgParser& args) : m_args(args) {}

    ArgParser& args() { return m_args; }
    const ArgParser& args() const { return m_args; }

    // --epoch (default: now). Returns JD(TT).
    double epoch() const;
    // Parses a date option (UTC string) or returns `fallback`.
    double dateOption(const std::string& name, double fallback) const;

    physics::BuildOptions buildOptions(double epochJD) const;
    physics::SimulationConfig simulationConfig() const;

    // Builds the system (prints warnings/notes depending on verbosity).
    physics::NBodySystem buildSystem(double epochJD, physics::BuildReport* report = nullptr) const;
    std::unique_ptr<physics::Simulation> makeSimulation(double epochJD) const;

    // Rendering configuration from the command line.
    void configureRenderer(render::SceneRenderer& renderer, int imageHeight) const;
    void configureCamera(render::Camera& cam, render::SceneRenderer& renderer, const physics::NBodySystem& sys,
                         int width, int height) const;
    render::PostSettings postSettings() const;
    // Index of --focus body (or -1).
    int focusIndex(const physics::NBodySystem& sys) const;

    bool verbose() const;
    bool quiet() const;

private:
    ArgParser& m_args;
};

// Registers every command and option with the parser.
void registerOptions(ArgParser& p);

// Commands (one translation unit each).
int cmdView(AppContext& ctx);
int cmdRender(AppContext& ctx);
int cmdInfo(AppContext& ctx);
int cmdEphemeris(AppContext& ctx);
int cmdEvents(AppContext& ctx);
int cmdMission(AppContext& ctx);
int cmdPrecession(AppContext& ctx);
int cmdLagrange(AppContext& ctx);
int cmdCompare(AppContext& ctx);
int cmdValidate(AppContext& ctx);
int cmdExport(AppContext& ctx);
int cmdSelfTest(AppContext& ctx);

// Shared console helpers.
void printBanner();
std::string progressBar(double fraction, int width = 30);
void printProgress(const std::string& label, double fraction);
void finishProgress();

} // namespace solsim::app
