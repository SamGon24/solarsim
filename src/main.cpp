// =============================================================================
//  SolarSim - main.cpp
//  Entry point: parses the command line and dispatches to a sub-command.
//
//      solarsim [command] [options]
//
//  With no command, an interactive terminal session starts the real-time
//  viewer; otherwise (output redirected) a summary table is printed.
// =============================================================================
#include "app/AppContext.hpp"
#include "core/ArgParser.hpp"
#include "core/Log.hpp"
#include "core/Table.hpp"
#include "platform/Terminal.hpp"

#include <cstdio>
#include <exception>
#include <functional>
#include <map>
#include <string>

using namespace solsim;

int main(int argc, char** argv) {
    const bool vt = platform::initConsole();
    platform::installInterruptHandler();

    ArgParser parser("solarsim");
    app::registerOptions(parser);
    try {
        parser.parse(argc, argv);
    } catch (const UsageError& e) {
        std::fprintf(stderr, "error: %s\nRun 'solarsim --help' for usage.\n", e.what());
        return 2;
    }

    if (parser.flag("quiet")) log::setLevel(log::Level::Warn);
    if (parser.flag("verbose")) log::setLevel(log::Level::Debug);
    log::enableColor(vt && platform::stdoutIsTerminal());
    if (parser.flag("ascii") || !vt) TextTable::setDefaultStyle(TextTable::Style::Ascii);
    if (parser.flag("markdown")) TextTable::setDefaultStyle(TextTable::Style::Markdown);

    std::string command = parser.command();
    if (command.empty()) {
        if (parser.flag("help")) command = "help";
        else command = (platform::stdoutIsTerminal() && platform::stdinIsTerminal()) ? "view" : "info";
    }

    if (command == "help" || parser.flag("help")) {
        if (command != "help") {
            std::printf("%s", parser.commandHelp(command).c_str());
        } else if (!parser.positionals().empty()) {
            std::printf("%s", parser.commandHelp(parser.positionals().front()).c_str());
        } else {
            app::printBanner();
            std::printf("\n%s", parser.helpText().c_str());
        }
        return 0;
    }
    if (command == "version" || command == "--version") {
        app::printBanner();
        std::printf("Built with %s, C++ standard %ld\n",
#if defined(__clang__)
                    "clang " __clang_version__,
#elif defined(__GNUC__)
                    "GCC " __VERSION__,
#elif defined(_MSC_VER)
                    "MSVC",
#else
                    "unknown compiler",
#endif
                    static_cast<long>(__cplusplus));
        return 0;
    }

    const std::map<std::string, std::function<int(app::AppContext&)>> commands = {
        {"view", app::cmdView},           {"render", app::cmdRender},     {"info", app::cmdInfo},
        {"ephem", app::cmdEphemeris},     {"events", app::cmdEvents},     {"mission", app::cmdMission},
        {"precession", app::cmdPrecession}, {"lagrange", app::cmdLagrange}, {"compare", app::cmdCompare},
        {"validate", app::cmdValidate},   {"export", app::cmdExport},     {"selftest", app::cmdSelfTest},
    };
    const auto it = commands.find(command);
    if (it == commands.end()) {
        std::fprintf(stderr, "error: unknown command '%s'\nRun 'solarsim --help' for the list of commands.\n",
                     command.c_str());
        return 2;
    }

    app::AppContext ctx(parser);
    try {
        return it->second(ctx);
    } catch (const UsageError& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "fatal: %s\n", e.what());
        return 1;
    }
}
