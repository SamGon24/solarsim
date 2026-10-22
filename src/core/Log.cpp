// =============================================================================
//  SolarSim - Log.cpp
// =============================================================================
#include "Log.hpp"

#include "StringUtil.hpp"

#include <cstdio>
#include <mutex>

namespace solsim::log {

namespace {
Level g_level = Level::Info;
bool g_color = false;
std::mutex g_mutex;

const char* tag(Level l) {
    switch (l) {
        case Level::Trace: return "trace";
        case Level::Debug: return "debug";
        case Level::Info: return "info";
        case Level::Warn: return "warn";
        case Level::Error: return "error";
        default: return "";
    }
}

const char* color(Level l) {
    switch (l) {
        case Level::Trace: return "\x1b[90m";
        case Level::Debug: return "\x1b[36m";
        case Level::Info: return "\x1b[32m";
        case Level::Warn: return "\x1b[33m";
        case Level::Error: return "\x1b[31m";
        default: return "";
    }
}
} // namespace

void setLevel(Level level) { g_level = level; }
Level level() { return g_level; }
void enableColor(bool on) { g_color = on; }

bool parseLevel(const std::string& name, Level& out) {
    const std::string n = str::toLower(name);
    if (n == "trace") out = Level::Trace;
    else if (n == "debug") out = Level::Debug;
    else if (n == "info") out = Level::Info;
    else if (n == "warn" || n == "warning") out = Level::Warn;
    else if (n == "error") out = Level::Error;
    else if (n == "silent" || n == "quiet" || n == "none") out = Level::Silent;
    else return false;
    return true;
}

void write(Level lvl, const std::string& message) {
    if (lvl < g_level || g_level == Level::Silent) return;
    std::lock_guard<std::mutex> lock(g_mutex);
    if (g_color)
        std::fprintf(stderr, "%s[%s]\x1b[0m %s\n", color(lvl), tag(lvl), message.c_str());
    else
        std::fprintf(stderr, "[%s] %s\n", tag(lvl), message.c_str());
    std::fflush(stderr);
}

} // namespace solsim::log
