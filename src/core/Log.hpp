// =============================================================================
//  SolarSim - Log.hpp
//  Minimal leveled logger writing to stderr (so stdout stays clean for tables
//  and CSV output that may be piped into other tools).
// =============================================================================
#pragma once

#include <string>

namespace solsim::log {

enum class Level { Trace = 0, Debug, Info, Warn, Error, Silent };

void setLevel(Level level);
Level level();
bool parseLevel(const std::string& name, Level& out);
void enableColor(bool on);

void write(Level level, const std::string& message);

inline void trace(const std::string& m) { write(Level::Trace, m); }
inline void debug(const std::string& m) { write(Level::Debug, m); }
inline void info(const std::string& m) { write(Level::Info, m); }
inline void warn(const std::string& m) { write(Level::Warn, m); }
inline void error(const std::string& m) { write(Level::Error, m); }

} // namespace solsim::log
