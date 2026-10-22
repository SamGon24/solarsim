// =============================================================================
//  SolarSim - StringUtil.hpp
//  String helpers: trimming, splitting, case conversion, strict numeric parsing
//  and printf-style formatting into std::string.
// =============================================================================
#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace solsim::str {

std::string trim(const std::string& s);
std::string toLower(std::string s);
std::string toUpper(std::string s);
bool startsWith(const std::string& s, const std::string& prefix);
bool endsWith(const std::string& s, const std::string& suffix);
bool iequals(const std::string& a, const std::string& b);
std::vector<std::string> split(const std::string& s, char delim, bool skipEmpty = true);
std::string join(const std::vector<std::string>& parts, const std::string& sep);
std::string replaceAll(std::string s, const std::string& from, const std::string& to);

// Strict parsers: the whole string must be consumed.
bool parseDouble(const std::string& s, double& out);
bool parseInt(const std::string& s, long long& out);
bool parseBool(const std::string& s, bool& out);

// printf-style formatting.
std::string format(const char* fmt, ...)
#if defined(__MINGW32__) && !defined(__clang__)
    __attribute__((format(gnu_printf, 1, 2)))
#elif defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;

std::string padRight(const std::string& s, std::size_t width);
std::string padLeft(const std::string& s, std::size_t width);

// Human readable helpers.
std::string formatDuration(double seconds);           // "1h 02m 03.4s"
std::string formatBytes(double bytes);                // "1.25 MiB"
std::string formatAngleDMS(double degrees, bool sign = true);   // "+12°34'56.7\""
std::string formatHoursHMS(double hours);             // "12h34m56.78s"
std::string formatEngineering(double v, int sig = 4); // "1.234e+05" or "123.4"

// Number of UTF-8 code points (for table alignment with non-ASCII glyphs).
std::size_t utf8Length(const std::string& s);

// Levenshtein edit distance (used for "did you mean ...?" suggestions).
std::size_t editDistance(const std::string& a, const std::string& b);

} // namespace solsim::str
