// =============================================================================
//  SolarSim - StringUtil.cpp
// =============================================================================
#include "StringUtil.hpp"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

namespace solsim::str {

std::string trim(const std::string& s) {
    std::size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::string toUpper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return s;
}

bool startsWith(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

bool endsWith(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool iequals(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (std::tolower(static_cast<unsigned char>(a[i])) != std::tolower(static_cast<unsigned char>(b[i])))
            return false;
    return true;
}

std::vector<std::string> split(const std::string& s, char delim, bool skipEmpty) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == delim) {
            if (!skipEmpty || !cur.empty()) out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!skipEmpty || !cur.empty()) out.push_back(cur);
    return out;
}

std::string join(const std::vector<std::string>& parts, const std::string& sep) {
    std::string out;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (i) out += sep;
        out += parts[i];
    }
    return out;
}

std::string replaceAll(std::string s, const std::string& from, const std::string& to) {
    if (from.empty()) return s;
    std::size_t pos = 0;
    while ((pos = s.find(from, pos)) != std::string::npos) {
        s.replace(pos, from.size(), to);
        pos += to.size();
    }
    return s;
}

bool parseDouble(const std::string& sIn, double& out) {
    const std::string s = trim(sIn);
    if (s.empty()) return false;
    errno = 0;
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    if (errno == ERANGE || end != s.c_str() + s.size()) return false;
    out = v;
    return true;
}

bool parseInt(const std::string& sIn, long long& out) {
    const std::string s = trim(sIn);
    if (s.empty()) return false;
    errno = 0;
    char* end = nullptr;
    const long long v = std::strtoll(s.c_str(), &end, 10);
    if (errno == ERANGE || end != s.c_str() + s.size()) return false;
    out = v;
    return true;
}

bool parseBool(const std::string& sIn, bool& out) {
    const std::string s = toLower(trim(sIn));
    if (s == "1" || s == "true" || s == "yes" || s == "on") { out = true; return true; }
    if (s == "0" || s == "false" || s == "no" || s == "off") { out = false; return true; }
    return false;
}

std::string format(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    va_list copy;
    va_copy(copy, args);
    const int n = std::vsnprintf(nullptr, 0, fmt, copy);
    va_end(copy);
    std::string out;
    if (n > 0) {
        out.resize(static_cast<std::size_t>(n) + 1);
        std::vsnprintf(&out[0], out.size(), fmt, args);
        out.resize(static_cast<std::size_t>(n));
    }
    va_end(args);
    return out;
}

std::size_t utf8Length(const std::string& s) {
    std::size_t n = 0;
    for (unsigned char c : s)
        if ((c & 0xC0) != 0x80) ++n;
    return n;
}

std::string padRight(const std::string& s, std::size_t width) {
    const std::size_t len = utf8Length(s);
    return len >= width ? s : s + std::string(width - len, ' ');
}

std::string padLeft(const std::string& s, std::size_t width) {
    const std::size_t len = utf8Length(s);
    return len >= width ? s : std::string(width - len, ' ') + s;
}

std::string formatDuration(double seconds) {
    if (!std::isfinite(seconds)) return "n/a";
    if (seconds < 1e-3) return format("%.1f us", seconds * 1e6);
    if (seconds < 1.0) return format("%.1f ms", seconds * 1e3);
    if (seconds < 60.0) return format("%.2f s", seconds);
    const long long total = static_cast<long long>(seconds);
    const long long h = total / 3600, m = (total % 3600) / 60;
    const double s = seconds - double(h * 3600 + m * 60);
    if (h > 0) return format("%lldh %02lldm %04.1fs", h, m, s);
    return format("%lldm %04.1fs", m, s);
}

std::string formatBytes(double bytes) {
    const char* units[] = {"B", "KiB", "MiB", "GiB", "TiB"};
    int u = 0;
    while (bytes >= 1024.0 && u < 4) { bytes /= 1024.0; ++u; }
    return u == 0 ? format("%.0f %s", bytes, units[u]) : format("%.2f %s", bytes, units[u]);
}

std::string formatAngleDMS(double degrees, bool sign) {
    const char s = degrees < 0 ? '-' : '+';
    double a = std::fabs(degrees);
    int d = static_cast<int>(a);
    a = (a - d) * 60.0;
    int m = static_cast<int>(a);
    double sec = (a - m) * 60.0;
    if (sec >= 59.95) { sec = 0.0; ++m; }
    if (m >= 60) { m = 0; ++d; }
    if (sign) return format("%c%02d\xC2\xB0%02d'%04.1f\"", s, d, m, sec);
    return format("%03d\xC2\xB0%02d'%04.1f\"", d, m, sec);
}

std::string formatHoursHMS(double hours) {
    double h = std::fmod(hours, 24.0);
    if (h < 0) h += 24.0;
    int hh = static_cast<int>(h);
    double rem = (h - hh) * 60.0;
    int mm = static_cast<int>(rem);
    double ss = (rem - mm) * 60.0;
    if (ss >= 59.995) { ss = 0.0; ++mm; }
    if (mm >= 60) { mm = 0; hh = (hh + 1) % 24; }
    return format("%02dh%02dm%05.2fs", hh, mm, ss);
}

std::string formatEngineering(double v, int sig) {
    if (v == 0.0) return "0";
    const double a = std::fabs(v);
    if (a >= 1e-3 && a < 1e6) return format("%.*g", sig, v);
    return format("%.*e", sig - 1, v);
}

std::size_t editDistance(const std::string& a, const std::string& b) {
    std::vector<std::size_t> prev(b.size() + 1), cur(b.size() + 1);
    for (std::size_t j = 0; j <= b.size(); ++j) prev[j] = j;
    for (std::size_t i = 1; i <= a.size(); ++i) {
        cur[0] = i;
        for (std::size_t j = 1; j <= b.size(); ++j) {
            const std::size_t cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            cur[j] = std::min({prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + cost});
        }
        std::swap(prev, cur);
    }
    return prev[b.size()];
}

} // namespace solsim::str
