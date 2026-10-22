// =============================================================================
//  SolarSim - JulianDate.cpp
// =============================================================================
#include "JulianDate.hpp"

#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"
#include "../core/StringUtil.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>

namespace solsim::astro {

using namespace solsim::constants;

double calendarToJD(int year, int month, double day) {
    int Y = year;
    int M = month;
    if (M <= 2) {
        Y -= 1;
        M += 12;
    }
    // Gregorian reform: 1582 October 15
    const bool gregorian = (year > 1582) || (year == 1582 && (month > 10 || (month == 10 && day >= 15.0)));
    int B = 0;
    if (gregorian) {
        const int A = static_cast<int>(std::floor(Y / 100.0));
        B = 2 - A + static_cast<int>(std::floor(A / 4.0));
    }
    return std::floor(365.25 * (Y + 4716)) + std::floor(30.6001 * (M + 1)) + day + B - 1524.5;
}

double calendarToJD(const CalendarDate& d) {
    return calendarToJD(d.year, d.month, d.day + d.dayFraction());
}

CalendarDate jdToCalendar(double jd) {
    const double jdp = jd + 0.5;
    const double Z = std::floor(jdp);
    const double F = jdp - Z;
    double A = Z;
    if (Z >= 2299161.0) {
        const double alpha = std::floor((Z - 1867216.25) / 36524.25);
        A = Z + 1.0 + alpha - std::floor(alpha / 4.0);
    }
    const double B = A + 1524.0;
    const double C = std::floor((B - 122.1) / 365.25);
    const double D = std::floor(365.25 * C);
    const double E = std::floor((B - D) / 30.6001);
    const double dayWithFrac = B - D - std::floor(30.6001 * E) + F;

    CalendarDate out;
    out.month = static_cast<int>(E < 14 ? E - 1 : E - 13);
    out.year = static_cast<int>(out.month > 2 ? C - 4716 : C - 4715);
    out.day = static_cast<int>(std::floor(dayWithFrac));
    double frac = dayWithFrac - out.day;
    // Round to the nearest millisecond to avoid 59.9999 s artifacts.
    double secondsOfDay = std::round(frac * DAY_S * 1000.0) / 1000.0;
    if (secondsOfDay >= DAY_S) {
        // Carry into the next day by recursion on the rounded value.
        return jdToCalendar(std::floor(jd + 0.5) + 0.5);
    }
    out.hour = static_cast<int>(secondsOfDay / 3600.0);
    secondsOfDay -= out.hour * 3600.0;
    out.minute = static_cast<int>(secondsOfDay / 60.0);
    out.second = secondsOfDay - out.minute * 60.0;
    return out;
}

bool isLeapYear(int y) {
    if (y > 1582) return (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0);
    return y % 4 == 0;
}

int dayOfYear(const CalendarDate& d) {
    const int k = isLeapYear(d.year) ? 1 : 2;
    return static_cast<int>(275 * d.month / 9) - k * static_cast<int>((d.month + 9) / 12) + d.day - 30;
}

int dayOfWeek(double jd) {
    return static_cast<int>(std::fmod(std::floor(jd + 1.5), 7.0));
}

const char* dayOfWeekName(int dow) {
    static const char* names[] = {"Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"};
    return (dow >= 0 && dow < 7) ? names[dow] : "?";
}

const char* monthName(int month) {
    static const char* names[] = {"January", "February", "March",     "April",   "May",      "June",
                                  "July",    "August",   "September", "October", "November", "December"};
    return (month >= 1 && month <= 12) ? names[month - 1] : "?";
}

double decimalYearFromJD(double jd) { return 2000.0 + (jd - J2000_JD) / JULIAN_YEAR_D; }

// Polynomial expressions for Delta T by Fred Espenak and Jean Meeus
// (NASA eclipse web site). Accuracy of a few seconds for 1800-2050.
double deltaTSeconds(double y) {
    double t;
    if (y < 1800.0 || y >= 2150.0) {
        const double u = (y - 1820.0) / 100.0;
        return -20.0 + 32.0 * u * u;
    }
    if (y < 1860.0) {
        t = y - 1800.0;
        return 13.72 - 0.332447 * t + 0.0068612 * t * t + 0.0041116 * t * t * t - 0.00037436 * std::pow(t, 4) +
               0.0000121272 * std::pow(t, 5) - 0.0000001699 * std::pow(t, 6) + 0.000000000875 * std::pow(t, 7);
    }
    if (y < 1900.0) {
        t = y - 1860.0;
        return 7.62 + 0.5737 * t - 0.251754 * t * t + 0.01680668 * t * t * t - 0.0004473624 * std::pow(t, 4) +
               std::pow(t, 5) / 233174.0;
    }
    if (y < 1920.0) {
        t = y - 1900.0;
        return -2.79 + 1.494119 * t - 0.0598939 * t * t + 0.0061966 * t * t * t - 0.000197 * std::pow(t, 4);
    }
    if (y < 1941.0) {
        t = y - 1920.0;
        return 21.20 + 0.84493 * t - 0.076100 * t * t + 0.0020936 * t * t * t;
    }
    if (y < 1961.0) {
        t = y - 1950.0;
        return 29.07 + 0.407 * t - t * t / 233.0 + t * t * t / 2547.0;
    }
    if (y < 1986.0) {
        t = y - 1975.0;
        return 45.45 + 1.067 * t - t * t / 260.0 - t * t * t / 718.0;
    }
    if (y < 2005.0) {
        t = y - 2000.0;
        return 63.86 + 0.3345 * t - 0.060374 * t * t + 0.0017275 * t * t * t + 0.000651814 * std::pow(t, 4) +
               0.00002373599 * std::pow(t, 5);
    }
    if (y < 2050.0) {
        // Espenak-Meeus' 2005-2050 polynomial overestimates the observed
        // value (69.2 s in 2024-2026) by several seconds; we blend it with the
        // observed IERS trend which is nearly flat since 2016.
        t = y - 2000.0;
        const double poly = 62.92 + 0.32217 * t + 0.005589 * t * t;
        const double observed = 69.2 + 0.05 * (y - 2025.0);
        // Weight 0 until 2030, rising to 1 in 2050 (continuous with the next branch).
        const double w = math::clamp((y - 2030.0) / 20.0, 0.0, 1.0);
        return (y < 2015.0) ? poly : (1.0 - w) * observed + w * poly;
    }
    const double u = (y - 1820.0) / 100.0;
    return -20.0 + 32.0 * u * u - 0.5628 * (2150.0 - y);
}

double utcToTT(double jdUtc) {
    return jdUtc + deltaTSeconds(decimalYearFromJD(jdUtc)) / DAY_S;
}

double ttToUTC(double jdTT) {
    double ut = jdTT;
    for (int i = 0; i < 3; ++i) ut = jdTT - deltaTSeconds(decimalYearFromJD(ut)) / DAY_S;
    return ut;
}

double julianCenturies(double jdTT) { return (jdTT - J2000_JD) / JULIAN_CENTURY_D; }

double gmstDegrees(double jdUt) {
    const double T = (jdUt - J2000_JD) / JULIAN_CENTURY_D;
    const double gmst = 280.46061837 + 360.98564736629 * (jdUt - J2000_JD) + 0.000387933 * T * T -
                        T * T * T / 38710000.0;
    return math::wrap360(gmst);
}

double nowJDUtc() {
    using namespace std::chrono;
    const auto now = system_clock::now();
    const double unixSeconds = duration<double>(now.time_since_epoch()).count();
    return 2440587.5 + unixSeconds / DAY_S;
}

namespace {

bool parseIntRange(const std::string& s, int lo, int hi, int& out) {
    long long v = 0;
    if (!str::parseInt(s, v) || v < lo || v > hi) return false;
    out = static_cast<int>(v);
    return true;
}

int daysInMonth(int year, int month) {
    static const int dm[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (month == 2 && isLeapYear(year)) return 29;
    return dm[month - 1];
}

} // namespace

bool parseDate(const std::string& textIn, double& jdTT, std::string& error) {
    const std::string text = str::trim(textIn);
    const std::string lower = str::toLower(text);
    if (lower.empty()) {
        error = "empty date";
        return false;
    }
    if (lower == "now" || lower == "today") {
        jdTT = utcToTT(nowJDUtc());
        return true;
    }
    if (lower == "j2000" || lower == "j2000.0") {
        jdTT = J2000_JD;
        return true;
    }
    double v = 0.0;
    if (str::startsWith(lower, "jd")) {
        if (!str::parseDouble(lower.substr(2), v)) {
            error = "invalid Julian Date '" + text + "'";
            return false;
        }
        jdTT = v;
        return true;
    }
    if (str::startsWith(lower, "mjd")) {
        if (!str::parseDouble(lower.substr(3), v)) {
            error = "invalid Modified Julian Date '" + text + "'";
            return false;
        }
        jdTT = v + MJD_OFFSET;
        return true;
    }
    // Plain number -> JD if large, decimal year otherwise (e.g. 2026.5).
    if (str::parseDouble(lower, v)) {
        if (v > 100000.0) jdTT = v;
        else jdTT = utcToTT(J2000_JD + (v - 2000.0) * JULIAN_YEAR_D);
        return true;
    }

    // ISO-like: YYYY-MM-DD[(T| )HH:MM[:SS[.sss]]][Z]
    std::string datePart = text, timePart;
    const auto sep = text.find_first_of("Tt ");
    if (sep != std::string::npos) {
        datePart = text.substr(0, sep);
        timePart = str::trim(text.substr(sep + 1));
    }
    if (!timePart.empty() && (timePart.back() == 'Z' || timePart.back() == 'z')) timePart.pop_back();
    if (str::endsWith(str::toLower(timePart), "utc")) timePart = str::trim(timePart.substr(0, timePart.size() - 3));

    // Allow a leading minus for BCE years: "-0500-03-01"
    bool negativeYear = false;
    if (!datePart.empty() && datePart[0] == '-') {
        negativeYear = true;
        datePart = datePart.substr(1);
    }
    const auto parts = str::split(datePart, '-', false);
    if (parts.size() != 3) {
        error = "unrecognized date '" + text + "' (expected YYYY-MM-DD[THH:MM[:SS]], JD<number>, J2000 or now)";
        return false;
    }
    CalendarDate cd;
    if (!parseIntRange(parts[0], 0, 99999, cd.year) || !parseIntRange(parts[1], 1, 12, cd.month)) {
        error = "invalid year/month in '" + text + "'";
        return false;
    }
    if (negativeYear) cd.year = -cd.year;
    if (!parseIntRange(parts[2], 1, daysInMonth(cd.year, cd.month), cd.day)) {
        error = "invalid day of month in '" + text + "'";
        return false;
    }
    if (!timePart.empty()) {
        const auto t = str::split(timePart, ':', false);
        if (t.size() < 2 || t.size() > 3 || !parseIntRange(t[0], 0, 23, cd.hour) ||
            !parseIntRange(t[1], 0, 59, cd.minute)) {
            error = "invalid time of day in '" + text + "'";
            return false;
        }
        if (t.size() == 3) {
            if (!str::parseDouble(t[2], cd.second) || cd.second < 0.0 || cd.second >= 61.0) {
                error = "invalid seconds in '" + text + "'";
                return false;
            }
        }
    }
    jdTT = utcToTT(calendarToJD(cd));
    return true;
}

std::string formatUTC(double jdTT, bool withSeconds) {
    const CalendarDate c = jdToCalendar(ttToUTC(jdTT));
    if (withSeconds)
        return str::format("%04d-%02d-%02d %02d:%02d:%02d UTC", c.year, c.month, c.day, c.hour, c.minute,
                           static_cast<int>(c.second));
    return str::format("%04d-%02d-%02d %02d:%02d UTC", c.year, c.month, c.day, c.hour, c.minute);
}

std::string formatDateOnly(double jdTT) {
    const CalendarDate c = jdToCalendar(ttToUTC(jdTT));
    return str::format("%04d-%02d-%02d", c.year, c.month, c.day);
}

std::string formatISO(double jdTT) {
    const CalendarDate c = jdToCalendar(ttToUTC(jdTT));
    return str::format("%04d-%02d-%02dT%02d:%02d:%02dZ", c.year, c.month, c.day, c.hour, c.minute,
                       static_cast<int>(c.second));
}

} // namespace solsim::astro
