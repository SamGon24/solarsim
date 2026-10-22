// =============================================================================
//  SolarSim - JulianDate.hpp
//  Calendar <-> Julian Date conversion (Meeus, "Astronomical Algorithms",
//  ch. 7), time-scale handling (UTC <-> TT via an approximation of Delta T)
//  and parsing/formatting of user supplied dates.
//
//  Convention: the simulation clock runs in TDB (~TT, differences < 2 ms).
//  User facing dates are UTC. The difference Delta T = TT - UT is ~69 s today.
// =============================================================================
#pragma once

#include <string>

namespace solsim::astro {

struct CalendarDate {
    int year = 2000;
    int month = 1;   // 1..12
    int day = 1;     // 1..31
    int hour = 0;
    int minute = 0;
    double second = 0.0;

    double dayFraction() const { return (hour + (minute + second / 60.0) / 60.0) / 24.0; }
};

// Gregorian calendar from 1582-10-15 onward, Julian calendar before.
double calendarToJD(const CalendarDate& date);
double calendarToJD(int year, int month, double day);
CalendarDate jdToCalendar(double jd);

// Delta T = TT - UT1 in seconds (Espenak & Meeus polynomial approximation).
double deltaTSeconds(double decimalYear);
double decimalYearFromJD(double jd);

double utcToTT(double jdUtc);   // returns JD(TT)
double ttToUTC(double jdTT);    // returns JD(UTC) (iterative inverse)

// Julian centuries since J2000.0
double julianCenturies(double jdTT);

int dayOfWeek(double jd);                     // 0 = Sunday
const char* dayOfWeekName(int dow);
const char* monthName(int month);
bool isLeapYear(int year);
int dayOfYear(const CalendarDate& d);

// Greenwich Mean Sidereal Time in degrees (IAU 1982 expression) for JD(UT1).
double gmstDegrees(double jdUt);

// Parses: "2026-10-21", "2026-10-21T13:45", "2026-10-21 13:45:30.5", "J2000",
// "JD2451545.0", "2451545.0", "MJD60000", "now". The result is JD(UTC)
// unless the string is an explicit JD/J2000 (taken as TT).
// Returns false and sets `error` on failure.
bool parseDate(const std::string& text, double& jdTT, std::string& error);

// "2026-10-21 13:45:30 UTC"  (input JD in TT; converted to UTC for display)
std::string formatUTC(double jdTT, bool withSeconds = true);
std::string formatDateOnly(double jdTT);
std::string formatISO(double jdTT); // 2026-10-21T13:45:30Z

// Current system time as JD(UTC).
double nowJDUtc();

} // namespace solsim::astro
