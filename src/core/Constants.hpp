// =============================================================================
//  SolarSim - Constants.hpp
//  Physical, astronomical and mathematical constants used across the project.
//
//  Internal unit system (unless stated otherwise):
//      length  : astronomical unit (AU)
//      time    : day (86400 SI seconds, TDB time scale)
//      mass    : expressed through the gravitational parameter GM [AU^3/day^2]
//
//  Using GM directly (instead of G and M separately) is standard practice in
//  celestial mechanics: GM is known to ~10 significant digits from spacecraft
//  tracking while G alone is only known to ~5.
// =============================================================================
#pragma once

namespace solsim::constants {

// ---- Mathematics ------------------------------------------------------------
inline constexpr double PI         = 3.14159265358979323846264338327950288;
inline constexpr double TWO_PI     = 2.0 * PI;
inline constexpr double HALF_PI    = 0.5 * PI;
inline constexpr double DEG2RAD    = PI / 180.0;
inline constexpr double RAD2DEG    = 180.0 / PI;
inline constexpr double ARCSEC2RAD = DEG2RAD / 3600.0;
inline constexpr double RAD2ARCSEC = 1.0 / ARCSEC2RAD;

// ---- Length -------------------------------------------------------------------
inline constexpr double AU_KM = 149597870.700;        // IAU 2012 Resolution B2 (exact)
inline constexpr double AU_M  = AU_KM * 1000.0;
inline constexpr double KM_PER_AU = AU_KM;
inline constexpr double AU_PER_KM = 1.0 / AU_KM;
inline constexpr double LIGHT_YEAR_AU = 63241.077084266;
inline constexpr double PARSEC_AU = 206264.80624709636;

// ---- Time ---------------------------------------------------------------------
inline constexpr double DAY_S            = 86400.0;
inline constexpr double HOUR_D           = 1.0 / 24.0;
inline constexpr double MINUTE_D         = 1.0 / 1440.0;
inline constexpr double JULIAN_YEAR_D    = 365.25;
inline constexpr double JULIAN_CENTURY_D = 36525.0;
inline constexpr double J2000_JD         = 2451545.0;        // 2000-01-01 12:00 TT
inline constexpr double MJD_OFFSET       = 2400000.5;

// ---- Speed of light -----------------------------------------------------------
inline constexpr double C_KM_S      = 299792.458;                       // exact
inline constexpr double C_AU_PER_DAY = C_KM_S * DAY_S / AU_KM;          // ~173.1446 AU/day
inline constexpr double LIGHT_TIME_AU_S = AU_KM / C_KM_S;               // ~499.005 s

// ---- Gravitation ------------------------------------------------------------
inline constexpr double G_SI = 6.67430e-11;   // CODATA 2018 [m^3 kg^-1 s^-2]
// Gaussian gravitational constant (historical definition of the AU).
inline constexpr double GAUSS_K = 0.01720209895;

// Conversion of GM from km^3/s^2 to AU^3/day^2.
inline constexpr double KM3S2_TO_AU3D2 = (DAY_S * DAY_S) / (AU_KM * AU_KM * AU_KM);
inline constexpr double AU3D2_TO_KM3S2 = 1.0 / KM3S2_TO_AU3D2;

// Heliocentric gravitational constant (DE440), km^3/s^2.
inline constexpr double GM_SUN_KM3S2 = 132712440041.279419;
inline constexpr double GM_SUN = GM_SUN_KM3S2 * KM3S2_TO_AU3D2;    // ~2.959122e-4

// Nominal solar mass (kg) derived from GM_sun / G.
inline constexpr double SOLAR_MASS_KG = GM_SUN_KM3S2 * 1.0e9 / G_SI;

// ---- Earth / Sun reference values -----------------------------------------
inline constexpr double EARTH_RADIUS_KM = 6371.0;          // mean radius
inline constexpr double EARTH_EQ_RADIUS_KM = 6378.137;     // WGS-84
inline constexpr double SUN_RADIUS_KM = 695700.0;          // IAU 2015 nominal
inline constexpr double SOLAR_LUMINOSITY_W = 3.828e26;     // IAU 2015 nominal
inline constexpr double SUN_ABS_MAG_V = 4.83;
inline constexpr double SUN_APP_MAG_V = -26.74;

// Mean obliquity of the ecliptic at J2000 (IAU 2006): 84381.406 arcsec.
inline constexpr double OBLIQUITY_J2000_DEG = 84381.406 / 3600.0;
inline constexpr double OBLIQUITY_J2000 = OBLIQUITY_J2000_DEG * DEG2RAD;

// Annual aberration constant (arcsec).
inline constexpr double ABERRATION_CONSTANT_ARCSEC = 20.49552;

// Stefan-Boltzmann constant (W m^-2 K^-4) - used for equilibrium temperatures.
inline constexpr double STEFAN_BOLTZMANN = 5.670374419e-8;

} // namespace solsim::constants
