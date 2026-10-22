// =============================================================================
//  SolarSim - LunarTheory.cpp
//  Coefficients: Meeus, Astronomical Algorithms, Tables 47.A and 47.B
//  (units: longitude/latitude 1e-6 degree, distance 1e-3 km).
// =============================================================================
#include "LunarTheory.hpp"

#include "Coordinates.hpp"
#include "JulianDate.hpp"

#include "../core/Constants.hpp"
#include "../core/MathUtil.hpp"

#include <cmath>
#include <cstdlib>

namespace solsim::astro {

using namespace solsim::constants;

namespace {

struct LRTerm {
    signed char D, M, Mp, F;
    int sl;  // longitude
    int sr;  // distance
};

struct BTerm {
    signed char D, M, Mp, F;
    int sb;  // latitude
};

const LRTerm kLR[] = {
    {0, 0, 1, 0, 6288774, -20905355}, {2, 0, -1, 0, 1274027, -3699111}, {2, 0, 0, 0, 658314, -2955968},
    {0, 0, 2, 0, 213618, -569925},    {0, 1, 0, 0, -185116, 48888},     {0, 0, 0, 2, -114332, -3149},
    {2, 0, -2, 0, 58793, 246158},     {2, -1, -1, 0, 57066, -152138},   {2, 0, 1, 0, 53322, -170733},
    {2, -1, 0, 0, 45758, -204586},    {0, 1, -1, 0, -40923, -129620},   {1, 0, 0, 0, -34720, 108743},
    {0, 1, 1, 0, -30383, 104755},     {2, 0, 0, -2, 15327, 10321},      {0, 0, 1, 2, -12528, 0},
    {0, 0, 1, -2, 10980, 79661},      {4, 0, -1, 0, 10675, -34782},     {0, 0, 3, 0, 10034, -23210},
    {4, 0, -2, 0, 8548, -21636},      {2, 1, -1, 0, -7888, 24208},      {2, 1, 0, 0, -6766, 30824},
    {1, 0, -1, 0, -5163, -8379},      {1, 1, 0, 0, 4987, -16675},       {2, -1, 1, 0, 4036, -12831},
    {2, 0, 2, 0, 3994, -10445},       {4, 0, 0, 0, 3861, -11650},       {2, 0, -3, 0, 3665, 14403},
    {0, 1, -2, 0, -2689, -7003},      {2, 0, -1, 2, -2602, 0},          {2, -1, -2, 0, 2390, 10056},
    {1, 0, 1, 0, -2348, 6322},        {2, -2, 0, 0, 2236, -9884},       {0, 1, 2, 0, -2120, 5751},
    {0, 2, 0, 0, -2069, 0},           {2, -2, -1, 0, 2048, -4950},      {2, 0, 1, -2, -1773, 4130},
    {2, 0, 0, 2, -1595, 0},           {4, -1, -1, 0, 1215, -3958},      {0, 0, 2, 2, -1110, 0},
    {3, 0, -1, 0, -892, 3258},        {2, 1, 1, 0, -810, 2616},         {4, -1, -2, 0, 759, -1897},
    {0, 2, -1, 0, -713, -2117},       {2, 2, -1, 0, -700, 2354},        {2, 1, -2, 0, 691, 0},
    {2, -1, 0, -2, 596, 0},           {4, 0, 1, 0, 549, -1423},         {0, 0, 4, 0, 537, -1117},
    {4, -1, 0, 0, 520, -1571},        {1, 0, -2, 0, -487, -1739},       {2, 1, 0, -2, -399, 0},
    {0, 0, 2, -2, -381, -4421},       {1, 1, 1, 0, 351, 0},             {3, 0, -2, 0, -340, 0},
    {4, 0, -3, 0, 330, 0},            {2, -1, 2, 0, 327, 0},            {0, 2, 1, 0, -323, 1165},
    {1, 1, -1, 0, 299, 0},            {2, 0, 3, 0, 294, 0},             {2, 0, -1, -2, 0, 8752},
};

const BTerm kB[] = {
    {0, 0, 0, 1, 5128122},  {0, 0, 1, 1, 280602},  {0, 0, 1, -1, 277693}, {2, 0, 0, -1, 173237},
    {2, 0, -1, 1, 55413},   {2, 0, -1, -1, 46271}, {2, 0, 0, 1, 32573},   {0, 0, 2, 1, 17198},
    {2, 0, 1, -1, 9266},    {0, 0, 2, -1, 8822},   {2, -1, 0, -1, 8216},  {2, 0, -2, -1, 4324},
    {2, 0, 1, 1, 4200},     {2, 1, 0, -1, -3359},  {2, -1, -1, 1, 2463},  {2, -1, 0, 1, 2211},
    {2, -1, -1, -1, 2065},  {0, 1, -1, -1, -1870}, {4, 0, -1, -1, 1828},  {0, 1, 0, 1, -1794},
    {0, 0, 0, 3, -1749},    {0, 1, -1, 1, -1565},  {1, 0, 0, 1, -1491},   {0, 1, 1, 1, -1475},
    {0, 1, 1, -1, -1410},   {0, 1, 0, -1, -1344},  {1, 0, 0, -1, -1335},  {0, 0, 3, 1, 1107},
    {4, 0, 0, -1, 1021},    {4, 0, -1, 1, 833},    {0, 0, 1, -3, 777},    {4, 0, -2, 1, 671},
    {2, 0, 0, -3, 607},     {2, 0, 2, -1, 596},    {2, -1, 1, -1, 491},   {2, 0, -2, 1, -451},
    {0, 0, 3, -1, 439},     {2, 0, 2, 1, 422},     {2, 0, -3, -1, 421},   {2, 1, -1, 1, -366},
    {2, 1, 0, 1, -351},     {4, 0, 0, 1, 331},     {2, -1, 1, 1, 315},    {2, -2, 0, -1, 302},
    {0, 0, 1, 3, -283},     {2, 1, 1, -1, -229},   {1, 1, 0, -1, 223},    {1, 1, 0, 1, 223},
    {0, 1, -2, -1, -220},   {2, 1, -1, -1, -220},  {1, 0, 1, 1, -185},    {2, -1, -2, -1, 181},
    {0, 1, 2, 1, -177},     {4, 0, -2, -1, 176},   {4, -1, -1, -1, 166},  {1, 0, 1, -1, -164},
    {4, 0, 1, -1, 132},     {1, 0, -1, -1, -119},  {4, -1, 0, -1, 115},   {2, -2, 0, 1, 107},
};

} // namespace

LunarPosition lunarPositionOfDate(double jdTT) {
    const double T = julianCenturies(jdTT);
    const double T2 = T * T, T3 = T2 * T, T4 = T3 * T;
    // Fundamental arguments (degrees), Meeus 47.1 - 47.5
    const double Lp = 218.3164477 + 481267.88123421 * T - 0.0015786 * T2 + T3 / 538841.0 - T4 / 65194000.0;
    const double D = 297.8501921 + 445267.1114034 * T - 0.0018819 * T2 + T3 / 545868.0 - T4 / 113065000.0;
    const double M = 357.5291092 + 35999.0502909 * T - 0.0001536 * T2 + T3 / 24490000.0;
    const double Mp = 134.9633964 + 477198.8675055 * T + 0.0087414 * T2 + T3 / 69699.0 - T4 / 14712000.0;
    const double F = 93.2720950 + 483202.0175233 * T - 0.0036539 * T2 - T3 / 3526000.0 + T4 / 863310000.0;
    const double A1 = 119.75 + 131.849 * T;
    const double A2 = 53.09 + 479264.290 * T;
    const double A3 = 313.45 + 481266.484 * T;
    // Eccentricity of Earth's orbit correction for terms containing M.
    const double E = 1.0 - 0.002516 * T - 0.0000074 * T2;

    const double d = math::wrap360(D) * DEG2RAD, m = math::wrap360(M) * DEG2RAD, mp = math::wrap360(Mp) * DEG2RAD,
                 f = math::wrap360(F) * DEG2RAD;
    auto eFactor = [E](int mCoeff) { return mCoeff == 0 ? 1.0 : (std::abs(mCoeff) == 1 ? E : E * E); };

    double sl = 0.0, sr = 0.0, sb = 0.0;
    for (const LRTerm& t : kLR) {
        const double arg = t.D * d + t.M * m + t.Mp * mp + t.F * f;
        const double ef = eFactor(t.M);
        sl += t.sl * ef * std::sin(arg);
        sr += t.sr * ef * std::cos(arg);
    }
    for (const BTerm& t : kB) {
        const double arg = t.D * d + t.M * m + t.Mp * mp + t.F * f;
        sb += t.sb * eFactor(t.M) * std::sin(arg);
    }
    // Additive terms: action of Venus (A1), Jupiter (A2) and Earth's flattening (L').
    const double lp = math::wrap360(Lp) * DEG2RAD;
    sl += 3958.0 * std::sin(A1 * DEG2RAD) + 1962.0 * std::sin(lp - f) + 318.0 * std::sin(A2 * DEG2RAD);
    sb += -2235.0 * std::sin(lp) + 382.0 * std::sin(A3 * DEG2RAD) + 175.0 * std::sin(A1 * DEG2RAD - f) +
          175.0 * std::sin(A1 * DEG2RAD + f) + 127.0 * std::sin(lp - mp) - 115.0 * std::sin(lp + mp);

    LunarPosition pos;
    pos.longitudeDeg = math::wrap360(Lp + sl * 1e-6);
    pos.latitudeDeg = sb * 1e-6;
    pos.distanceKm = 385000.56 + sr * 1e-3;
    return pos;
}

Vec3 lunarPositionJ2000(double jdTT) {
    const LunarPosition p = lunarPositionOfDate(jdTT);
    // Mean equinox of date -> J2000: remove the accumulated general precession
    // in longitude (the slow motion of the ecliptic itself is negligible here).
    const double lon = p.longitudeDeg * DEG2RAD - generalPrecessionLongitude(jdTT);
    return fromSpherical(lon, p.latitudeDeg * DEG2RAD, p.distanceKm / AU_KM);
}

void lunarStateJ2000(double jdTT, Vec3& r, Vec3& v) {
    const double h = 0.01; // days
    r = lunarPositionJ2000(jdTT);
    const Vec3 rp = lunarPositionJ2000(jdTT + h), rm = lunarPositionJ2000(jdTT - h);
    const Vec3 rp2 = lunarPositionJ2000(jdTT + 2 * h), rm2 = lunarPositionJ2000(jdTT - 2 * h);
    // 4th-order central difference.
    v = (rm2 - rm * 8.0 + rp * 8.0 - rp2) / (12.0 * h);
}

} // namespace solsim::astro
