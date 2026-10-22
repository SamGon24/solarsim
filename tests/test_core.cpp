// =============================================================================
//  SolarSim - test_core.cpp : math, strings, parsers, time scales
// =============================================================================
#include "TestFramework.hpp"

#include "astro/JulianDate.hpp"
#include "core/ArgParser.hpp"
#include "core/Constants.hpp"
#include "core/IniParser.hpp"
#include "core/Mat3.hpp"
#include "core/MathUtil.hpp"
#include "core/Random.hpp"
#include "core/StringUtil.hpp"
#include "core/ThreadPool.hpp"

#include <atomic>

using namespace solsim;

TEST(Math, VectorAlgebra) {
    const Vec3 a{1, 2, 3}, b{-2, 0.5, 4};
    CHECK_NEAR(a.dot(b), 11.0, 1e-15);
    const Vec3 c = a.cross(b);
    CHECK_NEAR(c.dot(a), 0.0, 1e-12);
    CHECK_NEAR(c.dot(b), 0.0, 1e-12);
    CHECK_NEAR(Vec3(3, 4, 0).norm(), 5.0, 1e-15);
    CHECK_NEAR(angleBetween(Vec3::unitX(), Vec3::unitY()), constants::HALF_PI, 1e-15);
    CHECK_NEAR(angleBetween(Vec3{1, 1e-9, 0}, Vec3::unitX()), 1e-9, 1e-20);
}

TEST(Math, MatrixInverseAndRotation) {
    const Mat3 r = Mat3::rotZ(0.3) * Mat3::rotX(-1.1) * Mat3::rotY(2.0);
    Mat3 inv;
    CHECK(r.inverse(inv));
    const Mat3 id = r * inv;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) CHECK_NEAR(id.m[i][j], i == j ? 1.0 : 0.0, 1e-14);
    CHECK_NEAR(r.determinant(), 1.0, 1e-14);
    const Mat3 aa = Mat3::axisAngle(Vec3::unitZ(), 0.7);
    const Vec3 v = aa * Vec3::unitX();
    CHECK_NEAR(v.x, std::cos(0.7), 1e-15);
    CHECK_NEAR(v.y, std::sin(0.7), 1e-15);
}

TEST(Math, AngleWrapping) {
    CHECK_NEAR(math::wrap360(-30.0), 330.0, 1e-12);
    CHECK_NEAR(math::wrap180(190.0), -170.0, 1e-12);
    CHECK_NEAR(math::wrapPi(3.5 * constants::PI), -0.5 * constants::PI, 1e-12);
    std::vector<double> a = {3.0, -3.0, -2.8};
    math::unwrapAngles(a);
    CHECK(a[1] > 3.0);
    CHECK(a[2] > a[1]);
}

TEST(Math, RootFinders) {
    auto f = [](double x) { return std::cos(x) - x; };
    const auto b = math::brent(f, 0.0, 1.0, 1e-14);
    CHECK(b.converged);
    CHECK_NEAR(b.x, 0.7390851332151607, 1e-12);
    const auto bi = math::bisect(f, 0.0, 1.0, 1e-12);
    CHECK_NEAR(bi.x, 0.7390851332151607, 1e-10);
    const auto n = math::newton(f, [](double x) { return -std::sin(x) - 1.0; }, 0.5);
    CHECK_NEAR(n.x, 0.7390851332151607, 1e-13);
    const double m = math::goldenSectionMax([](double x) { return -(x - 1.234) * (x - 1.234); }, 0, 3, 1e-10);
    CHECK_NEAR(m, 1.234, 1e-8);
}

TEST(Math, LeastSquares) {
    std::vector<double> x, y;
    for (int i = 0; i < 50; ++i) {
        x.push_back(i);
        y.push_back(3.0 - 0.5 * i);
    }
    const auto fit = math::fitLine(x, y);
    CHECK_NEAR(fit.slope, -0.5, 1e-12);
    CHECK_NEAR(fit.intercept, 3.0, 1e-10);
    CHECK_NEAR(fit.r2, 1.0, 1e-12);
}

TEST(Math, RandomDeterministicAndUniform) {
    Random a(42), b(42);
    for (int i = 0; i < 100; ++i) CHECK_EQ(a.nextU64(), b.nextU64());
    Random r(7);
    double sum = 0, sum2 = 0;
    const int n = 200000;
    for (int i = 0; i < n; ++i) {
        const double g = r.normal();
        sum += g;
        sum2 += g * g;
    }
    CHECK_NEAR(sum / n, 0.0, 0.01);
    CHECK_NEAR(sum2 / n, 1.0, 0.02);
}

TEST(Core, StringHelpers) {
    CHECK_EQ(str::trim("  abc \t"), std::string("abc"));
    CHECK_EQ(str::split("a,,b,c", ',').size(), static_cast<std::size_t>(3));
    double d = 0;
    CHECK(str::parseDouble("1.5e3", d));
    CHECK_NEAR(d, 1500.0, 0);
    CHECK(!str::parseDouble("1.5x", d));
    CHECK_EQ(str::editDistance("integrator", "integratr"), static_cast<std::size_t>(1));
    CHECK_EQ(str::utf8Length("\xC2\xB0" "C"), static_cast<std::size_t>(2));
    CHECK_EQ(str::formatHoursHMS(12.5), std::string("12h30m00.00s"));
}

TEST(Core, ArgParser) {
    ArgParser p("prog");
    p.addCommand("run", "run it");
    p.addOption("dt", "DAYS", "step", "g");
    p.addFlag("gr", "relativity", "g");
    p.addFlag("verbose", "v", "g", 'v');
    const char* argv[] = {"prog", "run", "--dt=0.5", "--gr", "-v", "extra"};
    p.parse(6, argv);
    CHECK_EQ(p.command(), std::string("run"));
    CHECK_NEAR(p.getDouble("dt", 0), 0.5, 0);
    CHECK(p.flag("gr"));
    CHECK(p.flag("verbose"));
    CHECK_EQ(p.positionals().size(), static_cast<std::size_t>(1));
    bool threw = false;
    try {
        const char* bad[] = {"prog", "--dtt", "1"};
        p.parse(3, bad);
    } catch (const UsageError& e) {
        threw = std::string(e.what()).find("did you mean --dt") != std::string::npos;
    }
    CHECK(threw);
}

TEST(Core, IniParser) {
    IniDocument doc;
    std::string err;
    CHECK(doc.loadString("; comment\n[body]\nname = X ; trailing\ne=0.5\n[body]\nname=Y\n", err));
    CHECK_EQ(doc.sectionsNamed("body").size(), static_cast<std::size_t>(2));
    CHECK_EQ(doc.sections()[0].get("name"), std::string("X"));
    CHECK_NEAR(doc.sections()[0].getDouble("e", 0), 0.5, 0);
    CHECK(!doc.loadString("[body]\nnoequals\n", err));
}

TEST(Core, ThreadPoolCoversRange) {
    ThreadPool pool(4);
    std::vector<int> hits(10000, 0);
    std::atomic<long long> sum{0};
    for (int rep = 0; rep < 20; ++rep) {
        pool.parallelFor(hits.size(), [&](std::size_t b, std::size_t e) {
            long long local = 0;
            for (std::size_t i = b; i < e; ++i) {
                ++hits[i];
                local += static_cast<long long>(i);
            }
            sum += local;
        });
    }
    for (int h : hits) CHECK_EQ(h, 20);
    CHECK_EQ(sum.load(), 20LL * 9999LL * 10000LL / 2LL);
}

TEST(Time, JulianDates) {
    CHECK_NEAR(astro::calendarToJD(2000, 1, 1.5), 2451545.0, 1e-9);
    CHECK_NEAR(astro::calendarToJD(1957, 10, 4.81), 2436116.31, 1e-6);   // Meeus 7.a
    CHECK_NEAR(astro::calendarToJD(333, 1, 27.5), 1842713.0, 1e-9);      // Meeus 7.b (Julian calendar)
    CHECK_NEAR(astro::calendarToJD(1582, 10, 15.0), 2299160.5, 1e-9);    // Gregorian reform
    const auto c = astro::jdToCalendar(2436116.31);
    CHECK_EQ(c.year, 1957);
    CHECK_EQ(c.month, 10);
    CHECK_EQ(c.day, 4);
    CHECK_EQ(c.hour, 19);
    CHECK_EQ(astro::dayOfWeek(2434923.5), 3); // 1954 June 30 was a Wednesday
}

TEST(Time, ParsingAndTimeScales) {
    double jd = 0;
    std::string err;
    CHECK(astro::parseDate("J2000", jd, err));
    CHECK_NEAR(jd, 2451545.0, 0);
    CHECK(astro::parseDate("2024-03-20T03:06", jd, err));
    const double ut = astro::ttToUTC(jd);
    const auto c = astro::jdToCalendar(ut);
    CHECK_EQ(c.hour, 3);
    CHECK_EQ(c.minute, 6);
    CHECK(!astro::parseDate("2024-02-30", jd, err));
    CHECK(!astro::parseDate("yesterday", jd, err));
    CHECK_NEAR(astro::deltaTSeconds(2024.0), 69.2, 1.5);
    CHECK_NEAR(astro::deltaTSeconds(1990.0), 56.9, 1.0);
    CHECK_NEAR(astro::gmstDegrees(2451545.0), 280.46061837, 1e-6);
}
