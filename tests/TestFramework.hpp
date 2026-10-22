// =============================================================================
//  SolarSim - TestFramework.hpp
//  A dependency-free unit test harness:
//      TEST(suite, name) { CHECK(cond); CHECK_NEAR(a, b, tol); CHECK_EQ(a, b); }
//  Tests self-register through static initialisation; run them all with
//  `solarsim_tests` or filter with `solarsim_tests <substring>`.
// =============================================================================
#pragma once

#include <cmath>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace testing {

struct TestCase {
    std::string suite;
    std::string name;
    std::function<void()> fn;
};

std::vector<TestCase>& registry();
void reportFailure(const char* file, int line, const std::string& message);

struct Registrar {
    Registrar(const char* suite, const char* name, std::function<void()> fn) {
        registry().push_back({suite, name, std::move(fn)});
    }
};

} // namespace testing

#define TEST_CONCAT_INNER(a, b) a##b
#define TEST_CONCAT(a, b) TEST_CONCAT_INNER(a, b)

#define TEST(suite, name)                                                                     \
    static void TEST_CONCAT(test_##suite##_, name)();                                         \
    static ::testing::Registrar TEST_CONCAT(registrar_##suite##_, name)(#suite, #name,         \
                                                                       &TEST_CONCAT(test_##suite##_, name)); \
    static void TEST_CONCAT(test_##suite##_, name)()

#define CHECK(cond)                                                                    \
    do {                                                                               \
        if (!(cond)) ::testing::reportFailure(__FILE__, __LINE__, "CHECK(" #cond ")"); \
    } while (0)

#define CHECK_EQ(a, b)                                                                              \
    do {                                                                                            \
        const auto va_ = (a);                                                                       \
        const auto vb_ = (b);                                                                       \
        if (!(va_ == vb_)) {                                                                        \
            std::ostringstream os_;                                                                 \
            os_ << "CHECK_EQ(" #a ", " #b ") values: " << va_ << " vs " << vb_;                     \
            ::testing::reportFailure(__FILE__, __LINE__, os_.str());                                \
        }                                                                                           \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                                        \
    do {                                                                                             \
        const double va_ = static_cast<double>(a);                                                   \
        const double vb_ = static_cast<double>(b);                                                   \
        if (!(std::fabs(va_ - vb_) <= (tol))) {                                                      \
            std::ostringstream os_;                                                                  \
            os_.precision(15);                                                                       \
            os_ << "CHECK_NEAR(" #a ", " #b ", " #tol ") values: " << va_ << " vs " << vb_            \
                << " (diff " << std::fabs(va_ - vb_) << ")";                                         \
            ::testing::reportFailure(__FILE__, __LINE__, os_.str());                                 \
        }                                                                                            \
    } while (0)
