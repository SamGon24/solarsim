// =============================================================================
//  SolarSim - test_main.cpp : test runner
// =============================================================================
#include "TestFramework.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>

namespace testing {

namespace {
int g_failures = 0;
bool g_currentFailed = false;
} // namespace

std::vector<TestCase>& registry() {
    static std::vector<TestCase> tests;
    return tests;
}

void reportFailure(const char* file, int line, const std::string& message) {
    ++g_failures;
    g_currentFailed = true;
    std::printf("    FAILED %s:%d: %s\n", file, line, message.c_str());
}

} // namespace testing

int main(int argc, char** argv) {
    using namespace testing;
    const char* filter = argc > 1 ? argv[1] : nullptr;
    int run = 0, failedTests = 0;
    const auto start = std::chrono::steady_clock::now();
    for (const auto& t : registry()) {
        const std::string full = t.suite + "." + t.name;
        if (filter && !std::strstr(full.c_str(), filter)) continue;
        g_currentFailed = false;
        const auto t0 = std::chrono::steady_clock::now();
        try {
            t.fn();
        } catch (const std::exception& e) {
            reportFailure(__FILE__, __LINE__, std::string("unexpected exception: ") + e.what());
        }
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        std::printf("[%s] %-48s %8.1f ms\n", g_currentFailed ? "FAIL" : " OK ", full.c_str(), ms);
        ++run;
        failedTests += g_currentFailed;
    }
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::printf("\n%d tests, %d failed, %d failed checks (%.2f s)\n", run, failedTests, g_failures, s);
    return failedTests == 0 ? 0 : 1;
}
