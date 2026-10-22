// =============================================================================
//  SolarSim - Forces.cpp
// =============================================================================
#include "Forces.hpp"

#include "../core/Constants.hpp"
#include "../core/StringUtil.hpp"
#include "../core/ThreadPool.hpp"

#include <algorithm>
#include <cmath>

namespace solsim::physics {

// -----------------------------------------------------------------------------
// NewtonianGravity
// -----------------------------------------------------------------------------
NewtonianGravity::NewtonianGravity(double softeningAU, ThreadPool* pool)
    : m_soft2(softeningAU * softeningAU), m_pool(pool) {}

double NewtonianGravity::softening() const { return std::sqrt(m_soft2); }

std::string NewtonianGravity::description() const {
    std::string d = "Newtonian point-mass gravity (direct O(N^2) summation)";
    if (m_soft2 > 0) d += str::format(", Plummer softening %.3g AU", softening());
    if (!m_excluded.empty()) d += str::format(", %zu excluded pair(s)", m_excluded.size());
    return d;
}

std::unique_ptr<ForceModel> NewtonianGravity::clone() const { return std::make_unique<NewtonianGravity>(*this); }

void NewtonianGravity::excludePair(int source, int target) {
    if (!excluded(static_cast<std::size_t>(source), static_cast<std::size_t>(target)))
        m_excluded.emplace_back(source, target);
}

bool NewtonianGravity::excluded(std::size_t source, std::size_t target) const {
    for (const auto& p : m_excluded)
        if (static_cast<std::size_t>(p.first) == source && static_cast<std::size_t>(p.second) == target) return true;
    return false;
}

void NewtonianGravity::accumulate(const NBodySystem& sys, const Vec3* pos, const Vec3* /*vel*/, Vec3* acc) const {
    const std::size_t nM = sys.massiveCount();
    const std::size_t n = sys.size();
    const double* gm = sys.gm().data();
    const double eps2 = m_soft2;
    const bool hasExclusions = !m_excluded.empty();
    const std::size_t skip = m_skipSource >= 0 ? static_cast<std::size_t>(m_skipSource) : static_cast<std::size_t>(-1);

    // 1) Massive <-> massive: symmetric pair loop (Newton's third law halves the work).
    for (std::size_t i = 0; i < nM; ++i) {
        const Vec3 pi = pos[i];
        Vec3 ai{};
        for (std::size_t j = i + 1; j < nM; ++j) {
            const Vec3 d = pos[j] - pi;
            const double r2 = d.norm2() + eps2;
            const double inv = 1.0 / (r2 * std::sqrt(r2));
            const bool jOnI = j != skip && !(hasExclusions && excluded(j, i));
            const bool iOnJ = i != skip && !(hasExclusions && excluded(i, j));
            if (jOnI) ai += d * (gm[j] * inv);
            if (iOnJ) acc[j] -= d * (gm[i] * inv);
        }
        acc[i] += ai;
    }

    // 2) Massive -> massless: embarrassingly parallel over test bodies.
    if (n <= nM) return;
    auto testKernel = [&](std::size_t begin, std::size_t end) {
        for (std::size_t k = begin; k < end; ++k) {
            const std::size_t t = nM + k;
            const Vec3 pt = pos[t];
            Vec3 at{};
            for (std::size_t j = 0; j < nM; ++j) {
                if (j == skip) continue;
                if (hasExclusions && excluded(j, t)) continue;
                const Vec3 d = pos[j] - pt;
                const double r2 = d.norm2() + eps2;
                at += d * (gm[j] / (r2 * std::sqrt(r2)));
            }
            acc[t] += at;
        }
    };
    const std::size_t nTest = n - nM;
    if (m_pool && nTest >= 256) m_pool->parallelFor(nTest, testKernel);
    else testKernel(0, nTest);
}

// -----------------------------------------------------------------------------
// PostNewtonianGR
// -----------------------------------------------------------------------------
void PostNewtonianGR::accumulate(const NBodySystem& sys, const Vec3* pos, const Vec3* vel, Vec3* acc) const {
    if (m_central < 0 || static_cast<std::size_t>(m_central) >= sys.size()) return;
    const std::size_t c = static_cast<std::size_t>(m_central);
    const double mu = sys.gm()[c];
    constexpr double c2 = constants::C_AU_PER_DAY * constants::C_AU_PER_DAY;
    const std::size_t n = sys.size();
    for (std::size_t i = 0; i < n; ++i) {
        if (i == c) continue;
        const Vec3 r = pos[i] - pos[c];
        const Vec3 v = vel[i] - vel[c];
        const double r2 = r.norm2();
        const double rm = std::sqrt(r2);
        if (rm < 1e-6) continue; // inside the star - ignore
        const double v2 = v.norm2();
        const double rv = r.dot(v);
        const double f = mu / (c2 * r2 * rm);
        acc[i] += (r * (4.0 * mu / rm - v2) + v * (4.0 * rv)) * f;
    }
}

// -----------------------------------------------------------------------------
// ForceSet
// -----------------------------------------------------------------------------
ForceSet::ForceSet(const ForceSet& other) : m_evaluations(other.m_evaluations) {
    for (const auto& m : other.m_models) m_models.push_back(m->clone());
}

ForceSet& ForceSet::operator=(const ForceSet& other) {
    if (this == &other) return *this;
    m_models.clear();
    for (const auto& m : other.m_models) m_models.push_back(m->clone());
    m_evaluations = other.m_evaluations;
    return *this;
}

void ForceSet::add(std::unique_ptr<ForceModel> model) { m_models.push_back(std::move(model)); }

bool ForceSet::remove(const std::string& name) {
    const auto it = std::find_if(m_models.begin(), m_models.end(),
                                 [&](const std::unique_ptr<ForceModel>& m) { return m->name() == name; });
    if (it == m_models.end()) return false;
    m_models.erase(it);
    return true;
}

ForceModel* ForceSet::find(const std::string& name) const {
    for (const auto& m : m_models)
        if (m->name() == name) return m.get();
    return nullptr;
}

NewtonianGravity* ForceSet::gravity() const {
    for (const auto& m : m_models)
        if (auto* g = dynamic_cast<NewtonianGravity*>(m.get())) return g;
    return nullptr;
}

void ForceSet::compute(const NBodySystem& sys, const Vec3* pos, const Vec3* vel, Vec3* acc) const {
    const std::size_t n = sys.size();
    for (std::size_t i = 0; i < n; ++i) acc[i] = Vec3{};
    for (const auto& m : m_models) m->accumulate(sys, pos, vel, acc);
    ++m_evaluations;
}

} // namespace solsim::physics
