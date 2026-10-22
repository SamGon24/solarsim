// =============================================================================
//  SolarSim - Forces.hpp
//  Force (acceleration) models. Each model ADDS its contribution for the given
//  candidate positions/velocities, so integrators can evaluate accelerations at
//  intermediate stages without touching the system state.
//
//   * NewtonianGravity      a_i = sum_j GM_j (r_j - r_i) / (|r_j - r_i|^2 + eps^2)^(3/2)
//   * PostNewtonianGR       1PN Schwarzschild correction from the central star
//                           (Einstein-Infeld-Hoffmann test-particle limit):
//                             a = GM/(c^2 r^3) [ (4GM/r - v^2) r + 4 (r.v) v ]
//                           reproduces Mercury's 43"/century perihelion advance
//   * ForceSet              ordered collection; zeroes and accumulates
// =============================================================================
#pragma once

#include "NBodySystem.hpp"

#include "../core/Vec3.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace solsim {
class ThreadPool;
}

namespace solsim::physics {

class ForceModel {
public:
    virtual ~ForceModel() = default;
    virtual std::string name() const = 0;
    virtual std::string description() const = 0;
    virtual bool velocityDependent() const { return false; }
    virtual void accumulate(const NBodySystem& sys, const Vec3* pos, const Vec3* vel, Vec3* acc) const = 0;
    virtual std::unique_ptr<ForceModel> clone() const = 0;
};

class NewtonianGravity final : public ForceModel {
public:
    explicit NewtonianGravity(double softeningAU = 0.0, ThreadPool* pool = nullptr);

    std::string name() const override { return "newton"; }
    std::string description() const override;
    void accumulate(const NBodySystem& sys, const Vec3* pos, const Vec3* vel, Vec3* acc) const override;
    std::unique_ptr<ForceModel> clone() const override;

    // Exclude the gravity of body `source` on body `target` (used by the
    // patched-conic mission planner).
    void excludePair(int source, int target);
    void clearExclusions() { m_excluded.clear(); }
    const std::vector<std::pair<int, int>>& exclusions() const { return m_excluded; }

    // Skip a body as a SOURCE entirely (Wisdom-Holman interaction Hamiltonian
    // omits the central star).
    void setSkipSource(int index) { m_skipSource = index; }
    void setThreadPool(ThreadPool* pool) { m_pool = pool; }
    void setSoftening(double eps) { m_soft2 = eps * eps; }
    double softening() const;

private:
    bool excluded(std::size_t source, std::size_t target) const;

    double m_soft2 = 0.0;
    ThreadPool* m_pool = nullptr;
    int m_skipSource = -1;
    std::vector<std::pair<int, int>> m_excluded;
};

class PostNewtonianGR final : public ForceModel {
public:
    explicit PostNewtonianGR(int centralIndex = 0) : m_central(centralIndex) {}
    std::string name() const override { return "gr1pn"; }
    std::string description() const override {
        return "1PN Schwarzschild correction (general relativity) around the central star";
    }
    bool velocityDependent() const override { return true; }
    void accumulate(const NBodySystem& sys, const Vec3* pos, const Vec3* vel, Vec3* acc) const override;
    std::unique_ptr<ForceModel> clone() const override { return std::make_unique<PostNewtonianGR>(*this); }

private:
    int m_central;
};

class ForceSet {
public:
    ForceSet() = default;
    ForceSet(const ForceSet& other);
    ForceSet& operator=(const ForceSet& other);
    ForceSet(ForceSet&&) noexcept = default;
    ForceSet& operator=(ForceSet&&) noexcept = default;

    void add(std::unique_ptr<ForceModel> model);
    bool remove(const std::string& name);
    ForceModel* find(const std::string& name) const;
    NewtonianGravity* gravity() const; // first NewtonianGravity, if any
    bool has(const std::string& name) const { return find(name) != nullptr; }
    std::size_t count() const { return m_models.size(); }
    const std::vector<std::unique_ptr<ForceModel>>& models() const { return m_models; }

    // acc[i] = sum of all models. Increments the evaluation counter.
    void compute(const NBodySystem& sys, const Vec3* pos, const Vec3* vel, Vec3* acc) const;
    std::uint64_t evaluations() const { return m_evaluations; }
    // For integrators that evaluate force models directly (Wisdom-Holman).
    void countEvaluation() const { ++m_evaluations; }
    void resetCounter() const { m_evaluations = 0; }

private:
    std::vector<std::unique_ptr<ForceModel>> m_models;
    mutable std::uint64_t m_evaluations = 0;
};

} // namespace solsim::physics
