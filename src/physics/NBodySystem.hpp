// =============================================================================
//  SolarSim - NBodySystem.hpp
//  Storage for the gravitating system in structure-of-arrays layout:
//
//      index range                     content
//      [0, massiveCount)               named bodies with GM > 0 (Sun first)
//      [massiveCount, namedCount)      named massless bodies (comets, probes)
//      [namedCount, size)              anonymous test particles (belts)
//
//  Test particles feel the massive bodies but do not attract anything, which
//  reduces the cost from O(N^2) to O(N_massive * N). All positions and
//  velocities are barycentric, ecliptic J2000, in AU and AU/day.
// =============================================================================
#pragma once

#include "../astro/BodyCatalog.hpp"
#include "../core/Vec3.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace solsim::physics {

struct BodyInfo {
    std::string name;
    astro::BodyKind kind = astro::BodyKind::Planet;
    int parent = -1;                 // index of the parent body (-1 for the root)
    double radius = 0.0;             // AU
    double rotationPeriod = 0.0;     // days, signed (negative = retrograde)
    double rotationPhase0 = 0.0;     // rad at epoch (cosmetic)
    Vec3 pole{0, 0, 1};              // spin axis, ecliptic J2000 unit vector
    double albedo = 0.3;
    double typicalDistance = 0.0;    // AU from parent (display scaling, dt choice)
    astro::SurfaceStyle style = astro::SurfaceStyle::Rocky;
    std::uint32_t color = 0xAAAAAA;
    std::uint32_t color2 = 0x666666;
    astro::RingSystem rings;
    bool approximate = false;
    std::string notes;
};

struct ParticleGroup {
    std::string name;
    std::uint32_t color = 0x888888;
    std::size_t begin = 0;  // absolute index range [begin, end)
    std::size_t end = 0;
};

class NBodySystem {
public:
    // ---- structure ----------------------------------------------------------
    // Inserts a named body keeping the layout invariant; returns its index.
    int addBody(const BodyInfo& info, double gm, const Vec3& r, const Vec3& v);
    // Creates (or extends) a particle group and returns its index.
    int addParticleGroup(const std::string& name, std::uint32_t color);
    void addParticle(int group, const Vec3& r, const Vec3& v);
    void removeBody(int index); // named bodies only

    std::size_t size() const { return m_pos.size(); }
    std::size_t massiveCount() const { return m_massive; }
    std::size_t namedCount() const { return m_info.size(); }
    std::size_t particleCount() const { return size() - namedCount(); }
    std::uint64_t revision() const { return m_revision; }

    int indexOf(const std::string& name) const; // case-insensitive, -1 if absent
    const BodyInfo& info(std::size_t i) const { return m_info[i]; }
    BodyInfo& info(std::size_t i) { return m_info[i]; }
    const std::vector<ParticleGroup>& groups() const { return m_groups; }

    // ---- state access -------------------------------------------------------
    std::vector<Vec3>& positions() { return m_pos; }
    std::vector<Vec3>& velocities() { return m_vel; }
    std::vector<Vec3>& accelerations() { return m_acc; }
    const std::vector<Vec3>& positions() const { return m_pos; }
    const std::vector<Vec3>& velocities() const { return m_vel; }
    const std::vector<Vec3>& accelerations() const { return m_acc; }
    const std::vector<double>& gm() const { return m_gm; }
    void setGM(std::size_t i, double gm); // keeps massive/massless classification

    const Vec3& pos(std::size_t i) const { return m_pos[i]; }
    const Vec3& vel(std::size_t i) const { return m_vel[i]; }

    // State relative to a body's parent (or to `relativeTo` if >= 0).
    Vec3 relativePosition(std::size_t i, int relativeTo) const;
    Vec3 relativeVelocity(std::size_t i, int relativeTo) const;
    // GM of the two-body problem body+parent (for osculating elements).
    double twoBodyMu(std::size_t i) const;

    // ---- global quantities --------------------------------------------------
    double totalGM() const;
    Vec3 centerOfMass() const;
    Vec3 centerOfMassVelocity() const;
    // Removes the barycentre offset and drift (sets total momentum to zero).
    void moveToBarycentricFrame();

    // Sun = index 0 by construction; returns -1 if empty.
    int centralBody() const { return m_info.empty() ? -1 : 0; }

private:
    void shiftIndicesFrom(std::size_t index, int delta);

    std::vector<Vec3> m_pos, m_vel, m_acc;
    std::vector<double> m_gm;
    std::vector<BodyInfo> m_info;
    std::vector<ParticleGroup> m_groups;
    std::size_t m_massive = 0;
    std::uint64_t m_revision = 0;
};

} // namespace solsim::physics
