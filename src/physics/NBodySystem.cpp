// =============================================================================
//  SolarSim - NBodySystem.cpp
// =============================================================================
#include "NBodySystem.hpp"

#include "../core/StringUtil.hpp"

#include <stdexcept>

namespace solsim::physics {

void NBodySystem::shiftIndicesFrom(std::size_t index, int delta) {
    for (auto& inf : m_info)
        if (inf.parent >= 0 && static_cast<std::size_t>(inf.parent) >= index) inf.parent += delta;
    for (auto& g : m_groups) {
        if (g.begin >= index) g.begin = static_cast<std::size_t>(static_cast<long long>(g.begin) + delta);
        if (g.end >= index) g.end = static_cast<std::size_t>(static_cast<long long>(g.end) + delta);
    }
}

int NBodySystem::addBody(const BodyInfo& info, double gm, const Vec3& r, const Vec3& v) {
    if (gm < 0.0) throw std::invalid_argument("negative GM for body " + info.name);
    const bool massive = gm > 0.0;
    const std::size_t index = massive ? m_massive : m_info.size();

    // Parent indices referring to positions >= index move by one.
    BodyInfo copy = info;
    shiftIndicesFrom(index, +1);
    if (copy.parent >= 0 && static_cast<std::size_t>(copy.parent) >= index) copy.parent += 1;

    m_pos.insert(m_pos.begin() + static_cast<std::ptrdiff_t>(index), r);
    m_vel.insert(m_vel.begin() + static_cast<std::ptrdiff_t>(index), v);
    m_acc.insert(m_acc.begin() + static_cast<std::ptrdiff_t>(index), Vec3{});
    m_gm.insert(m_gm.begin() + static_cast<std::ptrdiff_t>(index), gm);
    m_info.insert(m_info.begin() + static_cast<std::ptrdiff_t>(index), copy);
    if (massive) ++m_massive;
    ++m_revision;
    return static_cast<int>(index);
}

void NBodySystem::removeBody(int index) {
    if (index < 0 || static_cast<std::size_t>(index) >= m_info.size())
        throw std::out_of_range("removeBody: index out of range");
    const std::size_t i = static_cast<std::size_t>(index);
    for (auto& inf : m_info)
        if (inf.parent == index) inf.parent = m_info[i].parent;
    if (m_gm[i] > 0.0) --m_massive;
    m_pos.erase(m_pos.begin() + index);
    m_vel.erase(m_vel.begin() + index);
    m_acc.erase(m_acc.begin() + index);
    m_gm.erase(m_gm.begin() + index);
    m_info.erase(m_info.begin() + index);
    for (auto& inf : m_info)
        if (inf.parent > index) inf.parent -= 1;
    for (auto& g : m_groups) {
        if (g.begin > i) --g.begin;
        if (g.end > i) --g.end;
    }
    ++m_revision;
}

int NBodySystem::addParticleGroup(const std::string& name, std::uint32_t color) {
    ParticleGroup g;
    g.name = name;
    g.color = color;
    g.begin = g.end = size();
    m_groups.push_back(g);
    ++m_revision;
    return static_cast<int>(m_groups.size()) - 1;
}

void NBodySystem::addParticle(int group, const Vec3& r, const Vec3& v) {
    if (group < 0 || static_cast<std::size_t>(group) >= m_groups.size())
        throw std::out_of_range("addParticle: bad group");
    ParticleGroup& g = m_groups[static_cast<std::size_t>(group)];
    if (g.end != size()) throw std::logic_error("particles must be appended to the most recent group");
    m_pos.push_back(r);
    m_vel.push_back(v);
    m_acc.emplace_back();
    m_gm.push_back(0.0);
    ++g.end;
    ++m_revision;
}

void NBodySystem::setGM(std::size_t i, double gm) {
    const bool wasMassive = m_gm[i] > 0.0;
    if (wasMassive != (gm > 0.0)) throw std::logic_error("setGM cannot change massive/massless class");
    m_gm[i] = gm;
    ++m_revision;
}

int NBodySystem::indexOf(const std::string& name) const {
    for (std::size_t i = 0; i < m_info.size(); ++i)
        if (str::iequals(m_info[i].name, name)) return static_cast<int>(i);
    return -1;
}

Vec3 NBodySystem::relativePosition(std::size_t i, int relativeTo) const {
    int ref = relativeTo;
    if (ref < 0 && i < m_info.size()) ref = m_info[i].parent;
    if (ref < 0 && i >= m_info.size()) ref = centralBody();
    return ref < 0 ? m_pos[i] : m_pos[i] - m_pos[static_cast<std::size_t>(ref)];
}

Vec3 NBodySystem::relativeVelocity(std::size_t i, int relativeTo) const {
    int ref = relativeTo;
    if (ref < 0 && i < m_info.size()) ref = m_info[i].parent;
    if (ref < 0 && i >= m_info.size()) ref = centralBody();
    return ref < 0 ? m_vel[i] : m_vel[i] - m_vel[static_cast<std::size_t>(ref)];
}

double NBodySystem::twoBodyMu(std::size_t i) const {
    int p = (i < m_info.size()) ? m_info[i].parent : centralBody();
    if (p < 0) return m_gm[i];
    return m_gm[static_cast<std::size_t>(p)] + m_gm[i];
}

double NBodySystem::totalGM() const {
    double s = 0.0;
    for (std::size_t i = 0; i < m_massive; ++i) s += m_gm[i];
    return s;
}

Vec3 NBodySystem::centerOfMass() const {
    Vec3 c;
    double m = 0.0;
    for (std::size_t i = 0; i < m_massive; ++i) {
        c += m_pos[i] * m_gm[i];
        m += m_gm[i];
    }
    return m > 0.0 ? c / m : c;
}

Vec3 NBodySystem::centerOfMassVelocity() const {
    Vec3 c;
    double m = 0.0;
    for (std::size_t i = 0; i < m_massive; ++i) {
        c += m_vel[i] * m_gm[i];
        m += m_gm[i];
    }
    return m > 0.0 ? c / m : c;
}

void NBodySystem::moveToBarycentricFrame() {
    const Vec3 c = centerOfMass();
    const Vec3 cv = centerOfMassVelocity();
    for (std::size_t i = 0; i < size(); ++i) {
        m_pos[i] -= c;
        m_vel[i] -= cv;
    }
    ++m_revision;
}

} // namespace solsim::physics
