// =============================================================================
//  SolarSim - Camera.cpp
// =============================================================================
#include "Camera.hpp"

#include "../core/MathUtil.hpp"

#include <algorithm>
#include <cmath>

namespace solsim::render {

Camera::Camera() { update(); }

void Camera::setViewport(int width, int height, double pixelAspect) {
    m_w = std::max(1, width);
    m_h = std::max(1, height);
    m_pixelAspect = pixelAspect > 0 ? pixelAspect : 1.0;
    update();
}

void Camera::setDistance(double d) {
    m_distance = std::clamp(d, 1e-7, 1e7);
    update();
}

void Camera::setAngles(double yawRad, double pitchRad) {
    m_yaw = math::wrapPi(yawRad);
    m_pitch = std::clamp(pitchRad, -1.5690, 1.5690); // avoid the pole singularity
    update();
}

void Camera::update() {
    const double cp = std::cos(m_pitch), sp = std::sin(m_pitch);
    const Vec3 offset{cp * std::cos(m_yaw), cp * std::sin(m_yaw), sp};
    m_eye = m_target + offset * m_distance;
    m_forward = (-offset).normalized();
    Vec3 r = m_forward.cross(Vec3::unitZ());
    if (r.norm() < 1e-9) r = Vec3::unitX();
    m_right = r.normalized();
    m_up = m_right.cross(m_forward).normalized();
    m_focal = 0.5 * m_h / std::tan(0.5 * m_fovY);
    m_near = std::max(1e-9, m_distance * 1e-4);
}

bool Camera::project(const Vec3& world, Projected& out) const {
    const Vec3 d = world - m_eye;
    const double z = d.dot(m_forward);
    if (z <= m_near) return false;
    out.depth = z;
    out.x = 0.5 * m_w + m_focal * d.dot(m_right) / z / m_pixelAspect;
    out.y = 0.5 * m_h - m_focal * d.dot(m_up) / z;
    return true;
}

bool Camera::projectDirection(const Vec3& dir, Projected& out) const {
    const double z = dir.dot(m_forward);
    if (z <= 1e-6) return false;
    out.depth = 1e30;
    out.x = 0.5 * m_w + m_focal * dir.dot(m_right) / z / m_pixelAspect;
    out.y = 0.5 * m_h - m_focal * dir.dot(m_up) / z;
    return true;
}

Vec3 Camera::rayDirection(double sx, double sy) const {
    const double cx = (sx - 0.5 * m_w) * m_pixelAspect / m_focal;
    const double cy = (0.5 * m_h - sy) / m_focal;
    return (m_forward + m_right * cx + m_up * cy).normalized();
}

} // namespace solsim::render
