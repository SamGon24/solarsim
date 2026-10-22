// =============================================================================
//  SolarSim - Camera.hpp
//  Orbit camera with perspective projection. The camera looks at `target`
//  from a point on a sphere of radius `distance` (spherical angles yaw/pitch
//  around the ecliptic z axis).
// =============================================================================
#pragma once

#include "../core/Vec3.hpp"

namespace solsim::render {

struct Projected {
    double x = 0.0;      // pixels
    double y = 0.0;      // pixels (down)
    double depth = 0.0;  // camera-space z (distance along view direction)
};

class Camera {
public:
    Camera();

    void setViewport(int width, int height, double pixelAspect = 1.0);
    int width() const { return m_w; }
    int height() const { return m_h; }

    void setTarget(const Vec3& t) { m_target = t; update(); }
    void setDistance(double d);
    void setAngles(double yawRad, double pitchRad);
    void setFieldOfView(double fovYRad) { m_fovY = fovYRad; update(); }
    void orbit(double dYaw, double dPitch) { setAngles(m_yaw + dYaw, m_pitch + dPitch); }
    void zoom(double factor) { setDistance(m_distance * factor); }

    const Vec3& target() const { return m_target; }
    double distance() const { return m_distance; }
    double yaw() const { return m_yaw; }
    double pitch() const { return m_pitch; }
    double fieldOfView() const { return m_fovY; }

    const Vec3& eye() const { return m_eye; }
    const Vec3& forward() const { return m_forward; }
    const Vec3& right() const { return m_right; }
    const Vec3& up() const { return m_up; }
    double focalPixels() const { return m_focal; }
    double nearPlane() const { return m_near; }

    // World point -> screen. Returns false when behind the near plane.
    bool project(const Vec3& world, Projected& out) const;
    // Direction at infinity (stars) -> screen.
    bool projectDirection(const Vec3& dir, Projected& out) const;
    // Primary ray through pixel coordinates (sx, sy).
    Vec3 rayDirection(double sx, double sy) const;

    // Projected radius (pixels) of a sphere of radius R at camera depth z.
    double projectedRadius(double R, double depth) const { return depth > 0 ? m_focal * R / depth : 0.0; }

private:
    void update();

    Vec3 m_target{};
    double m_distance = 10.0;
    double m_yaw = -1.5707963267948966;
    double m_pitch = 0.5;
    double m_fovY = 0.7853981633974483; // 45 deg
    int m_w = 640, m_h = 480;
    double m_pixelAspect = 1.0;
    double m_focal = 1.0;
    double m_near = 1e-6;
    Vec3 m_eye, m_forward, m_right, m_up;
};

} // namespace solsim::render
