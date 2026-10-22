// =============================================================================
//  SolarSim - Vec3.hpp
//  Small, header-only 2D/3D vector templates used for positions, velocities,
//  accelerations, colors and screen coordinates.
// =============================================================================
#pragma once

#include <cmath>
#include <cstddef>
#include <ostream>

namespace solsim {

template <typename T>
struct Vector3 {
    T x{}, y{}, z{};

    constexpr Vector3() = default;
    constexpr Vector3(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

    template <typename U>
    constexpr explicit Vector3(const Vector3<U>& o)
        : x(static_cast<T>(o.x)), y(static_cast<T>(o.y)), z(static_cast<T>(o.z)) {}

    static constexpr Vector3 zero() { return {T(0), T(0), T(0)}; }
    static constexpr Vector3 unitX() { return {T(1), T(0), T(0)}; }
    static constexpr Vector3 unitY() { return {T(0), T(1), T(0)}; }
    static constexpr Vector3 unitZ() { return {T(0), T(0), T(1)}; }

    constexpr T& operator[](std::size_t i) { return i == 0 ? x : (i == 1 ? y : z); }
    constexpr const T& operator[](std::size_t i) const { return i == 0 ? x : (i == 1 ? y : z); }

    constexpr Vector3 operator+(const Vector3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr Vector3 operator-(const Vector3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr Vector3 operator-() const { return {-x, -y, -z}; }
    constexpr Vector3 operator*(T s) const { return {x * s, y * s, z * s}; }
    constexpr Vector3 operator/(T s) const { return {x / s, y / s, z / s}; }

    constexpr Vector3& operator+=(const Vector3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    constexpr Vector3& operator-=(const Vector3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    constexpr Vector3& operator*=(T s) { x *= s; y *= s; z *= s; return *this; }
    constexpr Vector3& operator/=(T s) { x /= s; y /= s; z /= s; return *this; }

    constexpr bool operator==(const Vector3& o) const { return x == o.x && y == o.y && z == o.z; }
    constexpr bool operator!=(const Vector3& o) const { return !(*this == o); }

    constexpr T dot(const Vector3& o) const { return x * o.x + y * o.y + z * o.z; }
    constexpr Vector3 cross(const Vector3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    constexpr Vector3 hadamard(const Vector3& o) const { return {x * o.x, y * o.y, z * o.z}; }

    constexpr T norm2() const { return x * x + y * y + z * z; }
    T norm() const { return std::sqrt(norm2()); }

    Vector3 normalized() const {
        const T n = norm();
        return n > T(0) ? (*this) / n : Vector3{};
    }

    T maxAbs() const {
        const T ax = std::fabs(x), ay = std::fabs(y), az = std::fabs(z);
        return ax > ay ? (ax > az ? ax : az) : (ay > az ? ay : az);
    }

    bool isFinite() const { return std::isfinite(x) && std::isfinite(y) && std::isfinite(z); }
};

template <typename T>
constexpr Vector3<T> operator*(T s, const Vector3<T>& v) { return v * s; }

template <typename T>
inline T distance(const Vector3<T>& a, const Vector3<T>& b) { return (a - b).norm(); }

template <typename T>
inline T angleBetween(const Vector3<T>& a, const Vector3<T>& b) {
    // atan2 formulation is accurate for both small and large angles.
    return std::atan2(a.cross(b).norm(), a.dot(b));
}

template <typename T>
constexpr Vector3<T> lerp(const Vector3<T>& a, const Vector3<T>& b, T t) { return a + (b - a) * t; }

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vector3<T>& v) {
    return os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
}

using Vec3  = Vector3<double>;
using Vec3f = Vector3<float>;

// -----------------------------------------------------------------------------
template <typename T>
struct Vector2 {
    T x{}, y{};
    constexpr Vector2() = default;
    constexpr Vector2(T x_, T y_) : x(x_), y(y_) {}
    constexpr Vector2 operator+(const Vector2& o) const { return {x + o.x, y + o.y}; }
    constexpr Vector2 operator-(const Vector2& o) const { return {x - o.x, y - o.y}; }
    constexpr Vector2 operator*(T s) const { return {x * s, y * s}; }
    constexpr T dot(const Vector2& o) const { return x * o.x + y * o.y; }
    constexpr T cross(const Vector2& o) const { return x * o.y - y * o.x; }
    T norm() const { return std::sqrt(x * x + y * y); }
};

using Vec2 = Vector2<double>;

} // namespace solsim
