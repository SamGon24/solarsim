// =============================================================================
//  SolarSim - Mat3.hpp
//  3x3 matrices: rotations between reference frames (ecliptic <-> equatorial,
//  perifocal -> inertial, precession), camera bases and small linear solves.
// =============================================================================
#pragma once

#include "Vec3.hpp"

#include <cmath>

namespace solsim {

struct Mat3 {
    // Row-major storage: m[row][col]
    double m[3][3]{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    static Mat3 identity() { return Mat3{}; }

    static Mat3 zero() {
        Mat3 r;
        for (auto& row : r.m)
            for (double& v : row) v = 0.0;
        return r;
    }

    static Mat3 fromRows(const Vec3& r0, const Vec3& r1, const Vec3& r2) {
        Mat3 r;
        r.m[0][0] = r0.x; r.m[0][1] = r0.y; r.m[0][2] = r0.z;
        r.m[1][0] = r1.x; r.m[1][1] = r1.y; r.m[1][2] = r1.z;
        r.m[2][0] = r2.x; r.m[2][1] = r2.y; r.m[2][2] = r2.z;
        return r;
    }

    static Mat3 fromColumns(const Vec3& c0, const Vec3& c1, const Vec3& c2) {
        return fromRows(c0, c1, c2).transposed();
    }

    // Active rotations (rotate a vector by +angle about the axis).
    static Mat3 rotX(double a) {
        const double c = std::cos(a), s = std::sin(a);
        return fromRows({1, 0, 0}, {0, c, -s}, {0, s, c});
    }
    static Mat3 rotY(double a) {
        const double c = std::cos(a), s = std::sin(a);
        return fromRows({c, 0, s}, {0, 1, 0}, {-s, 0, c});
    }
    static Mat3 rotZ(double a) {
        const double c = std::cos(a), s = std::sin(a);
        return fromRows({c, -s, 0}, {s, c, 0}, {0, 0, 1});
    }

    // Rodrigues' rotation formula.
    static Mat3 axisAngle(const Vec3& axisIn, double angle) {
        const Vec3 k = axisIn.normalized();
        const double c = std::cos(angle), s = std::sin(angle), t = 1.0 - c;
        return fromRows({t * k.x * k.x + c,       t * k.x * k.y - s * k.z, t * k.x * k.z + s * k.y},
                        {t * k.x * k.y + s * k.z, t * k.y * k.y + c,       t * k.y * k.z - s * k.x},
                        {t * k.x * k.z - s * k.y, t * k.y * k.z + s * k.x, t * k.z * k.z + c});
    }

    Vec3 row(int r) const { return {m[r][0], m[r][1], m[r][2]}; }
    Vec3 col(int c) const { return {m[0][c], m[1][c], m[2][c]}; }

    Vec3 operator*(const Vec3& v) const {
        return {m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
                m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
                m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
    }

    Mat3 operator*(const Mat3& o) const {
        Mat3 r = Mat3::zero();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k) r.m[i][j] += m[i][k] * o.m[k][j];
        return r;
    }

    Mat3 operator*(double s) const {
        Mat3 r = *this;
        for (auto& row : r.m)
            for (double& v : row) v *= s;
        return r;
    }

    Mat3 transposed() const {
        Mat3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = m[j][i];
        return r;
    }

    double determinant() const {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
               m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
               m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }

    // Returns false if the matrix is (numerically) singular.
    bool inverse(Mat3& out) const {
        const double det = determinant();
        if (std::fabs(det) < 1e-300) return false;
        const double inv = 1.0 / det;
        out.m[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv;
        out.m[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv;
        out.m[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv;
        out.m[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv;
        out.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv;
        out.m[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv;
        out.m[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv;
        out.m[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv;
        out.m[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv;
        return true;
    }
};

} // namespace solsim
