#pragma once

#include <algorithm>
#include <cmath>

namespace phys {

#ifdef PHYS_USE_FLOAT
using Real = float;
#else
using Real = double;
#endif

constexpr Real kPi = Real(3.14159265358979323846);
constexpr Real kEpsilon = Real(1e-9);

struct Vec3 {
    Real x = 0, y = 0, z = 0;

    constexpr Vec3() = default;
    constexpr Vec3(Real x_, Real y_, Real z_) : x(x_), y(y_), z(z_) {}

    Real& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
    Real operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }

    Vec3 operator-() const { return {-x, -y, -z}; }
    Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    Vec3 operator*(Real s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(Real s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    Vec3& operator-=(const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    Vec3& operator*=(Real s) { x *= s; y *= s; z *= s; return *this; }

    Real lengthSq() const { return x * x + y * y + z * z; }
    Real length() const { return std::sqrt(lengthSq()); }
    Vec3 normalized() const {
        Real len = length();
        return len > kEpsilon ? *this / len : Vec3{};
    }
};

inline Vec3 operator*(Real s, const Vec3& v) { return v * s; }
inline Real dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline Vec3 mul(const Vec3& a, const Vec3& b) { return {a.x * b.x, a.y * b.y, a.z * b.z}; }
inline Vec3 vmin(const Vec3& a, const Vec3& b) {
    return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)};
}
inline Vec3 vmax(const Vec3& a, const Vec3& b) {
    return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
}
inline Vec3 vabs(const Vec3& a) { return {std::abs(a.x), std::abs(a.y), std::abs(a.z)}; }

// Builds two unit tangents t1, t2 so that (n, t1, t2) is an orthonormal basis.
// Deterministic for a given n, which keeps friction warm starting consistent.
inline void orthonormalBasis(const Vec3& n, Vec3& t1, Vec3& t2) {
    if (std::abs(n.x) >= Real(0.57735)) {
        t1 = Vec3{n.y, -n.x, 0}.normalized();
    } else {
        t1 = Vec3{0, n.z, -n.y}.normalized();
    }
    t2 = cross(n, t1);
}

// 3x3 matrix stored as rows.
struct Mat3 {
    Vec3 r[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    Mat3() = default;
    Mat3(const Vec3& r0, const Vec3& r1, const Vec3& r2) : r{r0, r1, r2} {}

    static Mat3 identity() { return {}; }
    static Mat3 zero() { return {Vec3{}, Vec3{}, Vec3{}}; }
    static Mat3 diagonal(const Vec3& d) { return {{d.x, 0, 0}, {0, d.y, 0}, {0, 0, d.z}}; }

    Vec3 col(int i) const { return {r[0][i], r[1][i], r[2][i]}; }

    Vec3 operator*(const Vec3& v) const { return {dot(r[0], v), dot(r[1], v), dot(r[2], v)}; }

    Mat3 operator*(const Mat3& m) const {
        Mat3 out;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                out.r[i][j] = r[i][0] * m.r[0][j] + r[i][1] * m.r[1][j] + r[i][2] * m.r[2][j];
        return out;
    }

    Mat3 transposed() const { return {col(0), col(1), col(2)}; }

    Real determinant() const { return dot(r[0], cross(r[1], r[2])); }

    Mat3 inverse() const {
        Real det = determinant();
        if (std::abs(det) < kEpsilon) return zero();
        // Columns of the inverse are cross products of rows, divided by det.
        Vec3 c0 = cross(r[1], r[2]) / det;
        Vec3 c1 = cross(r[2], r[0]) / det;
        Vec3 c2 = cross(r[0], r[1]) / det;
        return Mat3{c0, c1, c2}.transposed();
    }
};

// Unit quaternion for orientation. w is the scalar part.
struct Quat {
    Real w = 1, x = 0, y = 0, z = 0;

    constexpr Quat() = default;
    constexpr Quat(Real w_, Real x_, Real y_, Real z_) : w(w_), x(x_), y(y_), z(z_) {}

    static Quat identity() { return {}; }

    static Quat fromAxisAngle(const Vec3& axis, Real angle) {
        Vec3 a = axis.normalized();
        Real s = std::sin(angle * Real(0.5));
        return {std::cos(angle * Real(0.5)), a.x * s, a.y * s, a.z * s};
    }

    Quat operator*(const Quat& q) const {
        return {w * q.w - x * q.x - y * q.y - z * q.z,
                w * q.x + x * q.w + y * q.z - z * q.y,
                w * q.y - x * q.z + y * q.w + z * q.x,
                w * q.z + x * q.y - y * q.x + z * q.w};
    }

    Quat conjugate() const { return {w, -x, -y, -z}; }

    Real lengthSq() const { return w * w + x * x + y * y + z * z; }

    Quat normalized() const {
        Real len = std::sqrt(lengthSq());
        if (len < kEpsilon) return identity();
        Real inv = Real(1) / len;
        return {w * inv, x * inv, y * inv, z * inv};
    }

    Vec3 rotate(const Vec3& v) const {
        // v' = v + 2w(u x v) + 2u x (u x v), with u the vector part.
        Vec3 u{x, y, z};
        Vec3 t = cross(u, v) * Real(2);
        return v + t * w + cross(u, t);
    }

    Vec3 inverseRotate(const Vec3& v) const { return conjugate().rotate(v); }

    Mat3 toMat3() const {
        Real xx = x * x, yy = y * y, zz = z * z;
        Real xy = x * y, xz = x * z, yz = y * z;
        Real wx = w * x, wy = w * y, wz = w * z;
        return {{1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy)},
                {2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx)},
                {2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy)}};
    }

    // Normalized linear blend along the shortest path. Good for small steps (rendering).
    static Quat nlerp(const Quat& a, Quat b, Real t) {
        if (a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z < 0) b = {-b.w, -b.x, -b.y, -b.z};
        return Quat{a.w + (b.w - a.w) * t, a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t,
                    a.z + (b.z - a.z) * t}
            .normalized();
    }

    // Advances the orientation by angular velocity w (world space) over dt.
    Quat integrated(const Vec3& omega, Real dt) const {
        Quat spin{0, omega.x, omega.y, omega.z};
        Quat dq = spin * (*this);
        Real h = Real(0.5) * dt;
        return Quat{w + dq.w * h, x + dq.x * h, y + dq.y * h, z + dq.z * h}.normalized();
    }
};

}  // namespace phys
