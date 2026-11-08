#pragma once

#include "phys/math.h"

namespace phys {

// Order matters: the collision dispatcher sorts pairs by this value.
enum class ShapeType { Sphere = 0, Box = 1, Plane = 2 };

struct Shape {
    ShapeType type = ShapeType::Sphere;
    Real radius = Real(0.5);              // Sphere
    Vec3 halfExtents{0.5, 0.5, 0.5};      // Box
    // Plane: the body's local +Y axis is the normal, passing through the body's position.

    static Shape sphere(Real r) {
        Shape s;
        s.type = ShapeType::Sphere;
        s.radius = r;
        return s;
    }
    static Shape box(const Vec3& halfExtents) {
        Shape s;
        s.type = ShapeType::Box;
        s.halfExtents = halfExtents;
        return s;
    }
    static Shape plane() {
        Shape s;
        s.type = ShapeType::Plane;
        return s;
    }

    Real volume() const {
        switch (type) {
            case ShapeType::Sphere: return Real(4) / 3 * kPi * radius * radius * radius;
            case ShapeType::Box: return 8 * halfExtents.x * halfExtents.y * halfExtents.z;
            case ShapeType::Plane: return 0;
        }
        return 0;
    }

    // Diagonal of the inertia tensor for the given mass, in the shape's local frame.
    Vec3 inertiaDiagonal(Real mass) const {
        switch (type) {
            case ShapeType::Sphere: {
                Real i = Real(0.4) * mass * radius * radius;
                return {i, i, i};
            }
            case ShapeType::Box: {
                Vec3 s = halfExtents * Real(2);
                Real k = mass / 12;
                return {k * (s.y * s.y + s.z * s.z), k * (s.x * s.x + s.z * s.z),
                        k * (s.x * s.x + s.y * s.y)};
            }
            case ShapeType::Plane: return {};
        }
        return {};
    }
};

}  // namespace phys
