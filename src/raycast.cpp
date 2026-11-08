#include "phys/raycast.h"

#include <limits>

#include "phys/body.h"

namespace phys {

std::optional<RayHit> raycastBody(RigidBody& body, const Ray& ray) {
    Vec3 dir = ray.direction.normalized();
    if (dir.lengthSq() == 0) return std::nullopt;
    const Shape& shape = body.shape();

    switch (shape.type) {
        case ShapeType::Sphere: {
            // Solve |o + t d - c|^2 = r^2 for the smallest t >= 0.
            Vec3 m = ray.origin - body.position;
            Real b = dot(m, dir);
            Real c = m.lengthSq() - shape.radius * shape.radius;
            if (c <= 0) return std::nullopt;   // Starts inside.
            if (b > 0) return std::nullopt;    // Pointing away.
            Real disc = b * b - c;
            if (disc < 0) return std::nullopt;
            Real t = -b - std::sqrt(disc);
            if (t > ray.maxDistance) return std::nullopt;
            Vec3 p = ray.origin + dir * t;
            return RayHit{&body, p, (p - body.position).normalized(), t};
        }
        case ShapeType::Box: {
            // Slab test in the box's local frame.
            Vec3 o = body.toLocal(ray.origin);
            Vec3 d = body.orientation.inverseRotate(dir);
            Vec3 e = shape.halfExtents;
            Real tMin = 0, tMax = ray.maxDistance;
            int hitAxis = -1;
            Real hitSign = 0;
            for (int i = 0; i < 3; ++i) {
                if (std::abs(d[i]) < kEpsilon) {
                    if (o[i] < -e[i] || o[i] > e[i]) return std::nullopt;
                    continue;
                }
                Real inv = 1 / d[i];
                Real t1 = (-e[i] - o[i]) * inv;
                Real t2 = (e[i] - o[i]) * inv;
                Real sign = -1;  // Entering through the -e face.
                if (t1 > t2) { std::swap(t1, t2); sign = 1; }
                if (t1 > tMin) { tMin = t1; hitAxis = i; hitSign = sign; }
                tMax = std::min(tMax, t2);
                if (tMin > tMax) return std::nullopt;
            }
            if (hitAxis < 0) return std::nullopt;  // Started inside.
            Vec3 localN;
            localN[hitAxis] = hitSign;
            return RayHit{&body, ray.origin + dir * tMin, body.orientation.rotate(localN), tMin};
        }
        case ShapeType::Plane: {
            Vec3 n = body.orientation.rotate({0, 1, 0});
            Real denom = dot(n, dir);
            Real dist = dot(ray.origin - body.position, n);
            if (dist <= 0 || denom >= 0) return std::nullopt;  // Behind, or not heading toward it.
            Real t = -dist / denom;
            if (t > ray.maxDistance) return std::nullopt;
            return RayHit{&body, ray.origin + dir * t, n, t};
        }
    }
    return std::nullopt;
}

}  // namespace phys
