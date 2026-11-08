#pragma once

#include <optional>

#include "phys/math.h"

namespace phys {

class RigidBody;

struct Ray {
    Vec3 origin;
    Vec3 direction;  // Does not need to be normalized.
    Real maxDistance = Real(1e30);
};

struct RayHit {
    RigidBody* body = nullptr;
    Vec3 point;
    Vec3 normal;     // Surface normal at the hit, pointing toward the ray origin side.
    Real distance = 0;
};

// Casts a ray against one body. Rays starting inside a shape do not hit it.
std::optional<RayHit> raycastBody(RigidBody& body, const Ray& ray);

}  // namespace phys
