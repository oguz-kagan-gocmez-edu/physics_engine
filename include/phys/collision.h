#pragma once

#include "phys/math.h"

namespace phys {

class RigidBody;

constexpr int kMaxContacts = 4;

struct ContactPoint {
    Vec3 pointA;      // World point on A's surface.
    Vec3 pointB;      // World point on B's surface.
    Real depth = 0;   // Penetration depth, > 0 when overlapping.
    Vec3 localA;      // pointA in A's local frame.
    Vec3 localB;      // pointB in B's local frame.

    // Accumulated normal impulse, carried between steps for warm starting.
    Real normalImpulse = 0;

    // Solver scratch data, rebuilt every step.
    Vec3 rA, rB;
    Real normalMass = 0;
    Real velocityBias = 0;
};

// Contact between two bodies. The normal points from A to B.
//
// Friction is solved once per manifold at the center of the contact patch
// (two tangent directions plus twist around the normal) instead of at every
// point. This converges better for resting boxes and costs less.
struct Manifold {
    RigidBody* a = nullptr;
    RigidBody* b = nullptr;
    Vec3 normal;
    Vec3 tangent1, tangent2;
    Real friction = 0;
    Real restitution = 0;
    int count = 0;
    ContactPoint points[kMaxContacts];

    // Accumulated friction impulses, carried between steps.
    Real tangentImpulse1 = 0;
    Real tangentImpulse2 = 0;
    Real twistImpulse = 0;

    // Friction scratch data, rebuilt every step.
    Vec3 rA, rB;              // Patch center relative to each body.
    Real tangentMass1 = 0;
    Real tangentMass2 = 0;
    Real twistMass = 0;
    Real frictionRadius = 0;  // Average distance of points from the center.

    // Total normal impulse of the last step. Divide by dt for the contact force.
    Real totalNormalImpulse() const {
        Real sum = 0;
        for (int i = 0; i < count; ++i) sum += points[i].normalImpulse;
        return sum;
    }
};

// Runs the exact shape test for a pair. Returns true and fills `m` when the shapes
// overlap or are closer than `margin`. Points within the margin have negative depth;
// the solver treats them as speculative contacts, which keeps resting contacts from
// flickering on and off.
bool collide(RigidBody* a, RigidBody* b, Manifold& m, Real margin = 0);

}  // namespace phys
