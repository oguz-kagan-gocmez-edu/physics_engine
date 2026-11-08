#include "phys/solver.h"

#include "phys/body.h"

namespace phys {
namespace solver {
namespace {

Real effectiveMass(const RigidBody& a, const RigidBody& b, const Vec3& rA, const Vec3& rB,
                   const Vec3& dir) {
    Vec3 ra = cross(rA, dir), rb = cross(rB, dir);
    Real k = a.invMass() + b.invMass() + dot(ra, a.invInertiaWorld() * ra) +
             dot(rb, b.invInertiaWorld() * rb);
    return k > kEpsilon ? 1 / k : 0;
}

void applyImpulse(RigidBody& a, RigidBody& b, const Vec3& rA, const Vec3& rB, const Vec3& P) {
    a.linearVelocity -= P * a.invMass();
    a.angularVelocity -= a.invInertiaWorld() * cross(rA, P);
    b.linearVelocity += P * b.invMass();
    b.angularVelocity += b.invInertiaWorld() * cross(rB, P);
}

void applyAngularImpulse(RigidBody& a, RigidBody& b, const Vec3& L) {
    a.angularVelocity -= a.invInertiaWorld() * L;
    b.angularVelocity += b.invInertiaWorld() * L;
}

Vec3 relativeVelocity(const RigidBody& a, const RigidBody& b, const Vec3& rA, const Vec3& rB) {
    return b.linearVelocity + cross(b.angularVelocity, rB) - a.linearVelocity -
           cross(a.angularVelocity, rA);
}

}  // namespace

void prepare(std::vector<Manifold>& manifolds, const SolverSettings& s, Real dt) {
    for (Manifold& m : manifolds) {
        RigidBody& a = *m.a;
        RigidBody& b = *m.b;
        m.friction = std::sqrt(a.friction * b.friction);
        m.restitution = std::max(a.restitution, b.restitution);
        orthonormalBasis(m.normal, m.tangent1, m.tangent2);

        Vec3 center;
        for (int i = 0; i < m.count; ++i) {
            ContactPoint& c = m.points[i];
            Vec3 p = (c.pointA + c.pointB) * Real(0.5);
            center += p;
            c.rA = p - a.position;
            c.rB = p - b.position;
            c.normalMass = effectiveMass(a, b, c.rA, c.rB, m.normal);

            // Speculative contacts (depth < 0) may close the gap this step but not cross it.
            Real vn = dot(relativeVelocity(a, b, c.rA, c.rB), m.normal);
            Real gapBias = c.depth < 0 ? c.depth / dt : 0;
            bool reachesThisStep = vn <= gapBias;
            if (reachesThisStep && vn < -s.restitutionThreshold && m.restitution > 0) {
                // Bounce target, measured before any impulse is applied this step.
                c.velocityBias = -m.restitution * vn;
            } else {
                c.velocityBias = gapBias;
            }
        }

        // Friction at the patch center.
        center = center / Real(m.count);
        m.rA = center - a.position;
        m.rB = center - b.position;
        m.tangentMass1 = effectiveMass(a, b, m.rA, m.rB, m.tangent1);
        m.tangentMass2 = effectiveMass(a, b, m.rA, m.rB, m.tangent2);
        Real kTwist = dot(m.normal, a.invInertiaWorld() * m.normal) +
                      dot(m.normal, b.invInertiaWorld() * m.normal);
        m.twistMass = kTwist > kEpsilon ? 1 / kTwist : 0;
        m.frictionRadius = 0;
        for (int i = 0; i < m.count; ++i) {
            Vec3 d = (m.points[i].pointA + m.points[i].pointB) * Real(0.5) - center;
            d -= m.normal * dot(d, m.normal);
            m.frictionRadius += d.length();
        }
        m.frictionRadius /= Real(m.count);

        if (!s.warmStarting) {
            for (int i = 0; i < m.count; ++i) m.points[i].normalImpulse = 0;
            m.tangentImpulse1 = m.tangentImpulse2 = m.twistImpulse = 0;
        }
    }
}

void warmStart(std::vector<Manifold>& manifolds) {
    for (Manifold& m : manifolds) {
        for (int i = 0; i < m.count; ++i) {
            const ContactPoint& c = m.points[i];
            applyImpulse(*m.a, *m.b, c.rA, c.rB, m.normal * c.normalImpulse);
        }
        applyImpulse(*m.a, *m.b, m.rA, m.rB,
                     m.tangent1 * m.tangentImpulse1 + m.tangent2 * m.tangentImpulse2);
        applyAngularImpulse(*m.a, *m.b, m.normal * m.twistImpulse);
    }
}

void solveVelocities(std::vector<Manifold>& manifolds) {
    for (Manifold& m : manifolds) {
        RigidBody& a = *m.a;
        RigidBody& b = *m.b;

        // Friction first, limited by the normal impulses from the previous pass.
        Real totalNormal = 0;
        for (int i = 0; i < m.count; ++i) totalNormal += m.points[i].normalImpulse;
        Real maxF = m.friction * totalNormal;

        Vec3 dv = relativeVelocity(a, b, m.rA, m.rB);
        Real old1 = m.tangentImpulse1, old2 = m.tangentImpulse2;
        Real t1 = old1 - dot(dv, m.tangent1) * m.tangentMass1;
        Real t2 = old2 - dot(dv, m.tangent2) * m.tangentMass2;
        // Coulomb cone: limit the combined tangent impulse, not each axis separately.
        Real len = std::sqrt(t1 * t1 + t2 * t2);
        if (len > maxF && len > 0) {
            t1 *= maxF / len;
            t2 *= maxF / len;
        }
        m.tangentImpulse1 = t1;
        m.tangentImpulse2 = t2;
        applyImpulse(a, b, m.rA, m.rB, m.tangent1 * (t1 - old1) + m.tangent2 * (t2 - old2));

        // Twist friction resists spinning around the normal.
        Real maxTwist = maxF * m.frictionRadius;
        Real wn = dot(b.angularVelocity - a.angularVelocity, m.normal);
        Real oldTwist = m.twistImpulse;
        m.twistImpulse = std::clamp(oldTwist - wn * m.twistMass, -maxTwist, maxTwist);
        applyAngularImpulse(a, b, m.normal * (m.twistImpulse - oldTwist));

        // Non-penetration: contacts can push but never pull.
        for (int i = 0; i < m.count; ++i) {
            ContactPoint& c = m.points[i];
            Real vn = dot(relativeVelocity(a, b, c.rA, c.rB), m.normal);
            Real lambda = -c.normalMass * (vn - c.velocityBias);
            Real old = c.normalImpulse;
            c.normalImpulse = std::max(old + lambda, Real(0));
            applyImpulse(a, b, c.rA, c.rB, m.normal * (c.normalImpulse - old));
        }
    }
}

bool solvePositions(std::vector<Manifold>& manifolds, const SolverSettings& s) {
    Real minSeparation = 0;
    for (Manifold& m : manifolds) {
        RigidBody& a = *m.a;
        RigidBody& b = *m.b;
        for (int i = 0; i < m.count; ++i) {
            const ContactPoint& c = m.points[i];
            // Current separation from the bodies' updated poses.
            Vec3 pA = a.toWorld(c.localA);
            Vec3 pB = b.toWorld(c.localB);
            Real separation = dot(pB - pA, m.normal);
            minSeparation = std::min(minSeparation, separation);

            Real C = std::clamp(s.baumgarte * (separation + s.linearSlop), -s.maxCorrection, Real(0));
            if (C >= 0) continue;

            Vec3 rA = pA - a.position;
            Vec3 rB = pB - b.position;
            Real mass = effectiveMass(a, b, rA, rB, m.normal);
            Vec3 P = m.normal * (-C * mass);

            // Move the bodies directly (no velocity change, so no energy is added).
            a.position -= P * a.invMass();
            a.orientation = a.orientation.integrated(-(a.invInertiaWorld() * cross(rA, P)), 1);
            b.position += P * b.invMass();
            b.orientation = b.orientation.integrated(b.invInertiaWorld() * cross(rB, P), 1);
        }
    }
    return minSeparation >= -3 * s.linearSlop;
}

}  // namespace solver
}  // namespace phys
