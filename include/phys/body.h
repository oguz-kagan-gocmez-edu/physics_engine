#pragma once

#include <cstdint>

#include "phys/aabb.h"
#include "phys/math.h"
#include "phys/shape.h"

namespace phys {

enum class BodyType { Static, Kinematic, Dynamic };

// Everything needed to create a body. Fill in what you need, keep the rest.
struct BodyDesc {
    BodyType type = BodyType::Dynamic;
    Shape shape;
    Vec3 position;
    Quat orientation;
    Vec3 linearVelocity;
    Vec3 angularVelocity;
    Real density = 1000;        // kg/m^3 (water). Ignored when mass > 0.
    Real mass = 0;              // If > 0, overrides density.
    Real friction = Real(0.5);
    Real restitution = 0;
    Real linearDamping = Real(0.01);
    Real angularDamping = Real(0.05);
    Real gravityScale = 1;
    uint32_t category = 1;      // Which group this body is in.
    uint32_t mask = 0xFFFFFFFF; // Which groups it collides with.
    void* userData = nullptr;
};

class RigidBody {
public:
    explicit RigidBody(const BodyDesc& desc, uint32_t id);

    uint32_t id() const { return id_; }
    BodyType type() const { return type_; }
    const Shape& shape() const { return shape_; }
    bool isStatic() const { return type_ == BodyType::Static; }
    bool isDynamic() const { return type_ == BodyType::Dynamic; }

    // State
    Vec3 position;
    Quat orientation;
    Vec3 linearVelocity;
    Vec3 angularVelocity;  // World space, rad/s.

    // Material and behavior
    Real friction;
    Real restitution;
    Real linearDamping;
    Real angularDamping;
    Real gravityScale;
    uint32_t category;
    uint32_t mask;
    void* userData;

    Real mass() const { return invMass_ > 0 ? 1 / invMass_ : 0; }
    Real invMass() const { return invMass_; }
    const Mat3& invInertiaWorld() const { return invInertiaWorld_; }
    void setMass(Real mass);

    // Forces are cleared after every step.
    void applyForce(const Vec3& force) { force_ += force; }
    void applyForceAtPoint(const Vec3& force, const Vec3& worldPoint);
    void applyTorque(const Vec3& torque) { torque_ += torque; }
    // Impulses change velocity right away.
    void applyImpulse(const Vec3& impulse) { linearVelocity += impulse * invMass_; }
    void applyImpulseAtPoint(const Vec3& impulse, const Vec3& worldPoint);

    Vec3 toWorld(const Vec3& localPoint) const { return position + orientation.rotate(localPoint); }
    Vec3 toLocal(const Vec3& worldPoint) const { return orientation.inverseRotate(worldPoint - position); }
    Vec3 velocityAt(const Vec3& worldPoint) const {
        return linearVelocity + cross(angularVelocity, worldPoint - position);
    }

    // Pose at the start of the last step, for smooth rendering between steps.
    Vec3 previousPosition;
    Quat previousOrientation;
    Vec3 interpolatedPosition(Real alpha) const {
        return previousPosition + (position - previousPosition) * alpha;
    }
    Quat interpolatedOrientation(Real alpha) const {
        return Quat::nlerp(previousOrientation, orientation, alpha);
    }

    Real kineticEnergy() const;
    AABB computeAABB() const;

    // Used by the world each step.
    void integrateVelocity(const Vec3& gravity, Real dt);
    void integratePosition(Real dt);
    void updateInertia();
    void clearForces() { force_ = {}; torque_ = {}; }

private:
    uint32_t id_;
    BodyType type_;
    Shape shape_;
    Real invMass_ = 0;
    Vec3 invInertiaLocal_;
    Mat3 invInertiaWorld_ = Mat3::zero();
    Vec3 force_;
    Vec3 torque_;
};

}  // namespace phys
