#include "phys/body.h"

namespace phys {

RigidBody::RigidBody(const BodyDesc& d, uint32_t id)
    : position(d.position),
      orientation(d.orientation.normalized()),
      linearVelocity(d.linearVelocity),
      angularVelocity(d.angularVelocity),
      friction(d.friction),
      restitution(d.restitution),
      linearDamping(d.linearDamping),
      angularDamping(d.angularDamping),
      gravityScale(d.gravityScale),
      category(d.category),
      mask(d.mask),
      userData(d.userData),
      id_(id),
      type_(d.type),
      shape_(d.shape) {
    previousPosition = position;
    previousOrientation = orientation;
    // Planes are infinite, so they can only be static.
    if (shape_.type == ShapeType::Plane) type_ = BodyType::Static;

    if (type_ == BodyType::Dynamic) {
        setMass(d.mass > 0 ? d.mass : d.density * shape_.volume());
    } else {
        linearVelocity = type_ == BodyType::Static ? Vec3{} : linearVelocity;
        angularVelocity = type_ == BodyType::Static ? Vec3{} : angularVelocity;
    }
    updateInertia();
}

void RigidBody::setMass(Real mass) {
    if (type_ != BodyType::Dynamic || mass <= 0) {
        invMass_ = 0;
        invInertiaLocal_ = {};
        return;
    }
    invMass_ = 1 / mass;
    Vec3 inertia = shape_.inertiaDiagonal(mass);
    invInertiaLocal_ = {inertia.x > 0 ? 1 / inertia.x : 0, inertia.y > 0 ? 1 / inertia.y : 0,
                        inertia.z > 0 ? 1 / inertia.z : 0};
    updateInertia();
}

void RigidBody::applyForceAtPoint(const Vec3& force, const Vec3& worldPoint) {
    force_ += force;
    torque_ += cross(worldPoint - position, force);
}

void RigidBody::applyImpulseAtPoint(const Vec3& impulse, const Vec3& worldPoint) {
    linearVelocity += impulse * invMass_;
    angularVelocity += invInertiaWorld_ * cross(worldPoint - position, impulse);
}

Real RigidBody::kineticEnergy() const {
    if (invMass_ == 0) return 0;
    Real linear = Real(0.5) * mass() * linearVelocity.lengthSq();
    // Angular energy: 0.5 * w^T I w, computed in the local frame where I is diagonal.
    Vec3 wl = orientation.inverseRotate(angularVelocity);
    Real angular = 0;
    for (int i = 0; i < 3; ++i)
        if (invInertiaLocal_[i] > 0) angular += wl[i] * wl[i] / invInertiaLocal_[i];
    return linear + Real(0.5) * angular;
}

AABB RigidBody::computeAABB() const {
    switch (shape_.type) {
        case ShapeType::Sphere: {
            Vec3 r{shape_.radius, shape_.radius, shape_.radius};
            return {position - r, position + r};
        }
        case ShapeType::Box: {
            // Extent along each world axis = |R| * halfExtents.
            Mat3 R = orientation.toMat3();
            Vec3 e = shape_.halfExtents;
            Vec3 ext{dot(vabs(R.r[0]), e), dot(vabs(R.r[1]), e), dot(vabs(R.r[2]), e)};
            return {position - ext, position + ext};
        }
        case ShapeType::Plane: {
            const Real big = Real(1e30);
            return {{-big, -big, -big}, {big, big, big}};
        }
    }
    return {position, position};
}

void RigidBody::updateInertia() {
    if (invMass_ == 0) {
        invInertiaWorld_ = Mat3::zero();
        return;
    }
    Mat3 R = orientation.toMat3();
    invInertiaWorld_ = R * Mat3::diagonal(invInertiaLocal_) * R.transposed();
}

void RigidBody::integrateVelocity(const Vec3& gravity, Real dt) {
    if (type_ != BodyType::Dynamic) return;
    linearVelocity += (gravity * gravityScale + force_ * invMass_) * dt;
    angularVelocity += (invInertiaWorld_ * torque_) * dt;
    // Implicit damping: stable for any damping value.
    linearVelocity *= 1 / (1 + dt * linearDamping);
    angularVelocity *= 1 / (1 + dt * angularDamping);
}

void RigidBody::integratePosition(Real dt) {
    if (type_ == BodyType::Static) return;
    position += linearVelocity * dt;
    orientation = orientation.integrated(angularVelocity, dt);
}

}  // namespace phys
