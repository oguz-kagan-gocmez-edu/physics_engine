#include "phys/world.h"

#include <algorithm>
#include <cmath>

namespace phys {

namespace {
bool pairLess(const Manifold& m, uint32_t a, uint32_t b) {
    return m.a->id() < a || (m.a->id() == a && m.b->id() < b);
}
}  // namespace

World::World(const WorldSettings& s) : settings(s) {}

RigidBody* World::createBody(const BodyDesc& desc) {
    owned_.push_back(std::make_unique<RigidBody>(desc, nextId_++));
    RigidBody* body = owned_.back().get();
    bodies_.push_back(body);
    return body;
}

void World::removeBody(RigidBody* body) {
    auto touches = [body](const Manifold& m) { return m.a == body || m.b == body; };
    manifolds_.erase(std::remove_if(manifolds_.begin(), manifolds_.end(), touches), manifolds_.end());
    previous_.erase(std::remove_if(previous_.begin(), previous_.end(), touches), previous_.end());
    bodies_.erase(std::remove(bodies_.begin(), bodies_.end(), body), bodies_.end());
    owned_.erase(std::remove_if(owned_.begin(), owned_.end(),
                                [body](const std::unique_ptr<RigidBody>& p) { return p.get() == body; }),
                 owned_.end());
}

// Copies accumulated impulses from last step's matching contacts (warm starting).
// Both lists are sorted by body pair, so a single forward cursor finds matches.
void World::carryImpulses(Manifold& m, size_t& cursor) {
    uint32_t ia = m.a->id(), ib = m.b->id();
    while (cursor < previous_.size() && pairLess(previous_[cursor], ia, ib)) ++cursor;
    if (cursor >= previous_.size()) return;
    const Manifold& old = previous_[cursor];
    if (old.a != m.a || old.b != m.b) return;
    if (dot(old.normal, m.normal) < Real(0.95)) return;  // Contact changed too much.

    // Friction: carry the old impulse as a world vector and project it onto the new
    // tangents, since the tangent basis turns slightly when the normal does.
    Vec3 oldFriction = old.tangent1 * old.tangentImpulse1 + old.tangent2 * old.tangentImpulse2;
    orthonormalBasis(m.normal, m.tangent1, m.tangent2);
    m.tangentImpulse1 = dot(oldFriction, m.tangent1);
    m.tangentImpulse2 = dot(oldFriction, m.tangent2);
    m.twistImpulse = old.twistImpulse;

    Real tolSq = settings.contactMatchDistance * settings.contactMatchDistance;
    for (int i = 0; i < m.count; ++i) {
        ContactPoint& c = m.points[i];
        int bestIndex = -1;
        Real bestDist = tolSq;
        for (int j = 0; j < old.count; ++j) {
            Real d = (old.points[j].localA - c.localA).lengthSq();
            if (d < bestDist) { bestDist = d; bestIndex = j; }
        }
        if (bestIndex >= 0) {
            c.normalImpulse = old.points[bestIndex].normalImpulse;
            ++stats_.warmStartedContacts;
        }
    }
}

void World::step(Real dt) {
    if (dt <= 0) return;
    for (RigidBody* b : bodies_) {
        b->previousPosition = b->position;
        b->previousOrientation = b->orientation;
    }

    // Several short sub-steps converge much better than one long step with many
    // solver iterations, for about the same cost (tall stacks stay still).
    int n = std::max(1, settings.subSteps);
    Real h = dt / n;
    for (int i = 0; i < n; ++i) subStep(h);

    for (RigidBody* b : bodies_) b->clearForces();
    time_ += dt;
}

void World::subStep(Real dt) {
    // 1. Forces -> velocities.
    for (RigidBody* b : bodies_) {
        b->updateInertia();
        b->integrateVelocity(settings.gravity, dt);
    }

    // 2. Broad phase.
    broadPhase_.findPairs(bodies_, pairs_);

    // 3. Narrow phase, with impulses carried over from the last step.
    previous_.swap(manifolds_);
    manifolds_.clear();
    size_t cursor = 0;
    stats_.warmStartedContacts = 0;
    for (const BodyPair& pair : pairs_) {
        Manifold m;
        if (!collide(pair.first, pair.second, m, settings.contactMargin)) continue;
        carryImpulses(m, cursor);
        manifolds_.push_back(m);
    }

    // 4. Contact velocities.
    solver::prepare(manifolds_, settings.solver, dt);
    if (settings.solver.warmStarting) solver::warmStart(manifolds_);
    for (int i = 0; i < settings.solver.velocityIterations; ++i) solver::solveVelocities(manifolds_);

    // 5. Velocities -> positions, then remove leftover penetration.
    for (RigidBody* b : bodies_) b->integratePosition(dt);
    for (int i = 0; i < settings.solver.positionIterations; ++i)
        if (solver::solvePositions(manifolds_, settings.solver)) break;

    stats_.bodies = bodies_.size();
    stats_.broadPhasePairs = pairs_.size();
    stats_.manifolds = manifolds_.size();
    stats_.contacts = 0;
    for (const Manifold& m : manifolds_) stats_.contacts += m.count;
}

Real World::update(Real frameTime) {
    accumulator_ += frameTime;
    int steps = 0;
    while (accumulator_ >= settings.fixedDt && steps < settings.maxStepsPerUpdate) {
        step(settings.fixedDt);
        accumulator_ -= settings.fixedDt;
        ++steps;
    }
    // If we fell too far behind, drop the backlog instead of trying to catch up forever.
    if (accumulator_ >= settings.fixedDt) accumulator_ = std::fmod(accumulator_, settings.fixedDt);
    return accumulator_ / settings.fixedDt;
}

std::optional<RayHit> World::raycast(const Ray& ray, uint32_t mask) const {
    std::optional<RayHit> best;
    Ray r = ray;
    for (RigidBody* b : bodies_) {
        if ((b->category & mask) == 0) continue;
        if (auto hit = raycastBody(*b, r)) {
            best = hit;
            r.maxDistance = hit->distance;  // Only closer hits from now on.
        }
    }
    return best;
}

}  // namespace phys
