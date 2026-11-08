#pragma once

#include <memory>
#include <optional>
#include <vector>

#include "phys/body.h"
#include "phys/broadphase.h"
#include "phys/collision.h"
#include "phys/raycast.h"
#include "phys/solver.h"

namespace phys {

struct WorldSettings {
    Vec3 gravity{0, Real(-9.81), 0};
    Real fixedDt = Real(1) / 60;
    int maxStepsPerUpdate = 8;        // Max fixed steps per update() call.
    int subSteps = 4;                 // Each step is split into this many sub-steps.
    Real contactMatchDistance = Real(0.05);  // For carrying impulses between steps (m).
    Real contactMargin = Real(0.02);         // Keep contacts this close before touching (m).
    SolverSettings solver;
};

struct WorldStats {
    size_t bodies = 0;
    size_t broadPhasePairs = 0;
    size_t manifolds = 0;
    size_t contacts = 0;
    size_t warmStartedContacts = 0;  // Contacts that reused last step's impulses.
};

class World {
public:
    explicit World(const WorldSettings& settings = {});

    World(const World&) = delete;
    World& operator=(const World&) = delete;

    // The returned pointer stays valid until removeBody() or the world is destroyed.
    RigidBody* createBody(const BodyDesc& desc);
    void removeBody(RigidBody* body);
    const std::vector<RigidBody*>& bodies() const { return bodies_; }

    // Advances exactly one step of length dt.
    void step(Real dt);
    // Advances by frameTime using fixed steps. Returns the interpolation factor
    // in [0, 1) for rendering between the last two steps.
    Real update(Real frameTime);

    // Closest hit along the ray, among bodies whose category matches `mask`.
    std::optional<RayHit> raycast(const Ray& ray, uint32_t mask = 0xFFFFFFFF) const;

    const std::vector<Manifold>& manifolds() const { return manifolds_; }
    const WorldStats& stats() const { return stats_; }
    Real time() const { return time_; }

    WorldSettings settings;

private:
    void subStep(Real dt);
    void carryImpulses(Manifold& m, size_t& cursor);

    std::vector<std::unique_ptr<RigidBody>> owned_;
    std::vector<RigidBody*> bodies_;
    uint32_t nextId_ = 0;
    SweepAndPrune broadPhase_;
    std::vector<BodyPair> pairs_;
    std::vector<Manifold> manifolds_;
    std::vector<Manifold> previous_;
    WorldStats stats_;
    Real accumulator_ = 0;
    Real time_ = 0;
};

}  // namespace phys
