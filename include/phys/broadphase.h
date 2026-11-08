#pragma once

#include <utility>
#include <vector>

#include "phys/aabb.h"

namespace phys {

class RigidBody;

using BodyPair = std::pair<RigidBody*, RigidBody*>;

// Finds pairs of bodies whose bounding boxes overlap.
class BroadPhase {
public:
    virtual ~BroadPhase() = default;
    // Pairs come back with the lower body id first, sorted by (idA, idB).
    virtual void findPairs(const std::vector<RigidBody*>& bodies, std::vector<BodyPair>& out) = 0;
};

// Sorts bounding boxes along X and sweeps through them.
class SweepAndPrune : public BroadPhase {
public:
    explicit SweepAndPrune(Real margin = Real(0.02)) : margin_(margin) {}
    void findPairs(const std::vector<RigidBody*>& bodies, std::vector<BodyPair>& out) override;

private:
    struct Entry {
        AABB box;
        RigidBody* body;
    };
    Real margin_;
    std::vector<Entry> entries_;
};

// True if the two bodies are allowed to collide (filter bits, body types).
bool shouldCollide(const RigidBody& a, const RigidBody& b);

}  // namespace phys
