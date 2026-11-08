#include "phys/broadphase.h"

#include <algorithm>

#include "phys/body.h"

namespace phys {

bool shouldCollide(const RigidBody& a, const RigidBody& b) {
    // At least one body must be dynamic for a contact to matter.
    if (!a.isDynamic() && !b.isDynamic()) return false;
    return (a.category & b.mask) != 0 && (b.category & a.mask) != 0;
}

void SweepAndPrune::findPairs(const std::vector<RigidBody*>& bodies, std::vector<BodyPair>& out) {
    out.clear();
    entries_.clear();
    entries_.reserve(bodies.size());
    for (RigidBody* b : bodies) entries_.push_back({b->computeAABB().expanded(margin_), b});

    // Sort by min X; ties broken by id so the order is deterministic.
    std::sort(entries_.begin(), entries_.end(), [](const Entry& a, const Entry& b) {
        if (a.box.min.x != b.box.min.x) return a.box.min.x < b.box.min.x;
        return a.body->id() < b.body->id();
    });

    for (size_t i = 0; i < entries_.size(); ++i) {
        const Entry& a = entries_[i];
        for (size_t j = i + 1; j < entries_.size(); ++j) {
            const Entry& b = entries_[j];
            if (b.box.min.x > a.box.max.x) break;  // Nothing further along X can overlap.
            if (!a.box.overlaps(b.box)) continue;
            if (!shouldCollide(*a.body, *b.body)) continue;
            if (a.body->id() < b.body->id())
                out.emplace_back(a.body, b.body);
            else
                out.emplace_back(b.body, a.body);
        }
    }

    std::sort(out.begin(), out.end(), [](const BodyPair& p, const BodyPair& q) {
        if (p.first->id() != q.first->id()) return p.first->id() < q.first->id();
        return p.second->id() < q.second->id();
    });
}

}  // namespace phys
