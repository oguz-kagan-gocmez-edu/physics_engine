#pragma once

#include <vector>

#include "phys/collision.h"

namespace phys {

struct SolverSettings {
    int velocityIterations = 4;            // Per sub-step.
    int positionIterations = 2;            // Per sub-step.
    Real linearSlop = Real(0.005);         // Allowed penetration (m).
    Real baumgarte = Real(0.2);            // Fraction of penetration fixed per pass.
    Real maxCorrection = Real(0.2);        // Max position fix per contact per pass (m).
    Real restitutionThreshold = Real(1.0); // Below this approach speed, no bounce (m/s).
    bool warmStarting = true;
};

// Sequential-impulse contact solver. Call in this order every step:
// prepare -> warmStart -> solveVelocities (N times) -> [integrate positions] -> solvePositions (N times)
namespace solver {

void prepare(std::vector<Manifold>& manifolds, const SolverSettings& s, Real dt);
void warmStart(std::vector<Manifold>& manifolds);
void solveVelocities(std::vector<Manifold>& manifolds);
// Returns true when every contact is within the allowed penetration.
bool solvePositions(std::vector<Manifold>& manifolds, const SolverSettings& s);

}  // namespace solver
}  // namespace phys
