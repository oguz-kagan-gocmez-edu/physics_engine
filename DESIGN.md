# Physics Engine Design

A small, embeddable 3D rigid-body physics engine in C++17, meant to be the
simulation core of simulator programs (vehicles, robots, mechanisms, games,
teaching tools). No dependencies beyond the standard library.

## Goals

- **Correct and stable first, fast second.** Stacks of boxes must stand still,
  resting bodies must not jitter, and energy must not grow on its own.
- **Deterministic.** The same inputs on the same build give the same results,
  so simulations can be replayed and tested.
- **Fixed time step.** Physics always advances in fixed steps (default 1/60 s).
  The host program can render at any frame rate and interpolate.
- **Easy to embed.** Plain C++ classes, raw pointers to bodies that stay valid
  until removed, no global state, no rendering code inside the engine.
- **Readable.** Each stage of the pipeline lives in its own file so it can be
  understood, tested and replaced on its own.

Non-goals for v1: soft bodies, fluids, joints, GPU acceleration, multithreading.
These are on the roadmap below and the architecture leaves room for them.

## Units and conventions

- SI units: meters, kilograms, seconds, radians.
- Right-handed coordinates, **+Y is up**; default gravity is (0, -9.81, 0).
- `Real` is `double` by default (simulators often run with large worlds and
  long durations). Build with `-DPHYS_USE_FLOAT` to use `float`.
- Orientation is a unit quaternion. Angular velocity is in world space.
- Contact normals point **from body A to body B**.

## Architecture

```
            ┌───────────────────────────────────────────────┐
            │                    World                      │
            │  bodies · settings · contact cache · step()   │
            └───────────────────────────────────────────────┘
                                   │ step(dt)
   ┌──────────────┬────────────────┼────────────────┬──────────────────┐
   ▼              ▼                ▼                ▼                  ▼
1. Integrate   2. Broad phase   3. Narrow phase   4. Solve          5. Integrate
   forces         (sweep and       (shape-pair       contacts          positions,
   → velocity     prune on AABBs)  tests → contact   (sequential       then fix
                  → candidate      manifolds, warm   impulses with     penetration
                  pairs            start from cache) friction)         (position solver)
```

### Modules

| File | Responsibility |
|------|----------------|
| `math.h` | `Vec3`, `Quat`, `Mat3`, helpers. Header-only. |
| `shape.h` | Shape descriptions (sphere, box, plane) and mass properties. |
| `body.h/.cpp` | `RigidBody`: state, mass, forces, transforms. |
| `aabb.h` | Axis-aligned bounding boxes. |
| `broadphase.h/.cpp` | Finds pairs whose AABBs overlap. |
| `collision.h/.cpp` | Exact shape-vs-shape tests producing contact manifolds. |
| `solver.h/.cpp` | Contact constraint solver (velocity and position). |
| `raycast.h/.cpp` | Ray queries against bodies (sensors, picking, lidar). |
| `world.h/.cpp` | Owns bodies and runs the pipeline. |

### Rigid body

Each body stores:

- **State:** position, orientation, linear velocity, angular velocity.
- **Mass:** inverse mass and inverse inertia (diagonal, in the body's local
  frame). The world-space inverse inertia `R · I⁻¹ · Rᵀ` is recomputed once per
  step. Static bodies have zero inverse mass and inertia.
- **Material:** friction coefficient, restitution (bounciness), density.
- **Damping:** linear and angular, to model air drag and help settling.
- **Accumulators:** force and torque, cleared after each step.
- **Filtering:** `category` and `mask` bits decide which bodies can collide.

Body types:

- `Static`: never moves (ground, walls). Infinite mass.
- `Kinematic`: moved by the program through its velocity; pushes dynamic bodies
  but is not pushed back (elevators, conveyor parts, scripted platforms).
- `Dynamic`: fully simulated.

### Shapes (v1)

| Shape | Parameters | Notes |
|-------|------------|-------|
| Sphere | radius | Cheapest shape. |
| Box | half extents | Oriented box; most simulator props. |
| Plane | (local +Y normal through the body's position) | Static only; infinite ground or walls. |

Mass and inertia come from shape volume × density, using the standard solid
sphere and solid box formulas.

### Time stepping

`World::step(dt)` advances exactly one fixed step. Each step is split into
**sub-steps** (default 4), and the whole pipeline (collision detection
included) runs once per sub-step. In testing, 4 sub-steps × 4 solver iterations
left an 8-box stack 10× steadier than 1 step × 10 iterations, at the same CPU
cost. This is the approach used by modern engines (Box2D v3, Jolt, PhysX TGS).

`World::update(frameTime)`
keeps a time accumulator, runs as many fixed steps as needed (capped to avoid a
"spiral of death"), and returns an interpolation factor `alpha` in [0, 1) that
the renderer can use to blend between the previous and current state.

Integration is **semi-implicit (symplectic) Euler**: velocities are updated
from forces first, then positions are updated from the new velocities. It is
cheap, stable for stiff contacts, and keeps orbits and oscillators from
gaining energy the way explicit Euler does.

Orientation update: `q ← normalize(q + ½·dt·(0, ω)·q)`.

### Broad phase

**Sweep and prune** on the X axis: sort AABBs by their minimum X, sweep through
the sorted list, and test Y and Z overlap only for boxes that overlap in X.
Static-vs-static pairs and pairs rejected by the collision filter are skipped.
This is simple, deterministic, and fast for a few thousand bodies. The
broad phase sits behind a small interface so a dynamic AABB tree can replace it
for larger or very spread-out worlds.

### Narrow phase

Each shape pair has a dedicated test that returns a **contact manifold**: a
shared normal and up to 4 contact points, each with a penetration depth and a
point on each body's surface.

Contacts are kept while the shapes are still up to `contactMargin` (2 cm)
apart. These **speculative contacts** have negative depth: the solver lets the
bodies close the gap within the step but not pass it. Without the margin, a
box corner lifting by a hair drops its contact, the box tips, and stacks rock.

| Pair | Method |
|------|--------|
| sphere–sphere | Center distance. |
| sphere–plane | Signed distance of the center. |
| sphere–box | Clamp the sphere center into the box (special case when the center is inside). |
| box–plane | Test the 8 corners against the plane, keep up to 4 deepest. |
| box–box | **Separating axis test** on 15 axes (3 + 3 face normals, 9 edge cross products). Face contacts use **Sutherland–Hodgman clipping** of the incident face against the reference face; edge contacts use closest points between the two edges. Face axes are preferred over edge axes with a small tolerance so resting boxes produce stable face contacts. |

The pair dispatcher orders shapes so each test is written once
(sphere < box < plane) and flips the normal when it swaps the bodies.

### Contact persistence and warm starting

Manifolds are cached between steps by body-pair key. When a new contact is
within a small distance of an old one (compared in body-local coordinates), it
inherits the old contact's accumulated impulses. Starting the solver from last
step's answer ("warm starting") is what lets tall stacks settle in a few
iterations instead of slowly sinking. Friction impulses are carried as a world
vector and projected onto the new tangents, since the tangent directions turn
slightly when the normal does. `WorldStats::warmStartedContacts` reports how
many contacts were matched.

### Contact solver

**Sequential impulses** (projected Gauss–Seidel), as in Box2D and Bullet.

- **Normal (per contact point):** relative velocity along the normal must be
  ≥ 0 (or, for speculative contacts, may close at most the remaining gap).
  Accumulated impulse is clamped to ≥ 0 (contacts push, never pull).
- **Friction (per manifold, at the center of the contact patch):** two tangent
  directions limited together by the Coulomb cone `|λt| ≤ μ·Σλn`, plus a
  **twist** term that resists spinning around the normal, limited by
  `μ·Σλn·r` with `r` the average distance of the points from the center.
  Friction is `sqrt(μA·μB)`. Solving friction once per patch instead of at
  each of the 4 corners halved the residual jitter of a resting stack: four
  per-corner friction impulses can drift against each other without ever
  converging.
- **Restitution:** if the bodies approach faster than a threshold
  (default 1 m/s), the normal constraint targets a separating speed of
  `e · v_approach`, with `e = max(eA, eB)`. The threshold stops resting bodies
  from bouncing forever.

Velocity iterations (default 4 per sub-step) loop over all contacts.

**Position correction** runs after positions are integrated, as a few
iterations (default 2 per sub-step) of non-linear Gauss–Seidel: it recomputes each contact's
current separation from the bodies' new poses and moves the bodies directly
apart. Unlike Baumgarte velocity bias, this removes penetration **without
adding kinetic energy**, so stacks do not pop or jitter. A small allowed
penetration ("slop", default 5 mm) keeps contacts alive between steps.

### Queries

- `World::raycast(origin, direction, maxDistance)` returns the closest hit
  (body, point, normal, distance). Useful for lidar and distance sensors,
  wheel ground checks, and mouse picking.
- `World::manifolds()` exposes the current contacts so the program can read
  impact forces (accumulated normal impulse ÷ dt) or play sounds.

### Determinism

- Bodies get increasing integer IDs; pairs are always ordered (lower ID first)
  and processed in sorted order.
- No hash-map iteration order leaks into the solver order.
- No randomness, no threads, no time-dependent behavior.

## Default settings

| Setting | Default | Meaning |
|---------|---------|---------|
| `gravity` | (0, -9.81, 0) | m/s² |
| `fixedDt` | 1/60 s | Physics step length |
| `subSteps` | 4 | Sub-steps per step |
| `maxStepsPerUpdate` | 8 | Max steps per `update()` call |
| `velocityIterations` | 4 | Solver passes for velocities, per sub-step |
| `positionIterations` | 2 | Solver passes for penetration, per sub-step |
| `contactMargin` | 0.02 m | Distance at which speculative contacts start |
| `linearSlop` | 0.005 m | Allowed penetration |
| `baumgarte` | 0.2 | Fraction of penetration fixed per position pass |
| `maxCorrection` | 0.2 m | Max position fix per contact per pass |
| `restitutionThreshold` | 1.0 m/s | Below this approach speed, no bounce |

## Testing

`tests/tests.cpp` is a self-contained test runner (no framework needed):

- Math: quaternion rotation, matrix inverse and transpose.
- Free fall matches the exact discrete semi-implicit Euler result.
- Sphere and box come to rest on the ground at the right height.
- A stack of boxes stays standing for 10 simulated seconds.
- Elastic head-on collision swaps velocities and conserves momentum.
- A sliding box stops after roughly `v² / (2μg)`; a frictionless one does not.
- A rotated box lands on another box without exploding (energy stays bounded).
- Edge-on-edge box contact gives one point with the exact depth.
- Collision filter bits, kinematic bodies, body removal, fixed-step `update()`.
- Raycasts hit the right body at the right distance.
- Two identical worlds produce bit-identical results (determinism).

## Roadmap

1. **Sleeping and islands:** group touching bodies into islands and put
   islands to sleep when they stop moving. Big speed-up for large scenes, and
   removes the last ~2 cm/s of solver jitter at the top of tall stacks.
2. **More shapes:** capsule, cylinder, convex hull (GJK + EPA), triangle mesh
   for static terrain, heightfield.
3. **Joints:** ball-and-socket, hinge (with motor and limits), slider, fixed,
   distance/spring. They plug into the same sequential-impulse solver.
4. **Continuous collision detection** for fast, small bodies (bullets).
5. **Dynamic AABB tree** broad phase.
6. **Gyroscopic torque** for fast-spinning bodies (flywheels, tops).
7. **Vehicle model:** raycast wheels with suspension and tire friction.
8. **Multithreading:** parallel narrow phase and island solving.
