# phys — a 3D rigid-body physics engine

A small C++17 physics engine for simulator programs. No dependencies.

- Spheres, boxes and infinite planes
- Static, kinematic and dynamic bodies
- Exact collision detection (box–box uses the separating axis test with face clipping)
- Sequential-impulse solver with warm starting, sub-stepping, Coulomb and twist friction, restitution
- Fixed time step with render interpolation
- Raycasts, collision filtering, deterministic results
- `double` precision by default, `float` with `-DPHYS_USE_FLOAT`

See [DESIGN.md](DESIGN.md) for how it works and the roadmap.

## Build

With Make:

```bash
make test
```

```bash
make demo
```

With CMake:

```bash
cmake -B build-cmake && cmake --build build-cmake && ctest --test-dir build-cmake
```

## Usage

```cpp
#include "phys/phys.h"
using namespace phys;

World world;  // Gravity -9.81 on Y, 60 Hz steps.

BodyDesc ground;
ground.type = BodyType::Static;
ground.shape = Shape::plane();          // Infinite plane, normal +Y.
world.createBody(ground);

BodyDesc crate;
crate.shape = Shape::box({0.5, 0.5, 0.5});
crate.position = {0, 5, 0};
crate.density = 600;                    // kg/m^3 (or set crate.mass).
crate.friction = 0.6;
RigidBody* box = world.createBody(crate);

// In your program's loop:
Real alpha = world.update(frameSeconds);   // Runs as many fixed steps as needed.
Vec3 drawPos = box->interpolatedPosition(alpha);
Quat drawRot = box->interpolatedOrientation(alpha);

// Push things around.
box->applyForce({0, 0, 100});                    // Newtons, for the next step.
box->applyImpulseAtPoint({1, 0, 0}, box->position + Vec3{0, 0.5, 0});

// Sensors.
if (auto hit = world.raycast({{0, 10, 0}, {0, -1, 0}, 50})) {
    // hit->body, hit->point, hit->normal, hit->distance
}

// Contact forces.
for (const Manifold& m : world.manifolds()) {
    Real force = m.totalNormalImpulse() / world.settings.fixedDt;
}
```

## Layout

```
include/phys/   public headers (phys.h includes everything)
src/            engine source
tests/          self-contained test runner
demos/          stack_demo: tower + wrecking ball, optional CSV output
```
