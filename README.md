# phys — a 3D rigid-body physics engine

A small C++17 physics engine for simulator programs. No dependencies.

- Spheres, boxes and infinite planes
- Static, kinematic and dynamic bodies
- Exact collision detection (box–box uses the separating axis test with face clipping)
- Sequential-impulse solver with warm starting, sub-stepping, Coulomb and twist friction, restitution
- Fixed time step with render interpolation
- Raycasts, collision filtering, deterministic results
- Recording to JSON and a 3D web viewer for playback
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

## 3D viewer

`viewer/index.html` plays back recorded simulations in the browser (three.js):
orbit camera, timeline scrubbing, slow motion, contact points and normals, and
click-to-inspect for any body (position, rotation, velocity).

```bash
make view
```

This records the three demo scenes (`tower`, `dominoes`, `pyramid`) into
`build/` and serves the project at http://localhost:8765/viewer/. Pick a scene
from the menu, or drag any recording onto the page.

Keys: **Space** play/pause · **← →** step one frame (Shift: 10) · **R** restart ·
**C** contacts · **Esc** deselect.

Record your own simulation with `Recorder`:

```cpp
Recorder rec;                        // Pass false to skip contact points.
for (int i = 0; i < 600; ++i) {
    world.step(world.settings.fixedDt);
    rec.capture(world);
}
rec.save("build/my_sim.json");
```

Then open http://localhost:8765/viewer/?file=../build/my_sim.json.

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
demos/          stack_demo: tower, dominoes and pyramid scenes; JSON or CSV output
viewer/         index.html: 3D playback of recordings
```
