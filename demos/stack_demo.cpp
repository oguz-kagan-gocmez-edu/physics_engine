// Demo scenes for the engine. Prints a summary while it runs, and can save the
// run for the web viewer (viewer/index.html) or as CSV for plotting.
//
//   ./build/stack_demo [tower|dominoes|pyramid] [--record out.json] [--csv out.csv] [--seconds N]
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <string>

#include "phys/phys.h"

using namespace phys;

namespace {

void addGround(World& world) {
    BodyDesc ground;
    ground.type = BodyType::Static;
    ground.shape = Shape::plane();
    world.createBody(ground);
}

// Called once per step with the current time; lets a scene add things later.
using SceneEvents = std::function<void(World&)>;

// A tower of boxes and some bouncy balls, then a steel ball knocks it over.
SceneEvents buildTower(World& world) {
    addGround(world);
    for (int i = 0; i < 10; ++i) {
        BodyDesc box;
        box.shape = Shape::box({0.5, 0.25, 0.5});
        box.position = {0, Real(0.25 + i * 0.5), 0};
        box.orientation = Quat::fromAxisAngle({0, 1, 0}, Real(0.15 * i));
        world.createBody(box);
    }
    for (int i = 0; i < 5; ++i) {
        BodyDesc ball;
        ball.shape = Shape::sphere(0.3);
        ball.position = {Real(2.0 + i * 0.1), Real(3.0 + i * 1.5), Real(0.2 * i)};
        ball.restitution = 0.6;
        ball.density = 500;
        world.createBody(ball);
    }
    bool thrown = false;
    return [thrown](World& w) mutable {
        if (thrown || w.time() < 3) return;
        thrown = true;
        BodyDesc wrecker;
        wrecker.shape = Shape::sphere(0.5);
        wrecker.position = {-8, 1.5, 0};
        wrecker.density = 7800;  // Steel.
        wrecker.linearVelocity = {12, 0, 0};
        w.createBody(wrecker);
    };
}

// A curved line of dominoes, started by a small push.
SceneEvents buildDominoes(World& world) {
    addGround(world);
    const int count = 55;
    const Real radius = 6;
    RigidBody* first = nullptr;
    for (int i = 0; i < count; ++i) {
        Real a = Real(i) / count * Real(1.4) * kPi;  // Three quarters of a circle.
        BodyDesc d;
        // Thin along local Z, which the rotation below turns along the curve.
        d.shape = Shape::box({0.2, 0.4, 0.05});
        d.position = {radius * std::cos(a), 0.4, radius * std::sin(a)};
        // Local X points outward from the center, local Z along the curve.
        d.orientation = Quat::fromAxisAngle({0, 1, 0}, -a);
        d.friction = 0.4;
        d.density = 600;  // Wood.
        RigidBody* b = world.createBody(d);
        if (!first) first = b;
    }
    bool pushed = false;
    return [first, pushed](World& w) mutable {
        if (pushed || w.time() < 0.5) return;
        pushed = true;
        // Tip the first domino toward the second one (+Z at angle 0).
        first->applyImpulseAtPoint({0, 0, Real(0.6) * first->mass()},
                                   first->position + Vec3{0, 0.35, 0});
    };
}

// A pyramid of boxes, then a few boxes rain down on it.
SceneEvents buildPyramid(World& world) {
    addGround(world);
    const int base = 7;
    for (int row = 0; row < base; ++row) {
        int n = base - row;
        for (int i = 0; i < n; ++i) {
            BodyDesc d;
            d.shape = Shape::box({0.4, 0.4, 0.4});
            d.position = {Real((i - (n - 1) * 0.5) * 0.85), Real(0.4 + row * 0.8), 0};
            world.createBody(d);
        }
    }
    int dropped = 0;
    return [dropped](World& w) mutable {
        if (dropped >= 6 || w.time() < 2 + dropped * 0.6) return;
        BodyDesc d;
        d.shape = dropped % 2 ? Shape::sphere(0.35) : Shape::box({0.3, 0.3, 0.3});
        d.position = {Real(-1.5 + dropped * 0.6), 9, Real(0.3 * (dropped % 3) - 0.3)};
        d.orientation = Quat::fromAxisAngle({1, 0.5, 0.2}, Real(0.7 * dropped));
        d.density = 2500;
        w.createBody(d);
        ++dropped;
    };
}

void usage() {
    std::printf("usage: stack_demo [tower|dominoes|pyramid] [--record out.json] [--csv out.csv] "
                "[--seconds N]\n");
}

}  // namespace

int main(int argc, char** argv) {
    std::string scene = "tower";
    const char* recordPath = nullptr;
    const char* csvPath = nullptr;
    Real duration = 0;  // 0 = scene default.
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--record") && i + 1 < argc) recordPath = argv[++i];
        else if (!std::strcmp(argv[i], "--csv") && i + 1 < argc) csvPath = argv[++i];
        else if (!std::strcmp(argv[i], "--seconds") && i + 1 < argc) duration = Real(std::atof(argv[++i]));
        else if (argv[i][0] != '-') scene = argv[i];
        else { usage(); return 1; }
    }

    World world;
    SceneEvents events;
    if (scene == "tower") events = buildTower(world);
    else if (scene == "dominoes") events = buildDominoes(world);
    else if (scene == "pyramid") events = buildPyramid(world);
    else { usage(); return 1; }
    if (duration <= 0) duration = scene == "dominoes" ? 11 : 8;  // The wave needs ~10 s.

    Recorder recorder;
    std::FILE* csv = csvPath ? std::fopen(csvPath, "w") : nullptr;
    if (csv) std::fprintf(csv, "time,id,shape,px,py,pz,qw,qx,qy,qz\n");

    const Real dt = world.settings.fixedDt;
    int frame = 0;
    if (recordPath) recorder.capture(world);
    std::printf("scene: %s\n%6s %7s %9s %12s\n", scene.c_str(), "time", "bodies", "contacts",
                "energy (J)");
    while (world.time() < duration - dt * Real(0.5)) {
        events(world);
        world.step(dt);
        if (recordPath) recorder.capture(world);

        if (csv) {
            for (const RigidBody* b : world.bodies()) {
                std::fprintf(csv, "%.4f,%u,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", world.time(),
                             b->id(), int(b->shape().type), b->position.x, b->position.y,
                             b->position.z, b->orientation.w, b->orientation.x, b->orientation.y,
                             b->orientation.z);
            }
        }

        if (++frame % 60 == 0) {
            Real energy = 0;
            for (const RigidBody* b : world.bodies()) energy += b->kineticEnergy();
            std::printf("%6.2f %7zu %9zu %12.2f\n", world.time(), world.bodies().size(),
                        world.stats().contacts, energy);
        }
    }

    if (csv) {
        std::fclose(csv);
        std::printf("Wrote CSV to %s\n", csvPath);
    }
    if (recordPath) {
        if (!recorder.save(recordPath)) {
            std::fprintf(stderr, "Could not write %s\n", recordPath);
            return 1;
        }
        std::printf("Wrote %zu frames to %s\n", recorder.frameCount(), recordPath);
    }
    return 0;
}
