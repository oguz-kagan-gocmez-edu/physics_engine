// Drops a box tower and a few balls onto the ground and prints what happens.
// Pass a file name to also write every body's pose per frame as CSV:
//   ./build/stack_demo frames.csv
#include <cstdio>

#include "phys/phys.h"

using namespace phys;

int main(int argc, char** argv) {
    World world;

    BodyDesc ground;
    ground.type = BodyType::Static;
    ground.shape = Shape::plane();
    world.createBody(ground);

    // A tower of boxes, each turned a little.
    for (int i = 0; i < 10; ++i) {
        BodyDesc box;
        box.shape = Shape::box({0.5, 0.25, 0.5});
        box.position = {0, 0.25 + i * 0.5, 0};
        box.orientation = Quat::fromAxisAngle({0, 1, 0}, 0.15 * i);
        world.createBody(box);
    }

    // Bouncy balls dropped next to the tower.
    for (int i = 0; i < 5; ++i) {
        BodyDesc ball;
        ball.shape = Shape::sphere(0.3);
        ball.position = {2.0 + i * 0.1, 3.0 + i * 1.5, 0.2 * i};
        ball.restitution = 0.6;
        ball.density = 500;
        world.createBody(ball);
    }

    // A heavy ball thrown at the tower after 3 seconds.
    BodyDesc wrecker;
    wrecker.shape = Shape::sphere(0.5);
    wrecker.position = {-8, 1.5, 0};
    wrecker.density = 7800;  // Steel.
    RigidBody* ball = nullptr;

    std::FILE* csv = argc > 1 ? std::fopen(argv[1], "w") : nullptr;
    if (csv) std::fprintf(csv, "time,id,shape,px,py,pz,qw,qx,qy,qz\n");

    const Real duration = 8;
    const Real dt = world.settings.fixedDt;
    int frame = 0;
    std::printf("%6s %8s %10s %10s %12s\n", "time", "contacts", "top box y", "ball x", "energy (J)");
    while (world.time() < duration) {
        if (!ball && world.time() >= 3) {
            ball = world.createBody(wrecker);
            ball->linearVelocity = {12, 0, 0};
        }
        world.step(dt);

        if (csv) {
            for (const RigidBody* b : world.bodies()) {
                std::fprintf(csv, "%.4f,%u,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", world.time(),
                             b->id(), int(b->shape().type), b->position.x, b->position.y,
                             b->position.z, b->orientation.w, b->orientation.x, b->orientation.y,
                             b->orientation.z);
            }
        }

        if (++frame % 30 == 0) {
            Real energy = 0;
            for (const RigidBody* b : world.bodies()) energy += b->kineticEnergy();
            std::printf("%6.2f %8zu %10.3f %10.3f %12.2f\n", world.time(), world.stats().contacts,
                        world.bodies()[10]->position.y, ball ? ball->position.x : 0.0, energy);
        }
    }
    if (csv) {
        std::fclose(csv);
        std::printf("\nWrote frames to %s\n", argv[1]);
    }
    return 0;
}
