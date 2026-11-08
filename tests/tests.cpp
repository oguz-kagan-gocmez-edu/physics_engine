// Self-contained test runner. Build with `make test`.
#include <cstdio>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

#include "phys/phys.h"

using namespace phys;

namespace {

struct TestCase {
    const char* name;
    std::function<void()> fn;
};
std::vector<TestCase>& registry() {
    static std::vector<TestCase> tests;
    return tests;
}
struct Register {
    Register(const char* name, std::function<void()> fn) { registry().push_back({name, std::move(fn)}); }
};

int gFailures = 0;
bool gCurrentFailed = false;

#define TEST(name)                                   \
    static void name();                              \
    static Register reg_##name(#name, name);         \
    static void name()

#define CHECK(cond)                                                               \
    do {                                                                          \
        if (!(cond)) {                                                            \
            std::printf("    FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);       \
            gCurrentFailed = true;                                                \
        }                                                                         \
    } while (0)

// Single precision cannot meet the tightest tolerances, so they are relaxed there.
#ifdef PHYS_USE_FLOAT
constexpr double kMinRelTol = 1e-4;
#else
constexpr double kMinRelTol = 0;
#endif

#define CHECK_NEAR(a, b, tol)                                                          \
    do {                                                                               \
        double va_ = (a), vb_ = (b);                                                   \
        double tol_ = std::max(double(tol), kMinRelTol * (1 + std::abs(vb_)));         \
        if (!(std::abs(va_ - vb_) <= tol_)) {                                          \
            std::printf("    FAIL %s:%d: %s = %.6f, expected %.6f (tol %.6f)\n", __FILE__, \
                        __LINE__, #a, va_, vb_, tol_);                                 \
            gCurrentFailed = true;                                                     \
        }                                                                              \
    } while (0)

RigidBody* addGround(World& w, Real friction = Real(0.5)) {
    BodyDesc d;
    d.type = BodyType::Static;
    d.shape = Shape::plane();
    d.friction = friction;
    return w.createBody(d);
}

RigidBody* addBox(World& w, Vec3 pos, Vec3 half, Quat q = {}) {
    BodyDesc d;
    d.shape = Shape::box(half);
    d.position = pos;
    d.orientation = q;
    return w.createBody(d);
}

RigidBody* addSphere(World& w, Vec3 pos, Real r) {
    BodyDesc d;
    d.shape = Shape::sphere(r);
    d.position = pos;
    return w.createBody(d);
}

void run(World& w, Real seconds) {
    int steps = int(std::lround(seconds / w.settings.fixedDt));
    for (int i = 0; i < steps; ++i) w.step(w.settings.fixedDt);
}

}  // namespace

// ---------------------------------------------------------------------------
// Math

TEST(quaternion_rotates_vectors) {
    Quat q = Quat::fromAxisAngle({0, 0, 1}, kPi / 2);
    Vec3 v = q.rotate({1, 0, 0});
    CHECK_NEAR(v.x, 0, 1e-9);
    CHECK_NEAR(v.y, 1, 1e-9);
    Vec3 back = q.inverseRotate(v);
    CHECK_NEAR(back.x, 1, 1e-9);
    // Matrix form agrees with the quaternion.
    Vec3 m = q.toMat3() * Vec3{1, 2, 3};
    Vec3 r = q.rotate({1, 2, 3});
    CHECK_NEAR(m.x, r.x, 1e-9);
    CHECK_NEAR(m.y, r.y, 1e-9);
    CHECK_NEAR(m.z, r.z, 1e-9);
}

TEST(matrix_inverse) {
    Mat3 a{{2, 1, 0}, {0, 3, 1}, {1, 0, 4}};
    Mat3 p = a * a.inverse();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) CHECK_NEAR(p.r[i][j], i == j ? 1 : 0, 1e-9);
}

TEST(box_inertia) {
    Shape s = Shape::box({1, 0.5, 0.25});
    Vec3 i = s.inertiaDiagonal(12);
    // I_x = m/12 (h^2 + d^2) with full sizes 2, 1, 0.5.
    CHECK_NEAR(i.x, 1 * (1 + 0.25), 1e-9);
    CHECK_NEAR(i.y, 1 * (4 + 0.25), 1e-9);
    CHECK_NEAR(i.z, 1 * (4 + 1), 1e-9);
}

// ---------------------------------------------------------------------------
// Integration

TEST(free_fall_matches_semi_implicit_euler) {
    World w;
    RigidBody* s = addSphere(w, {0, 100, 0}, 0.5);
    s->linearDamping = 0;
    const int steps = 60;
    run(w, steps * w.settings.fixedDt);
    // Semi-implicit Euler with n sub-steps of length h: y_n = y0 - g h^2 n(n+1)/2.
    int n = steps * w.settings.subSteps;
    Real h = w.settings.fixedDt / w.settings.subSteps, g = 9.81;
    CHECK_NEAR(s->position.y, 100 - g * h * h * n * (n + 1) / 2, 1e-9);
    CHECK_NEAR(s->linearVelocity.y, -g * h * n, 1e-9);
}

TEST(torque_spins_body_up) {
    World w;
    w.settings.gravity = {};
    RigidBody* b = addBox(w, {0, 0, 0}, {0.5, 0.5, 0.5});
    b->angularDamping = 0;
    b->applyTorque({0, 1, 0});
    w.step(1);
    // I = m/12 * (1 + 1) with m = 1000 kg (1 m^3 of water).
    Real I = 1000.0 / 12 * 2;
    CHECK_NEAR(b->angularVelocity.y, 1 / I, 1e-9);
}

// ---------------------------------------------------------------------------
// Resting contact

TEST(sphere_rests_on_ground) {
    World w;
    addGround(w);
    RigidBody* s = addSphere(w, {0, 2, 0}, 0.5);
    run(w, 3);
    CHECK_NEAR(s->position.y, 0.5, 2 * w.settings.solver.linearSlop);
    CHECK(s->linearVelocity.length() < 0.01);
}

TEST(box_rests_flat_on_ground) {
    World w;
    addGround(w);
    RigidBody* b = addBox(w, {0, 1, 0}, {0.5, 0.25, 0.5});
    run(w, 3);
    CHECK_NEAR(b->position.y, 0.25, 2 * w.settings.solver.linearSlop);
    CHECK(b->linearVelocity.length() < 0.01);
    CHECK(b->angularVelocity.length() < 0.01);
    CHECK_NEAR(b->orientation.w, 1, 1e-3);
}

TEST(tilted_box_settles_flat) {
    World w;
    addGround(w);
    RigidBody* b = addBox(w, {0, 2, 0}, {0.5, 0.5, 0.5}, Quat::fromAxisAngle({1, 0.3, 0.2}, 0.7));
    run(w, 6);
    CHECK_NEAR(b->position.y, 0.5, 0.02);
    CHECK(b->linearVelocity.length() < 0.02);
}

TEST(box_stack_stays_standing) {
    World w;
    addGround(w);
    std::vector<RigidBody*> boxes;
    for (int i = 0; i < 8; ++i) boxes.push_back(addBox(w, {0, Real(0.5 + i * 1.0), 0}, {0.5, 0.5, 0.5}));
    run(w, 10);
    for (int i = 0; i < 8; ++i) {
        CHECK(std::abs(boxes[i]->position.x) < 0.05);
        CHECK(std::abs(boxes[i]->position.z) < 0.05);
        CHECK_NEAR(boxes[i]->position.y, 0.5 + i * 1.0, 0.1);
    }
    CHECK(boxes.back()->linearVelocity.length() < 0.05);
}

TEST(sphere_rests_on_box) {
    World w;
    addGround(w);
    addBox(w, {0, 0.5, 0}, {1, 0.5, 1});
    RigidBody* s = addSphere(w, {0.2, 3, 0.1}, 0.3);
    run(w, 4);
    CHECK_NEAR(s->position.y, 1.3, 0.02);
}

// ---------------------------------------------------------------------------
// Collisions

TEST(elastic_collision_swaps_velocities) {
    World w;
    w.settings.gravity = {};
    RigidBody* a = addSphere(w, {-2, 0, 0}, 0.5);
    RigidBody* b = addSphere(w, {2, 0, 0}, 0.5);
    for (RigidBody* s : {a, b}) {
        s->restitution = 1;
        s->linearDamping = 0;
        s->friction = 0;
    }
    a->linearVelocity = {3, 0, 0};
    Real p0 = a->mass() * a->linearVelocity.x + b->mass() * b->linearVelocity.x;
    run(w, 2);
    Real p1 = a->mass() * a->linearVelocity.x + b->mass() * b->linearVelocity.x;
    CHECK_NEAR(p1, p0, 1e-6 * std::abs(p0));
    CHECK_NEAR(a->linearVelocity.x, 0, 0.01);
    CHECK_NEAR(b->linearVelocity.x, 3, 0.01);
}

TEST(bouncing_ball_loses_height_with_restitution) {
    World w;
    addGround(w);
    RigidBody* s = addSphere(w, {0, 5, 0}, 0.5);
    s->restitution = 0.5;
    s->linearDamping = 0;
    Real maxAfterBounce = 0;
    bool bounced = false;
    for (int i = 0; i < 300; ++i) {
        w.step(w.settings.fixedDt);
        if (s->linearVelocity.y > 0) bounced = true;
        if (bounced) maxAfterBounce = std::max(maxAfterBounce, s->position.y);
    }
    CHECK(bounced);
    // Ideal: rebound height = e^2 * drop height = 0.25 * 4.5 m above the rest height.
    CHECK_NEAR(maxAfterBounce - 0.5, 0.25 * 4.5, 0.15);
}

TEST(friction_stops_sliding_box) {
    World w;
    addGround(w, 0.5);
    RigidBody* b = addBox(w, {0, 0.25, 0}, {0.5, 0.25, 0.5});
    b->friction = 0.5;
    b->linearDamping = 0;
    run(w, 0.5);  // Settle.
    b->linearVelocity = {4, 0, 0};
    Real x0 = b->position.x;
    run(w, 3);
    Real mu = 0.5, g = 9.81, v = 4;
    Real expected = v * v / (2 * mu * g);
    CHECK(std::abs(b->linearVelocity.x) < 0.01);
    CHECK_NEAR(b->position.x - x0, expected, 0.15 * expected);
}

TEST(frictionless_box_keeps_sliding) {
    World w;
    addGround(w, 0);
    RigidBody* b = addBox(w, {0, 0.25, 0}, {0.5, 0.25, 0.5});
    b->friction = 0;
    b->linearDamping = 0;
    run(w, 0.5);
    b->linearVelocity = {4, 0, 0};
    run(w, 2);
    CHECK_NEAR(b->linearVelocity.x, 4, 0.01);
}

TEST(rotated_box_lands_on_box_without_exploding) {
    World w;
    addGround(w);
    addBox(w, {0, 0.5, 0}, {0.5, 0.5, 0.5});
    RigidBody* top =
        addBox(w, {0.1, 2.5, 0}, {0.5, 0.5, 0.5}, Quat::fromAxisAngle({0.3, 1, 0.5}, 0.9));
    Real e0 = top->mass() * 9.81 * top->position.y;
    Real maxKE = 0;
    for (int i = 0; i < 600; ++i) {
        w.step(w.settings.fixedDt);
        maxKE = std::max(maxKE, top->kineticEnergy());
    }
    // Kinetic energy can never exceed the potential energy it started with.
    CHECK(maxKE < e0);
    CHECK(top->linearVelocity.length() < 0.05);
    // It either stays on the lower box or tumbles off onto the ground; both are valid,
    // but it must end up resting flat on one of them.
    bool onBox = std::abs(top->position.y - 1.5) < 0.02;
    bool onGround = std::abs(top->position.y - 0.5) < 0.02;
    CHECK(onBox || onGround);
}

TEST(edge_on_edge_contact_is_detected) {
    World w;
    w.settings.gravity = {};
    RigidBody* a = addBox(w, {0, 0, 0}, {0.5, 0.5, 0.5}, Quat::fromAxisAngle({1, 0, 0}, kPi / 4));
    RigidBody* b = addBox(w, {0, 1.4, 0}, {0.5, 0.5, 0.5}, Quat::fromAxisAngle({0, 0, 1}, kPi / 4));
    Manifold m;
    // Edge tips reach to ~0.707 each, so the boxes overlap by ~0.014.
    CHECK(collide(a, b, m));
    CHECK(m.count == 1);
    CHECK_NEAR(m.normal.y, 1, 1e-6);
    CHECK_NEAR(m.points[0].depth, 2 * std::sqrt(0.5) - 1.4, 1e-6);
}

TEST(filter_bits_prevent_collision) {
    World w;
    addGround(w);
    RigidBody* s = addSphere(w, {0, 2, 0}, 0.5);
    s->mask = 0;  // Collides with nothing.
    run(w, 2);
    CHECK(s->position.y < -5);
}

TEST(kinematic_body_pushes_dynamic_body) {
    World w;
    w.settings.gravity = {};
    BodyDesc d;
    d.type = BodyType::Kinematic;
    d.shape = Shape::box({0.5, 0.5, 0.5});
    d.position = {-2, 0, 0};
    d.linearVelocity = {2, 0, 0};
    RigidBody* pusher = w.createBody(d);
    RigidBody* s = addSphere(w, {0, 0, 0}, 0.5);
    run(w, 2);
    CHECK(s->position.x > pusher->position.x);
    CHECK_NEAR(pusher->linearVelocity.x, 2, 1e-12);
}

// ---------------------------------------------------------------------------
// Queries and world

TEST(raycast_hits_closest_body) {
    World w;
    addGround(w);
    RigidBody* box = addBox(w, {0, 0.5, 5}, {0.5, 0.5, 0.5});
    RigidBody* ball = addSphere(w, {0, 0.5, 10}, 0.5);
    auto hit = w.raycast({{0, 0.5, 0}, {0, 0, 1}, 100});
    CHECK(hit && hit->body == box);
    if (hit) {
        CHECK_NEAR(hit->distance, 4.5, 1e-9);
        CHECK_NEAR(hit->normal.z, -1, 1e-9);
    }
    auto down = w.raycast({{0, 3, 0}, {0, -1, 0}, 100});
    CHECK(down && down->body->shape().type == ShapeType::Plane);
    if (down) CHECK_NEAR(down->distance, 3, 1e-9);
    auto behind = w.raycast({{0, 0.5, 8}, {0, 0, 1}, 100});
    CHECK(behind && behind->body == ball);
    auto miss = w.raycast({{0, 3, 0}, {0, 1, 0}, 100});
    CHECK(!miss);
}

TEST(update_runs_fixed_steps_and_interpolates) {
    World w;
    Real alpha = w.update(2.5 * w.settings.fixedDt);
    CHECK_NEAR(w.time(), 2 * w.settings.fixedDt, 1e-12);
    CHECK_NEAR(alpha, 0.5, 1e-9);
}

TEST(remove_body_is_safe) {
    World w;
    addGround(w);
    RigidBody* a = addBox(w, {0, 0.5, 0}, {0.5, 0.5, 0.5});
    addBox(w, {0, 1.5, 0}, {0.5, 0.5, 0.5});
    run(w, 1);
    w.removeBody(a);
    run(w, 2);
    CHECK(w.bodies().size() == 2);
    CHECK_NEAR(w.bodies()[1]->position.y, 0.5, 0.02);
}

TEST(simulation_is_deterministic) {
    auto build = [](World& w) {
        addGround(w);
        for (int i = 0; i < 20; ++i)
            addBox(w, {Real((i % 4) * 0.3), Real(1 + i * 1.1), Real((i % 3) * 0.2)}, {0.4, 0.3, 0.5},
                   Quat::fromAxisAngle({1, 1, 0}, 0.1 * i));
    };
    World w1, w2;
    build(w1);
    build(w2);
    run(w1, 5);
    run(w2, 5);
    bool same = true;
    for (size_t i = 0; i < w1.bodies().size(); ++i) {
        const RigidBody* a = w1.bodies()[i];
        const RigidBody* b = w2.bodies()[i];
        same = same && std::memcmp(&a->position, &b->position, sizeof(Vec3)) == 0 &&
               std::memcmp(&a->orientation, &b->orientation, sizeof(Quat)) == 0;
    }
    CHECK(same);
}

int main(int argc, char** argv) {
    const char* filter = argc > 1 ? argv[1] : nullptr;
    int run_ = 0;
    for (const TestCase& t : registry()) {
        if (filter && !std::strstr(t.name, filter)) continue;
        gCurrentFailed = false;
        t.fn();
        ++run_;
        std::printf("%s %s\n", gCurrentFailed ? "[FAIL]" : "[ ok ]", t.name);
        if (gCurrentFailed) ++gFailures;
    }
    std::printf("\n%d/%d tests passed\n", run_ - gFailures, run_);
    return gFailures == 0 ? 0 : 1;
}
