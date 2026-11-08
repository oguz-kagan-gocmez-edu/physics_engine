#include "phys/recorder.h"

#include <cstdio>

namespace phys {
namespace {

void appendNumber(std::string& out, double v) {
    char buf[32];
    // 5 decimals is 0.01 mm: plenty for viewing, and keeps files small.
    int n = std::snprintf(buf, sizeof(buf), "%.5f", v);
    // Trim trailing zeros ("1.50000" -> "1.5", "2.00000" -> "2").
    while (n > 1 && buf[n - 1] == '0') --n;
    if (n > 1 && buf[n - 1] == '.') --n;
    if (n == 2 && buf[0] == '-' && buf[1] == '0') { buf[0] = '0'; n = 1; }
    out.append(buf, n);
}

const char* shapeName(ShapeType t) {
    switch (t) {
        case ShapeType::Sphere: return "sphere";
        case ShapeType::Box: return "box";
        case ShapeType::Plane: return "plane";
    }
    return "unknown";
}

const char* typeName(BodyType t) {
    switch (t) {
        case BodyType::Static: return "static";
        case BodyType::Kinematic: return "kinematic";
        case BodyType::Dynamic: return "dynamic";
    }
    return "unknown";
}

}  // namespace

void Recorder::capture(const World& world) {
    if (firstTime_ < 0) firstTime_ = world.time();
    lastTime_ = world.time();

    // Describe bodies the first time they show up.
    for (const RigidBody* b : world.bodies()) {
        if (!known_.insert(b->id()).second) continue;
        if (!bodies_.empty()) bodies_ += ",\n";
        bodies_ += "{\"id\":" + std::to_string(b->id()) + ",\"shape\":\"" + shapeName(b->shape().type) +
                   "\",\"type\":\"" + typeName(b->type()) + "\",\"size\":[";
        const Shape& s = b->shape();
        if (s.type == ShapeType::Sphere) {
            appendNumber(bodies_, s.radius);
        } else if (s.type == ShapeType::Box) {
            appendNumber(bodies_, s.halfExtents.x);
            bodies_ += ',';
            appendNumber(bodies_, s.halfExtents.y);
            bodies_ += ',';
            appendNumber(bodies_, s.halfExtents.z);
        }
        bodies_ += "]}";
    }

    if (frames_ > 0) frameJson_ += ",\n";
    frameJson_ += "{\"t\":";
    appendNumber(frameJson_, world.time());
    frameJson_ += ",\"p\":[";
    bool first = true;
    for (const RigidBody* b : world.bodies()) {
        if (!first) frameJson_ += ',';
        first = false;
        frameJson_ += std::to_string(b->id());
        const double v[7] = {b->position.x,    b->position.y,    b->position.z,   b->orientation.w,
                             b->orientation.x, b->orientation.y, b->orientation.z};
        for (double x : v) {
            frameJson_ += ',';
            appendNumber(frameJson_, x);
        }
    }
    frameJson_ += ']';

    if (recordContacts_) {
        frameJson_ += ",\"c\":[";
        first = true;
        for (const Manifold& m : world.manifolds()) {
            for (int i = 0; i < m.count; ++i) {
                const ContactPoint& c = m.points[i];
                if (c.depth < 0) continue;  // Skip speculative (not yet touching) points.
                Vec3 p = (c.pointA + c.pointB) * Real(0.5);
                const double v[6] = {p.x, p.y, p.z, m.normal.x, m.normal.y, m.normal.z};
                for (double x : v) {
                    if (!first) frameJson_ += ',';
                    first = false;
                    appendNumber(frameJson_, x);
                }
            }
        }
        frameJson_ += ']';
    }
    frameJson_ += '}';
    ++frames_;
}

bool Recorder::save(const std::string& path) const {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    double dt = frames_ > 1 ? (lastTime_ - firstTime_) / double(frames_ - 1) : 0;
    std::fprintf(f, "{\"version\":1,\"dt\":%.6f,\n\"bodies\":[\n%s],\n\"frames\":[\n%s]}\n", dt,
                 bodies_.c_str(), frameJson_.c_str());
    return std::fclose(f) == 0;
}

}  // namespace phys
