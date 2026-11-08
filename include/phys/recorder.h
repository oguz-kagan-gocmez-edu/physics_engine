#pragma once

#include <string>
#include <unordered_set>

#include "phys/world.h"

namespace phys {

// Records a simulation to a JSON file that viewer/index.html can play back.
//
//   Recorder rec;
//   while (...) { world.step(dt); rec.capture(world); }
//   rec.save("recording.json");
//
// Format (version 1):
//   { "version": 1, "dt": <seconds between frames>,
//     "bodies": [ { "id", "shape": "sphere"|"box"|"plane", "type", "size": [..] } ],
//     "frames": [ { "t": time,
//                   "p": [id, px, py, pz, qw, qx, qy, qz, id, ...],
//                   "c": [px, py, pz, nx, ny, nz, ...] } ] }
// A body that is not listed in a frame does not exist at that time.
class Recorder {
public:
    // contacts: also record contact points and normals for each frame.
    explicit Recorder(bool contacts = true) : recordContacts_(contacts) {}

    void capture(const World& world);
    bool save(const std::string& path) const;

    size_t frameCount() const { return frames_; }

private:
    bool recordContacts_;
    std::unordered_set<uint32_t> known_;
    std::string bodies_;  // JSON fragments, built as we go.
    std::string frameJson_;
    size_t frames_ = 0;
    Real firstTime_ = -1;
    Real lastTime_ = 0;
};

}  // namespace phys
