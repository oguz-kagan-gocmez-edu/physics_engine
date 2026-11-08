#pragma once

#include "phys/math.h"

namespace phys {

struct AABB {
    Vec3 min, max;

    bool overlaps(const AABB& o) const {
        return min.x <= o.max.x && max.x >= o.min.x && min.y <= o.max.y && max.y >= o.min.y &&
               min.z <= o.max.z && max.z >= o.min.z;
    }

    AABB expanded(Real margin) const {
        Vec3 m{margin, margin, margin};
        return {min - m, max + m};
    }
};

}  // namespace phys
