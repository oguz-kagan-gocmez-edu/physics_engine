#include "phys/collision.h"

#include <limits>

#include "phys/body.h"

namespace phys {
namespace {

const Vec3 kUp{0, 1, 0};

// Raw contact before it is stored in a manifold.
struct RawContact {
    Vec3 pointA, pointB;
    Real depth;
};

struct RawManifold {
    Vec3 normal;  // A -> B
    int count = 0;
    RawContact points[8];

    void add(const Vec3& pa, const Vec3& pb) {
        if (count >= 8) return;
        points[count++] = {pa, pb, dot(pa - pb, normal)};
    }
};

Vec3 planeNormal(const RigidBody& plane) { return plane.orientation.rotate(kUp); }

// ---------------------------------------------------------------------------
// Sphere vs sphere
bool sphereSphere(const RigidBody& A, const RigidBody& B, Real margin, RawManifold& m) {
    Real ra = A.shape().radius, rb = B.shape().radius;
    Vec3 d = B.position - A.position;
    Real distSq = d.lengthSq();
    Real r = ra + rb + margin;
    if (distSq > r * r) return false;
    Real dist = std::sqrt(distSq);
    m.normal = dist > kEpsilon ? d / dist : kUp;
    m.add(A.position + m.normal * ra, B.position - m.normal * rb);
    return true;
}

// Sphere vs plane
bool spherePlane(const RigidBody& A, const RigidBody& B, Real margin, RawManifold& m) {
    Vec3 np = planeNormal(B);
    Real r = A.shape().radius;
    Real dist = dot(A.position - B.position, np);
    if (dist > r + margin) return false;
    m.normal = -np;
    m.add(A.position - np * r, A.position - np * dist);
    return true;
}

// Sphere vs box
bool sphereBox(const RigidBody& A, const RigidBody& B, Real margin, RawManifold& m) {
    Real r = A.shape().radius;
    Vec3 e = B.shape().halfExtents;
    Vec3 c = B.toLocal(A.position);
    Vec3 clamped{std::clamp(c.x, -e.x, e.x), std::clamp(c.y, -e.y, e.y), std::clamp(c.z, -e.z, e.z)};
    Vec3 diff = c - clamped;
    Real distSq = diff.lengthSq();

    if (distSq > kEpsilon) {
        // Center is outside the box.
        if (distSq > (r + margin) * (r + margin)) return false;
        Real dist = std::sqrt(distSq);
        m.normal = B.orientation.rotate(-diff / dist);
        m.add(A.position + m.normal * r, B.toWorld(clamped));
        return true;
    }

    // Center is inside the box: push out through the nearest face.
    int axis = 0;
    Real best = e.x - std::abs(c.x);
    for (int i = 1; i < 3; ++i) {
        Real d = e[i] - std::abs(c[i]);
        if (d < best) { best = d; axis = i; }
    }
    Real sign = c[axis] >= 0 ? Real(1) : Real(-1);
    Vec3 localOut;
    localOut[axis] = sign;
    Vec3 outward = B.orientation.rotate(localOut);
    Vec3 facePoint = c;
    facePoint[axis] = sign * e[axis];
    m.normal = -outward;
    m.add(A.position + m.normal * r, B.toWorld(facePoint));
    return true;
}

// Box vs plane: test the eight corners.
bool boxPlane(const RigidBody& A, const RigidBody& B, Real margin, RawManifold& m) {
    Vec3 np = planeNormal(B);
    Vec3 e = A.shape().halfExtents;
    m.normal = -np;

    struct Corner { Vec3 p; Real d; };
    Corner corners[8];
    int n = 0;
    for (int i = 0; i < 8; ++i) {
        Vec3 local{(i & 1) ? e.x : -e.x, (i & 2) ? e.y : -e.y, (i & 4) ? e.z : -e.z};
        Vec3 p = A.toWorld(local);
        Real d = dot(p - B.position, np);
        if (d <= margin) corners[n++] = {p, d};
    }
    if (n == 0) return false;

    // Keep the deepest four (stable insertion sort keeps the result deterministic).
    for (int i = 1; i < n; ++i)
        for (int j = i; j > 0 && corners[j].d < corners[j - 1].d; --j) std::swap(corners[j], corners[j - 1]);
    n = std::min(n, kMaxContacts);
    for (int i = 0; i < n; ++i) m.add(corners[i].p, corners[i].p - np * corners[i].d);
    return true;
}

// ---------------------------------------------------------------------------
// Box vs box: separating axis test plus face clipping.

struct BoxFrame {
    Vec3 center;
    Vec3 axis[3];
    Vec3 e;
};

BoxFrame frameOf(const RigidBody& b) {
    Mat3 R = b.orientation.toMat3();
    return {b.position, {R.col(0), R.col(1), R.col(2)}, b.shape().halfExtents};
}

Real projectedRadius(const BoxFrame& b, const Vec3& L) {
    return b.e.x * std::abs(dot(b.axis[0], L)) + b.e.y * std::abs(dot(b.axis[1], L)) +
           b.e.z * std::abs(dot(b.axis[2], L));
}

// Clips a polygon against the half space dot(n, p) <= offset.
int clipPolygon(const Vec3* in, int count, const Vec3& n, Real offset, Vec3* out) {
    int outCount = 0;
    if (count == 0) return 0;
    Vec3 prev = in[count - 1];
    Real prevDist = dot(n, prev) - offset;
    for (int i = 0; i < count; ++i) {
        Vec3 cur = in[i];
        Real curDist = dot(n, cur) - offset;
        if (prevDist <= 0 && curDist <= 0) {
            out[outCount++] = cur;
        } else if (prevDist <= 0 && curDist > 0) {
            out[outCount++] = prev + (cur - prev) * (prevDist / (prevDist - curDist));
        } else if (prevDist > 0 && curDist <= 0) {
            out[outCount++] = prev + (cur - prev) * (prevDist / (prevDist - curDist));
            out[outCount++] = cur;
        }
        prev = cur;
        prevDist = curDist;
    }
    return outCount;
}

// Face contact: clip the incident box's face against the reference face.
// `nRI` is the reference face normal, pointing from the reference box toward the incident box.
// Writes points as (point on reference, point on incident).
int clipFaces(const BoxFrame& R, int refAxis, const Vec3& nRI, const BoxFrame& I, Real margin,
              Vec3* onRef, Vec3* onInc) {
    Vec3 refCenter = R.center + nRI * R.e[refAxis];
    int ua = (refAxis + 1) % 3, va = (refAxis + 2) % 3;
    Vec3 u = R.axis[ua], v = R.axis[va];
    Real eu = R.e[ua], ev = R.e[va];

    // Incident face: the face of I most opposed to the reference normal.
    int incAxis = 0;
    Real bestDot = -1;
    for (int i = 0; i < 3; ++i) {
        Real d = std::abs(dot(I.axis[i], nRI));
        if (d > bestDot) { bestDot = d; incAxis = i; }
    }
    Real s = dot(I.axis[incAxis], nRI) > 0 ? Real(-1) : Real(1);
    Vec3 incCenter = I.center + I.axis[incAxis] * (s * I.e[incAxis]);
    Vec3 p = I.axis[(incAxis + 1) % 3] * I.e[(incAxis + 1) % 3];
    Vec3 q = I.axis[(incAxis + 2) % 3] * I.e[(incAxis + 2) % 3];

    Vec3 bufA[8] = {incCenter + p + q, incCenter - p + q, incCenter - p - q, incCenter + p - q};
    Vec3 bufB[8];
    int n = 4;
    Real cu = dot(u, refCenter), cv = dot(v, refCenter);
    n = clipPolygon(bufA, n, u, cu + eu, bufB);
    n = clipPolygon(bufB, n, -u, -cu + eu, bufA);
    n = clipPolygon(bufA, n, v, cv + ev, bufB);
    n = clipPolygon(bufB, n, -v, -cv + ev, bufA);

    int count = 0;
    for (int i = 0; i < n; ++i) {
        Real sep = dot(bufA[i] - refCenter, nRI);
        if (sep <= margin) {
            onInc[count] = bufA[i];
            onRef[count] = bufA[i] - nRI * sep;
            ++count;
        }
    }
    return count;
}

// Picks up to four points that cover the largest area, starting with the deepest.
void reduceContacts(RawManifold& m) {
    if (m.count <= kMaxContacts) return;
    RawContact* pts = m.points;
    auto mid = [](const RawContact& c) { return (c.pointA + c.pointB) * Real(0.5); };

    int i0 = 0;
    for (int i = 1; i < m.count; ++i)
        if (pts[i].depth > pts[i0].depth) i0 = i;
    Vec3 p0 = mid(pts[i0]);

    int i1 = -1;
    Real best = -1;
    for (int i = 0; i < m.count; ++i) {
        if (i == i0) continue;
        Real d = (mid(pts[i]) - p0).lengthSq();
        if (d > best) { best = d; i1 = i; }
    }
    Vec3 p1 = mid(pts[i1]);

    int i2 = -1;
    Real bestArea = -1;
    Real side2 = 0;
    for (int i = 0; i < m.count; ++i) {
        if (i == i0 || i == i1) continue;
        Real area = dot(cross(p1 - p0, mid(pts[i]) - p0), m.normal);
        if (std::abs(area) > bestArea) { bestArea = std::abs(area); i2 = i; side2 = area; }
    }

    int i3 = -1;
    bestArea = std::numeric_limits<Real>::lowest();
    for (int i = 0; i < m.count; ++i) {
        if (i == i0 || i == i1 || i == i2) continue;
        Real area = dot(cross(p1 - p0, mid(pts[i]) - p0), m.normal);
        // Prefer the point on the opposite side of the p0-p1 line from p2.
        Real score = side2 >= 0 ? -area : area;
        if (score > bestArea) { bestArea = score; i3 = i; }
    }

    RawContact keep[4] = {pts[i0], pts[i1], pts[i2], pts[i3]};
    for (int i = 0; i < 4; ++i) pts[i] = keep[i];
    m.count = 4;
}

bool boxBox(const RigidBody& bodyA, const RigidBody& bodyB, Real margin, RawManifold& m) {
    BoxFrame A = frameOf(bodyA), B = frameOf(bodyB);
    Vec3 t = B.center - A.center;

    struct Axis {
        Real overlap = std::numeric_limits<Real>::max();
        Vec3 normal;  // A -> B
        int i = -1, j = -1;
    };
    Axis faceA, faceB, edge;

    auto test = [&](const Vec3& L, Axis& best, int i, int j) {
        Real d = dot(t, L);
        Real overlap = projectedRadius(A, L) + projectedRadius(B, L) - std::abs(d);
        if (overlap < -margin) return false;
        if (overlap < best.overlap) {
            best.overlap = overlap;
            best.normal = d < 0 ? -L : L;
            best.i = i;
            best.j = j;
        }
        return true;
    };

    for (int i = 0; i < 3; ++i)
        if (!test(A.axis[i], faceA, i, -1)) return false;
    for (int j = 0; j < 3; ++j)
        if (!test(B.axis[j], faceB, -1, j)) return false;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Vec3 L = cross(A.axis[i], B.axis[j]);
            Real lenSq = L.lengthSq();
            if (lenSq < Real(1e-6)) continue;  // Parallel edges: covered by face axes.
            if (!test(L / std::sqrt(lenSq), edge, i, j)) return false;
        }
    }

    // Prefer face contacts: they give stable multi-point manifolds for resting boxes.
    const Real kRelTol = Real(0.95), kAbsTol = Real(0.0005);
    const Axis* best = &faceA;
    bool refIsA = true, isEdge = false;
    if (faceB.overlap < kRelTol * faceA.overlap - kAbsTol) { best = &faceB; refIsA = false; }
    if (edge.i >= 0 && edge.overlap < kRelTol * best->overlap - kAbsTol) { best = &edge; isEdge = true; }

    m.normal = best->normal;

    if (isEdge) {
        const Vec3& n = m.normal;
        Vec3 cA = A.center, cB = B.center;
        for (int k = 0; k < 3; ++k) {
            if (k != best->i) cA += A.axis[k] * (dot(A.axis[k], n) > 0 ? A.e[k] : -A.e[k]);
            if (k != best->j) cB += B.axis[k] * (dot(B.axis[k], n) < 0 ? B.e[k] : -B.e[k]);
        }
        Vec3 d1 = A.axis[best->i], d2 = B.axis[best->j];
        Vec3 r = cA - cB;
        Real b = dot(d1, d2), c = dot(d1, r), f = dot(d2, r);
        Real denom = 1 - b * b;
        Real s = denom > kEpsilon ? (b * f - c) / denom : 0;
        s = std::clamp(s, -A.e[best->i], A.e[best->i]);
        Real u = std::clamp(f + b * s, -B.e[best->j], B.e[best->j]);
        m.add(cA + d1 * s, cB + d2 * u);
        return true;
    }

    Vec3 onRef[8], onInc[8];
    int n;
    if (refIsA) {
        n = clipFaces(A, best->i, m.normal, B, margin, onRef, onInc);
        for (int k = 0; k < n; ++k) m.add(onRef[k], onInc[k]);
    } else {
        n = clipFaces(B, best->j, -m.normal, A, margin, onRef, onInc);
        for (int k = 0; k < n; ++k) m.add(onInc[k], onRef[k]);
    }
    reduceContacts(m);
    return m.count > 0;
}

}  // namespace

bool collide(RigidBody* a, RigidBody* b, Manifold& out, Real margin) {
    // Order the pair so each test only handles typeA <= typeB.
    bool flip = a->shape().type > b->shape().type;
    const RigidBody& A = flip ? *b : *a;
    const RigidBody& B = flip ? *a : *b;

    RawManifold m;
    bool hit = false;
    ShapeType ta = A.shape().type, tb = B.shape().type;
    if (ta == ShapeType::Sphere && tb == ShapeType::Sphere) hit = sphereSphere(A, B, margin, m);
    else if (ta == ShapeType::Sphere && tb == ShapeType::Box) hit = sphereBox(A, B, margin, m);
    else if (ta == ShapeType::Sphere && tb == ShapeType::Plane) hit = spherePlane(A, B, margin, m);
    else if (ta == ShapeType::Box && tb == ShapeType::Box) hit = boxBox(A, B, margin, m);
    else if (ta == ShapeType::Box && tb == ShapeType::Plane) hit = boxPlane(A, B, margin, m);
    if (!hit || m.count == 0) return false;

    out.a = a;
    out.b = b;
    out.normal = flip ? -m.normal : m.normal;
    out.count = std::min(m.count, kMaxContacts);
    for (int i = 0; i < out.count; ++i) {
        const RawContact& rc = m.points[i];
        ContactPoint& cp = out.points[i];
        cp = ContactPoint{};
        cp.pointA = flip ? rc.pointB : rc.pointA;
        cp.pointB = flip ? rc.pointA : rc.pointB;
        cp.depth = rc.depth;  // Unchanged by the flip.
        cp.localA = a->toLocal(cp.pointA);
        cp.localB = b->toLocal(cp.pointB);
    }
    return true;
}

}  // namespace phys
