#include <SDL.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kEps = 1e-8;

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
};

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
    Vec3 operator-(const Vec3& other) const { return {x - other.x, y - other.y, z - other.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
};

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    };
}

double length(const Vec3& v) {
    return std::sqrt(dot(v, v));
}

Vec3 normalize(const Vec3& v, const Vec3& fallback = {0.0, 1.0, 0.0}) {
    const double len = length(v);
    if (len < kEps) {
        return fallback;
    }
    return v / len;
}

double radiansToDegrees(double radians) {
    return radians * 180.0 / kPi;
}

double planeAngleDegrees(const Vec3& a, const Vec3& b) {
    const double c = std::clamp(std::abs(dot(normalize(a), normalize(b))), 0.0, 1.0);
    return radiansToDegrees(std::acos(c));
}

struct Plane {
    Vec3 n = {0.0, 1.0, 0.0};
    double d = 0.0;

    double signedDistance(const Vec3& p) const { return dot(n, p) - d; }
};

Plane makePlane(const Vec3& normal, const Vec3& point) {
    const Vec3 n = normalize(normal);
    return {n, dot(n, point)};
}

std::optional<Vec3> intersectPlanes(const Plane& a, const Plane& b, const Plane& c) {
    const Vec3 bxc = cross(b.n, c.n);
    const Vec3 cxa = cross(c.n, a.n);
    const Vec3 axb = cross(a.n, b.n);
    const double denom = dot(a.n, bxc);
    if (std::abs(denom) < kEps) {
        return std::nullopt;
    }
    return (bxc * a.d + cxa * b.d + axb * c.d) / denom;
}

struct Face {
    std::vector<int> v;
    SDL_Color color = {180, 185, 190, 255};
};

struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<Face> faces;
};

Vec3 faceNormal(const Mesh& mesh, const Face& face) {
    Vec3 n{};
    const int count = static_cast<int>(face.v.size());
    for (int i = 0; i < count; ++i) {
        const Vec3& a = mesh.vertices[face.v[i]];
        const Vec3& b = mesh.vertices[face.v[(i + 1) % count]];
        n.x += (a.y - b.y) * (a.z + b.z);
        n.y += (a.z - b.z) * (a.x + b.x);
        n.z += (a.x - b.x) * (a.y + b.y);
    }
    return normalize(n);
}

Vec3 faceCentroid(const Mesh& mesh, const Face& face) {
    Vec3 c{};
    for (int id : face.v) {
        c += mesh.vertices[id];
    }
    return c / static_cast<double>(face.v.size());
}

Plane facePlane(const Mesh& mesh, const Face& face) {
    return makePlane(faceNormal(mesh, face), mesh.vertices[face.v.front()]);
}

bool faceHasEdge(const Face& face, int a, int b) {
    for (size_t i = 0; i < face.v.size(); ++i) {
        const int u = face.v[i];
        const int v = face.v[(i + 1) % face.v.size()];
        if ((u == a && v == b) || (u == b && v == a)) {
            return true;
        }
    }
    return false;
}

std::optional<int> adjacentFaceOnEdge(const Mesh& mesh, int faceId, int a, int b) {
    const Vec3 fmNormal = faceNormal(mesh, mesh.faces[faceId]);
    std::optional<int> best;
    double bestOrthogonality = std::numeric_limits<double>::infinity();

    for (int i = 0; i < static_cast<int>(mesh.faces.size()); ++i) {
        if (i == faceId) {
            continue;
        }
        if (!faceHasEdge(mesh.faces[i], a, b)) {
            continue;
        }
        const double score = std::abs(dot(fmNormal, faceNormal(mesh, mesh.faces[i])));
        if (!best || score < bestOrthogonality) {
            best = i;
            bestOrthogonality = score;
        }
    }
    return best;
}

std::vector<int> facesOnEdge(const Mesh& mesh, int a, int b) {
    std::vector<int> out;
    for (int i = 0; i < static_cast<int>(mesh.faces.size()); ++i) {
        if (faceHasEdge(mesh.faces[i], a, b)) {
            out.push_back(i);
        }
    }
    return out;
}

std::vector<int> facesOnVertex(const Mesh& mesh, int vertexId) {
    std::vector<int> out;
    for (int i = 0; i < static_cast<int>(mesh.faces.size()); ++i) {
        if (std::find(mesh.faces[i].v.begin(), mesh.faces[i].v.end(), vertexId) != mesh.faces[i].v.end()) {
            out.push_back(i);
        }
    }
    return out;
}

int vertexIndexInFace(const Face& face, int vertexId) {
    for (int i = 0; i < static_cast<int>(face.v.size()); ++i) {
        if (face.v[i] == vertexId) {
            return i;
        }
    }
    return -1;
}

double distancePointSegment2D(const Vec2& p, const Vec2& a, const Vec2& b) {
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double denom = dx * dx + dy * dy;
    if (denom < kEps) {
        const double px = p.x - a.x;
        const double py = p.y - a.y;
        return std::sqrt(px * px + py * py);
    }
    const double t = std::clamp(((p.x - a.x) * dx + (p.y - a.y) * dy) / denom, 0.0, 1.0);
    const double qx = a.x + t * dx;
    const double qy = a.y + t * dy;
    const double px = p.x - qx;
    const double py = p.y - qy;
    return std::sqrt(px * px + py * py);
}

void replaceVertex(Face& face, int oldId, int newId) {
    for (int& id : face.v) {
        if (id == oldId) {
            id = newId;
        }
    }
}

void removeUnusedVertices(Mesh& mesh) {
    std::vector<int> remap(mesh.vertices.size(), -1);
    std::vector<Vec3> compacted;

    for (const Face& face : mesh.faces) {
        for (int id : face.v) {
            if (id >= 0 && id < static_cast<int>(mesh.vertices.size()) && remap[id] < 0) {
                remap[id] = static_cast<int>(compacted.size());
                compacted.push_back(mesh.vertices[id]);
            }
        }
    }

    for (Face& face : mesh.faces) {
        for (int& id : face.v) {
            id = remap[id];
        }
    }
    mesh.vertices = std::move(compacted);
}

double polygonAreaMagnitude(const Mesh& mesh, const Face& face) {
    Vec3 sum{};
    for (size_t i = 0; i < face.v.size(); ++i) {
        sum += cross(mesh.vertices[face.v[i]], mesh.vertices[face.v[(i + 1) % face.v.size()]]);
    }
    return 0.5 * length(sum);
}

void cleanupMesh(Mesh& mesh) {
    for (Face& face : mesh.faces) {
        std::vector<int> cleaned;
        for (int id : face.v) {
            if (cleaned.empty() || cleaned.back() != id) {
                cleaned.push_back(id);
            }
        }
        if (cleaned.size() > 1 && cleaned.front() == cleaned.back()) {
            cleaned.pop_back();
        }
        face.v = std::move(cleaned);
    }

    mesh.faces.erase(
        std::remove_if(mesh.faces.begin(), mesh.faces.end(), [&](const Face& face) {
            if (face.v.size() < 3) {
                return true;
            }
            std::unordered_set<int> unique;
            for (int id : face.v) {
                unique.insert(id);
            }
            return unique.size() < 3 || polygonAreaMagnitude(mesh, face) < 1e-5;
        }),
        mesh.faces.end());
    removeUnusedVertices(mesh);
}

void weldCloseVertices(Mesh& mesh, double epsilon = 1e-5) {
    std::vector<int> remap(mesh.vertices.size());
    for (int i = 0; i < static_cast<int>(mesh.vertices.size()); ++i) {
        remap[i] = i;
    }

    for (int i = 0; i < static_cast<int>(mesh.vertices.size()); ++i) {
        for (int j = 0; j < i; ++j) {
            if (length(mesh.vertices[i] - mesh.vertices[j]) < epsilon) {
                remap[i] = remap[j];
                break;
            }
        }
    }

    for (Face& face : mesh.faces) {
        for (int& id : face.v) {
            id = remap[id];
        }
    }
    cleanupMesh(mesh);
}

struct PushPullResult {
    bool changed = false;
    bool clamped = false;
    int insertedFaces = 0;
    int adjustedFaces = 0;
    int modifiedFaces = 0;
    std::string message;
};

struct FaceModification {
    int faceId = -1;
    Plane target;
    Vec3 direction = {0.0, 1.0, 0.0};
    double thetaDegrees = 30.0;
};

struct EdgeSupport {
    Plane plane;
    std::optional<int> adjacent;
    bool insert = false;
};

std::optional<int> modifiedIndexForFace(const std::vector<FaceModification>& mods, int faceId) {
    for (int i = 0; i < static_cast<int>(mods.size()); ++i) {
        if (mods[i].faceId == faceId) {
            return i;
        }
    }
    return std::nullopt;
}

double clampAtLocalEvents(const Mesh& mesh, const std::vector<int>& faceIds, double amount) {
    if (std::abs(amount) < kEps) {
        return amount;
    }

    std::unordered_set<int> selected;
    std::vector<Plane> starts;
    for (int faceId : faceIds) {
        if (faceId < 0 || faceId >= static_cast<int>(mesh.faces.size())) {
            continue;
        }
        const Face& face = mesh.faces[faceId];
        starts.push_back(facePlane(mesh, face));
        selected.insert(face.v.begin(), face.v.end());
    }
    if (starts.empty()) {
        return amount;
    }

    double best = amount;
    const double sign = amount > 0.0 ? 1.0 : -1.0;
    for (int i = 0; i < static_cast<int>(mesh.vertices.size()); ++i) {
        if (selected.count(i) != 0) {
            continue;
        }
        for (const Plane& start : starts) {
            const double dist = start.signedDistance(mesh.vertices[i]);
            if (sign > 0.0 && dist > 1e-4 && dist < best) {
                best = std::max(0.0, dist - 1e-3);
            } else if (sign < 0.0 && dist < -1e-4 && dist > best) {
                best = std::min(0.0, dist + 1e-3);
            }
        }
    }
    return best;
}

std::vector<EdgeSupport> buildEdgeSupports(
    const Mesh& mesh,
    const std::vector<FaceModification>& mods,
    const FaceModification& mod,
    const Face& originalFace) {
    const int count = static_cast<int>(originalFace.v.size());
    std::vector<EdgeSupport> supports(count);
    const Vec3 targetNormal = mod.target.n;

    for (int i = 0; i < count; ++i) {
        const int a = originalFace.v[i];
        const int b = originalFace.v[(i + 1) % count];
        const Vec3 edge = mesh.vertices[b] - mesh.vertices[a];
        const std::optional<int> adj = adjacentFaceOnEdge(mesh, mod.faceId, a, b);
        bool insert = !adj.has_value();

        Plane support;
        if (adj) {
            const std::optional<int> adjacentMod = modifiedIndexForFace(mods, *adj);
            const Vec3 adjacentNormal = adjacentMod ? mods[*adjacentMod].target.n : faceNormal(mesh, mesh.faces[*adj]);
            const double angle = planeAngleDegrees(targetNormal, adjacentNormal);
            insert = angle < (90.0 - mod.thetaDegrees);

            if (!insert) {
                support = adjacentMod ? mods[*adjacentMod].target : facePlane(mesh, mesh.faces[*adj]);
            }
        }

        if (insert) {
            Vec3 sideNormal = cross(normalize(edge), normalize(mod.direction, targetNormal));
            if (length(sideNormal) < kEps) {
                sideNormal = cross(normalize(edge), targetNormal);
            }
            support = makePlane(sideNormal, mesh.vertices[a]);
        }

        supports[i] = {support, adj, insert};
    }

    return supports;
}

PushPullResult modifyFaces(Mesh& mesh, const std::vector<FaceModification>& rawMods) {
    PushPullResult result;

    std::vector<FaceModification> mods;
    std::unordered_set<int> seen;
    for (const FaceModification& mod : rawMods) {
        if (mod.faceId >= 0 && mod.faceId < static_cast<int>(mesh.faces.size()) &&
            mesh.faces[mod.faceId].v.size() >= 3 && seen.insert(mod.faceId).second) {
            mods.push_back(mod);
        }
    }

    if (mods.empty()) {
        result.message = "No modifiable faces";
        return result;
    }

    struct FaceUpdate {
        Face original;
        std::vector<EdgeSupport> supports;
        std::vector<int> newIds;
    };

    std::vector<FaceUpdate> updates;
    updates.reserve(mods.size());
    for (const FaceModification& mod : mods) {
        FaceUpdate update;
        update.original = mesh.faces[mod.faceId];
        update.supports = buildEdgeSupports(mesh, mods, mod, update.original);
        update.newIds.resize(update.original.v.size());
        updates.push_back(update);
    }

    for (int m = 0; m < static_cast<int>(mods.size()); ++m) {
        const FaceModification& mod = mods[m];
        FaceUpdate& update = updates[m];
        const int count = static_cast<int>(update.original.v.size());

        for (int i = 0; i < count; ++i) {
            const Plane& prevSupport = update.supports[(i + count - 1) % count].plane;
            const Plane& nextSupport = update.supports[i].plane;
            const Vec3 fallback = mesh.vertices[update.original.v[i]] + mod.direction;
            const std::optional<Vec3> intersection = intersectPlanes(mod.target, prevSupport, nextSupport);
            update.newIds[i] = static_cast<int>(mesh.vertices.size());
            mesh.vertices.push_back(intersection.value_or(fallback));
        }
    }

    for (int m = 0; m < static_cast<int>(mods.size()); ++m) {
        mesh.faces[mods[m].faceId].v = updates[m].newIds;
        mesh.faces[mods[m].faceId].color = {220, 170, 70, 255};
    }

    for (int m = 0; m < static_cast<int>(mods.size()); ++m) {
        const FaceUpdate& update = updates[m];
        const int count = static_cast<int>(update.original.v.size());

        for (int i = 0; i < count; ++i) {
            const int a = update.original.v[i];
            const int b = update.original.v[(i + 1) % count];
            const int na = update.newIds[i];
            const int nb = update.newIds[(i + 1) % count];
            const EdgeSupport& support = update.supports[i];

            if (support.insert) {
                Face side;
                side.v = {a, b, nb, na};
                side.color = {100, 150, 220, 255};
                mesh.faces.push_back(side);
                ++result.insertedFaces;
            } else if (support.adjacent && !modifiedIndexForFace(mods, *support.adjacent)) {
                Face& adjacent = mesh.faces[*support.adjacent];
                replaceVertex(adjacent, a, na);
                replaceVertex(adjacent, b, nb);
                adjacent.color = {105, 185, 130, 255};
                ++result.adjustedFaces;
            }
        }
    }

    weldCloseVertices(mesh);
    result.changed = true;
    result.modifiedFaces = static_cast<int>(mods.size());
    std::ostringstream out;
    out << "PushPull++: " << result.modifiedFaces << " modified, " << result.adjustedFaces << " adjusted, "
        << result.insertedFaces << " inserted";
    result.message = out.str();
    return result;
}

PushPullResult pushPullFaceAlong(Mesh& mesh, int faceId, const Vec3& direction, double requestedAmount, double thetaDegrees);

PushPullResult pushPullFace(Mesh& mesh, int faceId, double requestedAmount, double thetaDegrees) {
    if (faceId < 0 || faceId >= static_cast<int>(mesh.faces.size())) {
        PushPullResult result;
        result.message = "No face selected";
        return result;
    }
    return pushPullFaceAlong(mesh, faceId, faceNormal(mesh, mesh.faces[faceId]), requestedAmount, thetaDegrees);
}

PushPullResult pushPullFaceAlong(Mesh& mesh, int faceId, const Vec3& direction, double requestedAmount, double thetaDegrees) {
    PushPullResult result;
    if (faceId < 0 || faceId >= static_cast<int>(mesh.faces.size())) {
        result.message = "No face selected";
        return result;
    }

    const double amount = clampAtLocalEvents(mesh, {faceId}, requestedAmount);
    if (std::abs(amount) < 1e-5) {
        result.message = "Stopped at local collapse event";
        return result;
    }

    const Plane startPlane = facePlane(mesh, mesh.faces[faceId]);
    const Vec3 delta = normalize(direction, startPlane.n) * amount;
    FaceModification mod;
    mod.faceId = faceId;
    mod.target = {startPlane.n, startPlane.d + dot(startPlane.n, delta)};
    mod.direction = delta;
    mod.thetaDegrees = thetaDegrees;
    result = modifyFaces(mesh, {mod});
    result.clamped = std::abs(amount - requestedAmount) > 1e-6;
    if (result.clamped) {
        result.message += " (clamped at event)";
    }
    return result;
}

std::optional<Plane> targetPlaneForDraggedEdge(const Mesh& mesh, int faceId, int a, int b, const Vec3& delta) {
    const Face& face = mesh.faces[faceId];
    int farthest = -1;
    double best = -1.0;
    const Vec3 pa = mesh.vertices[a];
    const Vec3 pb = mesh.vertices[b];
    const Vec3 edge = pb - pa;
    const double edgeLen2 = dot(edge, edge);

    for (int id : face.v) {
        if (id == a || id == b) {
            continue;
        }
        const Vec3 rel = mesh.vertices[id] - pa;
        const Vec3 proj = pa + edge * (dot(rel, edge) / std::max(edgeLen2, kEps));
        const double dist = length(mesh.vertices[id] - proj);
        if (dist > best) {
            best = dist;
            farthest = id;
        }
    }
    if (farthest < 0) {
        return std::nullopt;
    }

    const Vec3 vf = mesh.vertices[farthest];
    const Vec3 e = normalize(edge);
    Vec3 normal = cross(normalize((pa + delta) - vf, faceNormal(mesh, face)), e);
    if (dot(normal, faceNormal(mesh, face)) < 0.0) {
        normal = normal * -1.0;
    }
    return makePlane(normal, vf);
}

std::optional<Plane> targetPlaneForDraggedVertex(const Mesh& mesh, int faceId, int vertexId, const Vec3& delta) {
    const Face& face = mesh.faces[faceId];
    const int vertexPos = vertexIndexInFace(face, vertexId);
    if (vertexPos < 0) {
        return std::nullopt;
    }

    int farthest = -1;
    double best = -1.0;
    const Vec3 v = mesh.vertices[vertexId];
    for (int id : face.v) {
        if (id == vertexId) {
            continue;
        }
        const double dist = length(mesh.vertices[id] - v);
        if (dist > best) {
            best = dist;
            farthest = id;
        }
    }
    if (farthest < 0) {
        return std::nullopt;
    }

    const Vec3 vf = mesh.vertices[farthest];
    const Vec3 oldRay = normalize(v - vf);
    const Vec3 hinge = normalize(cross(oldRay, faceNormal(mesh, face)));
    Vec3 normal = cross(normalize((v + delta) - vf, faceNormal(mesh, face)), hinge);
    if (dot(normal, faceNormal(mesh, face)) < 0.0) {
        normal = normal * -1.0;
    }
    return makePlane(normal, vf);
}

PushPullResult pushPullEdge(Mesh& mesh, int a, int b, const Vec3& direction, double requestedAmount, double thetaDegrees) {
    PushPullResult result;
    const double amount = clampAtLocalEvents(mesh, facesOnEdge(mesh, a, b), requestedAmount);
    const Vec3 delta = normalize(direction) * amount;
    std::vector<FaceModification> mods;
    for (int faceId : facesOnEdge(mesh, a, b)) {
        std::optional<Plane> target = targetPlaneForDraggedEdge(mesh, faceId, a, b, delta);
        if (target) {
            mods.push_back({faceId, *target, delta, thetaDegrees});
        }
    }
    result = modifyFaces(mesh, mods);
    result.clamped = std::abs(amount - requestedAmount) > 1e-6;
    if (result.clamped) {
        result.message += " (clamped at event)";
    }
    return result;
}

PushPullResult pushPullVertex(Mesh& mesh, int vertexId, const Vec3& direction, double requestedAmount, double thetaDegrees) {
    PushPullResult result;
    const std::vector<int> faces = facesOnVertex(mesh, vertexId);
    const double amount = clampAtLocalEvents(mesh, faces, requestedAmount);
    const Vec3 delta = normalize(direction) * amount;
    std::vector<FaceModification> mods;
    for (int faceId : faces) {
        std::optional<Plane> target = targetPlaneForDraggedVertex(mesh, faceId, vertexId, delta);
        if (target) {
            mods.push_back({faceId, *target, delta, thetaDegrees});
        }
    }
    result = modifyFaces(mesh, mods);
    result.clamped = std::abs(amount - requestedAmount) > 1e-6;
    if (result.clamped) {
        result.message += " (clamped at event)";
    }
    return result;
}

Mesh makeCube() {
    Mesh mesh;
    mesh.vertices = {
        {-1.0, -1.0, -1.0},
        {1.0, -1.0, -1.0},
        {1.0, 1.0, -1.0},
        {-1.0, 1.0, -1.0},
        {-1.0, -1.0, 1.0},
        {1.0, -1.0, 1.0},
        {1.0, 1.0, 1.0},
        {-1.0, 1.0, 1.0},
    };
    mesh.faces = {
        {{4, 5, 6, 7}, {180, 185, 190, 255}},
        {{1, 0, 3, 2}, {180, 185, 190, 255}},
        {{5, 1, 2, 6}, {180, 185, 190, 255}},
        {{0, 4, 7, 3}, {180, 185, 190, 255}},
        {{7, 6, 2, 3}, {180, 185, 190, 255}},
        {{0, 1, 5, 4}, {180, 185, 190, 255}},
    };
    return mesh;
}

Mesh makeSlantedBlock() {
    Mesh mesh;
    mesh.vertices = {
        {-1.2, -1.0, -1.0},
        {1.2, -1.0, -1.0},
        {1.2, -1.0, 1.0},
        {-1.2, -1.0, 1.0},
        {-1.2, 0.25, -1.0},
        {1.2, 0.25, -1.0},
        {1.2, 1.15, 1.0},
        {-1.2, 1.15, 1.0},
    };
    mesh.faces = {
        {{0, 1, 2, 3}, {180, 185, 190, 255}},
        {{4, 7, 6, 5}, {180, 185, 190, 255}},
        {{3, 2, 6, 7}, {180, 185, 190, 255}},
        {{1, 0, 4, 5}, {180, 185, 190, 255}},
        {{2, 1, 5, 6}, {180, 185, 190, 255}},
        {{0, 3, 7, 4}, {180, 185, 190, 255}},
    };
    return mesh;
}

struct Camera {
    double yaw = -0.62;
    double pitch = 0.55;
    double scale = 170.0;
    int width = 1120;
    int height = 780;
};

enum class ToolMode {
    Face,
    Edge,
    Vertex,
    Split,
};

struct EdgeSelection {
    int a = -1;
    int b = -1;
};

struct DirectionChoice {
    Vec3 direction;
    double theta = 30.0;
    std::string name;
};

Vec3 cameraSpace(const Camera& camera, const Vec3& p) {
    const double cy = std::cos(camera.yaw);
    const double sy = std::sin(camera.yaw);
    const double cp = std::cos(camera.pitch);
    const double sp = std::sin(camera.pitch);

    const double x1 = cy * p.x + sy * p.z;
    const double z1 = -sy * p.x + cy * p.z;
    const double y2 = cp * p.y - sp * z1;
    const double z2 = sp * p.y + cp * z1;
    return {x1, y2, z2};
}

Vec2 project(const Camera& camera, const Vec3& p) {
    const Vec3 q = cameraSpace(camera, p);
    return {
        camera.width * 0.5 + q.x * camera.scale,
        camera.height * 0.54 - q.y * camera.scale,
    };
}

Vec3 inverseCameraSpace(const Camera& camera, const Vec3& q) {
    const double cy = std::cos(camera.yaw);
    const double sy = std::sin(camera.yaw);
    const double cp = std::cos(camera.pitch);
    const double sp = std::sin(camera.pitch);

    const double y1 = cp * q.y + sp * q.z;
    const double z1 = -sp * q.y + cp * q.z;
    const double x = cy * q.x - sy * z1;
    const double z = sy * q.x + cy * z1;
    return {x, y1, z};
}

std::optional<Vec3> screenPointOnPlane(const Camera& camera, int mouseX, int mouseY, const Plane& plane) {
    const double qx = (mouseX - camera.width * 0.5) / camera.scale;
    const double qy = -(mouseY - camera.height * 0.54) / camera.scale;
    const Vec3 rayOrigin = inverseCameraSpace(camera, {qx, qy, 0.0});
    const Vec3 rayDirection = normalize(inverseCameraSpace(camera, {0.0, 0.0, 1.0}) - inverseCameraSpace(camera, {0.0, 0.0, 0.0}));
    const double denom = dot(plane.n, rayDirection);
    if (std::abs(denom) < kEps) {
        return std::nullopt;
    }
    const double t = (plane.d - dot(plane.n, rayOrigin)) / denom;
    return rayOrigin + rayDirection * t;
}

double faceDepth(const Mesh& mesh, const Face& face, const Camera& camera) {
    double sum = 0.0;
    for (int id : face.v) {
        sum += cameraSpace(camera, mesh.vertices[id]).z;
    }
    return sum / static_cast<double>(face.v.size());
}

bool pointInPolygon(const std::vector<Vec2>& polygon, double x, double y) {
    bool inside = false;
    for (size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
        const Vec2& a = polygon[i];
        const Vec2& b = polygon[j];
        const bool crosses = ((a.y > y) != (b.y > y)) &&
                             (x < (b.x - a.x) * (y - a.y) / ((b.y - a.y) + kEps) + a.x);
        if (crosses) {
            inside = !inside;
        }
    }
    return inside;
}

void setColor(SDL_Renderer* renderer, SDL_Color color) {
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
}

std::array<std::uint8_t, 7> glyph(char input) {
    const char c = static_cast<char>(std::toupper(static_cast<unsigned char>(input)));
    switch (c) {
        case 'A': return {0x0E, 0x11, 0x11, 0x1F, 0x11, 0x11, 0x11};
        case 'B': return {0x1E, 0x11, 0x11, 0x1E, 0x11, 0x11, 0x1E};
        case 'C': return {0x0E, 0x11, 0x10, 0x10, 0x10, 0x11, 0x0E};
        case 'D': return {0x1E, 0x11, 0x11, 0x11, 0x11, 0x11, 0x1E};
        case 'E': return {0x1F, 0x10, 0x10, 0x1E, 0x10, 0x10, 0x1F};
        case 'F': return {0x1F, 0x10, 0x10, 0x1E, 0x10, 0x10, 0x10};
        case 'G': return {0x0E, 0x11, 0x10, 0x17, 0x11, 0x11, 0x0F};
        case 'H': return {0x11, 0x11, 0x11, 0x1F, 0x11, 0x11, 0x11};
        case 'I': return {0x0E, 0x04, 0x04, 0x04, 0x04, 0x04, 0x0E};
        case 'J': return {0x07, 0x02, 0x02, 0x02, 0x12, 0x12, 0x0C};
        case 'K': return {0x11, 0x12, 0x14, 0x18, 0x14, 0x12, 0x11};
        case 'L': return {0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x1F};
        case 'M': return {0x11, 0x1B, 0x15, 0x15, 0x11, 0x11, 0x11};
        case 'N': return {0x11, 0x19, 0x15, 0x13, 0x11, 0x11, 0x11};
        case 'O': return {0x0E, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0E};
        case 'P': return {0x1E, 0x11, 0x11, 0x1E, 0x10, 0x10, 0x10};
        case 'Q': return {0x0E, 0x11, 0x11, 0x11, 0x15, 0x12, 0x0D};
        case 'R': return {0x1E, 0x11, 0x11, 0x1E, 0x14, 0x12, 0x11};
        case 'S': return {0x0F, 0x10, 0x10, 0x0E, 0x01, 0x01, 0x1E};
        case 'T': return {0x1F, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04};
        case 'U': return {0x11, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0E};
        case 'V': return {0x11, 0x11, 0x11, 0x11, 0x11, 0x0A, 0x04};
        case 'W': return {0x11, 0x11, 0x11, 0x15, 0x15, 0x15, 0x0A};
        case 'X': return {0x11, 0x11, 0x0A, 0x04, 0x0A, 0x11, 0x11};
        case 'Y': return {0x11, 0x11, 0x0A, 0x04, 0x04, 0x04, 0x04};
        case 'Z': return {0x1F, 0x01, 0x02, 0x04, 0x08, 0x10, 0x1F};
        case '0': return {0x0E, 0x11, 0x13, 0x15, 0x19, 0x11, 0x0E};
        case '1': return {0x04, 0x0C, 0x04, 0x04, 0x04, 0x04, 0x0E};
        case '2': return {0x0E, 0x11, 0x01, 0x02, 0x04, 0x08, 0x1F};
        case '3': return {0x1E, 0x01, 0x01, 0x0E, 0x01, 0x01, 0x1E};
        case '4': return {0x02, 0x06, 0x0A, 0x12, 0x1F, 0x02, 0x02};
        case '5': return {0x1F, 0x10, 0x10, 0x1E, 0x01, 0x01, 0x1E};
        case '6': return {0x06, 0x08, 0x10, 0x1E, 0x11, 0x11, 0x0E};
        case '7': return {0x1F, 0x01, 0x02, 0x04, 0x08, 0x08, 0x08};
        case '8': return {0x0E, 0x11, 0x11, 0x0E, 0x11, 0x11, 0x0E};
        case '9': return {0x0E, 0x11, 0x11, 0x0F, 0x01, 0x02, 0x0C};
        case '+': return {0x00, 0x04, 0x04, 0x1F, 0x04, 0x04, 0x00};
        case '-': return {0x00, 0x00, 0x00, 0x1F, 0x00, 0x00, 0x00};
        case '/': return {0x01, 0x01, 0x02, 0x04, 0x08, 0x10, 0x10};
        case ':': return {0x00, 0x04, 0x04, 0x00, 0x04, 0x04, 0x00};
        case '.': return {0x00, 0x00, 0x00, 0x00, 0x00, 0x0C, 0x0C};
        case ',': return {0x00, 0x00, 0x00, 0x00, 0x04, 0x04, 0x08};
        case '[': return {0x0E, 0x08, 0x08, 0x08, 0x08, 0x08, 0x0E};
        case ']': return {0x0E, 0x02, 0x02, 0x02, 0x02, 0x02, 0x0E};
        case ' ': return {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
        default: return {0x1F, 0x11, 0x02, 0x04, 0x04, 0x00, 0x04};
    }
}

int textWidth(const std::string& text, int scale) {
    return static_cast<int>(text.size()) * 6 * scale;
}

void drawText(SDL_Renderer* renderer, int x, int y, const std::string& text, SDL_Color color, int scale = 2) {
    setColor(renderer, color);
    int cursor = x;
    for (char c : text) {
        const std::array<std::uint8_t, 7> rows = glyph(c);
        for (int row = 0; row < 7; ++row) {
            for (int col = 0; col < 5; ++col) {
                if ((rows[row] & (1 << (4 - col))) != 0) {
                    SDL_Rect pixel{cursor + col * scale, y + row * scale, scale, scale};
                    SDL_RenderFillRect(renderer, &pixel);
                }
            }
        }
        cursor += 6 * scale;
    }
}

void drawFilledRect(SDL_Renderer* renderer, SDL_Rect rect, SDL_Color color) {
    setColor(renderer, color);
    SDL_RenderFillRect(renderer, &rect);
}

void drawOutlinedRect(SDL_Renderer* renderer, SDL_Rect rect, SDL_Color color) {
    setColor(renderer, color);
    SDL_RenderDrawRect(renderer, &rect);
}

void drawLine(SDL_Renderer* renderer, const Vec2& a, const Vec2& b) {
    SDL_RenderDrawLine(renderer, static_cast<int>(std::lround(a.x)), static_cast<int>(std::lround(a.y)),
                       static_cast<int>(std::lround(b.x)), static_cast<int>(std::lround(b.y)));
}

void fillPolygon(SDL_Renderer* renderer, const std::vector<Vec2>& points, SDL_Color color) {
    if (points.size() < 3) {
        return;
    }

    double minY = points.front().y;
    double maxY = points.front().y;
    for (const Vec2& p : points) {
        minY = std::min(minY, p.y);
        maxY = std::max(maxY, p.y);
    }

    setColor(renderer, color);
    const int y0 = static_cast<int>(std::floor(minY));
    const int y1 = static_cast<int>(std::ceil(maxY));
    std::vector<double> nodes;

    for (int y = y0; y <= y1; ++y) {
        nodes.clear();
        const double scanY = y + 0.5;
        for (size_t i = 0, j = points.size() - 1; i < points.size(); j = i++) {
            const Vec2& a = points[i];
            const Vec2& b = points[j];
            if ((a.y < scanY && b.y >= scanY) || (b.y < scanY && a.y >= scanY)) {
                nodes.push_back(a.x + (scanY - a.y) / (b.y - a.y) * (b.x - a.x));
            }
        }
        std::sort(nodes.begin(), nodes.end());
        for (size_t i = 0; i + 1 < nodes.size(); i += 2) {
            SDL_RenderDrawLine(renderer, static_cast<int>(std::lround(nodes[i])), y,
                               static_cast<int>(std::lround(nodes[i + 1])), y);
        }
    }
}

std::vector<int> sortedFacesByDepth(const Mesh& mesh, const Camera& camera, bool frontToBack) {
    std::vector<int> order(mesh.faces.size());
    for (int i = 0; i < static_cast<int>(order.size()); ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return frontToBack ? faceDepth(mesh, mesh.faces[a], camera) > faceDepth(mesh, mesh.faces[b], camera)
                           : faceDepth(mesh, mesh.faces[a], camera) < faceDepth(mesh, mesh.faces[b], camera);
    });
    return order;
}

SDL_Color shadedColor(const Mesh& mesh, const Face& face, SDL_Color base, bool selected, bool hover) {
    const Vec3 n = faceNormal(mesh, face);
    const Vec3 light = normalize(Vec3{-0.25, 0.9, 0.45});
    const double shade = std::clamp(0.48 + 0.52 * std::max(0.0, dot(n, light)), 0.35, 1.0);
    SDL_Color out = {
        static_cast<Uint8>(std::clamp<int>(static_cast<int>(base.r * shade), 0, 255)),
        static_cast<Uint8>(std::clamp<int>(static_cast<int>(base.g * shade), 0, 255)),
        static_cast<Uint8>(std::clamp<int>(static_cast<int>(base.b * shade), 0, 255)),
        255,
    };
    if (selected) {
        out = {235, 175, 70, 255};
    } else if (hover) {
        out = {210, 210, 150, 255};
    }
    return out;
}

std::optional<int> pickFace(const Mesh& mesh, const Camera& camera, int mouseX, int mouseY) {
    for (int faceId : sortedFacesByDepth(mesh, camera, true)) {
        std::vector<Vec2> polygon;
        for (int id : mesh.faces[faceId].v) {
            polygon.push_back(project(camera, mesh.vertices[id]));
        }
        if (pointInPolygon(polygon, mouseX, mouseY)) {
            return faceId;
        }
    }
    return std::nullopt;
}

std::optional<int> pickVertex(const Mesh& mesh, const Camera& camera, int mouseX, int mouseY) {
    const Vec2 p{static_cast<double>(mouseX), static_cast<double>(mouseY)};
    std::optional<int> best;
    double bestDist = 12.0;
    for (int i = 0; i < static_cast<int>(mesh.vertices.size()); ++i) {
        const Vec2 q = project(camera, mesh.vertices[i]);
        const double dx = p.x - q.x;
        const double dy = p.y - q.y;
        const double dist = std::sqrt(dx * dx + dy * dy);
        if (dist < bestDist) {
            best = i;
            bestDist = dist;
        }
    }
    return best;
}

std::optional<EdgeSelection> pickEdge(const Mesh& mesh, const Camera& camera, int mouseX, int mouseY) {
    const Vec2 p{static_cast<double>(mouseX), static_cast<double>(mouseY)};
    std::optional<EdgeSelection> best;
    double bestDist = 10.0;
    std::unordered_set<std::uint64_t> visited;

    auto edgeKey = [](int a, int b) {
        const int lo = std::min(a, b);
        const int hi = std::max(a, b);
        return (static_cast<std::uint64_t>(lo) << 32U) | static_cast<std::uint32_t>(hi);
    };

    for (int faceId : sortedFacesByDepth(mesh, camera, true)) {
        const Face& face = mesh.faces[faceId];
        for (size_t i = 0; i < face.v.size(); ++i) {
            const int a = face.v[i];
            const int b = face.v[(i + 1) % face.v.size()];
            const std::uint64_t key = edgeKey(a, b);
            if (!visited.insert(key).second) {
                continue;
            }
            const double dist = distancePointSegment2D(p, project(camera, mesh.vertices[a]), project(camera, mesh.vertices[b]));
            if (dist < bestDist) {
                best = EdgeSelection{a, b};
                bestDist = dist;
            }
        }
    }
    return best;
}

std::vector<DirectionChoice> faceDirections(const Mesh& mesh, int faceId, double overrideTheta) {
    std::vector<DirectionChoice> out;
    if (faceId < 0 || faceId >= static_cast<int>(mesh.faces.size())) {
        return out;
    }

    const Face& face = mesh.faces[faceId];
    const Vec3 n = faceNormal(mesh, face);
    out.push_back({n, overrideTheta, "normal"});
    const Vec3 projected = normalize(Vec3{n.x, 0.0, n.z}, n);
    if (length(projected) > 0.5) {
        out.push_back({projected, 15.0, "projected normal"});
    }
    out.push_back({{0.0, 1.0, 0.0}, 15.0, "world Y"});

    for (size_t i = 0; i < face.v.size(); ++i) {
        const int a = face.v[i];
        const int b = face.v[(i + 1) % face.v.size()];
        const std::optional<int> adj = adjacentFaceOnEdge(mesh, faceId, a, b);
        if (!adj) {
            continue;
        }
        const Vec3 edge = normalize(mesh.vertices[b] - mesh.vertices[a]);
        Vec3 dir = normalize(cross(faceNormal(mesh, mesh.faces[*adj]), edge), n);
        if (dot(dir, n) < 0.0) {
            dir = dir * -1.0;
        }
        out.push_back({dir, 60.0, "adjacent"});
    }
    return out;
}

std::vector<DirectionChoice> edgeDirections(const Mesh& mesh, const EdgeSelection& edge) {
    std::vector<DirectionChoice> out = {
        {{1.0, 0.0, 0.0}, 60.0, "world X"},
        {{0.0, 1.0, 0.0}, 60.0, "world Y"},
        {{0.0, 0.0, 1.0}, 60.0, "world Z"},
    };
    const Vec3 e = normalize(mesh.vertices[edge.b] - mesh.vertices[edge.a]);
    for (int faceId : facesOnEdge(mesh, edge.a, edge.b)) {
        Vec3 dir = normalize(cross(faceNormal(mesh, mesh.faces[faceId]), e), {0.0, 1.0, 0.0});
        out.push_back({dir, 60.0, "adjacent"});
    }
    return out;
}

std::vector<DirectionChoice> vertexDirections(const Mesh& mesh, int vertexId) {
    std::vector<DirectionChoice> out = {
        {{1.0, 0.0, 0.0}, 60.0, "world X"},
        {{0.0, 1.0, 0.0}, 60.0, "world Y"},
        {{0.0, 0.0, 1.0}, 60.0, "world Z"},
    };
    for (int faceId : facesOnVertex(mesh, vertexId)) {
        const Face& face = mesh.faces[faceId];
        const int pos = vertexIndexInFace(face, vertexId);
        if (pos < 0) {
            continue;
        }
        const int prev = face.v[(pos + static_cast<int>(face.v.size()) - 1) % face.v.size()];
        const int next = face.v[(pos + 1) % face.v.size()];
        out.push_back({normalize(mesh.vertices[prev] - mesh.vertices[vertexId]), 60.0, "edge"});
        out.push_back({normalize(mesh.vertices[next] - mesh.vertices[vertexId]), 60.0, "edge"});
    }
    return out;
}

std::vector<DirectionChoice> filterDirections(std::vector<DirectionChoice> directions) {
    std::vector<DirectionChoice> out;
    for (const DirectionChoice& choice : directions) {
        bool duplicate = false;
        for (const DirectionChoice& kept : out) {
            if (planeAngleDegrees(choice.direction, kept.direction) < 10.0) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            out.push_back(choice);
        }
        if (out.size() >= 5) {
            break;
        }
    }
    return out;
}

void drawDirections(SDL_Renderer* renderer, const Camera& camera, const Vec3& origin,
                    const std::vector<DirectionChoice>& directions, int active) {
    for (int i = 0; i < static_cast<int>(directions.size()); ++i) {
        const Vec2 a = project(camera, origin);
        const Vec2 b = project(camera, origin + normalize(directions[i].direction) * 0.65);
        setColor(renderer, i == active ? SDL_Color{245, 115, 55, 255} : SDL_Color{60, 125, 190, 255});
        drawLine(renderer, a, b);
        SDL_Rect tip{
            static_cast<int>(std::lround(b.x)) - 3,
            static_cast<int>(std::lround(b.y)) - 3,
            6,
            6,
        };
        SDL_RenderFillRect(renderer, &tip);
    }
}

double cross2(const Vec2& a, const Vec2& b) {
    return a.x * b.y - a.y * b.x;
}

std::optional<std::pair<double, double>> lineSegmentIntersection2D(const Vec2& p, const Vec2& r, const Vec2& q, const Vec2& s) {
    const double denom = cross2(r, s);
    if (std::abs(denom) < 1e-9) {
        return std::nullopt;
    }
    const Vec2 qp{q.x - p.x, q.y - p.y};
    const double t = cross2(qp, s) / denom;
    const double u = cross2(qp, r) / denom;
    if (u < -1e-7 || u > 1.0 + 1e-7) {
        return std::nullopt;
    }
    return std::make_pair(t, std::clamp(u, 0.0, 1.0));
}

bool splitFaceByLine(Mesh& mesh, int faceId, const Vec3& worldA, const Vec3& worldB) {
    if (faceId < 0 || faceId >= static_cast<int>(mesh.faces.size()) || mesh.faces[faceId].v.size() < 4) {
        return false;
    }

    const Face face = mesh.faces[faceId];
    const Plane plane = facePlane(mesh, face);
    const Vec3 u = normalize(worldB - worldA);
    if (length(u) < kEps) {
        return false;
    }
    const Vec3 v = normalize(cross(plane.n, u));
    const Vec2 lineP{0.0, 0.0};
    const Vec2 lineR{length(worldB - worldA), 0.0};

    struct Hit {
        int edge = -1;
        double edgeT = 0.0;
        Vec3 point;
    };
    std::vector<Hit> hits;

    for (int i = 0; i < static_cast<int>(face.v.size()); ++i) {
        const Vec3 a3 = mesh.vertices[face.v[i]];
        const Vec3 b3 = mesh.vertices[face.v[(i + 1) % face.v.size()]];
        const Vec3 ar = a3 - worldA;
        const Vec3 br = b3 - worldA;
        const Vec2 a2{dot(ar, u), dot(ar, v)};
        const Vec2 b2{dot(br, u), dot(br, v)};
        const Vec2 seg{b2.x - a2.x, b2.y - a2.y};
        const std::optional<std::pair<double, double>> hit = lineSegmentIntersection2D(lineP, lineR, a2, seg);
        if (!hit) {
            continue;
        }
        if (!hits.empty() && length((a3 + (b3 - a3) * hit->second) - hits.back().point) < 1e-5) {
            continue;
        }
        hits.push_back({i, hit->second, a3 + (b3 - a3) * hit->second});
    }

    if (hits.size() < 2) {
        return false;
    }

    std::sort(hits.begin(), hits.end(), [&](const Hit& a, const Hit& b) {
        return length(a.point - worldA) < length(b.point - worldA);
    });
    Hit h0 = hits.front();
    Hit h1 = hits.back();
    if (h0.edge == h1.edge || length(h0.point - h1.point) < 1e-5) {
        return false;
    }

    const int i0 = static_cast<int>(mesh.vertices.size());
    mesh.vertices.push_back(h0.point);
    const int i1 = static_cast<int>(mesh.vertices.size());
    mesh.vertices.push_back(h1.point);

    std::vector<int> loop;
    for (int i = 0; i < static_cast<int>(face.v.size()); ++i) {
        loop.push_back(face.v[i]);
        if (i == h0.edge) {
            loop.push_back(i0);
        }
        if (i == h1.edge) {
            loop.push_back(i1);
        }
    }

    const int p0 = static_cast<int>(std::find(loop.begin(), loop.end(), i0) - loop.begin());
    const int p1 = static_cast<int>(std::find(loop.begin(), loop.end(), i1) - loop.begin());

    std::vector<int> aLoop;
    for (int i = p0;; i = (i + 1) % static_cast<int>(loop.size())) {
        aLoop.push_back(loop[i]);
        if (i == p1) {
            break;
        }
    }

    std::vector<int> bLoop;
    for (int i = p1;; i = (i + 1) % static_cast<int>(loop.size())) {
        bLoop.push_back(loop[i]);
        if (i == p0) {
            break;
        }
    }

    mesh.faces[faceId].v = aLoop;
    mesh.faces[faceId].color = {145, 190, 170, 255};
    Face newFace;
    newFace.v = bLoop;
    newFace.color = {145, 190, 170, 255};
    mesh.faces.push_back(newFace);
    cleanupMesh(mesh);
    return true;
}

void renderMesh(SDL_Renderer* renderer, const Mesh& mesh, const Camera& camera, int selectedFace, int hoverFace) {
    for (int faceId : sortedFacesByDepth(mesh, camera, false)) {
        const Face& face = mesh.faces[faceId];
        std::vector<Vec2> polygon;
        for (int id : face.v) {
            polygon.push_back(project(camera, mesh.vertices[id]));
        }

        const SDL_Color fill = shadedColor(mesh, face, face.color, faceId == selectedFace, faceId == hoverFace);
        fillPolygon(renderer, polygon, fill);

        setColor(renderer, {38, 43, 48, 255});
        for (size_t i = 0; i < polygon.size(); ++i) {
            drawLine(renderer, polygon[i], polygon[(i + 1) % polygon.size()]);
        }
    }
}

void drawFaceNormal(SDL_Renderer* renderer, const Mesh& mesh, const Camera& camera, int faceId) {
    if (faceId < 0 || faceId >= static_cast<int>(mesh.faces.size())) {
        return;
    }
    const Face& face = mesh.faces[faceId];
    const Vec3 c = faceCentroid(mesh, face);
    const Vec3 n = faceNormal(mesh, face);
    const Vec2 a = project(camera, c);
    const Vec2 b = project(camera, c + n * 0.55);
    setColor(renderer, {245, 115, 55, 255});
    drawLine(renderer, a, b);
    SDL_Rect tip{
        static_cast<int>(std::lround(b.x)) - 4,
        static_cast<int>(std::lround(b.y)) - 4,
        8,
        8,
    };
    SDL_RenderFillRect(renderer, &tip);
}

void drawVertexHandle(SDL_Renderer* renderer, const Mesh& mesh, const Camera& camera, int vertexId, SDL_Color color) {
    if (vertexId < 0 || vertexId >= static_cast<int>(mesh.vertices.size())) {
        return;
    }
    const Vec2 p = project(camera, mesh.vertices[vertexId]);
    setColor(renderer, color);
    SDL_Rect rect{
        static_cast<int>(std::lround(p.x)) - 5,
        static_cast<int>(std::lround(p.y)) - 5,
        10,
        10,
    };
    SDL_RenderFillRect(renderer, &rect);
}

void drawEdgeHandle(SDL_Renderer* renderer, const Mesh& mesh, const Camera& camera, const EdgeSelection& edge, SDL_Color color) {
    if (edge.a < 0 || edge.b < 0 || edge.a >= static_cast<int>(mesh.vertices.size()) ||
        edge.b >= static_cast<int>(mesh.vertices.size())) {
        return;
    }
    setColor(renderer, color);
    drawLine(renderer, project(camera, mesh.vertices[edge.a]), project(camera, mesh.vertices[edge.b]));
    drawLine(renderer, project(camera, mesh.vertices[edge.a] + Vec3{0.01, 0.01, 0.01}),
             project(camera, mesh.vertices[edge.b] + Vec3{0.01, 0.01, 0.01}));
}

void drawGrid(SDL_Renderer* renderer, const Camera& camera) {
    setColor(renderer, {216, 220, 218, 255});
    for (int i = -6; i <= 6; ++i) {
        const double f = static_cast<double>(i) * 0.5;
        drawLine(renderer, project(camera, {-3.0, -1.02, f}), project(camera, {3.0, -1.02, f}));
        drawLine(renderer, project(camera, {f, -1.02, -3.0}), project(camera, {f, -1.02, 3.0}));
    }

    setColor(renderer, {190, 90, 80, 255});
    drawLine(renderer, project(camera, {-3.3, -1.01, 0.0}), project(camera, {3.3, -1.01, 0.0}));
    setColor(renderer, {80, 130, 190, 255});
    drawLine(renderer, project(camera, {0.0, -1.01, -3.3}), project(camera, {0.0, -1.01, 3.3}));
}

struct RenderState {
    ToolMode mode = ToolMode::Face;
    int selectedFace = -1;
    EdgeSelection selectedEdge;
    int selectedVertex = -1;
    int hoverFace = -1;
    int hoverVertex = -1;
    std::optional<EdgeSelection> hoverEdge;
    std::vector<DirectionChoice> directions;
    int activeDirection = 0;
    std::optional<Vec3> splitStart;
    double theta = 30.0;
    std::string status;
};

std::string modeName(ToolMode mode);

void drawToolbarButton(SDL_Renderer* renderer, int x, int y, const std::string& label, bool active) {
    const int width = std::max(70, textWidth(label, 2) + 22);
    SDL_Rect rect{x, y, width, 28};
    drawFilledRect(renderer, rect, active ? SDL_Color{245, 115, 55, 255} : SDL_Color{232, 235, 232, 255});
    drawOutlinedRect(renderer, rect, active ? SDL_Color{120, 50, 30, 255} : SDL_Color{160, 166, 164, 255});
    drawText(renderer, x + 11, y + 7, label, active ? SDL_Color{255, 255, 248, 255} : SDL_Color{42, 48, 52, 255}, 2);
}

void drawOverlay(SDL_Renderer* renderer, const Mesh& mesh, const Camera& camera, const RenderState& state) {
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

    drawFilledRect(renderer, {0, 0, camera.width, 56}, {248, 249, 246, 235});
    drawFilledRect(renderer, {0, camera.height - 42, camera.width, 42}, {248, 249, 246, 235});
    drawText(renderer, 18, 17, "PUSHPULL++", {34, 40, 44, 255}, 3);

    int buttonX = 224;
    drawToolbarButton(renderer, buttonX, 14, "F FACE", state.mode == ToolMode::Face);
    buttonX += 92;
    drawToolbarButton(renderer, buttonX, 14, "E EDGE", state.mode == ToolMode::Edge);
    buttonX += 92;
    drawToolbarButton(renderer, buttonX, 14, "V VERTEX", state.mode == ToolMode::Vertex);
    buttonX += 116;
    drawToolbarButton(renderer, buttonX, 14, "S SPLIT", state.mode == ToolMode::Split);

    const int panelW = 276;
    const int panelX = camera.width - panelW - 18;
    drawFilledRect(renderer, {panelX, 76, panelW, 220}, {255, 255, 250, 230});
    drawOutlinedRect(renderer, {panelX, 76, panelW, 220}, {184, 190, 186, 255});
    drawText(renderer, panelX + 16, 94, "INSPECTOR", {48, 54, 58, 255}, 2);

    std::ostringstream stats;
    stats << "VERTS " << mesh.vertices.size() << "  FACES " << mesh.faces.size();
    drawText(renderer, panelX + 16, 124, stats.str(), {70, 78, 82, 255}, 2);

    std::ostringstream modeLine;
    modeLine << "MODE " << modeName(state.mode);
    drawText(renderer, panelX + 16, 148, modeLine.str(), {70, 78, 82, 255}, 2);

    std::ostringstream thetaLine;
    thetaLine << "THETA " << static_cast<int>(std::lround(state.theta)) << " DEG";
    drawText(renderer, panelX + 16, 172, thetaLine.str(), {70, 78, 82, 255}, 2);

    const std::string dirName = state.directions.empty() ? "NONE" : state.directions[state.activeDirection].name;
    drawText(renderer, panelX + 16, 196, "DIR " + dirName, {70, 78, 82, 255}, 2);
    drawText(renderer, panelX + 16, 230, "STATUS", {48, 54, 58, 255}, 2);
    drawText(renderer, panelX + 16, 254, state.status.substr(0, 28), {70, 78, 82, 255}, 2);

    drawText(renderer, 18, camera.height - 28,
             "LEFT DRAG EDIT   RIGHT DRAG ORBIT   WHEEL ZOOM   TAB DIR   [ ] THETA   1/2 SCENES   R RESET",
             {54, 60, 64, 255}, 2);
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_NONE);
}

void renderScene(SDL_Renderer* renderer, const Mesh& mesh, const Camera& camera, const RenderState& state) {
    SDL_SetRenderDrawColor(renderer, 246, 247, 244, 255);
    SDL_RenderClear(renderer);
    drawGrid(renderer, camera);
    renderMesh(renderer, mesh, camera, state.mode == ToolMode::Face ? state.selectedFace : -1, state.hoverFace);
    if (state.mode == ToolMode::Face) {
        drawFaceNormal(renderer, mesh, camera, state.selectedFace >= 0 ? state.selectedFace : state.hoverFace);
    } else if (state.mode == ToolMode::Edge) {
        if (state.hoverEdge && state.selectedEdge.a < 0) {
            drawEdgeHandle(renderer, mesh, camera, *state.hoverEdge, {120, 120, 90, 255});
        }
        drawEdgeHandle(renderer, mesh, camera, state.selectedEdge, {245, 115, 55, 255});
    } else if (state.mode == ToolMode::Vertex) {
        if (state.hoverVertex >= 0 && state.selectedVertex < 0) {
            drawVertexHandle(renderer, mesh, camera, state.hoverVertex, {120, 120, 90, 255});
        }
        drawVertexHandle(renderer, mesh, camera, state.selectedVertex, {245, 115, 55, 255});
    }
    if (!state.directions.empty() && state.mode != ToolMode::Split) {
        Vec3 origin{};
        if (state.mode == ToolMode::Face && state.selectedFace >= 0) {
            origin = faceCentroid(mesh, mesh.faces[state.selectedFace]);
        } else if (state.mode == ToolMode::Edge && state.selectedEdge.a >= 0) {
            origin = (mesh.vertices[state.selectedEdge.a] + mesh.vertices[state.selectedEdge.b]) * 0.5;
        } else if (state.mode == ToolMode::Vertex && state.selectedVertex >= 0) {
            origin = mesh.vertices[state.selectedVertex];
        }
        drawDirections(renderer, camera, origin, state.directions, state.activeDirection);
    }
    if (state.splitStart) {
        const Vec2 p = project(camera, *state.splitStart);
        setColor(renderer, {245, 115, 55, 255});
        SDL_Rect rect{static_cast<int>(std::lround(p.x)) - 4, static_cast<int>(std::lround(p.y)) - 4, 8, 8};
        SDL_RenderFillRect(renderer, &rect);
    }
    drawOverlay(renderer, mesh, camera, state);
}

std::string titleFor(double theta, const std::string& message) {
    std::ostringstream out;
    out << "PushPull++ SDL2 | drag selected face | theta=" << static_cast<int>(std::lround(theta))
        << " deg | " << message;
    return out.str();
}

std::string modeName(ToolMode mode) {
    switch (mode) {
        case ToolMode::Face:
            return "face";
        case ToolMode::Edge:
            return "edge";
        case ToolMode::Vertex:
            return "vertex";
        case ToolMode::Split:
            return "split";
    }
    return "unknown";
}

bool allFacesPlanar(const Mesh& mesh) {
    for (const Face& face : mesh.faces) {
        if (face.v.size() < 3) {
            return false;
        }
        const Plane plane = facePlane(mesh, face);
        for (int id : face.v) {
            if (std::abs(plane.signedDistance(mesh.vertices[id])) > 1e-4) {
                return false;
            }
        }
    }
    return true;
}

int runSelfTest() {
    {
        Mesh mesh = makeCube();
        const int initialFaces = static_cast<int>(mesh.faces.size());
        const PushPullResult result = pushPullFace(mesh, 4, 0.45, 60.0);
        if (!result.changed || !allFacesPlanar(mesh) || static_cast<int>(mesh.faces.size()) != initialFaces) {
            SDL_Log("Cube constrained move failed: %s", result.message.c_str());
            return 1;
        }
    }

    {
        Mesh mesh = makeSlantedBlock();
        const int initialFaces = static_cast<int>(mesh.faces.size());
        const PushPullResult result = pushPullFace(mesh, 1, 0.45, 15.0);
        if (!result.changed || !allFacesPlanar(mesh) || static_cast<int>(mesh.faces.size()) <= initialFaces) {
            SDL_Log("Adaptive insertion failed: %s", result.message.c_str());
            return 1;
        }
    }

    {
        Mesh mesh = makeSlantedBlock();
        const PushPullResult result = pushPullEdge(mesh, 4, 5, {0.0, 1.0, 0.0}, 0.25, 60.0);
        if (!result.changed || !allFacesPlanar(mesh)) {
            SDL_Log("Simultaneous edge modification failed: %s", result.message.c_str());
            return 1;
        }
    }

    {
        Mesh mesh = makeSlantedBlock();
        const PushPullResult result = pushPullVertex(mesh, 6, {0.0, 1.0, 0.0}, 0.2, 60.0);
        if (!result.changed || !allFacesPlanar(mesh)) {
            SDL_Log("Simultaneous vertex modification failed: %s", result.message.c_str());
            return 1;
        }
    }

    {
        Mesh mesh = makeCube();
        const int initialFaces = static_cast<int>(mesh.faces.size());
        const bool split = splitFaceByLine(mesh, 0, {-0.35, -1.0, 1.0}, {0.35, 1.0, 1.0});
        if (!split || !allFacesPlanar(mesh) || static_cast<int>(mesh.faces.size()) != initialFaces + 1) {
            SDL_Log("Face split failed");
            return 1;
        }
    }

    SDL_Log("Self-test passed");
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc > 1 && std::string(argv[1]) == "--self-test") {
        return runSelfTest();
    }

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        SDL_Log("SDL_Init failed: %s", SDL_GetError());
        return 1;
    }

    SDL_Window* window = SDL_CreateWindow("PushPull++ SDL2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                          1120, 780, SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    if (!window) {
        SDL_Log("SDL_CreateWindow failed: %s", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        SDL_Log("SDL_CreateRenderer failed: %s", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    Camera camera;
    Mesh mesh = makeSlantedBlock();
    Mesh resetMesh = mesh;
    Mesh dragBase = mesh;
    double theta = 30.0;
    ToolMode mode = ToolMode::Face;
    int selectedFace = -1;
    EdgeSelection selectedEdge;
    int selectedVertex = -1;
    int hoverFace = -1;
    int hoverVertex = -1;
    std::optional<EdgeSelection> hoverEdge;
    std::vector<DirectionChoice> directions;
    int activeDirection = 0;
    std::optional<Vec3> splitStart;
    int splitFace = -1;
    bool running = true;
    bool leftDragging = false;
    bool rightDragging = false;
    int dragStartY = 0;
    std::string status = "F/E/V select tools, S split, Tab cycles directions";

    auto resetSelection = [&]() {
        selectedFace = -1;
        selectedEdge = {};
        selectedVertex = -1;
        splitStart.reset();
        splitFace = -1;
        activeDirection = 0;
        directions.clear();
    };

    auto refreshDirections = [&]() {
        if (mode == ToolMode::Face) {
            directions = filterDirections(faceDirections(mesh, selectedFace, theta));
        } else if (mode == ToolMode::Edge && selectedEdge.a >= 0) {
            directions = filterDirections(edgeDirections(mesh, selectedEdge));
        } else if (mode == ToolMode::Vertex && selectedVertex >= 0) {
            directions = filterDirections(vertexDirections(mesh, selectedVertex));
        } else {
            directions.clear();
        }
        if (directions.empty()) {
            activeDirection = 0;
        } else {
            activeDirection = std::clamp(activeDirection, 0, static_cast<int>(directions.size()) - 1);
        }
    };

    SDL_SetWindowTitle(window, titleFor(theta, status).c_str());

    while (running) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
                case SDL_QUIT:
                    running = false;
                    break;
                case SDL_WINDOWEVENT:
                    if (event.window.event == SDL_WINDOWEVENT_SIZE_CHANGED) {
                        camera.width = event.window.data1;
                        camera.height = event.window.data2;
                    }
                    break;
                case SDL_MOUSEMOTION: {
                    if (rightDragging) {
                        camera.yaw += event.motion.xrel * 0.008;
                        camera.pitch = std::clamp(camera.pitch + event.motion.yrel * 0.008, -1.35, 1.35);
                    } else if (leftDragging) {
                        const double amount = (dragStartY - event.motion.y) / camera.scale;
                        mesh = dragBase;
                        PushPullResult result;
                        const DirectionChoice choice = directions.empty()
                                                           ? DirectionChoice{{0.0, 1.0, 0.0}, theta, "fallback"}
                                                           : directions[activeDirection];
                        if (mode == ToolMode::Face && selectedFace >= 0) {
                            result = pushPullFaceAlong(mesh, selectedFace, choice.direction, amount, choice.theta);
                        } else if (mode == ToolMode::Edge && selectedEdge.a >= 0) {
                            result = pushPullEdge(mesh, selectedEdge.a, selectedEdge.b, choice.direction, amount, choice.theta);
                        } else if (mode == ToolMode::Vertex && selectedVertex >= 0) {
                            result = pushPullVertex(mesh, selectedVertex, choice.direction, amount, choice.theta);
                        } else {
                            result.message = "Nothing selected";
                        }
                        status = result.message;
                        SDL_SetWindowTitle(window, titleFor(theta, status).c_str());
                    } else {
                        hoverFace = pickFace(mesh, camera, event.motion.x, event.motion.y).value_or(-1);
                        hoverVertex = pickVertex(mesh, camera, event.motion.x, event.motion.y).value_or(-1);
                        hoverEdge = pickEdge(mesh, camera, event.motion.x, event.motion.y);
                    }
                    break;
                }
                case SDL_MOUSEBUTTONDOWN:
                    if (event.button.button == SDL_BUTTON_RIGHT) {
                        rightDragging = true;
                    } else if (event.button.button == SDL_BUTTON_LEFT) {
                        if (mode == ToolMode::Split) {
                            const int faceId = pickFace(mesh, camera, event.button.x, event.button.y).value_or(-1);
                            if (faceId >= 0) {
                                const std::optional<Vec3> point =
                                    screenPointOnPlane(camera, event.button.x, event.button.y, facePlane(mesh, mesh.faces[faceId]));
                                if (point && !splitStart) {
                                    splitStart = *point;
                                    splitFace = faceId;
                                    status = "Split start placed; click a second point";
                                } else if (point && splitFace == faceId && splitFaceByLine(mesh, splitFace, *splitStart, *point)) {
                                    splitStart.reset();
                                    splitFace = -1;
                                    status = "Face split";
                                } else {
                                    splitStart.reset();
                                    splitFace = -1;
                                    status = "Split cancelled";
                                }
                            }
                        } else {
                            if (mode == ToolMode::Face) {
                                selectedFace = pickFace(mesh, camera, event.button.x, event.button.y).value_or(-1);
                            } else if (mode == ToolMode::Edge) {
                                selectedEdge = pickEdge(mesh, camera, event.button.x, event.button.y).value_or(EdgeSelection{});
                            } else if (mode == ToolMode::Vertex) {
                                selectedVertex = pickVertex(mesh, camera, event.button.x, event.button.y).value_or(-1);
                            }
                            refreshDirections();
                            const bool hasSelection =
                                (mode == ToolMode::Face && selectedFace >= 0) ||
                                (mode == ToolMode::Edge && selectedEdge.a >= 0) ||
                                (mode == ToolMode::Vertex && selectedVertex >= 0);
                            if (hasSelection) {
                                leftDragging = true;
                                dragStartY = event.button.y;
                                dragBase = mesh;
                                status = "Dragging " + modeName(mode);
                            } else {
                                status = "No " + modeName(mode) + " under cursor";
                            }
                        }
                        SDL_SetWindowTitle(window, titleFor(theta, status).c_str());
                    }
                    break;
                case SDL_MOUSEBUTTONUP:
                    if (event.button.button == SDL_BUTTON_RIGHT) {
                        rightDragging = false;
                    } else if (event.button.button == SDL_BUTTON_LEFT) {
                        leftDragging = false;
                        cleanupMesh(mesh);
                        hoverFace = pickFace(mesh, camera, event.button.x, event.button.y).value_or(-1);
                        hoverVertex = pickVertex(mesh, camera, event.button.x, event.button.y).value_or(-1);
                        hoverEdge = pickEdge(mesh, camera, event.button.x, event.button.y);
                        refreshDirections();
                    }
                    break;
                case SDL_MOUSEWHEEL:
                    camera.scale = std::clamp(camera.scale * (event.wheel.y > 0 ? 1.1 : 0.9), 50.0, 650.0);
                    break;
                case SDL_KEYDOWN:
                    if (event.key.keysym.sym == SDLK_ESCAPE) {
                        running = false;
                    } else if (event.key.keysym.sym == SDLK_1) {
                        mesh = makeCube();
                        resetMesh = mesh;
                        resetSelection();
                        status = "Loaded cube";
                    } else if (event.key.keysym.sym == SDLK_2) {
                        mesh = makeSlantedBlock();
                        resetMesh = mesh;
                        resetSelection();
                        status = "Loaded slanted block";
                    } else if (event.key.keysym.sym == SDLK_r) {
                        mesh = resetMesh;
                        resetSelection();
                        status = "Reset model";
                    } else if (event.key.keysym.sym == SDLK_f) {
                        mode = ToolMode::Face;
                        resetSelection();
                        status = "Face tool";
                    } else if (event.key.keysym.sym == SDLK_e) {
                        mode = ToolMode::Edge;
                        resetSelection();
                        status = "Edge tool";
                    } else if (event.key.keysym.sym == SDLK_v) {
                        mode = ToolMode::Vertex;
                        resetSelection();
                        status = "Vertex tool";
                    } else if (event.key.keysym.sym == SDLK_s) {
                        mode = ToolMode::Split;
                        resetSelection();
                        status = "Split tool";
                    } else if (event.key.keysym.sym == SDLK_TAB) {
                        if (!directions.empty()) {
                            activeDirection = (activeDirection + 1) % static_cast<int>(directions.size());
                            status = "Direction: " + directions[activeDirection].name;
                        }
                    } else if (event.key.keysym.sym == SDLK_LEFTBRACKET) {
                        theta = std::max(0.0, theta - 5.0);
                        refreshDirections();
                        status = "Lowered threshold";
                    } else if (event.key.keysym.sym == SDLK_RIGHTBRACKET) {
                        theta = std::min(85.0, theta + 5.0);
                        refreshDirections();
                        status = "Raised threshold";
                    }
                    SDL_SetWindowTitle(window, titleFor(theta, status).c_str());
                    break;
            }
        }

        RenderState renderState;
        renderState.mode = mode;
        renderState.selectedFace = selectedFace;
        renderState.selectedEdge = selectedEdge;
        renderState.selectedVertex = selectedVertex;
        renderState.hoverFace = hoverFace;
        renderState.hoverVertex = hoverVertex;
        renderState.hoverEdge = hoverEdge;
        renderState.directions = directions;
        renderState.activeDirection = activeDirection;
        renderState.splitStart = splitStart;
        renderState.theta = theta;
        renderState.status = status;
        renderScene(renderer, mesh, camera, renderState);
        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
