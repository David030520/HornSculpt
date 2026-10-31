#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <hornsculpt/gmsh_mesh.hpp>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>

namespace hornsculpt {
namespace {
using Id = std::uint64_t;
[[noreturn]] void fail(const std::string &s) {
    throw std::invalid_argument("Gmsh: " + s);
}
template <class T> void read(std::istream &in, T &x) {
    if constexpr (std::is_unsigned_v<T>) {
        std::string token;
        if (!(in >> token) || token.empty())
            fail("truncated unsigned value");
        T value{};
        const auto result = std::from_chars(token.data(), token.data() + token.size(), value);
        if (result.ec != std::errc{} || result.ptr != token.data() + token.size())
            fail("invalid unsigned value");
        x = value;
    } else if (!(in >> x))
        fail("truncated or invalid input");
}
void countValue(std::istream &in, std::size_t &value) {
    read(in, value);
    // CadSnapshot accepts at most two million surface triangles. Bound all
    // record lists to that same explicit project import contract.
    constexpr std::size_t maximumRecords = 2000000;
    if (value > maximumRecords)
        fail("record count exceeds the 2000000 item import limit");
}
void expect(std::istream &in, const char *expected) {
    std::string s;
    read(in, s);
    if (s != expected)
        fail(std::string("expected ") + expected + ", got " + s);
}
std::uint32_t tag32(Id tag) {
    if (tag == 0 || tag > std::numeric_limits<std::uint32_t>::max())
        fail("invalid entity or physical tag");
    return static_cast<std::uint32_t>(tag);
}
struct Parser {
    std::istream &in;
    double scale;
    std::stop_token stop;
    GmshSurfaceImport out;
    void cancelled() const {
        if (stop.stop_requested())
            throw std::runtime_error("Gmsh: Import abgebrochen");
    }
    std::unordered_map<Id, std::uint32_t> nodes;
    std::set<Id> elements;
    std::map<std::pair<int, std::uint32_t>, std::vector<std::uint32_t>> entities;
    bool hasEntities{}, hasNodes{}, hasElements{}, hasNames{};
    std::map<std::uint32_t, std::vector<std::uint32_t>> surfaceGroups;
    double version{};
    void vertex(Id id, double x, double y, double z) {
        cancelled();
        x *= scale;
        y *= scale;
        z *= scale;
        if (!id || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) ||
            out.mesh.vertices.size() >= std::numeric_limits<std::uint32_t>::max())
            fail("invalid node");
        auto idx = static_cast<std::uint32_t>(out.mesh.vertices.size());
        if (!nodes.emplace(id, idx).second)
            fail("duplicate node tag");
        out.mesh.vertices.push_back({x, y, z});
        out.nodeTags.push_back(id);
    }
    void element(Id id, int type, std::uint32_t entity, const std::vector<std::uint32_t> &physical,
                 const std::vector<Id> &tags) {
        cancelled();
        if (!id || !elements.insert(id).second)
            fail("duplicate or zero element tag");
        std::vector<std::uint32_t> indices;
        for (auto tag : tags) {
            auto it = nodes.find(tag);
            if (it == nodes.end())
                fail("unknown node reference");
            indices.push_back(it->second);
        }
        if (type == 15) {
            ++out.diagnostics.skippedPoints;
            return;
        }
        if (type == 1) {
            ++out.diagnostics.skippedLines;
            return;
        }
        auto canonical = physical;
        std::sort(canonical.begin(), canonical.end());
        const auto [group, inserted] = surfaceGroups.emplace(entity, canonical);
        if (!inserted && group->second != canonical)
            fail(
                "Eine Entity hat unterschiedliche Physical-Gruppen; getrennte Flächen exportieren");
        auto emit = [&](unsigned a, unsigned b, unsigned c) {
            if (out.mesh.triangles.size() >= 2000000)
                fail("surface exceeds the 2000000 triangle project limit");
            out.mesh.triangles.push_back(
                {{indices[a], indices[b], indices[c]}, entity, BoundaryCondition::wall});
            out.triangleMetadata.push_back({id, entity, physical});
        };
        if (type == 2) {
            emit(0, 1, 2);
            return;
        }
        // Reject warped, concave and self-crossing quads rather than change their contour.
        auto a = out.mesh.vertices[indices[0]], b = out.mesh.vertices[indices[1]],
             c = out.mesh.vertices[indices[2]], d = out.mesh.vertices[indices[3]];
        auto sub = [](Vec3 x, Vec3 y) { return Vec3{x.x - y.x, x.y - y.y, x.z - y.z}; };
        auto cross = [](Vec3 x, Vec3 y) {
            return Vec3{x.y * y.z - x.z * y.y, x.z * y.x - x.x * y.z, x.x * y.y - x.y * y.x};
        };
        auto dot = [](Vec3 x, Vec3 y) { return x.x * y.x + x.y * y.y + x.z * y.z; };
        auto length = [&](Vec3 x) { return std::sqrt(dot(x, x)); };
        const auto n = cross(sub(b, a), sub(c, a));
        const double size =
            std::max({length(sub(b, a)), length(sub(c, b)), length(sub(d, c)), length(sub(a, d))});
        const double normal = length(n);
        if (!(size > 0) || !std::isfinite(size) || !(normal > 1e-12 * size * size) ||
            std::abs(dot(sub(d, a), n)) > 1e-10 * size * normal)
            fail("quadrilateral must be planar and nondegenerate");
        const Vec3 v[4]{a, b, c, d};
        for (unsigned i = 0; i < 4; ++i) {
            const auto turn = cross(sub(v[(i + 1) % 4], v[i]), sub(v[(i + 2) % 4], v[(i + 1) % 4]));
            if (!(dot(turn, n) > 1e-12 * size * size * normal))
                fail("quadrilateral must be simple and strictly convex");
        }
        emit(0, 1, 2);
        emit(0, 2, 3);
    }
    void names() {
        if (hasNames)
            fail("repeated PhysicalNames section");
        hasNames = true;
        std::size_t n;
        countValue(in, n);
        std::set<std::pair<int, std::uint32_t>> seen;
        for (std::size_t i = 0; i < n; ++i) {
            cancelled();
            int dim;
            Id tag;
            std::string name;
            read(in, dim);
            read(in, tag);
            in >> std::ws;
            if (in.peek() != '"' || !(in >> std::quoted(name)))
                fail("physical name must be quoted");
            auto t = tag32(tag);
            if (dim < 0 || dim > 3 || !seen.emplace(dim, t).second)
                fail("invalid or duplicate physical name");
            out.physicalGroups.push_back({dim, t, name});
        }
        expect(in, "$EndPhysicalNames");
    }
    void entitySection() {
        if (hasEntities || hasNodes)
            fail("Entities must occur once before Nodes");
        hasEntities = true;
        std::size_t counts[4];
        for (auto &count : counts)
            countValue(in, count);
        for (int dim = 0; dim < 4; ++dim)
            for (std::size_t i = 0; i < counts[dim]; ++i) {
                cancelled();
                Id id;
                read(in, id);
                auto entity = tag32(id);
                for (int j = 0; j < (dim == 0 ? 3 : 6); ++j) {
                    double coordinate;
                    read(in, coordinate);
                    if (!std::isfinite(coordinate))
                        fail("invalid entity bounds");
                }
                std::size_t count;
                countValue(in, count);
                std::vector<std::uint32_t> physical;
                std::set<std::uint32_t> uniquePhysical;
                for (std::size_t j = 0; j < count; ++j) {
                    cancelled();
                    Id t;
                    read(in, t);
                    auto tag = tag32(t);
                    if (!uniquePhysical.insert(tag).second)
                        fail("duplicate entity physical tag");
                    physical.push_back(tag);
                }
                if (!entities.emplace(std::make_pair(dim, entity), std::move(physical)).second)
                    fail("duplicate entity");
                if (dim) {
                    countValue(in, count);
                    for (std::size_t j = 0; j < count; ++j) {
                        cancelled();
                        long long ref;
                        read(in, ref);
                        if (!ref)
                            fail("zero bounding entity");
                    }
                }
            }
        expect(in, "$EndEntities");
    }
    void nodeSection() {
        if (hasNodes || hasElements)
            fail("Nodes must occur once before Elements");
        hasNodes = true;
        if (version == 2.2) {
            std::size_t count;
            countValue(in, count);
            if (count > 1000000)
                fail("node count exceeds the 1000000 vertex project limit");
            for (std::size_t i = 0; i < count; ++i) {
                Id id;
                double x, y, z;
                read(in, id);
                read(in, x);
                read(in, y);
                read(in, z);
                vertex(id, x, y, z);
            }
        } else {
            std::size_t blocks, count;
            Id minimum, maximum;
            countValue(in, blocks);
            countValue(in, count);
            if (count > 1000000)
                fail("node count exceeds the 1000000 vertex project limit");
            read(in, minimum);
            read(in, maximum);
            std::size_t actual = 0;
            Id actualMin = std::numeric_limits<Id>::max(), actualMax = 0;
            for (std::size_t b = 0; b < blocks; ++b) {
                int dim, parametric;
                Id entity;
                std::size_t n;
                read(in, dim);
                read(in, entity);
                read(in, parametric);
                countValue(in, n);
                if (n > count - actual)
                    fail("node block count exceeds declared total");
                if (dim < 0 || dim > 3 || !entity || (parametric != 0 && parametric != 1))
                    fail("invalid node block");
                if (hasEntities && !entities.contains({dim, tag32(entity)}))
                    fail("unknown node entity");
                std::vector<Id> ids(n);
                for (auto &id : ids) {
                    cancelled();
                    read(in, id);
                    actualMin = std::min(actualMin, id);
                    actualMax = std::max(actualMax, id);
                }
                for (auto id : ids) {
                    double x, y, z;
                    read(in, x);
                    read(in, y);
                    read(in, z);
                    for (int j = 0; j < (parametric ? dim : 0); ++j) {
                        double u;
                        read(in, u);
                        if (!std::isfinite(u))
                            fail("invalid parameter");
                    }
                    vertex(id, x, y, z);
                }
                actual += n;
            }
            if (actual != count || (count && (actualMin != minimum || actualMax != maximum)))
                fail("node block count or range mismatch");
        }
        expect(in, "$EndNodes");
    }
    void elementSection() {
        if (!hasNodes || hasElements)
            fail("Elements require preceding Nodes");
        hasElements = true;
        if (version == 2.2) {
            std::size_t count;
            countValue(in, count);
            for (std::size_t i = 0; i < count; ++i) {
                Id id;
                int type;
                std::size_t n;
                read(in, id);
                read(in, type);
                countValue(in, n);
                const int nodeCount = type == 15  ? 1
                                      : type == 1 ? 2
                                      : type == 2 ? 3
                                      : type == 3 ? 4
                                                  : 0;
                if (!nodeCount)
                    fail("unsupported element type; only linear point, line, triangle and quad "
                         "supported");
                if (n != 2)
                    fail("MSH 2.2 requires physical and entity tags only; partition metadata is "
                         "unsupported");
                std::vector<Id> tags(n);
                for (auto &t : tags)
                    read(in, t);
                if (n != 2)
                    fail("MSH 2.2 requires physical and entity tags only; partition metadata is "
                         "unsupported");
                auto entity = tag32(tags[1]);
                std::vector<std::uint32_t> physical;
                if (tags[0])
                    physical.push_back(tag32(tags[0]));
                std::vector<Id> points(static_cast<std::size_t>(nodeCount));
                for (auto &v : points)
                    read(in, v);
                element(id, type, entity, physical, points);
            }
        } else {
            std::size_t blocks, count;
            Id minimum, maximum;
            countValue(in, blocks);
            countValue(in, count);
            read(in, minimum);
            read(in, maximum);
            std::size_t actual = 0;
            Id actualMin = std::numeric_limits<Id>::max(), actualMax = 0;
            for (std::size_t b = 0; b < blocks; ++b) {
                int dim, type;
                Id id;
                std::size_t n;
                read(in, dim);
                read(in, id);
                read(in, type);
                countValue(in, n);
                const int nodeCount = type == 15  ? 1
                                      : type == 1 ? 2
                                      : type == 2 ? 3
                                      : type == 3 ? 4
                                                  : 0;
                if (n > count - actual)
                    fail("element block count exceeds declared total");
                if (!nodeCount || (type == 15 && dim != 0) || (type == 1 && dim != 1) ||
                    ((type == 2 || type == 3) && dim != 2))
                    fail("unsupported or mismatched element block");
                auto entity = tag32(id);
                std::vector<std::uint32_t> physical;
                if (hasEntities) {
                    auto it = entities.find({dim, entity});
                    if (it == entities.end())
                        fail("unknown element entity");
                    physical = it->second;
                }
                for (std::size_t j = 0; j < n; ++j) {
                    Id elementTag;
                    read(in, elementTag);
                    std::vector<Id> points(static_cast<std::size_t>(nodeCount));
                    for (auto &v : points)
                        read(in, v);
                    actualMin = std::min(actualMin, elementTag);
                    actualMax = std::max(actualMax, elementTag);
                    element(elementTag, type, entity, physical, points);
                }
                actual += n;
            }
            if (actual != count || (count && (actualMin != minimum || actualMax != maximum)))
                fail("element block count or range mismatch");
        }
        expect(in, "$EndElements");
    }
    void compactSurfaceVertices() {
        std::vector<bool> used(out.mesh.vertices.size());
        for (const auto &triangle : out.mesh.triangles)
            for (auto index : triangle.vertices)
                used[index] = true;
        std::vector<std::uint32_t> remap(used.size());
        std::vector<Vec3> vertices;
        std::vector<std::uint64_t> tags;
        for (std::size_t i = 0; i < used.size(); ++i) {
            if ((i & 1023u) == 0)
                cancelled();
            if (used[i]) {
                remap[i] = static_cast<std::uint32_t>(vertices.size());
                vertices.push_back(out.mesh.vertices[i]);
                tags.push_back(out.nodeTags[i]);
            }
        }
        for (auto &triangle : out.mesh.triangles) {
            cancelled();
            for (auto &index : triangle.vertices)
                index = remap[index];
        }
        out.mesh.vertices = std::move(vertices);
        out.nodeTags = std::move(tags);
    }
    void diagnose() {
        struct Edge {
            unsigned count{};
            int direction{};
        };
        std::map<std::pair<std::uint32_t, std::uint32_t>, Edge> edges;
        for (const auto &t : out.mesh.triangles) {
            cancelled();
            auto a = out.mesh.vertices[t.vertices[0]], b = out.mesh.vertices[t.vertices[1]],
                 c = out.mesh.vertices[t.vertices[2]];
            const double x = (b.y - a.y) * (c.z - a.z) - (b.z - a.z) * (c.y - a.y);
            const double y = (b.z - a.z) * (c.x - a.x) - (b.x - a.x) * (c.z - a.z);
            const double z = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
            if (!(std::hypot(x, y, z) > 0))
                ++out.diagnostics.degenerateTriangles;
            for (unsigned i = 0; i < 3; ++i) {
                auto u = t.vertices[i], v = t.vertices[(i + 1) % 3];
                auto &edge = edges[std::minmax(u, v)];
                ++edge.count;
                edge.direction += (u < v ? 1 : -1);
            }
        }
        for (const auto &[key, e] : edges) {
            (void)key;
            if (e.count == 1)
                ++out.diagnostics.boundaryEdges;
            else if (e.count != 2)
                ++out.diagnostics.nonManifoldEdges;
            else if (e.direction != 0)
                ++out.diagnostics.inconsistentEdges;
        }
    }
    GmshSurfaceImport run() {
        cancelled();
        expect(in, "$MeshFormat");
        int binary, size;
        read(in, version);
        read(in, binary);
        read(in, size);
        if ((version != 2.2 && version != 4.1) || binary != 0 ||
            (size != 8 && !(version == 4.1 && size == 4)))
            fail("only ASCII MSH 2.2/4.1 supported");
        expect(in, "$EndMeshFormat");
        std::string section;
        while (in >> section) {
            if (section == "$PhysicalNames")
                names();
            else if (section == "$Entities" && version == 4.1)
                entitySection();
            else if (section == "$Nodes")
                nodeSection();
            else if (section == "$Elements")
                elementSection();
            else
                fail("unsupported or repeated section " + section);
        }
        if (!in.eof() || !hasNodes || !hasElements || out.mesh.triangles.empty())
            fail("missing mesh data");
        compactSurfaceVertices();
        diagnose();
        return std::move(out);
    }
};
} // namespace
GmshSurfaceImport importGmsh(std::istream &in, GmshImportOptions options) {
    if (!std::isfinite(options.lengthScaleToMetres) || options.lengthScaleToMetres <= 0)
        throw std::invalid_argument("Gmsh: explicit positive scale to metres required");
    return Parser{in, options.lengthScaleToMetres, options.stop}.run();
}
GmshSurfaceImport importGmshFile(const std::filesystem::path &path, GmshImportOptions options) {
    std::ifstream in(path, std::ios::binary);
    if (!in)
        throw std::runtime_error("Gmsh: cannot open file");
    return importGmsh(in, options);
}
} // namespace hornsculpt
