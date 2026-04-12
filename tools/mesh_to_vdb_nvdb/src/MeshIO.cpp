#include "MeshIO.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#define TINYOBJLOADER_IMPLEMENTATION
#include "chrono_thirdparty/tinyobjloader/tiny_obj_loader.h"
extern "C" {
#include "chrono_thirdparty/libstl/stlfile.h"
}

namespace sdfprep {
namespace {

std::string ToLower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return value;
}

FILE* OpenBinaryFileForRead(const std::filesystem::path& path) {
#ifdef _WIN32
    FILE* file = nullptr;
    return _wfopen_s(&file, path.wstring().c_str(), L"rb") == 0 ? file : nullptr;
#else
    return std::fopen(path.string().c_str(), "rb");
#endif
}

TriangleMeshData LoadObjMesh(const std::filesystem::path& path, double scale) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warning;
    std::string error;

    const auto base_dir = path.parent_path().string();
    const bool ok = tinyobj::LoadObj(&attrib, &shapes, &materials, &warning, &error, path.string().c_str(),
                                     base_dir.empty() ? nullptr : base_dir.c_str(), true);
    if (!ok) {
        throw std::runtime_error("Failed to load OBJ '" + path.string() + "': " + warning + " " + error);
    }

    TriangleMeshData mesh;
    mesh.points.reserve(attrib.vertices.size() / 3);
    for (std::size_t i = 0; i < attrib.vertices.size(); i += 3) {
        mesh.points.push_back(Float3{static_cast<float>(attrib.vertices[i + 0] * scale),
                                     static_cast<float>(attrib.vertices[i + 1] * scale),
                                     static_cast<float>(attrib.vertices[i + 2] * scale)});
    }

    for (const auto& shape : shapes) {
        std::size_t index_offset = 0;
        for (const unsigned char face_vertices : shape.mesh.num_face_vertices) {
            if (face_vertices != 3) {
                throw std::runtime_error("OBJ face is not triangular after triangulation in '" + path.string() + "'");
            }

            const auto& i0 = shape.mesh.indices[index_offset + 0];
            const auto& i1 = shape.mesh.indices[index_offset + 1];
            const auto& i2 = shape.mesh.indices[index_offset + 2];
            if (i0.vertex_index < 0 || i1.vertex_index < 0 || i2.vertex_index < 0) {
                throw std::runtime_error("OBJ contains invalid vertex indices in '" + path.string() + "'");
            }

            mesh.triangles.push_back(Triangle{i0.vertex_index, i1.vertex_index, i2.vertex_index});
            index_offset += face_vertices;
        }
    }

    return mesh;
}

TriangleMeshData LoadBinaryStlMesh(const std::filesystem::path& path, double scale) {
    FILE* file = OpenBinaryFileForRead(path);
    if (!file) {
        throw std::runtime_error("Failed to open STL '" + path.string() + "'");
    }

    char comment[80] = {};
    float* vertices = nullptr;
    vertex_t vertex_count = 0;
    vertex_t* triangles = nullptr;
    uint16_t* attributes = nullptr;
    triangle_t triangle_count = 0;

    const int rc = loadstl(file, comment, &vertices, &vertex_count, &triangles, &attributes, &triangle_count);
    std::fclose(file);
    if (rc != 0) {
        std::free(triangles);
        std::free(vertices);
        std::free(attributes);
        throw std::runtime_error("Failed to load STL '" + path.string() +
                                 "'. Only binary STL is supported by this tool.");
    }

    TriangleMeshData mesh;
    mesh.points.reserve(vertex_count);
    for (vertex_t i = 0, j = 0; i < vertex_count; ++i, j += 3) {
        mesh.points.push_back(Float3{static_cast<float>(vertices[j + 0] * scale),
                                     static_cast<float>(vertices[j + 1] * scale),
                                     static_cast<float>(vertices[j + 2] * scale)});
    }

    mesh.triangles.reserve(triangle_count);
    for (triangle_t i = 0, j = 0; i < triangle_count; ++i, j += 3) {
        mesh.triangles.push_back(Triangle{static_cast<int>(triangles[j + 0]), static_cast<int>(triangles[j + 1]),
                                          static_cast<int>(triangles[j + 2])});
    }

    std::free(triangles);
    std::free(vertices);
    std::free(attributes);
    return mesh;
}

double TriangleAreaSquared(const Float3& a, const Float3& b, const Float3& c) {
    const double abx = static_cast<double>(b.x) - static_cast<double>(a.x);
    const double aby = static_cast<double>(b.y) - static_cast<double>(a.y);
    const double abz = static_cast<double>(b.z) - static_cast<double>(a.z);

    const double acx = static_cast<double>(c.x) - static_cast<double>(a.x);
    const double acy = static_cast<double>(c.y) - static_cast<double>(a.y);
    const double acz = static_cast<double>(c.z) - static_cast<double>(a.z);

    const double cx = aby * acz - abz * acy;
    const double cy = abz * acx - abx * acz;
    const double cz = abx * acy - aby * acx;
    return cx * cx + cy * cy + cz * cz;
}

std::uint64_t MakeEdgeKey(int a, int b) {
    const auto lo = static_cast<std::uint32_t>(std::min(a, b));
    const auto hi = static_cast<std::uint32_t>(std::max(a, b));
    return (static_cast<std::uint64_t>(lo) << 32) | hi;
}

}  // namespace

TriangleMeshData LoadTriangleMesh(const std::filesystem::path& path, double scale) {
    if (scale <= 0.0) {
        throw std::runtime_error("Scale must be positive.");
    }
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("Input mesh does not exist: " + path.string());
    }

    const auto extension = ToLower(path.extension().string());
    if (extension == ".obj") {
        return LoadObjMesh(path, scale);
    }
    if (extension == ".stl") {
        return LoadBinaryStlMesh(path, scale);
    }

    throw std::runtime_error("Unsupported mesh format '" + extension +
                             "'. Supported inputs are OBJ and binary STL.");
}

BoundingBox ComputeBoundingBox(const TriangleMeshData& mesh) {
    BoundingBox bbox;
    if (mesh.points.empty()) {
        return bbox;
    }

    bbox.min = mesh.points.front();
    bbox.max = mesh.points.front();
    bbox.valid = true;

    for (const auto& p : mesh.points) {
        bbox.min.x = std::min(bbox.min.x, p.x);
        bbox.min.y = std::min(bbox.min.y, p.y);
        bbox.min.z = std::min(bbox.min.z, p.z);

        bbox.max.x = std::max(bbox.max.x, p.x);
        bbox.max.y = std::max(bbox.max.y, p.y);
        bbox.max.z = std::max(bbox.max.z, p.z);
    }

    return bbox;
}

MeshTopologyReport AnalyzeTopology(const TriangleMeshData& mesh) {
    MeshTopologyReport report;
    std::unordered_map<std::uint64_t, std::uint32_t> edge_use_count;
    edge_use_count.reserve(mesh.triangles.size() * 3);

    for (const auto& tri : mesh.triangles) {
        const std::array<int, 3> idx = {tri.v0, tri.v1, tri.v2};
        if (idx[0] == idx[1] || idx[1] == idx[2] || idx[0] == idx[2]) {
            ++report.degenerate_faces;
            continue;
        }

        if (tri.v0 < 0 || tri.v1 < 0 || tri.v2 < 0 || tri.v0 >= static_cast<int>(mesh.points.size()) ||
            tri.v1 >= static_cast<int>(mesh.points.size()) || tri.v2 >= static_cast<int>(mesh.points.size())) {
            ++report.degenerate_faces;
            continue;
        }

        const auto& a = mesh.points[tri.v0];
        const auto& b = mesh.points[tri.v1];
        const auto& c = mesh.points[tri.v2];
        if (TriangleAreaSquared(a, b, c) <= 1e-24) {
            ++report.degenerate_faces;
            continue;
        }

        ++edge_use_count[MakeEdgeKey(tri.v0, tri.v1)];
        ++edge_use_count[MakeEdgeKey(tri.v1, tri.v2)];
        ++edge_use_count[MakeEdgeKey(tri.v2, tri.v0)];
    }

    for (const auto& [edge, count] : edge_use_count) {
        (void)edge;
        if (count == 1) {
            ++report.boundary_edges;
        } else if (count > 2) {
            ++report.non_manifold_edges;
        }
    }

    return report;
}

}  // namespace sdfprep
