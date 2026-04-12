#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>

namespace sdfprep {

struct Float3 {
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
};

struct Triangle {
    int v0 = 0;
    int v1 = 0;
    int v2 = 0;
};

struct BoundingBox {
    Float3 min{};
    Float3 max{};
    bool valid = false;
};

struct TriangleMeshData {
    std::vector<Float3> points;
    std::vector<Triangle> triangles;
};

struct MeshTopologyReport {
    std::size_t boundary_edges = 0;
    std::size_t non_manifold_edges = 0;
    std::size_t degenerate_faces = 0;
};

TriangleMeshData LoadTriangleMesh(const std::filesystem::path& path, double scale = 1.0);
BoundingBox ComputeBoundingBox(const TriangleMeshData& mesh);
MeshTopologyReport AnalyzeTopology(const TriangleMeshData& mesh);

}  // namespace sdfprep
