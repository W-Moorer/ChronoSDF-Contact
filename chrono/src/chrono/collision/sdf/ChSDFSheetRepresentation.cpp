// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================

#include "chrono/collision/sdf/ChSDFSheetRepresentation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>
#include <unordered_map>

namespace chrono {
namespace {

struct FiberKey {
    int a = 0;
    int b = 0;

    bool operator==(const FiberKey& other) const { return a == other.a && b == other.b; }
};

struct FiberKeyHash {
    std::size_t operator()(const FiberKey& key) const {
        const std::size_t ha = static_cast<std::size_t>(static_cast<uint32_t>(key.a));
        const std::size_t hb = static_cast<std::size_t>(static_cast<uint32_t>(key.b));
        return ha * 73856093u ^ hb * 19349663u;
    }
};

struct SupportCellKey {
    int axis = -1;
    int u = 0;
    int v = 0;

    bool operator==(const SupportCellKey& other) const { return axis == other.axis && u == other.u && v == other.v; }
};

struct SupportCellKeyHash {
    std::size_t operator()(const SupportCellKey& key) const {
        const std::size_t ha = static_cast<std::size_t>(static_cast<uint32_t>(key.axis + 1));
        const std::size_t hu = static_cast<std::size_t>(static_cast<uint32_t>(key.u));
        const std::size_t hv = static_cast<std::size_t>(static_cast<uint32_t>(key.v));
        return ha * 73856093u ^ hu * 19349663u ^ hv * 83492791u;
    }
};

struct SheetCellTopologyStats {
    std::size_t layered_cell_count = 0;
    std::size_t max_layer_count_per_cell = 0;
    std::size_t largest_connected_sheet_cells = 0;
    double mean_layer_count_per_cell = 0;
    double fill_ratio = 0;
};

struct SpatialHashKey {
    int x = 0;
    int y = 0;
    int z = 0;

    bool operator==(const SpatialHashKey& other) const { return x == other.x && y == other.y && z == other.z; }
};

struct SpatialHashKeyHash {
    std::size_t operator()(const SpatialHashKey& key) const {
        const std::size_t hx = static_cast<std::size_t>(static_cast<uint32_t>(key.x));
        const std::size_t hy = static_cast<std::size_t>(static_cast<uint32_t>(key.y));
        const std::size_t hz = static_cast<std::size_t>(static_cast<uint32_t>(key.z));
        return hx * 73856093u ^ hy * 19349663u ^ hz * 83492791u;
    }
};

double AABBDistance(const ChAABB& a, const ChAABB& b);

struct LegacyFiberAccumulator {
    int dominant_axis = -1;
    ChVector2i lattice_coord = ChVector2i(0, 0);

    double measure_area = 0;
    double pressure_weight = 0;

    ChVector3d area_centroid_sum = VNULL;
    ChVector3d pressure_center_sum = VNULL;
    ChVector3d normal_sum = VNULL;
    ChVector3d force_sum = VNULL;

    ChAABB bounds_world;
    std::vector<std::size_t> source_sample_indices;
};

class UnionFind {
  public:
    explicit UnionFind(std::size_t count) : m_parent(count), m_rank(count, 0) {
        std::iota(m_parent.begin(), m_parent.end(), std::size_t{0});
    }

    std::size_t Find(std::size_t i) {
        if (m_parent[i] != i) {
            m_parent[i] = Find(m_parent[i]);
        }
        return m_parent[i];
    }

    void Unite(std::size_t a, std::size_t b) {
        a = Find(a);
        b = Find(b);
        if (a == b) {
            return;
        }

        if (m_rank[a] < m_rank[b]) {
            std::swap(a, b);
        }
        m_parent[b] = a;
        if (m_rank[a] == m_rank[b]) {
            ++m_rank[a];
        }
    }

  private:
    std::vector<std::size_t> m_parent;
    std::vector<int> m_rank;
};

ChVector3d SafeNormalized(const ChVector3d& v, const ChVector3d& fallback = VECT_Y) {
    const double length = v.Length();
    return length > 1.0e-12 ? v / length : fallback;
}

double NormalCosine(const ChVector3d& a, const ChVector3d& b) {
    const double la = a.Length();
    const double lb = b.Length();
    if (la <= 1.0e-12 || lb <= 1.0e-12) {
        return -1.0;
    }
    return Vdot(a / la, b / lb);
}

bool IsFiniteVector(const ChVector3d& v) {
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
}

void BuildOrthonormalBasis(const ChVector3d& normal, ChVector3d& tangent_u, ChVector3d& tangent_v) {
    const ChVector3d n = SafeNormalized(normal, VECT_Y);
    const ChVector3d reference = (std::abs(n.y()) < 0.9) ? VECT_Y : VECT_X;
    tangent_u = SafeNormalized(Vcross(reference, n), VECT_X);
    tangent_v = SafeNormalized(Vcross(n, tangent_u), VECT_Z);
}

double ResolveRegionSpacing(const ChSDFBrickPairWrenchResult& band_region) {
    if (band_region.region.sample_spacing > 0) {
        return band_region.region.sample_spacing;
    }

    for (const auto& sample : band_region.samples) {
        if (sample.resolution_length > 0) {
            return sample.resolution_length;
        }
    }

    return 0.0;
}

ChVector3d ResolveRegionNormal(const ChSDFBrickPairWrenchResult& band_region) {
    ChVector3d normal_sum = band_region.region.mean_normal_world;
    if (normal_sum.Length2() > 1.0e-24) {
        return SafeNormalized(normal_sum);
    }

    normal_sum = VNULL;
    for (const auto& sample : band_region.samples) {
        if (!sample.active || sample.quadrature_area <= 1.0e-16) {
            continue;
        }
        normal_sum += sample.region_sample.contact_normal_world * sample.quadrature_area;
    }

    return SafeNormalized(normal_sum, VECT_Y);
}

double ResolveSheetScale(const ChSDFBrickPairWrenchResult& band_region, const ChSDFSheetCollapseSettings& settings) {
    if (settings.fiber_lateral_tolerance > 0) {
        return settings.fiber_lateral_tolerance;
    }

    const double spacing = ResolveRegionSpacing(band_region);
    if (spacing > 0) {
        return spacing;
    }

    return 1.0;
}

double ResolveFiberTangentTolerance(const ChSDFBrickPairWrenchResult& band_region,
                                    const ChSDFSheetCollapseSettings& settings) {
    if (settings.fiber_tangent_tolerance > 0) {
        return settings.fiber_tangent_tolerance;
    }
    return ResolveSheetScale(band_region, settings);
}

double ResolveFiberPlaneTolerance(const ChSDFBrickPairWrenchResult& band_region,
                                  const ChSDFSheetCollapseSettings& settings) {
    if (settings.fiber_plane_tolerance > 0) {
        return settings.fiber_plane_tolerance;
    }
    return 0.5 * ResolveSheetScale(band_region, settings);
}

double ResolvePatchConnectionRadius(const ChSDFBrickPairWrenchResult& band_region,
                                    const ChSDFSheetCollapseSettings& settings) {
    if (settings.patch_connection_radius > 0) {
        return settings.patch_connection_radius;
    }
    return 2.0 * ResolveSheetScale(band_region, settings);
}

SpatialHashKey MakeSpatialHashKey(const ChVector3d& point_world, double cell_size) {
    if (!IsFiniteVector(point_world)) {
        return SpatialHashKey{};
    }

    const double inv_cell = cell_size > 1.0e-12 ? 1.0 / cell_size : 1.0;
    return SpatialHashKey{static_cast<int>(std::floor(point_world.x() * inv_cell)),
                          static_cast<int>(std::floor(point_world.y() * inv_cell)),
                          static_cast<int>(std::floor(point_world.z() * inv_cell))};
}

std::vector<ChVector2i> BuildPatchNeighborOffsets(int neighbor_mode) {
    std::vector<ChVector2i> offsets;
    offsets.reserve(8);

    for (int dv = -1; dv <= 1; ++dv) {
        for (int du = -1; du <= 1; ++du) {
            if (du == 0 && dv == 0) {
                continue;
            }

            const int l1 = std::abs(du) + std::abs(dv);
            const bool keep = (neighbor_mode == 4) ? (l1 == 1) : true;
            if (keep) {
                offsets.emplace_back(du, dv);
            }
        }
    }

    return offsets;
}

ChVector2i QuantizeProjectedCoord(const ChVector3d& point_world,
                                  const ChVector3d& projection_origin,
                                  const ChVector3d& region_normal,
                                  const ChVector3d& sample_normal,
                                  const ChVector3d& tangent_u,
                                  const ChVector3d& tangent_v,
                                  double lateral_tolerance) {
    const ChVector3d normal = SafeNormalized(region_normal, VECT_Y);
    const ChVector3d projection_dir = SafeNormalized(sample_normal, normal);

    ChVector3d projected = point_world;
    const double denom = Vdot(projection_dir, normal);
    if (std::abs(denom) > 1.0e-8) {
        const double alpha = Vdot(point_world - projection_origin, normal) / denom;
        projected = point_world - alpha * projection_dir;
    } else {
        projected = point_world - Vdot(point_world - projection_origin, normal) * normal;
    }

    const ChVector3d rel = projected - projection_origin;
    const double u = Vdot(rel, tangent_u);
    const double v = Vdot(rel, tangent_v);

    return ChVector2i(static_cast<int>(std::lround(u / lateral_tolerance)),
                      static_cast<int>(std::lround(v / lateral_tolerance)));
}

struct ProjectedHullMetrics;
ChSDFSheetLocalFootprint MakeSquareFootprint(const ChVector3d& origin_world,
                                             const ChVector3d& normal_world,
                                             double target_area);
void UpdateBBoxFromFootprint(const ChSDFSheetLocalFootprint& footprint, ChAABB& bbox);

void FinalizeLegacyPatch(ChSDFSheetPatch& patch,
                         const std::vector<ChSDFSheetFiberSample>& samples,
                         const ChVector3d& fallback_normal) {
    patch.support_columns = 0;
    patch.support_seed_count = 0;
    patch.measure_area = 0;
    patch.footprint_area = 0;
    patch.centroid_world = VNULL;
    patch.pressure_center_world = VNULL;
    patch.mean_normal_world = VNULL;
    patch.bounds_world = ChAABB();
    patch.support_bbox_world = ChAABB();
    patch.support_footprint = ChSDFSheetLocalFootprint();
    patch.layered_cell_count = 0;
    patch.max_layer_count_per_cell = 0;
    patch.largest_connected_sheet_cells = 0;
    patch.mean_layer_count_per_cell = 0;
    patch.sheet_fill_ratio = 0;
    patch.support_cells.clear();

    if (patch.sample_indices.empty()) {
        return;
    }

    ChVector3d centroid_sum = VNULL;
    ChVector3d normal_sum = VNULL;
    ChVector3d pressure_center_sum = VNULL;
    double pressure_weight_sum = 0;

    for (const auto sample_index : patch.sample_indices) {
        const auto& sample = samples[sample_index];
        patch.support_columns++;
        patch.support_seed_count += sample.support_seed_count;
        patch.measure_area += sample.measure_area;
        patch.footprint_area += sample.footprint_area;
        centroid_sum += sample.centroid_world * sample.measure_area;
        normal_sum += sample.normal_world * sample.measure_area;
        patch.bounds_world += sample.centroid_world;
        patch.support_bbox_world += sample.support_bbox_world;

        const double pressure_weight = sample.force_world.Length();
        pressure_weight_sum += pressure_weight;
        pressure_center_sum += sample.pressure_center_world * pressure_weight;
    }

    if (patch.measure_area > 0) {
        patch.centroid_world = centroid_sum / patch.measure_area;
    }
    patch.pressure_center_world =
        pressure_weight_sum > 0 ? pressure_center_sum / pressure_weight_sum : patch.centroid_world;
    patch.mean_normal_world = SafeNormalized(normal_sum, fallback_normal);
    patch.footprint_area = std::max(patch.footprint_area, patch.measure_area);
    patch.support_footprint = MakeSquareFootprint(patch.centroid_world, patch.mean_normal_world, patch.footprint_area);
    UpdateBBoxFromFootprint(patch.support_footprint, patch.support_bbox_world);
}

std::array<ChVector3d, 8> BuildAABBCorners(const ChAABB& box) {
    return {ChVector3d(box.min.x(), box.min.y(), box.min.z()), ChVector3d(box.min.x(), box.min.y(), box.max.z()),
            ChVector3d(box.min.x(), box.max.y(), box.min.z()), ChVector3d(box.min.x(), box.max.y(), box.max.z()),
            ChVector3d(box.max.x(), box.min.y(), box.min.z()), ChVector3d(box.max.x(), box.min.y(), box.max.z()),
            ChVector3d(box.max.x(), box.max.y(), box.min.z()), ChVector3d(box.max.x(), box.max.y(), box.max.z())};
}

double SignedPolygonArea(const std::vector<ChVector2d>& polygon) {
    if (polygon.size() < 3) {
        return 0;
    }

    double twice_area = 0;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        twice_area += a.x() * b.y() - a.y() * b.x();
    }
    return 0.5 * twice_area;
}

double PolygonArea(const std::vector<ChVector2d>& polygon) {
    return std::abs(SignedPolygonArea(polygon));
}

double Cross2D(const ChVector2d& a, const ChVector2d& b, const ChVector2d& c) {
    const ChVector2d ab = b - a;
    const ChVector2d ac = c - a;
    return ab.x() * ac.y() - ab.y() * ac.x();
}

ChVector2d PolygonCentroid(const std::vector<ChVector2d>& polygon) {
    if (polygon.empty()) {
        return ChVector2d(0, 0);
    }
    if (polygon.size() < 3) {
        ChVector2d centroid(0, 0);
        for (const auto& point : polygon) {
            centroid += point;
        }
        return centroid / static_cast<double>(polygon.size());
    }

    const double signed_area = SignedPolygonArea(polygon);
    if (std::abs(signed_area) <= 1.0e-16) {
        ChVector2d centroid(0, 0);
        for (const auto& point : polygon) {
            centroid += point;
        }
        return centroid / static_cast<double>(polygon.size());
    }

    ChVector2d centroid(0, 0);
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        const double cross = a.x() * b.y() - a.y() * b.x();
        centroid += (a + b) * cross;
    }
    return centroid / (6.0 * signed_area);
}

std::vector<ChVector2d> BuildConvexHull2D(std::vector<ChVector2d> points) {
    if (points.empty()) {
        return points;
    }

    std::sort(points.begin(), points.end(), [](const ChVector2d& a, const ChVector2d& b) {
        return a.x() != b.x() ? a.x() < b.x() : a.y() < b.y();
    });
    points.erase(std::unique(points.begin(), points.end(),
                             [](const ChVector2d& a, const ChVector2d& b) {
                                 return std::abs(a.x() - b.x()) <= 1.0e-12 && std::abs(a.y() - b.y()) <= 1.0e-12;
                             }),
                 points.end());

    if (points.size() <= 2) {
        return points;
    }

    std::vector<ChVector2d> hull;
    hull.reserve(points.size() * 2);

    for (const auto& point : points) {
        while (hull.size() >= 2 && Cross2D(hull[hull.size() - 2], hull[hull.size() - 1], point) <= 0) {
            hull.pop_back();
        }
        hull.push_back(point);
    }

    const std::size_t lower_size = hull.size();
    for (auto it = points.rbegin() + 1; it != points.rend(); ++it) {
        while (hull.size() > lower_size && Cross2D(hull[hull.size() - 2], hull[hull.size() - 1], *it) <= 0) {
            hull.pop_back();
        }
        hull.push_back(*it);
    }

    if (hull.size() > 1) {
        hull.pop_back();
    }
    return hull;
}

struct ProjectedHullMetrics {
    double area = 0;
    double perimeter = 0;
    double min_u = 0;
    double max_u = 0;
    double min_v = 0;
    double max_v = 0;
    ChVector2d centroid_uv = ChVector2d(0, 0);
    std::vector<ChVector2d> polygon_uv;
    bool valid = false;
};

ProjectedHullMetrics ComputeProjectedHullMetrics(const std::vector<ChVector3d>& support_points,
                                                 const ChVector3d& centroid_world,
                                                 const ChVector3d& normal_world) {
    ProjectedHullMetrics metrics;
    if (support_points.empty()) {
        return metrics;
    }

    ChVector3d tangent_u;
    ChVector3d tangent_v;
    BuildOrthonormalBasis(normal_world, tangent_u, tangent_v);

    std::vector<ChVector2d> projected;
    projected.reserve(support_points.size());

    metrics.min_u = std::numeric_limits<double>::infinity();
    metrics.max_u = -std::numeric_limits<double>::infinity();
    metrics.min_v = std::numeric_limits<double>::infinity();
    metrics.max_v = -std::numeric_limits<double>::infinity();

    for (const auto& point : support_points) {
        const ChVector3d rel = point - centroid_world;
        const double u = Vdot(rel, tangent_u);
        const double v = Vdot(rel, tangent_v);
        projected.emplace_back(u, v);
        metrics.min_u = std::min(metrics.min_u, u);
        metrics.max_u = std::max(metrics.max_u, u);
        metrics.min_v = std::min(metrics.min_v, v);
        metrics.max_v = std::max(metrics.max_v, v);
    }

    metrics.polygon_uv = BuildConvexHull2D(std::move(projected));
    if (metrics.polygon_uv.size() == 1) {
        metrics.valid = true;
        metrics.centroid_uv = metrics.polygon_uv.front();
        return metrics;
    }

    if (metrics.polygon_uv.empty()) {
        return metrics;
    }

    for (std::size_t i = 0; i < metrics.polygon_uv.size(); ++i) {
        const auto& a = metrics.polygon_uv[i];
        const auto& b = metrics.polygon_uv[(i + 1) % metrics.polygon_uv.size()];
        metrics.perimeter += (b - a).Length();
    }

    metrics.area = PolygonArea(metrics.polygon_uv);
    metrics.centroid_uv = PolygonCentroid(metrics.polygon_uv);
    metrics.valid = true;
    return metrics;
}

SupportCellKey MakeSupportCellKey(const ChSDFSheetSupportEvidence& evidence) {
    int iu = 0;
    int iv = 0;
    switch (evidence.source_carrier_axis) {
        case 0:
            iu = evidence.source_coord.y();
            iv = evidence.source_coord.z();
            break;
        case 1:
            iu = evidence.source_coord.x();
            iv = evidence.source_coord.z();
            break;
        case 2:
        default:
            iu = evidence.source_coord.x();
            iv = evidence.source_coord.y();
            break;
    }
    return SupportCellKey{evidence.source_carrier_axis, iu, iv};
}

ChVector2i QuantizePatchPlaneCell(const ChVector2d& uv, double cell_size) {
    const double inv_cell = cell_size > 1.0e-12 ? 1.0 / cell_size : 1.0;
    return ChVector2i(static_cast<int>(std::lround(uv.x() * inv_cell)),
                      static_cast<int>(std::lround(uv.y() * inv_cell)));
}

std::vector<ChSDFPatchPlaneLayeredCell> BuildPatchPlaneLayeredCells(const ChSDFSheetPatch& patch,
                                                                    const std::vector<ChSDFSheetFiberSample>& samples,
                                                                    const ChVector3d& patch_origin_world,
                                                                    const ChVector3d& patch_normal_world,
                                                                    double cell_size) {
    std::vector<ChSDFPatchPlaneLayeredCell> layered_cells;
    if (patch.sample_indices.empty() || cell_size <= 1.0e-12) {
        return layered_cells;
    }

    ChVector3d tangent_u;
    ChVector3d tangent_v;
    BuildOrthonormalBasis(patch_normal_world, tangent_u, tangent_v);

    std::unordered_map<SupportCellKey, ChSDFPatchPlaneLayeredCell, SupportCellKeyHash> cell_map;
    cell_map.reserve(patch.sample_indices.size() * 4);

    for (const auto sample_index : patch.sample_indices) {
        const auto& sample = samples[sample_index];
        for (const auto& evidence : sample.support_evidence) {
            const SupportCellKey key = MakeSupportCellKey(evidence);
            auto& cell = cell_map[key];
            cell.carrier_axis = evidence.source_carrier_axis;
            cell.cell_ij = ChVector2i(key.u, key.v);
            cell.half_extents_uv = ChVector2d(0.5 * cell_size, 0.5 * cell_size);
            cell.measure_area_sum += evidence.measure_area;

            const ChVector3d rel = evidence.seed_world - patch_origin_world;
            const ChVector2d uv(Vdot(rel, tangent_u), Vdot(rel, tangent_v));
            cell.center_uv += uv * evidence.measure_area;

            ChSDFPatchPlaneLayerSample layer_sample;
            layer_sample.source_carrier_axis = evidence.source_carrier_axis;
            layer_sample.source_coord = evidence.source_coord;
            layer_sample.point_world = evidence.point_world;
            layer_sample.seed_world = evidence.seed_world;
            layer_sample.normal_world = evidence.normal_world;
            layer_sample.measure_area = evidence.measure_area;
            layer_sample.normal_depth = Vdot(evidence.seed_world - patch_origin_world, patch_normal_world);
            layer_sample.source_sample_index = evidence.source_sample_index;
            cell.samples.push_back(std::move(layer_sample));
        }
    }

    layered_cells.reserve(cell_map.size());
    for (auto& item : cell_map) {
        auto& cell = item.second;
        if (cell.measure_area_sum > 1.0e-16) {
            cell.center_uv /= cell.measure_area_sum;
        }
        std::stable_sort(cell.samples.begin(), cell.samples.end(), [](const auto& a, const auto& b) {
            return a.normal_depth < b.normal_depth;
        });
        layered_cells.push_back(std::move(cell));
    }

    std::stable_sort(layered_cells.begin(), layered_cells.end(), [](const auto& a, const auto& b) {
        if (a.carrier_axis != b.carrier_axis) {
            return a.carrier_axis < b.carrier_axis;
        }
        return a.cell_ij.y() != b.cell_ij.y() ? a.cell_ij.y() < b.cell_ij.y() : a.cell_ij.x() < b.cell_ij.x();
    });
    return layered_cells;
}

std::vector<ChSDFPatchPlaneSupportCell> CollapseLayeredCellsToSheetCells(
    const std::vector<ChSDFPatchPlaneLayeredCell>& layered_cells,
    const ChVector3d& patch_origin_world,
    const ChVector3d& patch_normal_world,
    double cell_size) {
    std::vector<ChSDFPatchPlaneSupportCell> support_cells;
    if (layered_cells.empty() || cell_size <= 1.0e-12) {
        return support_cells;
    }

    ChVector3d tangent_u;
    ChVector3d tangent_v;
    BuildOrthonormalBasis(patch_normal_world, tangent_u, tangent_v);
    int min_u = std::numeric_limits<int>::max();
    int max_u = std::numeric_limits<int>::min();
    int min_v = std::numeric_limits<int>::max();
    int max_v = std::numeric_limits<int>::min();
    for (const auto& cell : layered_cells) {
        min_u = std::min(min_u, cell.cell_ij.x());
        max_u = std::max(max_u, cell.cell_ij.x());
        min_v = std::min(min_v, cell.cell_ij.y());
        max_v = std::max(max_v, cell.cell_ij.y());
    }

    const int span_u = max_u - min_u + 1;
    const int span_v = max_v - min_v + 1;
    const bool major_along_u = span_u >= span_v;
    const double area_per_sheet_cell = cell_size * cell_size;

    struct StripColumnAccumulator {
        int major_coord = 0;
        double measure_area = 0;
        double center_minor_sum = 0;
        double point_weight_sum = 0;
        ChVector3d normal_sum = VNULL;
        std::array<double, 3> axis_weight = {0.0, 0.0, 0.0};
        std::size_t strip_count = 0;
        std::vector<std::size_t> source_sample_indices;
    };

    std::unordered_map<int, StripColumnAccumulator> columns;
    columns.reserve(layered_cells.size());
    for (const auto& layered : layered_cells) {
        const int major_coord = major_along_u ? layered.cell_ij.x() : layered.cell_ij.y();
        const int minor_coord = major_along_u ? layered.cell_ij.y() : layered.cell_ij.x();
        auto& column = columns[major_coord];
        column.major_coord = major_coord;
        column.measure_area += layered.measure_area_sum;
        column.center_minor_sum += static_cast<double>(minor_coord) * layered.measure_area_sum;
        column.strip_count++;

        for (const auto& sample : layered.samples) {
            column.point_weight_sum += sample.measure_area;
            column.normal_sum += sample.normal_world * sample.measure_area;
            if (sample.source_carrier_axis >= 0 && sample.source_carrier_axis < 3) {
                column.axis_weight[static_cast<std::size_t>(sample.source_carrier_axis)] += sample.measure_area;
            }
            column.source_sample_indices.push_back(sample.source_sample_index);
        }
    }

    struct SheetBin {
        ChVector2i cell_ij = ChVector2i(0, 0);
        double measure_area = 0;
        ChVector3d normal_sum = VNULL;
        std::array<double, 3> axis_weight = {0.0, 0.0, 0.0};
        std::size_t contributing_strips = 0;
        std::vector<std::size_t> source_sample_indices;
    };

    std::unordered_map<FiberKey, SheetBin, FiberKeyHash> sheet_bins;
    sheet_bins.reserve(layered_cells.size() * 8);
    for (auto& item : columns) {
        auto& column = item.second;
        if (column.measure_area <= 1.0e-16) {
            continue;
        }

        const int width_cells =
            std::max(1, static_cast<int>(std::lround(column.measure_area / std::max(area_per_sheet_cell, 1.0e-16))));
        const double center_minor = column.center_minor_sum / column.measure_area;
        const int start_minor = static_cast<int>(std::lround(center_minor - 0.5 * static_cast<double>(width_cells - 1)));
        const double per_cell_measure = column.measure_area / static_cast<double>(width_cells);

        for (int offset = 0; offset < width_cells; ++offset) {
            const int minor_coord = start_minor + offset;
            const ChVector2i cell_ij =
                major_along_u ? ChVector2i(column.major_coord, minor_coord) : ChVector2i(minor_coord, column.major_coord);
            FiberKey key{cell_ij.x(), cell_ij.y()};
            auto& bin = sheet_bins[key];
            bin.cell_ij = cell_ij;
            bin.measure_area += per_cell_measure;
            bin.normal_sum += column.normal_sum / static_cast<double>(width_cells);
            for (std::size_t axis = 0; axis < bin.axis_weight.size(); ++axis) {
                bin.axis_weight[axis] += column.axis_weight[axis] / static_cast<double>(width_cells);
            }
            bin.contributing_strips = std::max(bin.contributing_strips, column.strip_count);
            bin.source_sample_indices.insert(bin.source_sample_indices.end(), column.source_sample_indices.begin(),
                                             column.source_sample_indices.end());
        }
    }

    support_cells.reserve(sheet_bins.size());
    for (auto& item : sheet_bins) {
        auto& bin = item.second;
        ChSDFPatchPlaneSupportCell cell;
        cell.cell_ij = bin.cell_ij;
        cell.half_extents_uv = ChVector2d(0.5 * cell_size, 0.5 * cell_size);
        cell.measure_area = bin.measure_area;
        cell.layer_count = std::max<std::size_t>(1, bin.contributing_strips);
        cell.occupied = true;
        cell.first_layer = true;
        cell.center_uv = ChVector2d(bin.cell_ij.x() * cell_size, bin.cell_ij.y() * cell_size);
        cell.representative_point_world = patch_origin_world + tangent_u * cell.center_uv.x() + tangent_v * cell.center_uv.y();
        cell.representative_normal_world = SafeNormalized(bin.normal_sum, patch_normal_world);
        cell.carrier_axis = static_cast<int>(std::distance(
            bin.axis_weight.begin(), std::max_element(bin.axis_weight.begin(), bin.axis_weight.end())));
        cell.source_sample_indices = std::move(bin.source_sample_indices);
        std::sort(cell.source_sample_indices.begin(), cell.source_sample_indices.end());
        cell.source_sample_indices.erase(
            std::unique(cell.source_sample_indices.begin(), cell.source_sample_indices.end()),
            cell.source_sample_indices.end());
        support_cells.push_back(std::move(cell));
    }

    std::stable_sort(support_cells.begin(), support_cells.end(), [](const auto& a, const auto& b) {
        if (a.carrier_axis != b.carrier_axis) {
            return a.carrier_axis < b.carrier_axis;
        }
        return a.cell_ij.y() != b.cell_ij.y() ? a.cell_ij.y() < b.cell_ij.y() : a.cell_ij.x() < b.cell_ij.x();
    });
    return support_cells;
}

SheetCellTopologyStats AnnotateSupportSheetTopology(std::size_t layered_cell_count,
                                                    std::vector<ChSDFPatchPlaneSupportCell>& support_cells) {
    SheetCellTopologyStats stats;
    stats.layered_cell_count = layered_cell_count;
    if (support_cells.empty()) {
        return stats;
    }

    std::unordered_map<SupportCellKey, std::size_t, SupportCellKeyHash> index_by_key;
    index_by_key.reserve(support_cells.size());

    int min_u = std::numeric_limits<int>::max();
    int max_u = std::numeric_limits<int>::min();
    int min_v = std::numeric_limits<int>::max();
    int max_v = std::numeric_limits<int>::min();
    std::size_t layer_sum = 0;

    for (std::size_t i = 0; i < support_cells.size(); ++i) {
        auto& cell = support_cells[i];
        const SupportCellKey key{cell.carrier_axis, cell.cell_ij.x(), cell.cell_ij.y()};
        index_by_key.emplace(key, i);
        min_u = std::min(min_u, cell.cell_ij.x());
        max_u = std::max(max_u, cell.cell_ij.x());
        min_v = std::min(min_v, cell.cell_ij.y());
        max_v = std::max(max_v, cell.cell_ij.y());
        layer_sum += cell.layer_count;
        stats.max_layer_count_per_cell = std::max(stats.max_layer_count_per_cell, cell.layer_count);
        cell.shell = false;
    }

    stats.mean_layer_count_per_cell =
        static_cast<double>(layer_sum) / static_cast<double>(support_cells.size());
    const std::size_t bbox_count = static_cast<std::size_t>(max_u - min_u + 1) * static_cast<std::size_t>(max_v - min_v + 1);
    if (bbox_count > 0) {
        stats.fill_ratio = static_cast<double>(support_cells.size()) / static_cast<double>(bbox_count);
    }

    std::vector<char> visited(support_cells.size(), 0);
    for (std::size_t i = 0; i < support_cells.size(); ++i) {
        const auto& cell = support_cells[i];
        const SupportCellKey key{cell.carrier_axis, cell.cell_ij.x(), cell.cell_ij.y()};
        const SupportCellKey neighbors[4] = {
            {key.axis, key.u - 1, key.v}, {key.axis, key.u + 1, key.v}, {key.axis, key.u, key.v - 1}, {key.axis, key.u, key.v + 1}};

        std::size_t neighbor_count = 0;
        for (const auto& neighbor : neighbors) {
            if (index_by_key.find(neighbor) != index_by_key.end()) {
                ++neighbor_count;
            }
        }
        support_cells[i].shell = neighbor_count < 4;

        if (visited[i]) {
            continue;
        }

        std::size_t component_size = 0;
        std::queue<std::size_t> frontier;
        frontier.push(i);
        visited[i] = 1;
        while (!frontier.empty()) {
            const auto current = frontier.front();
            frontier.pop();
            ++component_size;

            const auto& current_cell = support_cells[current];
            const SupportCellKey current_key{current_cell.carrier_axis, current_cell.cell_ij.x(), current_cell.cell_ij.y()};
            const SupportCellKey current_neighbors[4] = {{current_key.axis, current_key.u - 1, current_key.v},
                                                         {current_key.axis, current_key.u + 1, current_key.v},
                                                         {current_key.axis, current_key.u, current_key.v - 1},
                                                         {current_key.axis, current_key.u, current_key.v + 1}};
            for (const auto& neighbor : current_neighbors) {
                const auto it = index_by_key.find(neighbor);
                if (it == index_by_key.end() || visited[it->second]) {
                    continue;
                }
                visited[it->second] = 1;
                frontier.push(it->second);
            }
        }

        stats.largest_connected_sheet_cells = std::max(stats.largest_connected_sheet_cells, component_size);
    }

    return stats;
}

ChSDFSheetLocalFootprint MakeSquareFootprint(const ChVector3d& origin_world,
                                             const ChVector3d& normal_world,
                                             double target_area) {
    ChSDFSheetLocalFootprint footprint;
    footprint.origin_world = origin_world;
    BuildOrthonormalBasis(normal_world, footprint.tangent_u_world, footprint.tangent_v_world);
    if (target_area <= 1.0e-16) {
        return footprint;
    }

    const double half_extent = 0.5 * std::sqrt(target_area);
    footprint.polygon_uv = {ChVector2d(-half_extent, -half_extent), ChVector2d(half_extent, -half_extent),
                            ChVector2d(half_extent, half_extent), ChVector2d(-half_extent, half_extent)};
    footprint.area = PolygonArea(footprint.polygon_uv);
    footprint.centroid_uv = PolygonCentroid(footprint.polygon_uv);
    return footprint;
}

ChSDFSheetLocalFootprint BuildScaledProjectedFootprint(const ProjectedHullMetrics& hull,
                                                       const ChVector3d& origin_world,
                                                       const ChVector3d& normal_world,
                                                       double target_area) {
    const double min_reasonable_hull_area = std::max(1.0e-12, 0.10 * std::max(target_area, 0.0));
    if (!hull.valid || hull.polygon_uv.size() < 3 || target_area <= 1.0e-16 || hull.area < min_reasonable_hull_area) {
        return MakeSquareFootprint(origin_world, normal_world, target_area);
    }

    ChSDFSheetLocalFootprint footprint;
    footprint.origin_world = origin_world;
    BuildOrthonormalBasis(normal_world, footprint.tangent_u_world, footprint.tangent_v_world);
    footprint.polygon_uv = hull.polygon_uv;
    const double current_area = std::max(hull.area, 1.0e-16);
    const ChVector2d centroid_uv = PolygonCentroid(footprint.polygon_uv);
    const double scale = std::sqrt(target_area / current_area);
    for (auto& vertex : footprint.polygon_uv) {
        vertex = centroid_uv + (vertex - centroid_uv) * scale;
    }
    footprint.area = PolygonArea(footprint.polygon_uv);
    footprint.centroid_uv = PolygonCentroid(footprint.polygon_uv);
    return footprint;
}

void UpdateBBoxFromFootprint(const ChSDFSheetLocalFootprint& footprint, ChAABB& bbox) {
    bbox = ChAABB();
    if (!footprint.HasPolygon()) {
        return;
    }
    for (const auto& vertex_uv : footprint.polygon_uv) {
        const ChVector3d vertex_world = footprint.origin_world + footprint.tangent_u_world * vertex_uv.x() +
                                        footprint.tangent_v_world * vertex_uv.y();
        bbox += vertex_world;
    }
}

void AppendFootprintWorldVertices(const ChSDFSheetLocalFootprint& footprint, std::vector<ChVector3d>& points) {
    if (!footprint.HasPolygon()) {
        return;
    }

    points.reserve(points.size() + footprint.polygon_uv.size());
    for (const auto& vertex_uv : footprint.polygon_uv) {
        const ChVector3d vertex_world = footprint.origin_world + footprint.tangent_u_world * vertex_uv.x() +
                                        footprint.tangent_v_world * vertex_uv.y();
        points.push_back(vertex_world);
    }
}

double ResolvePatchRecoveryRadius(double patch_measure_area, std::size_t support_seed_count, double sheet_scale) {
    double radius = sheet_scale > 0 ? sheet_scale : 0.0;
    if (support_seed_count > 0 && patch_measure_area > 0) {
        radius = std::sqrt(patch_measure_area / static_cast<double>(support_seed_count));
    }
    if (sheet_scale > 0) {
        radius = std::max(0.5 * sheet_scale, std::min(radius, 1.5 * sheet_scale));
    }
    return std::max(radius, 0.0);
}

double EstimateAxisHalfStep(std::vector<double> coords, double fallback) {
    if (coords.size() < 2) {
        return std::max(fallback, 0.0);
    }

    std::sort(coords.begin(), coords.end());
    std::vector<double> unique_coords;
    unique_coords.reserve(coords.size());
    for (const double value : coords) {
        if (unique_coords.empty() || std::abs(value - unique_coords.back()) > 1.0e-8) {
            unique_coords.push_back(value);
        }
    }

    if (unique_coords.size() < 2) {
        return std::max(fallback, 0.0);
    }

    std::vector<double> diffs;
    diffs.reserve(unique_coords.size() - 1);
    for (std::size_t i = 1; i < unique_coords.size(); ++i) {
        const double diff = unique_coords[i] - unique_coords[i - 1];
        if (diff > 1.0e-8) {
            diffs.push_back(diff);
        }
    }

    if (diffs.empty()) {
        return std::max(fallback, 0.0);
    }

    std::sort(diffs.begin(), diffs.end());
    return 0.5 * diffs[diffs.size() / 2];
}

void ExpandProjectedSupportBBox(ChAABB& bbox,
                                const ChVector3d& centroid_world,
                                const ChVector3d& normal_world,
                                const std::vector<ChVector3d>& support_points,
                                double target_area) {
    if (support_points.empty() || target_area <= 1.0e-16) {
        return;
    }

    ChVector3d tangent_u;
    ChVector3d tangent_v;
    BuildOrthonormalBasis(normal_world, tangent_u, tangent_v);

    double min_u = std::numeric_limits<double>::infinity();
    double max_u = -std::numeric_limits<double>::infinity();
    double min_v = std::numeric_limits<double>::infinity();
    double max_v = -std::numeric_limits<double>::infinity();

    for (const auto& point : support_points) {
        const ChVector3d rel = point - centroid_world;
        const double u = Vdot(rel, tangent_u);
        const double v = Vdot(rel, tangent_v);
        min_u = std::min(min_u, u);
        max_u = std::max(max_u, u);
        min_v = std::min(min_v, v);
        max_v = std::max(max_v, v);
    }

    double half_u = 0.5 * (std::isfinite(min_u) ? std::max(max_u - min_u, 0.0) : 0.0);
    double half_v = 0.5 * (std::isfinite(min_v) ? std::max(max_v - min_v, 0.0) : 0.0);
    const double center_u = std::isfinite(min_u) && std::isfinite(max_u) ? 0.5 * (min_u + max_u) : 0.0;
    const double center_v = std::isfinite(min_v) && std::isfinite(max_v) ? 0.5 * (min_v + max_v) : 0.0;
    const double min_half = std::sqrt(target_area / CH_PI);

    if (half_u <= 1.0e-12 && half_v <= 1.0e-12) {
        half_u = min_half;
        half_v = min_half;
    } else {
        const double full_u = std::max(2.0 * half_u, 1.0e-12);
        const double full_v = std::max(2.0 * half_v, 1.0e-12);
        const double observed_rect_area = full_u * full_v;
        const double required_area = std::max(target_area, observed_rect_area);

        if (full_u >= full_v) {
            half_v = std::max(half_v, 0.5 * required_area / full_u);
        } else {
            half_u = std::max(half_u, 0.5 * required_area / full_v);
        }
    }

    const ChVector3d rect_center = centroid_world + tangent_u * center_u + tangent_v * center_v;
    const std::array<ChVector3d, 4> corners = {
        rect_center + tangent_u * half_u + tangent_v * half_v,
        rect_center + tangent_u * half_u - tangent_v * half_v,
        rect_center - tangent_u * half_u + tangent_v * half_v,
        rect_center - tangent_u * half_u - tangent_v * half_v};
    for (const auto& corner : corners) {
        bbox += corner;
    }
}

void FinalizeV2Patch(ChSDFSheetPatch& patch,
                     const std::vector<ChSDFSheetFiberSample>& samples,
                     double sheet_scale) {
    patch.support_columns = patch.sample_indices.size();
    patch.support_seed_count = 0;
    patch.measure_area = 0;
    patch.footprint_area = 0;
    patch.centroid_world = VNULL;
    patch.pressure_center_world = VNULL;
    patch.mean_normal_world = VNULL;
    patch.bounds_world = ChAABB();
    patch.support_bbox_world = ChAABB();
    patch.support_footprint = ChSDFSheetLocalFootprint();
    patch.layered_cell_count = 0;
    patch.max_layer_count_per_cell = 0;
    patch.largest_connected_sheet_cells = 0;
    patch.mean_layer_count_per_cell = 0;
    patch.sheet_fill_ratio = 0;
    patch.support_cells.clear();

    if (patch.sample_indices.empty()) {
        return;
    }

    ChVector3d centroid_sum = VNULL;
    ChVector3d normal_sum = VNULL;
    ChVector3d pressure_center_sum = VNULL;
    double pressure_weight_sum = 0;
    std::vector<ChVector3d> support_points;
    std::vector<ChVector3d> occupancy_points;
    std::vector<ChVector3d> column_centers;
    support_points.reserve(patch.sample_indices.size() * 8);
    occupancy_points.reserve(patch.sample_indices.size() * 8);
    column_centers.reserve(patch.sample_indices.size());

    for (const auto sample_index : patch.sample_indices) {
        const auto& sample = samples[sample_index];
        patch.support_seed_count += sample.support_seed_count;
        patch.measure_area += sample.measure_area;
        patch.footprint_area += sample.footprint_area;
        centroid_sum += sample.centroid_world * sample.measure_area;
        normal_sum += sample.normal_world * sample.measure_area;
        patch.bounds_world += sample.source_bounds_world;
        patch.support_bbox_world += sample.support_bbox_world;
        column_centers.push_back(sample.centroid_world);

        if (sample.support_footprint.HasPolygon()) {
            AppendFootprintWorldVertices(sample.support_footprint, support_points);
        } else if (!sample.support_bbox_world.IsInverted()) {
            const auto corners = BuildAABBCorners(sample.support_bbox_world);
            support_points.insert(support_points.end(), corners.begin(), corners.end());
        }
        if (!sample.source_bounds_world.IsInverted()) {
            const auto source_corners = BuildAABBCorners(sample.source_bounds_world);
            occupancy_points.insert(occupancy_points.end(), source_corners.begin(), source_corners.end());
        }
        if (sample.support_bbox_world.IsInverted() && sample.source_bounds_world.IsInverted()) {
            support_points.push_back(sample.centroid_world);
            occupancy_points.push_back(sample.centroid_world);
        }

        const double pressure_weight = sample.force_world.Length();
        pressure_weight_sum += pressure_weight;
        pressure_center_sum += sample.pressure_center_world * pressure_weight;
    }

    if (patch.measure_area > 0) {
        patch.centroid_world = centroid_sum / patch.measure_area;
    }
    patch.pressure_center_world =
        pressure_weight_sum > 0 ? pressure_center_sum / pressure_weight_sum : patch.centroid_world;
    patch.mean_normal_world = SafeNormalized(normal_sum, samples[patch.sample_indices.front()].normal_world);
    const auto layered_cells =
        BuildPatchPlaneLayeredCells(patch, samples, patch.centroid_world, patch.mean_normal_world, sheet_scale);
    patch.support_cells =
        CollapseLayeredCellsToSheetCells(layered_cells, patch.centroid_world, patch.mean_normal_world, sheet_scale);
    const auto topology_stats = AnnotateSupportSheetTopology(layered_cells.size(), patch.support_cells);
    patch.layered_cell_count = topology_stats.layered_cell_count;
    patch.max_layer_count_per_cell = topology_stats.max_layer_count_per_cell;
    patch.largest_connected_sheet_cells = topology_stats.largest_connected_sheet_cells;
    patch.mean_layer_count_per_cell = topology_stats.mean_layer_count_per_cell;
    patch.sheet_fill_ratio = topology_stats.fill_ratio;

    const auto support_metrics = ComputeProjectedHullMetrics(support_points, patch.centroid_world, patch.mean_normal_world);
    const auto occupancy_metrics =
        ComputeProjectedHullMetrics(occupancy_points, patch.centroid_world, patch.mean_normal_world);
    const double base_recovery_radius =
        ResolvePatchRecoveryRadius(patch.measure_area, patch.support_seed_count, sheet_scale);
    ChVector3d column_tangent_u;
    ChVector3d column_tangent_v;
    BuildOrthonormalBasis(patch.mean_normal_world, column_tangent_u, column_tangent_v);

    double recovered_area = patch.measure_area;
    double recovery_radius_u = base_recovery_radius;
    double recovery_radius_v = base_recovery_radius;
    double min_u = 0;
    double max_u = 0;
    double min_v = 0;
    double max_v = 0;
    if (support_metrics.valid || occupancy_metrics.valid) {
        const auto& primary_metrics = support_metrics.valid ? support_metrics : occupancy_metrics;
        min_u = primary_metrics.min_u;
        max_u = primary_metrics.max_u;
        min_v = primary_metrics.min_v;
        max_v = primary_metrics.max_v;

        const double primary_span_u = std::max(primary_metrics.max_u - primary_metrics.min_u, 0.0);
        const double primary_span_v = std::max(primary_metrics.max_v - primary_metrics.min_v, 0.0);
        const bool u_major_primary = primary_span_u >= primary_span_v;
        const double primary_minor = std::max(std::min(primary_span_u, primary_span_v), 1.0e-12);
        const double primary_major = std::max(primary_span_u, primary_span_v);
        const double primary_aspect = primary_major / primary_minor;

        if (support_metrics.valid && occupancy_metrics.valid && primary_aspect > 1.25) {
            if (u_major_primary) {
                min_v = std::min(min_v, occupancy_metrics.min_v);
                max_v = std::max(max_v, occupancy_metrics.max_v);
            } else {
                min_u = std::min(min_u, occupancy_metrics.min_u);
                max_u = std::max(max_u, occupancy_metrics.max_u);
            }
        }

        const double span_u = std::max(max_u - min_u, 0.0);
        const double span_v = std::max(max_v - min_v, 0.0);
        const double minor_span = std::max(std::min(span_u, span_v), 1.0e-12);
        const double major_span = std::max(span_u, span_v);
        const double aspect_ratio = major_span / minor_span;
        const double slenderness = std::max(0.0, std::min((aspect_ratio - 1.0) / 2.0, 1.0));
        const double major_scale = 0.125;
        const double minor_scale = 0.125 + 0.875 * slenderness;

        if (span_u >= span_v) {
            recovery_radius_u = major_scale * base_recovery_radius;
            recovery_radius_v = minor_scale * base_recovery_radius;
        } else {
            recovery_radius_u = minor_scale * base_recovery_radius;
            recovery_radius_v = major_scale * base_recovery_radius;
        }

        if (occupancy_metrics.valid) {
            min_u = std::max(min_u, occupancy_metrics.min_u - recovery_radius_u);
            max_u = std::min(max_u, occupancy_metrics.max_u + recovery_radius_u);
            min_v = std::max(min_v, occupancy_metrics.min_v - recovery_radius_v);
            max_v = std::min(max_v, occupancy_metrics.max_v + recovery_radius_v);

            if (max_u < min_u) {
                const double mid_u = 0.5 * (min_u + max_u);
                min_u = mid_u;
                max_u = mid_u;
            }
            if (max_v < min_v) {
                const double mid_v = 0.5 * (min_v + max_v);
                min_v = mid_v;
                max_v = mid_v;
            }
        }

        if (!column_centers.empty()) {
            std::vector<double> column_u;
            std::vector<double> column_v;
            column_u.reserve(column_centers.size());
            column_v.reserve(column_centers.size());

            double column_min_u = std::numeric_limits<double>::infinity();
            double column_max_u = -std::numeric_limits<double>::infinity();
            double column_min_v = std::numeric_limits<double>::infinity();
            double column_max_v = -std::numeric_limits<double>::infinity();

            for (const auto& center_world : column_centers) {
                const ChVector3d rel = center_world - patch.centroid_world;
                const double u = Vdot(rel, column_tangent_u);
                const double v = Vdot(rel, column_tangent_v);
                column_u.push_back(u);
                column_v.push_back(v);
                column_min_u = std::min(column_min_u, u);
                column_max_u = std::max(column_max_u, u);
                column_min_v = std::min(column_min_v, v);
                column_max_v = std::max(column_max_v, v);
            }

            const double half_step_u = EstimateAxisHalfStep(column_u, recovery_radius_u);
            const double half_step_v = EstimateAxisHalfStep(column_v, recovery_radius_v);
            const double centroid_min_u = column_min_u - half_step_u;
            const double centroid_max_u = column_max_u + half_step_u;
            const double centroid_min_v = column_min_v - half_step_v;
            const double centroid_max_v = column_max_v + half_step_v;

            min_u = std::max(min_u, centroid_min_u);
            max_u = std::min(max_u, centroid_max_u);
            min_v = std::max(min_v, centroid_min_v);
            max_v = std::min(max_v, centroid_max_v);

            if (max_u < min_u) {
                min_u = centroid_min_u;
                max_u = centroid_max_u;
            }
            if (max_v < min_v) {
                min_v = centroid_min_v;
                max_v = centroid_max_v;
            }
        }

        const double core_area = support_metrics.valid ? support_metrics.area : occupancy_metrics.area;
        recovered_area = core_area + 2.0 * recovery_radius_u * span_v + 2.0 * recovery_radius_v * span_u +
                         CH_PI * recovery_radius_u * recovery_radius_v;
    } else if (base_recovery_radius > 0) {
        recovered_area = CH_PI * base_recovery_radius * base_recovery_radius;
    }
    if (support_metrics.valid || occupancy_metrics.valid) {
        patch.support_footprint.origin_world = patch.centroid_world;
        BuildOrthonormalBasis(patch.mean_normal_world, patch.support_footprint.tangent_u_world,
                              patch.support_footprint.tangent_v_world);
        patch.support_footprint.polygon_uv = {ChVector2d(min_u - recovery_radius_u, min_v - recovery_radius_v),
                                              ChVector2d(max_u + recovery_radius_u, min_v - recovery_radius_v),
                                              ChVector2d(max_u + recovery_radius_u, max_v + recovery_radius_v),
                                              ChVector2d(min_u - recovery_radius_u, max_v + recovery_radius_v)};
        patch.support_footprint.area = PolygonArea(patch.support_footprint.polygon_uv);
        patch.support_footprint.centroid_uv = PolygonCentroid(patch.support_footprint.polygon_uv);
    }
    if (!patch.support_footprint.HasPolygon()) {
        const double fallback_area = recovered_area > 1.0e-16 ? recovered_area : patch.measure_area;
        patch.support_footprint = MakeSquareFootprint(patch.centroid_world, patch.mean_normal_world, fallback_area);
    }
    patch.footprint_area =
        patch.support_footprint.area > 1.0e-16 ? patch.support_footprint.area : std::max(recovered_area, 0.0);
    UpdateBBoxFromFootprint(patch.support_footprint, patch.support_bbox_world);
    if (patch.support_bbox_world.IsInverted()) {
        ExpandProjectedSupportBBox(patch.support_bbox_world, patch.centroid_world, patch.mean_normal_world, support_points,
                                   std::max(patch.footprint_area, std::max(recovered_area, 0.0)));
    }
}

bool PatchesBelongTogether(const ChSDFSheetPatch& a,
                           const ChSDFSheetPatch& b,
                           double patch_connection_radius,
                           double patch_normal_cosine) {
    const double cosine = NormalCosine(a.mean_normal_world, b.mean_normal_world);
    if (patch_normal_cosine > -1 && cosine < std::min(patch_normal_cosine, std::cos(CH_PI / 4.0))) {
        return false;
    }

    const double radius_a = a.footprint_area > 0 ? std::sqrt(a.footprint_area / CH_PI) : 0.0;
    const double radius_b = b.footprint_area > 0 ? std::sqrt(b.footprint_area / CH_PI) : 0.0;
    const double bbox_gap = AABBDistance(a.support_bbox_world, b.support_bbox_world);
    if (bbox_gap <= patch_connection_radius + 0.25 * (radius_a + radius_b) + 1.0e-12) {
        return true;
    }

    const ChVector3d n_bar = SafeNormalized(a.mean_normal_world + b.mean_normal_world, a.mean_normal_world);
    const ChVector3d delta = b.centroid_world - a.centroid_world;
    const double delta_n = std::abs(Vdot(n_bar, delta));
    const double delta_t = (delta - n_bar * Vdot(n_bar, delta)).Length();

    const double min_area = std::min(a.measure_area, b.measure_area);
    const double max_area = std::max(a.measure_area, b.measure_area);
    const bool small_fragment_bridge =
        max_area > 0 && min_area <= 0.35 * max_area &&
        bbox_gap <= 2.0 * patch_connection_radius + 0.5 * (radius_a + radius_b) + 1.0e-12 &&
        delta_n <= patch_connection_radius + 0.25 * (radius_a + radius_b) + 1.0e-12;
    if (small_fragment_bridge) {
        return true;
    }

    return delta_t <= patch_connection_radius + radius_a + radius_b + 1.0e-12 &&
           delta_n <= 0.5 * patch_connection_radius + 0.25 * (radius_a + radius_b) + 1.0e-12;
}

std::vector<ChSDFSheetPatch> MergeV2Patches(const std::vector<ChSDFSheetPatch>& input_patches,
                                            const std::vector<ChSDFSheetFiberSample>& samples,
                                            double sheet_scale,
                                            double patch_connection_radius,
                                            double patch_normal_cosine) {
    if (input_patches.size() <= 1) {
        return input_patches;
    }

    UnionFind uf(input_patches.size());
    for (std::size_t i = 0; i < input_patches.size(); ++i) {
        for (std::size_t j = i + 1; j < input_patches.size(); ++j) {
            if (PatchesBelongTogether(input_patches[i], input_patches[j], patch_connection_radius, patch_normal_cosine)) {
                uf.Unite(i, j);
            }
        }
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> merged_indices;
    merged_indices.reserve(input_patches.size());
    for (std::size_t i = 0; i < input_patches.size(); ++i) {
        auto& indices = merged_indices[uf.Find(i)];
        indices.insert(indices.end(), input_patches[i].sample_indices.begin(), input_patches[i].sample_indices.end());
    }

    std::vector<ChSDFSheetPatch> merged_patches;
    merged_patches.reserve(merged_indices.size());
    for (auto& item : merged_indices) {
        auto& indices = item.second;
        std::sort(indices.begin(), indices.end());
        indices.erase(std::unique(indices.begin(), indices.end()), indices.end());

        ChSDFSheetPatch patch;
        patch.patch_id = merged_patches.size();
        patch.sample_indices = std::move(indices);
        FinalizeV2Patch(patch, samples, sheet_scale);
        merged_patches.push_back(std::move(patch));
    }

    std::stable_sort(merged_patches.begin(), merged_patches.end(),
                     [](const ChSDFSheetPatch& a, const ChSDFSheetPatch& b) { return a.patch_id < b.patch_id; });
    return merged_patches;
}

std::vector<ChSDFSheetPatch> BuildLegacyPatches(const std::vector<ChSDFSheetFiberSample>& samples,
                                                int neighbor_mode,
                                                const ChVector3d& fallback_normal,
                                                double normal_cosine_threshold) {
    std::vector<ChSDFSheetPatch> patches;
    if (samples.empty()) {
        return patches;
    }

    std::unordered_map<FiberKey, std::vector<std::size_t>, FiberKeyHash> sample_by_coord;
    sample_by_coord.reserve(samples.size());
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto& sample = samples[i];
        sample_by_coord[FiberKey{sample.lattice_coord.x(), sample.lattice_coord.y()}].push_back(i);
    }

    const auto neighbor_offsets = BuildPatchNeighborOffsets(neighbor_mode);
    std::vector<char> visited(samples.size(), 0);

    for (std::size_t seed = 0; seed < samples.size(); ++seed) {
        if (visited[seed]) {
            continue;
        }

        ChSDFSheetPatch patch;
        patch.patch_id = patches.size();

        std::queue<std::size_t> frontier;
        frontier.push(seed);
        visited[seed] = 1;

        while (!frontier.empty()) {
            const auto current = frontier.front();
            frontier.pop();
            patch.sample_indices.push_back(current);

            const auto& coord = samples[current].lattice_coord;
            std::vector<FiberKey> candidate_coords;
            candidate_coords.reserve(neighbor_offsets.size() + 1);
            candidate_coords.push_back(FiberKey{coord.x(), coord.y()});
            for (const auto& offset : neighbor_offsets) {
                candidate_coords.push_back(FiberKey{coord.x() + offset.x(), coord.y() + offset.y()});
            }

            for (const auto& neighbor : candidate_coords) {
                auto it = sample_by_coord.find(neighbor);
                if (it == sample_by_coord.end()) {
                    continue;
                }

                for (const auto neighbor_index : it->second) {
                    if (visited[neighbor_index]) {
                        continue;
                    }

                    const double cosine =
                        NormalCosine(samples[current].normal_world, samples[neighbor_index].normal_world);
                    if (normal_cosine_threshold > -1 && cosine < normal_cosine_threshold) {
                        continue;
                    }

                    visited[neighbor_index] = 1;
                    frontier.push(neighbor_index);
                }
            }
        }

        FinalizeLegacyPatch(patch, samples, fallback_normal);
        patches.push_back(std::move(patch));
    }

    return patches;
}

ChSDFSheetSeed BuildSeed(const ChSDFBrickPairWrenchSample& sample,
                         std::size_t sample_index,
                         std::size_t region_id,
                         const ChSDFSheetCollapseSettings& settings,
                         const ChVector3d& fallback_normal) {
    ChSDFSheetSeed seed;
    if (!sample.active || sample.quadrature_area <= std::max(settings.min_sheet_sample_area, 0.0)) {
        return seed;
    }

    const ChVector3d seed_normal = SafeNormalized(sample.region_sample.contact_normal_world, fallback_normal);
    if (seed_normal.Length2() <= 1.0e-24 || !IsFiniteVector(sample.region_sample.point_world)) {
        return seed;
    }

    seed.region_id = region_id;
    seed.source_sample_index = sample_index;
    seed.source_carrier_axis = sample.region_sample.carrier_axis;
    seed.source_coord = sample.region_sample.coord;
    seed.band_point_world = sample.region_sample.point_world;
    seed.seed_world = seed.band_point_world;
    seed.seed_normal_world = seed_normal;
    seed.measure_area = sample.quadrature_area;
    seed.pressure = sample.pressure;
    seed.h_value = sample.quadrature_point.h_value;
    seed.force_world = sample.force_world;
    seed.source_bounds_world += seed.band_point_world;

    const double normal_gap = sample.region_sample.normal_gap;
    if (std::abs(normal_gap) > 1.0e-12) {
        const double grad_norm = std::abs(sample.quadrature_point.h_value / (-normal_gap));
        if (std::isfinite(grad_norm) && grad_norm >= settings.min_gradient_norm) {
            seed.grad_h_world = -seed_normal * grad_norm;
            const ChVector3d recovered =
                seed.band_point_world - seed.grad_h_world * (seed.h_value / (grad_norm * grad_norm));
            if (IsFiniteVector(recovered)) {
                seed.seed_world = recovered;
            }
        }
    }

    if (seed.grad_h_world.Length2() <= 1.0e-24 && sample.quadrature_point.normal_world.Length2() > 1.0e-24) {
        seed.seed_normal_world = SafeNormalized(sample.quadrature_point.normal_world, seed_normal);
    }

    if (!IsFiniteVector(seed.seed_world) || !IsFiniteVector(seed.seed_normal_world) || !std::isfinite(seed.measure_area) ||
        seed.measure_area <= 0) {
        return seed;
    }

    seed.source_bounds_world += seed.seed_world;
    seed.valid = true;
    return seed;
}

std::vector<ChSDFSheetSeed> BuildSeeds(const ChSDFBrickPairWrenchResult& band_region,
                                       const ChSDFSheetCollapseSettings& settings,
                                       const ChVector3d& fallback_normal) {
    std::vector<ChSDFSheetSeed> seeds;
    seeds.reserve(band_region.samples.size());

    for (std::size_t i = 0; i < band_region.samples.size(); ++i) {
        auto seed = BuildSeed(band_region.samples[i], i, band_region.region.region_id, settings, fallback_normal);
        if (seed.valid) {
            seeds.push_back(std::move(seed));
        }
    }

    return seeds;
}

bool SeedsBelongToSameFiber(const ChSDFSheetSeed& a,
                            const ChSDFSheetSeed& b,
                            double normal_cosine_threshold,
                            double tangent_tolerance,
                            double plane_tolerance) {
    const double cosine = NormalCosine(a.seed_normal_world, b.seed_normal_world);
    if (normal_cosine_threshold > -1 && cosine < normal_cosine_threshold) {
        return false;
    }

    const ChVector3d n_bar = SafeNormalized(a.seed_normal_world + b.seed_normal_world, a.seed_normal_world);
    const ChVector3d delta = b.seed_world - a.seed_world;
    const double delta_n = std::abs(Vdot(n_bar, delta));
    const ChVector3d delta_t_vec = delta - n_bar * Vdot(n_bar, delta);
    const double delta_t = delta_t_vec.Length();

    return delta_t <= tangent_tolerance + 1.0e-12 && delta_n <= plane_tolerance + 1.0e-12;
}

std::vector<std::vector<std::size_t>> BuildFiberClusters(const std::vector<ChSDFSheetSeed>& seeds,
                                                         double normal_cosine_threshold,
                                                         double tangent_tolerance,
                                                         double plane_tolerance) {
    std::vector<std::vector<std::size_t>> clusters;
    if (seeds.empty()) {
        return clusters;
    }

    const double search_radius =
        std::sqrt(tangent_tolerance * tangent_tolerance + plane_tolerance * plane_tolerance);
    const double cell_size = std::max(search_radius, 1.0e-6);

    std::unordered_map<SpatialHashKey, std::vector<std::size_t>, SpatialHashKeyHash> buckets;
    buckets.reserve(seeds.size());
    for (std::size_t i = 0; i < seeds.size(); ++i) {
        buckets[MakeSpatialHashKey(seeds[i].seed_world, cell_size)].push_back(i);
    }

    UnionFind uf(seeds.size());
    for (std::size_t i = 0; i < seeds.size(); ++i) {
        const auto base_key = MakeSpatialHashKey(seeds[i].seed_world, cell_size);
        for (int dz = -1; dz <= 1; ++dz) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    const SpatialHashKey neighbor_key{base_key.x + dx, base_key.y + dy, base_key.z + dz};
                    auto it = buckets.find(neighbor_key);
                    if (it == buckets.end()) {
                        continue;
                    }

                    for (const auto j : it->second) {
                        if (j <= i) {
                            continue;
                        }
                        if (SeedsBelongToSameFiber(seeds[i], seeds[j], normal_cosine_threshold, tangent_tolerance,
                                                   plane_tolerance)) {
                            uf.Unite(i, j);
                        }
                    }
                }
            }
        }
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> cluster_map;
    cluster_map.reserve(seeds.size());
    for (std::size_t i = 0; i < seeds.size(); ++i) {
        cluster_map[uf.Find(i)].push_back(i);
    }

    clusters.reserve(cluster_map.size());
    for (auto& item : cluster_map) {
        auto& cluster = item.second;
        std::sort(cluster.begin(), cluster.end());
        clusters.push_back(std::move(cluster));
    }

    std::stable_sort(clusters.begin(), clusters.end(),
                     [](const auto& a, const auto& b) { return a.front() < b.front(); });
    return clusters;
}

ChSDFSheetFiberSample CollapseFiberCluster(const std::vector<ChSDFSheetSeed>& seeds,
                                           const std::vector<std::size_t>& cluster,
                                           std::size_t region_id,
                                           int dominant_axis) {
    ChSDFSheetFiberSample sample;
    sample.region_id = region_id;
    sample.dominant_axis = dominant_axis;

    if (cluster.empty()) {
        return sample;
    }

    ChVector3d centroid_sum = VNULL;
    ChVector3d pressure_center_sum = VNULL;
    ChVector3d normal_sum = VNULL;
    ChVector3d force_sum = VNULL;
    double measure_area = 0;
    double pressure_weight_sum = 0;
    std::vector<ChVector3d> seed_points;
    seed_points.reserve(cluster.size());

    for (const auto seed_index : cluster) {
        const auto& seed = seeds[seed_index];
        measure_area += seed.measure_area;
        centroid_sum += seed.seed_world * seed.measure_area;
        normal_sum += seed.seed_normal_world * seed.measure_area;
        force_sum += seed.force_world;
        sample.source_bounds_world += seed.band_point_world;
        sample.source_sample_indices.push_back(seed.source_sample_index);
        sample.support_seed_count++;
        sample.support_evidence.push_back(
            ChSDFSheetSupportEvidence{seed.source_carrier_axis, seed.source_coord, seed.band_point_world, seed.seed_world,
                                      seed.seed_normal_world, seed.measure_area, seed.source_sample_index});
        seed_points.push_back(seed.seed_world);

        const double pressure_weight = seed.force_world.Length();
        pressure_weight_sum += pressure_weight;
        pressure_center_sum += seed.seed_world * pressure_weight;
    }

    if (measure_area <= 1.0e-16) {
        return sample;
    }

    sample.measure_area = measure_area;
    sample.centroid_world = centroid_sum / measure_area;
    sample.pressure_center_world =
        pressure_weight_sum > 0 ? pressure_center_sum / pressure_weight_sum : sample.centroid_world;
    sample.normal_world = SafeNormalized(normal_sum, seeds[cluster.front()].seed_normal_world);
    sample.force_world = force_sum;

    sample.support_bbox_world = ChAABB();
    for (const auto& point : seed_points) {
        sample.support_bbox_world += point;
    }

    ChVector3d tangent_u;
    ChVector3d tangent_v;
    BuildOrthonormalBasis(sample.normal_world, tangent_u, tangent_v);

    double min_u = std::numeric_limits<double>::infinity();
    double max_u = -std::numeric_limits<double>::infinity();
    double min_v = std::numeric_limits<double>::infinity();
    double max_v = -std::numeric_limits<double>::infinity();

    for (const auto& point : seed_points) {
        const ChVector3d rel = point - sample.centroid_world;
        const double u = Vdot(rel, tangent_u);
        const double v = Vdot(rel, tangent_v);
        min_u = std::min(min_u, u);
        max_u = std::max(max_u, u);
        min_v = std::min(min_v, v);
        max_v = std::max(max_v, v);
    }

    const double span_u = std::isfinite(min_u) ? (max_u - min_u) : 0.0;
    const double span_v = std::isfinite(min_v) ? (max_v - min_v) : 0.0;

    double half_u = 0.5 * std::max(span_u, 0.0);
    double half_v = 0.5 * std::max(span_v, 0.0);
    const double min_half = std::sqrt(std::max(sample.measure_area, 0.0) / CH_PI);

    if (half_u <= 1.0e-12 && half_v <= 1.0e-12) {
        half_u = min_half;
        half_v = min_half;
    } else {
        const double full_u = std::max(2.0 * half_u, 1.0e-12);
        const double full_v = std::max(2.0 * half_v, 1.0e-12);
        const double observed_rect_area = full_u * full_v;
        const double target_area = std::max(sample.measure_area, observed_rect_area);

        if (full_u >= full_v) {
            half_v = std::max(half_v, 0.5 * target_area / full_u);
        } else {
            half_u = std::max(half_u, 0.5 * target_area / full_v);
        }
    }

    sample.support_footprint.origin_world = sample.centroid_world;
    sample.support_footprint.tangent_u_world = tangent_u;
    sample.support_footprint.tangent_v_world = tangent_v;
    sample.support_footprint.polygon_uv = {ChVector2d(-half_u, -half_v), ChVector2d(half_u, -half_v),
                                           ChVector2d(half_u, half_v), ChVector2d(-half_u, half_v)};
    sample.support_footprint.area = PolygonArea(sample.support_footprint.polygon_uv);
    sample.support_footprint.centroid_uv = PolygonCentroid(sample.support_footprint.polygon_uv);
    sample.footprint_area =
        sample.support_footprint.area > 1.0e-16 ? sample.support_footprint.area : sample.measure_area;
    UpdateBBoxFromFootprint(sample.support_footprint, sample.support_bbox_world);
    for (const auto& point : seed_points) {
        sample.support_bbox_world += point;
    }

    return sample;
}

double EquivalentSupportRadius(const ChSDFSheetFiberSample& sample) {
    const double area = std::max(sample.footprint_area, sample.measure_area);
    return area > 0 ? std::sqrt(area / CH_PI) : 0.0;
}

double AABBDistance(const ChAABB& a, const ChAABB& b) {
    const auto axis_gap = [](double a_min, double a_max, double b_min, double b_max) {
        if (a_max < b_min) {
            return b_min - a_max;
        }
        if (b_max < a_min) {
            return a_min - b_max;
        }
        return 0.0;
    };

    const double dx = axis_gap(a.min.x(), a.max.x(), b.min.x(), b.max.x());
    const double dy = axis_gap(a.min.y(), a.max.y(), b.min.y(), b.max.y());
    const double dz = axis_gap(a.min.z(), a.max.z(), b.min.z(), b.max.z());
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

bool FibersBelongToSamePatch(const ChSDFSheetFiberSample& a,
                             const ChSDFSheetFiberSample& b,
                             double patch_connection_radius,
                             double patch_normal_cosine) {
    const double cosine = NormalCosine(a.normal_world, b.normal_world);
    const double radius_a = EquivalentSupportRadius(a);
    const double radius_b = EquivalentSupportRadius(b);
    const double adaptive_tangent_tol = patch_connection_radius + radius_a + radius_b;
    const double adaptive_plane_tol = 0.5 * patch_connection_radius + 0.25 * (radius_a + radius_b);

    const ChVector3d n_bar = SafeNormalized(a.normal_world + b.normal_world, a.normal_world);
    const ChVector3d delta = b.centroid_world - a.centroid_world;
    const double delta_n = std::abs(Vdot(n_bar, delta));
    const double delta_t = (delta - n_bar * Vdot(n_bar, delta)).Length();
    const double bbox_gap = AABBDistance(a.support_bbox_world, b.support_bbox_world);

    const bool centroid_adjacent = delta_t <= adaptive_tangent_tol + 1.0e-12 &&
                                   delta_n <= adaptive_plane_tol + 1.0e-12;
    const bool bbox_adjacent = bbox_gap <= (patch_connection_radius + 0.25 * (radius_a + radius_b) + 1.0e-12);

    if (!centroid_adjacent && !bbox_adjacent) {
        return false;
    }

    const double strict_cosine = patch_normal_cosine;
    const double overlap_cosine = std::min(patch_normal_cosine, std::cos(CH_PI / 4.0));
    const double required_cosine = bbox_adjacent ? overlap_cosine : strict_cosine;
    if (required_cosine > -1 && cosine < required_cosine) {
        return false;
    }

    return true;
}

std::vector<ChSDFSheetPatch> BuildPatchGraph(const std::vector<ChSDFSheetFiberSample>& samples,
                                             double sheet_scale,
                                             double patch_connection_radius,
                                             double patch_normal_cosine) {
    std::vector<ChSDFSheetPatch> patches;
    if (samples.empty()) {
        return patches;
    }

    double max_support_radius = 0;
    for (const auto& sample : samples) {
        max_support_radius = std::max(max_support_radius, EquivalentSupportRadius(sample));
    }

    const double cell_size = std::max(patch_connection_radius + 2.0 * max_support_radius, 1.0e-6);
    std::unordered_map<SpatialHashKey, std::vector<std::size_t>, SpatialHashKeyHash> buckets;
    buckets.reserve(samples.size());
    for (std::size_t i = 0; i < samples.size(); ++i) {
        buckets[MakeSpatialHashKey(samples[i].centroid_world, cell_size)].push_back(i);
    }

    UnionFind uf(samples.size());
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto base_key = MakeSpatialHashKey(samples[i].centroid_world, cell_size);
        for (int dz = -1; dz <= 1; ++dz) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    const SpatialHashKey neighbor_key{base_key.x + dx, base_key.y + dy, base_key.z + dz};
                    auto it = buckets.find(neighbor_key);
                    if (it == buckets.end()) {
                        continue;
                    }

                    for (const auto j : it->second) {
                        if (j <= i) {
                            continue;
                        }

                        if (FibersBelongToSamePatch(samples[i], samples[j], patch_connection_radius,
                                                    patch_normal_cosine)) {
                            uf.Unite(i, j);
                        }
                    }
                }
            }
        }
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> patch_map;
    patch_map.reserve(samples.size());
    for (std::size_t i = 0; i < samples.size(); ++i) {
        patch_map[uf.Find(i)].push_back(i);
    }

    patches.reserve(patch_map.size());
    for (auto& item : patch_map) {
        auto& indices = item.second;
        std::sort(indices.begin(), indices.end());

        ChSDFSheetPatch patch;
        patch.patch_id = patches.size();
        patch.sample_indices = std::move(indices);
        FinalizeV2Patch(patch, samples, sheet_scale);
        patches.push_back(std::move(patch));
    }

    patches = MergeV2Patches(patches, samples, sheet_scale, patch_connection_radius, patch_normal_cosine);
    std::stable_sort(patches.begin(), patches.end(),
                     [](const ChSDFSheetPatch& a, const ChSDFSheetPatch& b) { return a.patch_id < b.patch_id; });
    return patches;
}

double ComputeNormalSpread(const std::vector<ChSDFSheetFiberSample>& samples, const ChVector3d& mean_normal_world) {
    if (samples.empty()) {
        return 0.0;
    }

    double weighted_error_sum = 0;
    double weight_sum = 0;
    for (const auto& sample : samples) {
        const double cosine = std::max(-1.0, std::min(1.0, NormalCosine(sample.normal_world, mean_normal_world)));
        const double error = 1.0 - cosine;
        weighted_error_sum += sample.measure_area * error * error;
        weight_sum += sample.measure_area;
    }
    return weight_sum > 0 ? std::sqrt(weighted_error_sum / weight_sum) : 0.0;
}

bool ShouldFallbackToDominantAxis(const ChSDFSheetRegion& region, const ChSDFSheetCollapseSettings& settings) {
    if (!settings.allow_dominant_axis_fallback || !region.HasSamples()) {
        return false;
    }

    if (settings.min_sheet_area_ratio > 0 && region.area_ratio < settings.min_sheet_area_ratio) {
        return true;
    }
    if (settings.max_sheet_area_ratio > 0 && region.area_ratio > settings.max_sheet_area_ratio) {
        return true;
    }
    if (settings.max_patch_count_before_fallback > 0 && region.patch_count > settings.max_patch_count_before_fallback) {
        return true;
    }
    if (settings.max_fiber_count_before_fallback > 0 && region.fiber_count > settings.max_fiber_count_before_fallback) {
        return true;
    }
    return false;
}

ChSDFSheetRegion BuildRegionLegacyDominantAxis(const ChSDFBrickPairWrenchResult& band_region,
                                               const ChSDFSheetCollapseSettings& settings,
                                               bool use_projection) {
    ChSDFSheetRegion region;
    region.region_id = band_region.region.region_id;
    region.persistent_id = band_region.region.persistent_id;

    if (!settings.enable || !band_region.HasActiveContact() || band_region.samples.empty()) {
        return region;
    }

    const ChVector3d fallback_normal = ResolveRegionNormal(band_region);
    region.dominant_axis = ChSDFSheetBuilder::ComputeDominantAxis(fallback_normal);
    const double lateral_tolerance = ResolveSheetScale(band_region, settings);
    ChVector3d tangent_u;
    ChVector3d tangent_v;
    BuildOrthonormalBasis(fallback_normal, tangent_u, tangent_v);

    ChVector3d projection_origin = band_region.region.centroid_world;
    if (projection_origin.Length2() <= 1.0e-24) {
        double weight_sum = 0;
        projection_origin = VNULL;
        for (const auto& sample : band_region.samples) {
            if (!sample.active || sample.quadrature_area <= 1.0e-16) {
                continue;
            }
            projection_origin += sample.region_sample.point_world * sample.quadrature_area;
            weight_sum += sample.quadrature_area;
        }
        if (weight_sum > 0) {
            projection_origin /= weight_sum;
        }
    }

    std::unordered_map<FiberKey, std::vector<LegacyFiberAccumulator>, FiberKeyHash> fibers;
    fibers.reserve(band_region.samples.size());

    for (std::size_t sample_index = 0; sample_index < band_region.samples.size(); ++sample_index) {
        const auto& sample = band_region.samples[sample_index];
        if (!sample.active || sample.quadrature_area <= std::max(settings.min_sheet_sample_area, 0.0)) {
            continue;
        }

        const ChVector3d sample_normal = SafeNormalized(sample.region_sample.contact_normal_world, fallback_normal);
        ChVector2i lattice_coord =
            ChSDFSheetBuilder::ComputeLatticeCoord(sample.region_sample.coord, region.dominant_axis);
        if (use_projection) {
            lattice_coord = QuantizeProjectedCoord(sample.region_sample.point_world, projection_origin, fallback_normal,
                                                  sample_normal, tangent_u, tangent_v, lateral_tolerance);
        }

        const FiberKey key{lattice_coord.x(), lattice_coord.y()};
        auto& bucket = fibers[key];

        LegacyFiberAccumulator* fiber = nullptr;
        for (auto& candidate : bucket) {
            const ChVector3d candidate_normal =
                candidate.normal_sum.Length2() > 1.0e-24 ? SafeNormalized(candidate.normal_sum, fallback_normal)
                                                         : fallback_normal;
            const double cosine = NormalCosine(sample_normal, candidate_normal);
            if (settings.fiber_normal_cosine <= -1 || cosine >= settings.fiber_normal_cosine) {
                fiber = &candidate;
                break;
            }
        }
        if (!fiber) {
            bucket.emplace_back();
            fiber = &bucket.back();
            fiber->dominant_axis = region.dominant_axis;
            fiber->lattice_coord = lattice_coord;
        }

        fiber->measure_area += sample.quadrature_area;
        fiber->area_centroid_sum += sample.region_sample.point_world * sample.quadrature_area;
        fiber->normal_sum += sample.region_sample.contact_normal_world * sample.quadrature_area;
        fiber->force_sum += sample.force_world;
        fiber->bounds_world += sample.region_sample.point_world;
        fiber->source_sample_indices.push_back(sample_index);

        const double pressure_weight = std::max(sample.pressure * sample.quadrature_area, 0.0);
        fiber->pressure_weight += pressure_weight;
        fiber->pressure_center_sum += sample.region_sample.point_world * pressure_weight;
    }

    if (fibers.empty()) {
        return region;
    }

    const double footprint_area = lateral_tolerance * lateral_tolerance;
    region.samples.reserve(fibers.size());

    ChVector3d centroid_sum = VNULL;
    ChVector3d pressure_center_sum = VNULL;
    ChVector3d normal_sum = VNULL;
    double pressure_weight_sum = 0;
    std::size_t support_seed_sum = 0;

    std::vector<std::pair<FiberKey, LegacyFiberAccumulator>> ordered_fibers;
    ordered_fibers.reserve(fibers.size());
    for (auto& item : fibers) {
        for (auto& fiber : item.second) {
            ordered_fibers.push_back({item.first, std::move(fiber)});
        }
    }
    std::stable_sort(ordered_fibers.begin(), ordered_fibers.end(), [](const auto& a, const auto& b) {
        return a.first.a != b.first.a ? a.first.a < b.first.a : a.first.b < b.first.b;
    });

    for (const auto& item : ordered_fibers) {
        const auto& fiber_data = item.second;
        if (fiber_data.measure_area <= std::max(settings.min_sheet_sample_area, 0.0)) {
            continue;
        }

        ChSDFSheetFiberSample sample;
        sample.region_id = region.region_id;
        sample.fiber_id = region.samples.size();
        sample.dominant_axis = region.dominant_axis;
        sample.lattice_coord = fiber_data.lattice_coord;
        sample.measure_area = fiber_data.measure_area;
        sample.footprint_area = footprint_area;
        sample.support_seed_count = fiber_data.source_sample_indices.size();
        sample.centroid_world = fiber_data.area_centroid_sum / fiber_data.measure_area;
        sample.pressure_center_world =
            fiber_data.pressure_weight > 0 ? fiber_data.pressure_center_sum / fiber_data.pressure_weight
                                           : sample.centroid_world;
        sample.normal_world = SafeNormalized(fiber_data.normal_sum, fallback_normal);
        sample.force_world = fiber_data.force_sum;
        sample.source_bounds_world = fiber_data.bounds_world;
        sample.support_bbox_world += sample.centroid_world;
        sample.source_sample_indices = fiber_data.source_sample_indices;

        region.measure_area += sample.measure_area;
        region.footprint_area += sample.footprint_area;
        centroid_sum += sample.centroid_world * sample.measure_area;
        normal_sum += sample.normal_world * sample.measure_area;
        region.bounds_world += sample.source_bounds_world;
        region.support_bbox_world += sample.support_bbox_world;
        support_seed_sum += sample.support_seed_count;

        const double region_pressure_weight = sample.force_world.Length();
        pressure_weight_sum += region_pressure_weight;
        pressure_center_sum += sample.pressure_center_world * region_pressure_weight;

        region.samples.push_back(std::move(sample));
    }

    if (region.samples.empty()) {
        return region;
    }

    region.fiber_count = region.samples.size();
    region.centroid_world = centroid_sum / region.measure_area;
    region.mean_normal_world = SafeNormalized(normal_sum, fallback_normal);
    region.pressure_center_world =
        pressure_weight_sum > 0 ? pressure_center_sum / pressure_weight_sum : region.centroid_world;
    region.patches =
        BuildLegacyPatches(region.samples, settings.patch_neighbor_mode, region.mean_normal_world,
                           settings.fiber_normal_cosine);
    region.patch_count = region.patches.size();
    region.largest_patch_area = 0;
    for (const auto& patch : region.patches) {
        region.largest_patch_area = std::max(region.largest_patch_area, patch.measure_area);
    }
    region.area_ratio = band_region.active_area > 0 ? region.footprint_area / band_region.active_area : 0.0;
    region.mean_support_seed_count =
        region.samples.empty() ? 0.0 : static_cast<double>(support_seed_sum) / static_cast<double>(region.samples.size());
    region.normal_spread = ComputeNormalSpread(region.samples, region.mean_normal_world);
    return region;
}

ChSDFSheetRegion BuildRegionV2(const ChSDFBrickPairWrenchResult& band_region,
                               const ChSDFSheetCollapseSettings& settings) {
    ChSDFSheetRegion region;
    region.region_id = band_region.region.region_id;
    region.persistent_id = band_region.region.persistent_id;

    if (!settings.enable || !band_region.HasActiveContact() || band_region.samples.empty()) {
        return region;
    }

    const ChVector3d fallback_normal = ResolveRegionNormal(band_region);
    region.dominant_axis = ChSDFSheetBuilder::ComputeDominantAxis(fallback_normal);

    const auto seeds = BuildSeeds(band_region, settings, fallback_normal);
    if (seeds.empty()) {
        return region;
    }

    const double tangent_tolerance = ResolveFiberTangentTolerance(band_region, settings);
    const double plane_tolerance = ResolveFiberPlaneTolerance(band_region, settings);
    const auto clusters =
        BuildFiberClusters(seeds, settings.fiber_normal_cosine, tangent_tolerance, plane_tolerance);
    if (clusters.empty()) {
        return region;
    }

    ChVector3d centroid_sum = VNULL;
    ChVector3d pressure_center_sum = VNULL;
    ChVector3d normal_sum = VNULL;
    double pressure_weight_sum = 0;
    std::size_t support_seed_sum = 0;

    region.samples.reserve(clusters.size());
    for (const auto& cluster : clusters) {
        auto sample = CollapseFiberCluster(seeds, cluster, region.region_id, region.dominant_axis);
        if (!sample.HasSupport()) {
            continue;
        }

        sample.fiber_id = region.samples.size();
        region.measure_area += sample.measure_area;
        region.footprint_area += sample.footprint_area;
        centroid_sum += sample.centroid_world * sample.measure_area;
        normal_sum += sample.normal_world * sample.measure_area;
        region.bounds_world += sample.source_bounds_world;
        region.support_bbox_world += sample.support_bbox_world;
        support_seed_sum += sample.support_seed_count;

        const double region_pressure_weight = sample.force_world.Length();
        pressure_weight_sum += region_pressure_weight;
        pressure_center_sum += sample.pressure_center_world * region_pressure_weight;

        region.samples.push_back(std::move(sample));
    }

    if (region.samples.empty()) {
        return region;
    }

    region.fiber_count = region.samples.size();
    region.centroid_world = centroid_sum / region.measure_area;
    region.mean_normal_world = SafeNormalized(normal_sum, fallback_normal);
    region.pressure_center_world =
        pressure_weight_sum > 0 ? pressure_center_sum / pressure_weight_sum : region.centroid_world;
    region.area_ratio = band_region.active_area > 0 ? region.footprint_area / band_region.active_area : 0.0;
    region.mean_support_seed_count =
        static_cast<double>(support_seed_sum) / static_cast<double>(region.samples.size());
    region.normal_spread = ComputeNormalSpread(region.samples, region.mean_normal_world);
    const double sheet_scale = ResolveSheetScale(band_region, settings);
    region.patches =
        BuildPatchGraph(region.samples, sheet_scale, ResolvePatchConnectionRadius(band_region, settings),
                        settings.patch_normal_cosine);
    region.patch_count = region.patches.size();
    region.largest_patch_area = 0;
    region.support_bbox_world = ChAABB();
    region.footprint_area = 0;
    for (const auto& patch : region.patches) {
        region.largest_patch_area = std::max(region.largest_patch_area, patch.footprint_area);
        region.footprint_area += patch.footprint_area;
        region.support_bbox_world += patch.support_bbox_world;
    }

    return region;
}

}  // namespace

int ChSDFSheetBuilder::ComputeDominantAxis(const ChVector3d& normal_world) {
    const ChVector3d abs_normal(std::abs(normal_world.x()), std::abs(normal_world.y()), std::abs(normal_world.z()));
    if (abs_normal.x() >= abs_normal.y() && abs_normal.x() >= abs_normal.z()) {
        return 0;
    }
    if (abs_normal.y() >= abs_normal.x() && abs_normal.y() >= abs_normal.z()) {
        return 1;
    }
    return 2;
}

ChVector2i ChSDFSheetBuilder::ComputeLatticeCoord(const ChVector3i& coord, int dominant_axis) {
    switch (dominant_axis) {
        case 0:
            return ChVector2i(coord.y(), coord.z());
        case 1:
            return ChVector2i(coord.x(), coord.z());
        case 2:
        default:
            return ChVector2i(coord.x(), coord.y());
    }
}

ChSDFSheetRegion ChSDFSheetBuilder::BuildRegion(const ChSDFBrickPairWrenchResult& band_region,
                                                const ChSDFSheetCollapseSettings& settings) {
    if (!settings.use_step6_v2) {
        return BuildRegionLegacyDominantAxis(band_region, settings, settings.use_local_fiber_projection);
    }

    auto region = BuildRegionV2(band_region, settings);
    if (!region.HasSamples() && settings.allow_dominant_axis_fallback) {
        auto fallback = BuildRegionLegacyDominantAxis(band_region, settings, false);
        fallback.used_fallback = fallback.HasSamples();
        return fallback;
    }

    if (ShouldFallbackToDominantAxis(region, settings)) {
        auto fallback = BuildRegionLegacyDominantAxis(band_region, settings, false);
        fallback.used_fallback = fallback.HasSamples();
        return fallback;
    }

    return region;
}

std::vector<ChSDFSheetRegion> ChSDFSheetBuilder::BuildRegions(
    const std::vector<ChSDFBrickPairWrenchResult>& band_regions,
    const ChSDFSheetCollapseSettings& settings) {
    std::vector<ChSDFSheetRegion> regions;
    regions.reserve(band_regions.size());

    for (const auto& band_region : band_regions) {
        auto region = BuildRegion(band_region, settings);
        if (region.HasSamples()) {
            regions.push_back(std::move(region));
        }
    }

    return regions;
}

ChSDFSheetShapePairResult ChSDFSheetBuilder::BuildShapePair(const ChSDFShapePairContactResult& band_result,
                                                            const ChSDFSheetCollapseSettings& settings) {
    ChSDFSheetShapePairResult result;
    result.band_area = band_result.active_area;

    if (!settings.enable || !band_result.HasActiveContact()) {
        return result;
    }

    ChVector3d sheet_center_sum = VNULL;
    ChVector3d pressure_center_sum = VNULL;
    double pressure_weight_sum = 0;
    double support_seed_sum = 0;

    result.regions.reserve(band_result.regions.size());

    for (const auto& band_region : band_result.regions) {
        if (!band_region.HasActiveContact()) {
            continue;
        }

        const double spacing = ResolveRegionSpacing(band_region);
        if (spacing > 0) {
            result.occupied_area += static_cast<double>(band_region.active_samples) * spacing * spacing;
        }
        result.occupied_cells += band_region.active_samples;

        auto sheet_region = BuildRegion(band_region, settings);
        if (!sheet_region.HasSamples()) {
            continue;
        }

        result.sheet_area += sheet_region.measure_area;
        result.sheet_footprint_area += sheet_region.footprint_area;
        result.collapsed_samples += sheet_region.samples.size();
        result.patch_count += sheet_region.patch_count;
        result.fiber_count += sheet_region.fiber_count;
        result.largest_patch_area = std::max(result.largest_patch_area, sheet_region.largest_patch_area);
        result.fallback_regions += sheet_region.used_fallback ? 1u : 0u;
        result.used_fallback = result.used_fallback || sheet_region.used_fallback;
        sheet_center_sum += sheet_region.centroid_world * sheet_region.measure_area;
        support_seed_sum += sheet_region.mean_support_seed_count * static_cast<double>(sheet_region.samples.size());

        double region_pressure_weight = 0;
        result.support_bbox_world += sheet_region.support_bbox_world;
        for (const auto& sample : sheet_region.samples) {
            region_pressure_weight += sample.force_world.Length();
        }
        pressure_weight_sum += region_pressure_weight;
        pressure_center_sum += sheet_region.pressure_center_world * region_pressure_weight;

        result.regions.push_back(std::move(sheet_region));
    }

    if (result.sheet_area > 0) {
        result.sheet_center_world = sheet_center_sum / result.sheet_area;
    }
    result.pressure_center_world =
        pressure_weight_sum > 0 ? pressure_center_sum / pressure_weight_sum : result.sheet_center_world;
    result.area_ratio = result.band_area > 0 ? result.sheet_footprint_area / result.band_area : 0.0;
    result.mean_support_seed_count =
        result.collapsed_samples > 0 ? support_seed_sum / static_cast<double>(result.collapsed_samples) : 0.0;

    if (result.sheet_area > 0) {
        double weighted_error_sum = 0;
        for (const auto& region : result.regions) {
            weighted_error_sum += region.measure_area * region.normal_spread * region.normal_spread;
        }
        result.normal_spread = std::sqrt(weighted_error_sum / result.sheet_area);
    }

    return result;
}

}  // end namespace chrono
