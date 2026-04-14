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
        patch.support_bbox_world += sample.centroid_world;

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

    for (const auto seed_index : cluster) {
        const auto& seed = seeds[seed_index];
        measure_area += seed.measure_area;
        centroid_sum += seed.seed_world * seed.measure_area;
        normal_sum += seed.seed_normal_world * seed.measure_area;
        force_sum += seed.force_world;
        sample.source_bounds_world += seed.band_point_world;
        sample.support_bbox_world += seed.seed_world;
        sample.source_sample_indices.push_back(seed.source_sample_index);
        sample.support_seed_count++;

        const double pressure_weight = seed.force_world.Length();
        pressure_weight_sum += pressure_weight;
        pressure_center_sum += seed.seed_world * pressure_weight;
    }

    if (measure_area <= 1.0e-16) {
        return sample;
    }

    sample.measure_area = measure_area;
    sample.footprint_area = measure_area;
    sample.centroid_world = centroid_sum / measure_area;
    sample.pressure_center_world =
        pressure_weight_sum > 0 ? pressure_center_sum / pressure_weight_sum : sample.centroid_world;
    sample.normal_world = SafeNormalized(normal_sum, seeds[cluster.front()].seed_normal_world);
    sample.force_world = force_sum;
    return sample;
}

std::vector<ChSDFSheetPatch> BuildPatchGraph(const std::vector<ChSDFSheetFiberSample>& samples,
                                             double patch_connection_radius,
                                             double patch_normal_cosine) {
    std::vector<ChSDFSheetPatch> patches;
    if (samples.empty()) {
        return patches;
    }

    const double cell_size = std::max(patch_connection_radius, 1.0e-6);
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

                        const double cosine = NormalCosine(samples[i].normal_world, samples[j].normal_world);
                        if (patch_normal_cosine > -1 && cosine < patch_normal_cosine) {
                            continue;
                        }
                        if ((samples[j].centroid_world - samples[i].centroid_world).Length() >
                            patch_connection_radius + 1.0e-12) {
                            continue;
                        }
                        uf.Unite(i, j);
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
        patch.support_columns = patch.sample_indices.size();
        patch.support_seed_count = 0;

        ChVector3d centroid_sum = VNULL;
        ChVector3d pressure_center_sum = VNULL;
        ChVector3d normal_sum = VNULL;
        double pressure_weight_sum = 0;

        for (const auto sample_index : patch.sample_indices) {
            const auto& sample = samples[sample_index];
            patch.support_seed_count += sample.support_seed_count;
            patch.measure_area += sample.measure_area;
            patch.footprint_area += sample.footprint_area;
            centroid_sum += sample.centroid_world * sample.measure_area;
            normal_sum += sample.normal_world * sample.measure_area;
            patch.bounds_world += sample.centroid_world;
            patch.support_bbox_world += sample.centroid_world;

            const double pressure_weight = sample.force_world.Length();
            pressure_weight_sum += pressure_weight;
            pressure_center_sum += sample.pressure_center_world * pressure_weight;
        }

        if (patch.measure_area > 0) {
            patch.centroid_world = centroid_sum / patch.measure_area;
        }
        patch.pressure_center_world =
            pressure_weight_sum > 0 ? pressure_center_sum / pressure_weight_sum : patch.centroid_world;
        patch.mean_normal_world = SafeNormalized(normal_sum);
        patches.push_back(std::move(patch));
    }

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
        region.support_bbox_world += sample.centroid_world;
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
        region.support_bbox_world += sample.centroid_world;
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
    region.patches = BuildPatchGraph(region.samples, ResolvePatchConnectionRadius(band_region, settings),
                                     settings.patch_normal_cosine);
    region.patch_count = region.patches.size();
    region.largest_patch_area = 0;
    for (const auto& patch : region.patches) {
        region.largest_patch_area = std::max(region.largest_patch_area, patch.measure_area);
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
        for (const auto& sample : sheet_region.samples) {
            result.support_bbox_world += sample.centroid_world;
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
