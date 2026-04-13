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

#include "chrono/collision/sdf/ChSDFContactRegion.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>

#include "chrono/core/ChFrameMoving.h"
#include "chrono/core/ChMatrix33.h"

namespace chrono {
namespace {

struct CoordKey {
    int x = 0;
    int y = 0;
    int z = 0;

    bool operator==(const CoordKey& other) const { return x == other.x && y == other.y && z == other.z; }
};

struct CoordKeyHasher {
    std::size_t operator()(const CoordKey& key) const {
        const std::size_t hx = static_cast<std::size_t>(static_cast<uint32_t>(key.x) * 73856093u);
        const std::size_t hy = static_cast<std::size_t>(static_cast<uint32_t>(key.y) * 19349663u);
        const std::size_t hz = static_cast<std::size_t>(static_cast<uint32_t>(key.z) * 83492791u);
        return hx ^ hy ^ hz;
    }
};

CoordKey MakeCoordKey(const ChVector3i& coord) {
    return CoordKey{coord.x(), coord.y(), coord.z()};
}

double ResolvePositive(double value, double fallback) {
    return value > 0 ? value : fallback;
}

struct IndexRange {
    int min = 0;
    int max = -1;
};

double GridStep(std::size_t count, double half_length) {
    return count > 1 ? (2.0 * half_length / static_cast<double>(count - 1)) : (2.0 * half_length);
}

bool Overlap1D(double amin, double amax, double bmin, double bmax) {
    return !(amax < bmin || bmax < amin);
}

ChAABB TransformAABBToLocal(const ChAABB& aabb, const ChFrame<>& frame) {
    if (aabb.IsInverted()) {
        return ChAABB();
    }

    ChAABB transformed;
    for (int ix = 0; ix < 2; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int iz = 0; iz < 2; ++iz) {
                const ChVector3d corner(ix ? aabb.max.x() : aabb.min.x(), iy ? aabb.max.y() : aabb.min.y(),
                                        iz ? aabb.max.z() : aabb.min.z());
                transformed += frame.TransformPointParentToLocal(corner);
            }
        }
    }

    return transformed;
}

bool IsNeighborOffsetAccepted(int dx, int dy, int dz, int neighbor_mode) {
    const int manhattan = std::abs(dx) + std::abs(dy) + std::abs(dz);
    if (manhattan == 0) {
        return false;
    }

    if (neighbor_mode <= 6) {
        return manhattan == 1;
    }

    if (neighbor_mode <= 18) {
        return manhattan <= 2;
    }

    return true;
}

bool SamplesAreConnected(const ChSDFContactRegionSample& a,
                         const ChSDFContactRegionSample& b,
                         const ChSDFContactRegionBuilder::Settings& settings) {
    if (settings.max_distance_jump >= 0 &&
        std::abs(a.probe.distance - b.probe.distance) > settings.max_distance_jump) {
        return false;
    }

    if (settings.min_neighbor_normal_cosine > -1) {
        const double na = a.probe.normal_world.Length();
        const double nb = b.probe.normal_world.Length();
        if (na > 0 && nb > 0) {
            const double cosine = Vdot(a.probe.normal_world / na, b.probe.normal_world / nb);
            if (cosine < settings.min_neighbor_normal_cosine) {
                return false;
            }
        }
    }

    return true;
}

void FinalizeRegion(ChSDFContactRegion& region,
                    const ChFrame<>& seed_patch_frame_shape,
                    const ChSDFContactPatchSampler::Settings& patch_settings,
                    const ChNanoVDBGridInfo& grid_info) {
    if (region.samples.empty()) {
        return;
    }

    const double inv_count = 1.0 / static_cast<double>(region.samples.size());
    region.centroid_shape *= inv_count;
    region.centroid_patch *= inv_count;

    const double normal_length = region.mean_normal_shape.Length();
    if (normal_length > 0) {
        region.mean_normal_shape /= normal_length;
    }

    const ChVector3d origin_patch(region.centroid_patch.x(), region.centroid_patch.y(), 0.0);
    region.patch_frame_shape =
        ChFrame<>(seed_patch_frame_shape.TransformPointLocalToParent(origin_patch), seed_patch_frame_shape.GetRotMat());

    region.suggested_patch_settings = patch_settings;

    const double du = GridStep(patch_settings.samples_u, patch_settings.half_length_u);
    const double dv = GridStep(patch_settings.samples_v, patch_settings.half_length_v);
    const double min_du = du > 0 ? du : std::max(grid_info.voxel_size.x(), 1.0e-6);
    const double min_dv = dv > 0 ? dv : std::max(grid_info.voxel_size.y(), 1.0e-6);

    const double half_u =
        std::max(std::abs(region.patch_bounds.min.x() - region.centroid_patch.x()),
                 std::abs(region.patch_bounds.max.x() - region.centroid_patch.x())) +
        0.5 * min_du;
    const double half_v =
        std::max(std::abs(region.patch_bounds.min.y() - region.centroid_patch.y()),
                 std::abs(region.patch_bounds.max.y() - region.centroid_patch.y())) +
        0.5 * min_dv;

    region.suggested_patch_settings.half_length_u = std::max(0.5 * min_du, half_u);
    region.suggested_patch_settings.half_length_v = std::max(0.5 * min_dv, half_v);
    region.suggested_patch_settings.samples_u =
        std::max<std::size_t>(3, static_cast<std::size_t>(std::ceil(2.0 * region.suggested_patch_settings.half_length_u /
                                                                     std::max(min_du, 1.0e-12))) +
                                      1);
    region.suggested_patch_settings.samples_v =
        std::max<std::size_t>(3, static_cast<std::size_t>(std::ceil(2.0 * region.suggested_patch_settings.half_length_v /
                                                                     std::max(min_dv, 1.0e-12))) +
                                      1);
}

double MinComponent(const ChVector3d& v) {
    return std::min(v.x(), std::min(v.y(), v.z()));
}

double ResolveSampleSpacing(const ChNanoVDBGridInfo& info_a,
                            const ChNanoVDBGridInfo& info_b,
                            const ChSDFContactRegionBuilder::Settings& settings) {
    if (settings.sample_spacing > 0) {
        return settings.sample_spacing;
    }

    double spacing = std::numeric_limits<double>::infinity();
    if (info_a.valid) {
        spacing = std::min(spacing, MinComponent(info_a.voxel_size));
    }
    if (info_b.valid) {
        spacing = std::min(spacing, MinComponent(info_b.voxel_size));
    }

    return std::isfinite(spacing) && spacing > 0 ? spacing : 1.0e-3;
}

double ResolveSampleMaxAbsDistance(const ChNanoVDBGridInfo& info_a,
                                   const ChNanoVDBGridInfo& info_b,
                                   const ChSDFContactRegionBuilder::Settings& settings) {
    if (settings.sample_max_abs_distance > 0) {
        return settings.sample_max_abs_distance;
    }

    double spacing = ResolveSampleSpacing(info_a, info_b, settings);
    return 1.5 * spacing;
}

ChVector3i QuantizePoint(const ChVector3d& point, const ChVector3d& origin, double spacing) {
    return ChVector3i(static_cast<int>(std::llround((point.x() - origin.x()) / spacing)),
                      static_cast<int>(std::llround((point.y() - origin.y()) / spacing)),
                      static_cast<int>(std::llround((point.z() - origin.z()) / spacing)));
}

IndexRange ComputeGridRange(double min_value, double max_value, double origin, double spacing) {
    IndexRange range;
    range.min = static_cast<int>(std::ceil((min_value - origin) / spacing));
    range.max = static_cast<int>(std::floor((max_value - origin) / spacing));

    if (range.min > range.max) {
        const int center = static_cast<int>(std::llround((0.5 * (min_value + max_value) - origin) / spacing));
        range.min = center;
        range.max = center;
    }

    return range;
}

ChVector3d SurfacePoint(const ChVector3d& query_point_shape,
                        const ChSDFProbeResult& probe_shape,
                        const ChVector3d& normal_shape) {
    return query_point_shape - normal_shape * probe_shape.distance;
}

ChVector3d ProjectOntoTangent(const ChVector3d& direction, const ChVector3d& normal) {
    return direction - normal * Vdot(direction, normal);
}

ChVector3d ChoosePatchTangent(const ChSDFBrickPairRegion& region, const ChVector3d& normal) {
    ChVector3d best_tangent = VNULL;
    double best_length = 0;

    for (const auto& sample : region.samples) {
        const ChVector3d tangent = ProjectOntoTangent(sample.point_world - region.centroid_world, normal);
        const double length = tangent.Length();
        if (length > best_length) {
            best_tangent = tangent / length;
            best_length = length;
        }
    }

    if (best_length > 1.0e-12) {
        return best_tangent;
    }

    const ChVector3d fallback_axes[] = {VECT_X, VECT_Y, VECT_Z};
    for (const auto& axis : fallback_axes) {
        const ChVector3d tangent = ProjectOntoTangent(axis, normal);
        const double length = tangent.Length();
        if (length > 1.0e-12) {
            return tangent / length;
        }
    }

    return VECT_X;
}

bool PairSamplesAreConnected(const ChSDFBrickPairRegionSample& a,
                             const ChSDFBrickPairRegionSample& b,
                             const ChSDFContactRegionBuilder::Settings& settings) {
    if (settings.max_distance_jump >= 0 && std::abs(a.combined_gap - b.combined_gap) > settings.max_distance_jump) {
        return false;
    }

    if (settings.min_neighbor_normal_cosine > -1) {
        const double na = a.contact_normal_world.Length();
        const double nb = b.contact_normal_world.Length();
        if (na > 0 && nb > 0) {
            const double cosine = Vdot(a.contact_normal_world / na, b.contact_normal_world / nb);
            if (cosine < settings.min_neighbor_normal_cosine) {
                return false;
            }
        }
    }

    return true;
}

void FinalizeBrickPairRegion(ChSDFBrickPairRegion& region,
                             const ChFrame<>& shape_a_frame_abs,
                             const ChFrame<>& shape_b_frame_abs,
                             double sample_spacing) {
    if (region.samples.empty()) {
        return;
    }

    const double inv_count = 1.0 / static_cast<double>(region.samples.size());
    region.centroid_world *= inv_count;
    region.centroid_shape_a *= inv_count;
    region.centroid_shape_b *= inv_count;

    const double normal_length = region.mean_normal_world.Length();
    if (normal_length > 0) {
        region.mean_normal_world /= normal_length;
    } else {
        region.mean_normal_world = VECT_Z;
    }

    ChMatrix33<> rot;
    rot.SetFromAxisZ(region.mean_normal_world, ChoosePatchTangent(region, region.mean_normal_world));
    region.contact_frame_world = ChFrame<>(region.centroid_world, rot);
    region.contact_frame_shape_a =
        ChFrame<>(shape_a_frame_abs.TransformParentToLocal(ChFrameMoving<>(region.contact_frame_world)).GetCoordsys());
    region.contact_frame_shape_b =
        ChFrame<>(shape_b_frame_abs.TransformParentToLocal(ChFrameMoving<>(region.contact_frame_world)).GetCoordsys());

    region.sample_spacing = sample_spacing;
    region.suggested_patch_settings.max_abs_distance = sample_spacing;

    for (auto& sample : region.samples) {
        sample.point_patch = region.contact_frame_world.TransformPointParentToLocal(sample.point_world);
        region.patch_bounds += sample.point_patch;
    }

    const double half_u =
        std::max(std::abs(region.patch_bounds.min.x()), std::abs(region.patch_bounds.max.x())) + 0.5 * sample_spacing;
    const double half_v =
        std::max(std::abs(region.patch_bounds.min.y()), std::abs(region.patch_bounds.max.y())) + 0.5 * sample_spacing;

    region.suggested_patch_settings.half_length_u = std::max(0.5 * sample_spacing, half_u);
    region.suggested_patch_settings.half_length_v = std::max(0.5 * sample_spacing, half_v);
    region.suggested_patch_settings.samples_u =
        std::max<std::size_t>(3, static_cast<std::size_t>(
                                     std::ceil(2.0 * region.suggested_patch_settings.half_length_u /
                                               std::max(sample_spacing, 1.0e-12))) +
                                      1);
    region.suggested_patch_settings.samples_v =
        std::max<std::size_t>(3, static_cast<std::size_t>(
                                     std::ceil(2.0 * region.suggested_patch_settings.half_length_v /
                                               std::max(sample_spacing, 1.0e-12))) +
                                      1);
    region.suggested_patch_settings.require_zero_crossing = false;
    region.suggested_patch_settings.min_abs_normal_alignment = 0;
    region.suggested_patch_settings.max_abs_distance =
        std::max(region.suggested_patch_settings.max_abs_distance, 1.5 * sample_spacing);
    region.suggested_patch_settings.samples_u = std::max<std::size_t>(region.suggested_patch_settings.samples_u, 2);
    region.suggested_patch_settings.samples_v = std::max<std::size_t>(region.suggested_patch_settings.samples_v, 2);
    if (region.suggested_patch_settings.samples_u % 2 == 0) {
        region.suggested_patch_settings.samples_u += 1;
    }
    if (region.suggested_patch_settings.samples_v % 2 == 0) {
        region.suggested_patch_settings.samples_v += 1;
    }
}

}  // namespace

std::vector<ChSDFPatchBrickCandidate> ChSDFContactRegionBuilder::FindPatchBrickCandidates(
    const ChCollisionShapeSDF& shape,
    const ChFrame<>& patch_frame_shape,
    const ChSDFContactPatchSampler::Settings& patch_settings,
    const Settings& settings) {
    std::vector<ChSDFPatchBrickCandidate> candidates;

    if (!shape.IsLoaded()) {
        return candidates;
    }

    const auto& bricks = shape.GetLeafBricks();
    if (bricks.empty()) {
        return candidates;
    }

    const double half_u = patch_settings.half_length_u + settings.brick_margin;
    const double half_v = patch_settings.half_length_v + settings.brick_margin;
    const double half_w = ResolvePositive(settings.slab_half_thickness, patch_settings.max_abs_distance) + settings.brick_margin;

    for (std::size_t i = 0; i < bricks.size(); ++i) {
        const auto& brick = bricks[i];
        if (!brick.valid || brick.index_bounds.IsInverted() || brick.world_bounds.IsInverted()) {
            continue;
        }

        const ChAABB brick_patch_bounds = TransformAABBToLocal(brick.world_bounds, patch_frame_shape);
        if (brick_patch_bounds.IsInverted()) {
            continue;
        }

        if (!Overlap1D(brick_patch_bounds.min.x(), brick_patch_bounds.max.x(), -half_u, half_u)) {
            continue;
        }
        if (!Overlap1D(brick_patch_bounds.min.y(), brick_patch_bounds.max.y(), -half_v, half_v)) {
            continue;
        }
        if (!Overlap1D(brick_patch_bounds.min.z(), brick_patch_bounds.max.z(), -half_w, half_w)) {
            continue;
        }

        ChSDFPatchBrickCandidate candidate;
        candidate.brick_index = i;
        candidate.slab_distance = std::abs(brick_patch_bounds.Center().z());
        candidate.center_patch = brick_patch_bounds.Center();
        candidate.size_patch = brick_patch_bounds.Size();
        candidate.brick = brick;
        candidates.push_back(candidate);
    }

    std::stable_sort(candidates.begin(), candidates.end(), [](const ChSDFPatchBrickCandidate& a,
                                                              const ChSDFPatchBrickCandidate& b) {
        if (a.slab_distance != b.slab_distance) {
            return a.slab_distance < b.slab_distance;
        }
        return a.brick.min_abs_value < b.brick.min_abs_value;
    });

    return candidates;
}

std::vector<ChSDFContactRegion> ChSDFContactRegionBuilder::BuildPatchRegions(
    const ChCollisionShapeSDF& shape,
    const ChFrame<>& patch_frame_shape,
    const ChSDFContactPatchSampler::Settings& patch_settings,
    const Settings& settings) {
    std::vector<ChSDFContactRegion> regions;

    const auto level_set = shape.GetLevelSet();
    if (!shape.IsLoaded() || !level_set) {
        return regions;
    }

    const auto candidates = FindPatchBrickCandidates(shape, patch_frame_shape, patch_settings, settings);
    if (candidates.empty()) {
        return regions;
    }

    const ChNanoVDBGridInfo grid_info = shape.GetGridInfo();
    const ChVector3d plane_normal_shape =
        patch_frame_shape.TransformDirectionLocalToParent(VECT_Z).GetNormalized();
    const double slab_half_thickness = ResolvePositive(settings.slab_half_thickness, patch_settings.max_abs_distance);
    const double sample_max_abs_distance =
        ResolvePositive(settings.sample_max_abs_distance, patch_settings.max_abs_distance);

    std::vector<ChSDFContactRegionSample> accepted_samples;
    std::unordered_map<CoordKey, std::size_t, CoordKeyHasher> sample_lookup;
    accepted_samples.reserve(candidates.size() * 64);

    for (const auto& candidate : candidates) {
        const ChIntAABB& bbox = candidate.brick.index_bounds;
        for (int ix = bbox.min.x(); ix <= bbox.max.x(); ++ix) {
            for (int iy = bbox.min.y(); iy <= bbox.max.y(); ++iy) {
                for (int iz = bbox.min.z(); iz <= bbox.max.z(); ++iz) {
                    ChVector3i coord(ix, iy, iz);
                    CoordKey key = MakeCoordKey(coord);
                    if (sample_lookup.find(key) != sample_lookup.end()) {
                        continue;
                    }

                    const ChSDFProbeResult probe =
                        level_set->ProbeIndex(ChVector3d(ix + 0.5, iy + 0.5, iz + 0.5));
                    if (!probe.valid) {
                        continue;
                    }

                    const ChVector3d point_patch = patch_frame_shape.TransformPointParentToLocal(probe.point_world);
                    if (std::abs(point_patch.x()) > patch_settings.half_length_u + settings.brick_margin) {
                        continue;
                    }
                    if (std::abs(point_patch.y()) > patch_settings.half_length_v + settings.brick_margin) {
                        continue;
                    }
                    if (std::abs(point_patch.z()) > slab_half_thickness) {
                        continue;
                    }
                    if (std::abs(probe.distance) > sample_max_abs_distance) {
                        continue;
                    }

                    const double normal_alignment = std::abs(Vdot(probe.normal_world, plane_normal_shape));
                    if (normal_alignment < settings.min_abs_normal_alignment) {
                        continue;
                    }
                    if (settings.require_zero_crossing && !probe.zero_crossing) {
                        continue;
                    }

                    ChSDFContactRegionSample sample;
                    sample.brick_index = candidate.brick_index;
                    sample.coord = coord;
                    sample.point_shape = probe.point_world;
                    sample.point_patch = point_patch;
                    sample.normal_alignment = normal_alignment;
                    sample.probe = probe;

                    sample_lookup.emplace(key, accepted_samples.size());
                    accepted_samples.push_back(sample);
                }
            }
        }
    }

    if (accepted_samples.empty()) {
        return regions;
    }

    std::vector<char> visited(accepted_samples.size(), 0);

    for (std::size_t i = 0; i < accepted_samples.size(); ++i) {
        if (visited[i]) {
            continue;
        }

        std::queue<std::size_t> queue;
        queue.push(i);
        visited[i] = 1;

        ChSDFContactRegion region;
        region.region_id = regions.size();

        std::unordered_map<std::size_t, char> brick_flags;

        while (!queue.empty()) {
            const std::size_t current = queue.front();
            queue.pop();

            const auto& sample = accepted_samples[current];
            region.samples.push_back(sample);
            region.shape_bounds += sample.point_shape;
            region.patch_bounds += sample.point_patch;
            region.centroid_shape += sample.point_shape;
            region.centroid_patch += sample.point_patch;
            region.mean_normal_shape += sample.probe.normal_world;
            brick_flags[sample.brick_index] = 1;

            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        if (!IsNeighborOffsetAccepted(dx, dy, dz, settings.neighbor_mode)) {
                            continue;
                        }

                        const CoordKey neighbor_key{sample.coord.x() + dx, sample.coord.y() + dy, sample.coord.z() + dz};
                        const auto it = sample_lookup.find(neighbor_key);
                        if (it == sample_lookup.end()) {
                            continue;
                        }

                        const std::size_t neighbor_index = it->second;
                        if (visited[neighbor_index]) {
                            continue;
                        }

                        if (!SamplesAreConnected(sample, accepted_samples[neighbor_index], settings)) {
                            continue;
                        }

                        visited[neighbor_index] = 1;
                        queue.push(neighbor_index);
                    }
                }
            }
        }

        if (region.samples.size() < settings.min_region_samples) {
            continue;
        }

        region.brick_indices.reserve(brick_flags.size());
        for (const auto& item : brick_flags) {
            region.brick_indices.push_back(item.first);
        }
        std::sort(region.brick_indices.begin(), region.brick_indices.end());

        FinalizeRegion(region, patch_frame_shape, patch_settings, grid_info);
        regions.push_back(region);
    }

    std::stable_sort(regions.begin(), regions.end(),
                     [](const ChSDFContactRegion& a, const ChSDFContactRegion& b) {
                         return a.samples.size() > b.samples.size();
                     });

    for (std::size_t i = 0; i < regions.size(); ++i) {
        regions[i].region_id = i;
    }

    return regions;
}

std::vector<ChSDFBrickPairRegionSample> ChSDFContactRegionBuilder::BuildBrickPairSamples(
    const ChCollisionShapeSDF& shape_a,
    const ChFrame<>& shape_a_frame_abs,
    const ChCollisionShapeSDF& shape_b,
    const ChFrame<>& shape_b_frame_abs,
    const std::vector<ChSDFBrickPairCandidate>& brick_pairs,
    const Settings& settings) {
    std::vector<ChSDFBrickPairRegionSample> samples;

    if (!shape_a.IsLoaded() || !shape_b.IsLoaded() || brick_pairs.empty()) {
        return samples;
    }

    ChAABB union_bounds;
    bool has_overlap = false;
    for (const auto& brick_pair : brick_pairs) {
        if (brick_pair.overlap_world.IsInverted()) {
            continue;
        }
        union_bounds += brick_pair.overlap_world;
        has_overlap = true;
    }

    if (!has_overlap || union_bounds.IsInverted()) {
        return samples;
    }

    const ChNanoVDBGridInfo info_a = shape_a.GetGridInfo();
    const ChNanoVDBGridInfo info_b = shape_b.GetGridInfo();
    const double spacing = ResolveSampleSpacing(info_a, info_b, settings);
    const double max_abs_distance = ResolveSampleMaxAbsDistance(info_a, info_b, settings);

    const ChVector3d origin = union_bounds.min;
    std::unordered_map<CoordKey, std::size_t, CoordKeyHasher> sample_lookup;
    sample_lookup.reserve(brick_pairs.size() * 64);

    for (const auto& brick_pair : brick_pairs) {
        if (brick_pair.overlap_world.IsInverted()) {
            continue;
        }

        const IndexRange ix_range =
            ComputeGridRange(brick_pair.overlap_world.min.x(), brick_pair.overlap_world.max.x(), origin.x(), spacing);
        const IndexRange iy_range =
            ComputeGridRange(brick_pair.overlap_world.min.y(), brick_pair.overlap_world.max.y(), origin.y(), spacing);
        const IndexRange iz_range =
            ComputeGridRange(brick_pair.overlap_world.min.z(), brick_pair.overlap_world.max.z(), origin.z(), spacing);

        for (int ix = ix_range.min; ix <= ix_range.max; ++ix) {
            for (int iy = iy_range.min; iy <= iy_range.max; ++iy) {
                for (int iz = iz_range.min; iz <= iz_range.max; ++iz) {
                    const ChVector3i coord(ix, iy, iz);
                    const ChVector3d query_world(origin.x() + spacing * ix, origin.y() + spacing * iy,
                                                 origin.z() + spacing * iz);

                    const ChVector3d query_shape_a = shape_a_frame_abs.TransformPointParentToLocal(query_world);
                    const ChVector3d query_shape_b = shape_b_frame_abs.TransformPointParentToLocal(query_world);
                    const ChSDFProbeResult probe_a = shape_a.ProbeLocal(query_shape_a);
                    const ChSDFProbeResult probe_b = shape_b.ProbeLocal(query_shape_b);

                    if (!probe_a.valid || !probe_b.valid) {
                        continue;
                    }

                    if (std::abs(probe_a.distance) > max_abs_distance || std::abs(probe_b.distance) > max_abs_distance) {
                        continue;
                    }

                    if (settings.require_zero_crossing && !(probe_a.zero_crossing || probe_b.zero_crossing)) {
                        continue;
                    }

                    const double normal_a_length = probe_a.normal_world.Length();
                    const double normal_b_length = probe_b.normal_world.Length();
                    if (normal_a_length <= 1.0e-12 || normal_b_length <= 1.0e-12) {
                        continue;
                    }

                    const ChVector3d normal_shape_a = probe_a.normal_world / normal_a_length;
                    const ChVector3d normal_shape_b = probe_b.normal_world / normal_b_length;
                    const ChVector3d normal_world_a =
                        shape_a_frame_abs.TransformDirectionLocalToParent(normal_shape_a).GetNormalized();
                    const ChVector3d normal_world_b =
                        shape_b_frame_abs.TransformDirectionLocalToParent(normal_shape_b).GetNormalized();
                    const double normal_opposition = -Vdot(normal_world_a, normal_world_b);
                    if (normal_opposition < settings.min_opposed_normal_cosine) {
                        continue;
                    }

                    const ChVector3d surface_shape_a = SurfacePoint(query_shape_a, probe_a, normal_shape_a);
                    const ChVector3d surface_shape_b = SurfacePoint(query_shape_b, probe_b, normal_shape_b);
                    const ChVector3d surface_world_a = shape_a_frame_abs.TransformPointLocalToParent(surface_shape_a);
                    const ChVector3d surface_world_b = shape_b_frame_abs.TransformPointLocalToParent(surface_shape_b);

                    ChVector3d contact_normal_world = normal_world_a - normal_world_b;
                    const double contact_normal_length = contact_normal_world.Length();
                    if (contact_normal_length <= 1.0e-12) {
                        continue;
                    }
                    contact_normal_world /= contact_normal_length;

                    const double combined_gap = Vdot(surface_world_b - surface_world_a, contact_normal_world);
                    if (settings.max_combined_gap >= 0 && std::abs(combined_gap) > settings.max_combined_gap) {
                        continue;
                    }

                    ChSDFBrickPairRegionSample sample;
                    sample.brick_a_index = brick_pair.brick_a_index;
                    sample.brick_b_index = brick_pair.brick_b_index;
                    sample.coord = coord;
                    sample.point_world = 0.5 * (surface_world_a + surface_world_b);
                    sample.point_shape_a = query_shape_a;
                    sample.point_shape_b = query_shape_b;
                    sample.surface_world_a = surface_world_a;
                    sample.surface_world_b = surface_world_b;
                    sample.surface_shape_a = surface_shape_a;
                    sample.surface_shape_b = surface_shape_b;
                    sample.normal_world_a = normal_world_a;
                    sample.normal_world_b = normal_world_b;
                    sample.contact_normal_world = contact_normal_world;
                    sample.distance_a = probe_a.distance;
                    sample.distance_b = probe_b.distance;
                    sample.combined_gap = combined_gap;
                    sample.normal_opposition = normal_opposition;
                    sample.probe_a = probe_a;
                    sample.probe_b = probe_b;

                    const CoordKey key = MakeCoordKey(coord);
                    const auto it = sample_lookup.find(key);
                    if (it == sample_lookup.end()) {
                        sample_lookup.emplace(key, samples.size());
                        samples.push_back(sample);
                    } else {
                        const std::size_t old_index = it->second;
                        const double old_score =
                            std::max(std::abs(samples[old_index].distance_a), std::abs(samples[old_index].distance_b));
                        const double new_score = std::max(std::abs(sample.distance_a), std::abs(sample.distance_b));
                        if (new_score < old_score) {
                            samples[old_index] = sample;
                        }
                    }
                }
            }
        }
    }

    return samples;
}

std::vector<ChSDFBrickPairRegion> ChSDFContactRegionBuilder::BuildBrickPairRegions(
    const ChCollisionShapeSDF& shape_a,
    const ChFrame<>& shape_a_frame_abs,
    const ChCollisionShapeSDF& shape_b,
    const ChFrame<>& shape_b_frame_abs,
    const std::vector<ChSDFBrickPairCandidate>& brick_pairs,
    const Settings& settings) {
    std::vector<ChSDFBrickPairRegion> regions;

    const auto samples = BuildBrickPairSamples(shape_a, shape_a_frame_abs, shape_b, shape_b_frame_abs, brick_pairs, settings);
    if (samples.empty()) {
        return regions;
    }

    const ChNanoVDBGridInfo info_a = shape_a.GetGridInfo();
    const ChNanoVDBGridInfo info_b = shape_b.GetGridInfo();
    const double sample_spacing = ResolveSampleSpacing(info_a, info_b, settings);

    std::unordered_map<CoordKey, std::size_t, CoordKeyHasher> sample_lookup;
    sample_lookup.reserve(samples.size());
    for (std::size_t i = 0; i < samples.size(); ++i) {
        sample_lookup.emplace(MakeCoordKey(samples[i].coord), i);
    }

    std::vector<char> visited(samples.size(), 0);

    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (visited[i]) {
            continue;
        }

        std::queue<std::size_t> queue;
        queue.push(i);
        visited[i] = 1;

        ChSDFBrickPairRegion region;
        region.region_id = regions.size();

        std::unordered_map<std::size_t, char> brick_flags_a;
        std::unordered_map<std::size_t, char> brick_flags_b;

        while (!queue.empty()) {
            const std::size_t current = queue.front();
            queue.pop();

            const auto& sample = samples[current];
            region.samples.push_back(sample);
            region.world_bounds += sample.point_world;
            region.shape_a_bounds += sample.surface_shape_a;
            region.shape_b_bounds += sample.surface_shape_b;
            region.centroid_world += sample.point_world;
            region.centroid_shape_a += sample.surface_shape_a;
            region.centroid_shape_b += sample.surface_shape_b;
            region.mean_normal_world += sample.contact_normal_world;
            brick_flags_a[sample.brick_a_index] = 1;
            brick_flags_b[sample.brick_b_index] = 1;

            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        if (!IsNeighborOffsetAccepted(dx, dy, dz, settings.neighbor_mode)) {
                            continue;
                        }

                        const CoordKey neighbor_key{sample.coord.x() + dx, sample.coord.y() + dy, sample.coord.z() + dz};
                        const auto it = sample_lookup.find(neighbor_key);
                        if (it == sample_lookup.end()) {
                            continue;
                        }

                        const std::size_t neighbor_index = it->second;
                        if (visited[neighbor_index]) {
                            continue;
                        }

                        if (!PairSamplesAreConnected(sample, samples[neighbor_index], settings)) {
                            continue;
                        }

                        visited[neighbor_index] = 1;
                        queue.push(neighbor_index);
                    }
                }
            }
        }

        if (region.samples.size() < settings.min_region_samples) {
            continue;
        }

        region.brick_a_indices.reserve(brick_flags_a.size());
        for (const auto& item : brick_flags_a) {
            region.brick_a_indices.push_back(item.first);
        }
        std::sort(region.brick_a_indices.begin(), region.brick_a_indices.end());

        region.brick_b_indices.reserve(brick_flags_b.size());
        for (const auto& item : brick_flags_b) {
            region.brick_b_indices.push_back(item.first);
        }
        std::sort(region.brick_b_indices.begin(), region.brick_b_indices.end());

        FinalizeBrickPairRegion(region, shape_a_frame_abs, shape_b_frame_abs, sample_spacing);
        regions.push_back(region);
    }

    std::stable_sort(regions.begin(), regions.end(),
                     [](const ChSDFBrickPairRegion& a, const ChSDFBrickPairRegion& b) {
                         return a.samples.size() > b.samples.size();
                     });

    for (std::size_t i = 0; i < regions.size(); ++i) {
        regions[i].region_id = i;
    }

    return regions;
}

void ChSDFContactRegionBuilder::ReparameterizeBrickPairRegion(ChSDFBrickPairRegion& region,
                                                              const ChFrame<>& shape_a_frame_abs,
                                                              const ChFrame<>& shape_b_frame_abs) {
    if (region.samples.empty()) {
        return;
    }

    region.world_bounds = ChAABB();
    region.patch_bounds = ChAABB();
    region.shape_a_bounds = ChAABB();
    region.shape_b_bounds = ChAABB();
    region.centroid_world = VNULL;
    region.centroid_shape_a = VNULL;
    region.centroid_shape_b = VNULL;
    region.mean_normal_world = VNULL;

    for (const auto& sample : region.samples) {
        region.world_bounds += sample.point_world;
        region.shape_a_bounds += sample.surface_shape_a;
        region.shape_b_bounds += sample.surface_shape_b;
        region.centroid_world += sample.point_world;
        region.centroid_shape_a += sample.surface_shape_a;
        region.centroid_shape_b += sample.surface_shape_b;
        region.mean_normal_world += sample.contact_normal_world;
    }

    const double sample_spacing = region.sample_spacing > 0 ? region.sample_spacing : 1.0e-3;
    FinalizeBrickPairRegion(region, shape_a_frame_abs, shape_b_frame_abs, sample_spacing);
}

}  // end namespace chrono
