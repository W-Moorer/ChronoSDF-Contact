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
#include <queue>
#include <unordered_map>

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

}  // end namespace chrono
