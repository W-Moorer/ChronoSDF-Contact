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

struct CarrierKey {
    int axis = -1;
    int x = 0;
    int y = 0;
    int z = 0;

    bool operator==(const CarrierKey& other) const {
        return axis == other.axis && x == other.x && y == other.y && z == other.z;
    }
};

struct CarrierKeyHasher {
    std::size_t operator()(const CarrierKey& key) const {
        std::size_t seed = static_cast<std::size_t>(static_cast<uint32_t>(key.axis + 17) * 12582917u);
        seed ^= static_cast<std::size_t>(static_cast<uint32_t>(key.x) * 73856093u);
        seed ^= static_cast<std::size_t>(static_cast<uint32_t>(key.y) * 19349663u);
        seed ^= static_cast<std::size_t>(static_cast<uint32_t>(key.z) * 83492791u);
        return seed;
    }
};

CarrierKey MakeCarrierKey(int axis, const ChVector3d& surfel_index) {
    return CarrierKey{axis, static_cast<int>(std::llround(2.0 * surfel_index.x())),
                      static_cast<int>(std::llround(2.0 * surfel_index.y())),
                      static_cast<int>(std::llround(2.0 * surfel_index.z()))};
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

ChVector3d SafeNormalized(const ChVector3d& v, const ChVector3d& fallback = VECT_Z) {
    const double length = v.Length();
    return length > 1.0e-12 ? (v / length) : fallback;
}

ChVector3d AxisUnit(int axis) {
    switch (axis) {
        case 0:
            return VECT_X;
        case 1:
            return VECT_Y;
        default:
            return VECT_Z;
    }
}

double FaceAreaForAxis(int axis, const ChVector3d& voxel_size) {
    switch (axis) {
        case 0:
            return voxel_size.y() * voxel_size.z();
        case 1:
            return voxel_size.x() * voxel_size.z();
        default:
            return voxel_size.x() * voxel_size.y();
    }
}

int DominantAxis(const ChVector3d& v) {
    const double ax = std::abs(v.x());
    const double ay = std::abs(v.y());
    const double az = std::abs(v.z());

    if (ax >= ay && ax >= az) {
        return 0;
    }
    if (ay >= az) {
        return 1;
    }
    return 2;
}

bool HasSurfaceCrossing(double d0, double d1, double threshold) {
    const double eps = std::max(threshold * 1.0e-3, 1.0e-12);

    if (std::abs(d0) <= eps && std::abs(d1) <= eps) {
        return false;
    }
    if (std::abs(d0) <= eps || std::abs(d1) <= eps) {
        return true;
    }

    return d0 * d1 < 0.0;
}

struct NormalProjectionResult {
    bool valid = false;
    double signed_gap = 0;
    ChVector3d projected_world = VNULL;
    ChVector3d projected_shape_b = VNULL;
    ChSDFProbeResult projected_probe_b;
};

bool FindZeroAlongIndexSegment(const ChNanoVDBLevelSet& level_set,
                               const ChVector3d& index_begin,
                               const ChVector3d& index_end,
                               double f_begin,
                               double f_end,
                               ChSDFProbeResult& result) {
    const double distance_tolerance = 1.0e-8;
    const double index_tolerance = 1.0e-8;

    if (std::abs(f_begin) <= distance_tolerance) {
        result = level_set.ProbeIndex(index_begin);
        return result.valid;
    }
    if (std::abs(f_end) <= distance_tolerance) {
        result = level_set.ProbeIndex(index_end);
        return result.valid;
    }
    if (f_begin * f_end > 0.0) {
        return false;
    }

    ChVector3d low = index_begin;
    ChVector3d high = index_end;
    double flow = f_begin;

    for (int iter = 0; iter < 32; ++iter) {
        const ChVector3d mid = 0.5 * (low + high);
        const ChSDFProbeResult probe_mid = level_set.ProbeIndex(mid);
        if (!probe_mid.valid) {
            return false;
        }

        if (std::abs(probe_mid.distance) <= distance_tolerance || (high - low).Length() <= index_tolerance) {
            result = probe_mid;
            return true;
        }

        if (flow * probe_mid.distance <= 0.0) {
            high = mid;
        } else {
            low = mid;
            flow = probe_mid.distance;
        }
    }

    result = level_set.ProbeIndex(0.5 * (low + high));
    return result.valid;
}

bool FindZeroAlongRay(const ChCollisionShapeSDF& shape_b,
                      const ChFrame<>& shape_b_frame_abs,
                      const ChVector3d& origin_world,
                      const ChVector3d& ray_direction_world,
                      double s_begin,
                      double s_end,
                      double f_begin,
                      double f_end,
                      NormalProjectionResult& result) {
    if (f_begin * f_end > 0.0) {
        return false;
    }

    double low = s_begin;
    double high = s_end;
    double flow = f_begin;
    double fhigh = f_end;

    for (int iter = 0; iter < 32; ++iter) {
        const double mid = 0.5 * (low + high);
        const ChVector3d mid_world = origin_world + ray_direction_world * mid;
        const ChSDFProbeResult probe_mid =
            shape_b.ProbeLocal(shape_b_frame_abs.TransformPointParentToLocal(mid_world));
        if (!probe_mid.valid) {
            return false;
        }

        if (std::abs(probe_mid.distance) <= 1.0e-8 || std::abs(high - low) <= 1.0e-8) {
            result.valid = true;
            result.signed_gap = mid;
            result.projected_world = mid_world;
            result.projected_shape_b = shape_b_frame_abs.TransformPointParentToLocal(mid_world);
            result.projected_probe_b = probe_mid;
            return true;
        }

        if (flow * probe_mid.distance <= 0.0) {
            high = mid;
            fhigh = probe_mid.distance;
        } else {
            low = mid;
            flow = probe_mid.distance;
        }
    }

    const double root = 0.5 * (low + high);
    const ChVector3d root_world = origin_world + ray_direction_world * root;
    const ChSDFProbeResult probe_root =
        shape_b.ProbeLocal(shape_b_frame_abs.TransformPointParentToLocal(root_world));
    if (!probe_root.valid) {
        return false;
    }

    result.valid = true;
    result.signed_gap = root;
    result.projected_world = root_world;
    result.projected_shape_b = shape_b_frame_abs.TransformPointParentToLocal(root_world);
    result.projected_probe_b = probe_root;
    return true;
}

NormalProjectionResult ProjectPointToShapeAlongNormal(const ChCollisionShapeSDF& shape_b,
                                                      const ChFrame<>& shape_b_frame_abs,
                                                      const ChVector3d& carrier_world,
                                                      const ChVector3d& carrier_normal_world,
                                                      double max_search,
                                                      double step_length) {
    NormalProjectionResult best;

    const ChVector3d direction = SafeNormalized(carrier_normal_world, VECT_Y);
    const ChSDFProbeResult probe0 = shape_b.ProbeLocal(shape_b_frame_abs.TransformPointParentToLocal(carrier_world));
    if (!probe0.valid) {
        return best;
    }

    if (std::abs(probe0.distance) <= 1.0e-8) {
        best.valid = true;
        best.signed_gap = 0;
        best.projected_world = carrier_world;
        best.projected_shape_b = shape_b_frame_abs.TransformPointParentToLocal(carrier_world);
        best.projected_probe_b = probe0;
        return best;
    }

    const double safe_step = step_length > 1.0e-6 ? step_length : 1.0e-3;
    const int max_steps = std::max(2, static_cast<int>(std::ceil(max_search / safe_step)));

    for (const double sign : {1.0, -1.0}) {
        double prev_s = 0.0;
        double prev_f = probe0.distance;

        for (int step = 1; step <= max_steps; ++step) {
            const double current_s = sign * safe_step * static_cast<double>(step);
            const ChVector3d current_world = carrier_world + direction * current_s;
            const ChSDFProbeResult current_probe =
                shape_b.ProbeLocal(shape_b_frame_abs.TransformPointParentToLocal(current_world));
            if (!current_probe.valid) {
                break;
            }

            if (prev_f * current_probe.distance <= 0.0) {
                NormalProjectionResult candidate;
                if (FindZeroAlongRay(shape_b, shape_b_frame_abs, carrier_world, direction, prev_s, current_s, prev_f,
                                     current_probe.distance, candidate)) {
                    if (!best.valid || std::abs(candidate.signed_gap) < std::abs(best.signed_gap)) {
                        best = candidate;
                    }
                }
                break;
            }

            prev_s = current_s;
            prev_f = current_probe.distance;
        }
    }

    return best;
}

double NeighborDistanceThreshold(double sample_spacing, int neighbor_mode) {
    const double safe_spacing = sample_spacing > 1.0e-12 ? sample_spacing : 1.0;
    if (neighbor_mode <= 6) {
        return 1.05 * safe_spacing;
    }
    if (neighbor_mode <= 18) {
        return 1.05 * std::sqrt(2.0) * safe_spacing;
    }

    return 1.05 * std::sqrt(3.0) * safe_spacing;
}

ChVector3i QuantizePoint(const ChVector3d& point, const ChVector3d& origin, double spacing) {
    return ChVector3i(static_cast<int>(std::llround((point.x() - origin.x()) / spacing)),
                      static_cast<int>(std::llround((point.y() - origin.y()) / spacing)),
                      static_cast<int>(std::llround((point.z() - origin.z()) / spacing)));
}

bool PointInsideAABB(const ChVector3d& point, const ChAABB& aabb, double margin) {
    if (aabb.IsInverted()) {
        return false;
    }

    return point.x() >= aabb.min.x() - margin && point.x() <= aabb.max.x() + margin &&
           point.y() >= aabb.min.y() - margin && point.y() <= aabb.max.y() + margin &&
           point.z() >= aabb.min.z() - margin && point.z() <= aabb.max.z() + margin;
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

ChVector3d ProjectOntoTangent(const ChVector3d& direction, const ChVector3d& normal) {
    return direction - normal * Vdot(direction, normal);
}

ChVector3d ChoosePatchTangent(const ChSDFBrickPairRegion& region, const ChVector3d& normal) {
    const ChVector3d preferred_axes[] = {VECT_X, VECT_Y, VECT_Z};
    ChVector3d preferred_tangent = VNULL;
    double preferred_length = 0;

    for (const auto& axis : preferred_axes) {
        const ChVector3d tangent = ProjectOntoTangent(axis, normal);
        const double length = tangent.Length();
        if (length > preferred_length) {
            preferred_tangent = tangent / length;
            preferred_length = length;
        }
    }

    if (preferred_length > 1.0e-12) {
        return preferred_tangent;
    }

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

    return VECT_X;
}

bool PairSamplesAreConnected(const ChSDFBrickPairRegionSample& a,
                             const ChSDFBrickPairRegionSample& b,
                             const ChSDFContactRegionBuilder::Settings& settings) {
    if (settings.max_distance_jump >= 0 && std::abs(a.normal_gap - b.normal_gap) > settings.max_distance_jump) {
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

void FilterBrickPairRegionByDominantNormal(ChSDFBrickPairRegion& region, double min_region_normal_cosine) {
    if (region.samples.empty() || min_region_normal_cosine <= -1.0) {
        return;
    }

    ChVector3d dominant_normal = VNULL;
    for (const auto& sample : region.samples) {
        dominant_normal += sample.contact_normal_world;
    }
    dominant_normal = SafeNormalized(dominant_normal, VECT_Z);

    std::vector<ChSDFBrickPairRegionSample> filtered_samples;
    filtered_samples.reserve(region.samples.size());
    for (const auto& sample : region.samples) {
        const ChVector3d sample_normal = SafeNormalized(sample.contact_normal_world, dominant_normal);
        if (Vdot(sample_normal, dominant_normal) >= min_region_normal_cosine) {
            filtered_samples.push_back(sample);
        }
    }

    if (!filtered_samples.empty() && filtered_samples.size() < region.samples.size()) {
        region.samples = std::move(filtered_samples);
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

    const auto level_set_a = shape_a.GetLevelSet();
    if (!level_set_a) {
        return samples;
    }

    ChAABB union_bounds;
    bool has_overlap = false;
    std::unordered_map<std::size_t, std::vector<const ChSDFBrickPairCandidate*>> brick_pairs_by_a;
    brick_pairs_by_a.reserve(brick_pairs.size());
    for (const auto& brick_pair : brick_pairs) {
        if (brick_pair.overlap_world.IsInverted()) {
            continue;
        }
        union_bounds += brick_pair.overlap_world;
        has_overlap = true;
        brick_pairs_by_a[brick_pair.brick_a_index].push_back(&brick_pair);
    }

    if (!has_overlap || union_bounds.IsInverted() || brick_pairs_by_a.empty()) {
        return samples;
    }

    const auto& bricks_a = shape_a.GetLeafBricks();
    const ChNanoVDBGridInfo info_a = shape_a.GetGridInfo();
    const ChNanoVDBGridInfo info_b = shape_b.GetGridInfo();
    const double spacing = ResolveSampleSpacing(info_a, info_b, settings);
    const double max_abs_distance = ResolveSampleMaxAbsDistance(info_a, info_b, settings);
    const ChVector3d voxel_size_a =
        info_a.valid && info_a.voxel_size.x() > 0 && info_a.voxel_size.y() > 0 && info_a.voxel_size.z() > 0
            ? info_a.voxel_size
            : ChVector3d(spacing, spacing, spacing);
    const double near_surface_threshold = 0.5 * std::max(MinComponent(voxel_size_a), 1.0e-6);
    const double projection_search = std::max(2.0 * max_abs_distance, 2.0 * spacing);
    const double projection_step = 0.5 * std::max(MinComponent(voxel_size_a), 1.0e-6);
    std::unordered_map<CarrierKey, std::size_t, CarrierKeyHasher> sample_lookup;
    sample_lookup.reserve(brick_pairs_by_a.size() * 64);

    for (const auto& item : brick_pairs_by_a) {
        const std::size_t brick_a_index = item.first;
        if (brick_a_index >= bricks_a.size()) {
            continue;
        }

        const auto& brick_a = bricks_a[brick_a_index];
        if (!brick_a.valid || brick_a.index_bounds.IsInverted()) {
            continue;
        }

        const ChIntAABB& bbox = brick_a.index_bounds;
        for (int ix = bbox.min.x(); ix <= bbox.max.x(); ++ix) {
            for (int iy = bbox.min.y(); iy <= bbox.max.y(); ++iy) {
                for (int iz = bbox.min.z(); iz <= bbox.max.z(); ++iz) {
                    const ChVector3i coord(ix, iy, iz);
                    const ChVector3d cell_index(ix + 0.5, iy + 0.5, iz + 0.5);
                    const ChSDFProbeResult probe_cell_a = level_set_a->ProbeIndex(cell_index);
                    if (!probe_cell_a.valid) {
                        continue;
                    }

                    for (int axis = 0; axis < 3; ++axis) {
                        ChVector3i neighbor_coord = coord;
                        if (axis == 0) {
                            neighbor_coord.x() += 1;
                        } else if (axis == 1) {
                            neighbor_coord.y() += 1;
                        } else {
                            neighbor_coord.z() += 1;
                        }

                        const ChVector3d neighbor_index(neighbor_coord.x() + 0.5, neighbor_coord.y() + 0.5,
                                                        neighbor_coord.z() + 0.5);
                        const ChSDFProbeResult probe_neighbor_a = level_set_a->ProbeIndex(neighbor_index);
                        if (!probe_neighbor_a.valid) {
                            continue;
                        }

                        if (!HasSurfaceCrossing(probe_cell_a.distance, probe_neighbor_a.distance, near_surface_threshold)) {
                            continue;
                        }

                        ChSDFProbeResult surface_probe_a;
                        if (!FindZeroAlongIndexSegment(*level_set_a, cell_index, neighbor_index, probe_cell_a.distance,
                                                       probe_neighbor_a.distance, surface_probe_a)) {
                            continue;
                        }

                        const ChVector3d normal_shape_a = SafeNormalized(surface_probe_a.normal_world, AxisUnit(axis));
                        if (DominantAxis(normal_shape_a) != axis) {
                            continue;
                        }
                        const ChVector3d surface_shape_a = surface_probe_a.point_world;
                        const ChVector3d surface_world_a =
                            shape_a_frame_abs.TransformPointLocalToParent(surface_shape_a);

                        const ChVector3d query_shape_b = shape_b_frame_abs.TransformPointParentToLocal(surface_world_a);
                        const ChSDFProbeResult probe_b_at_carrier = shape_b.ProbeLocal(query_shape_b);
                        if (!probe_b_at_carrier.valid || probe_b_at_carrier.distance > max_abs_distance) {
                            continue;
                        }

                        const ChVector3d normal_world_a =
                            shape_a_frame_abs.TransformDirectionLocalToParent(normal_shape_a).GetNormalized();
                        const NormalProjectionResult projection = ProjectPointToShapeAlongNormal(
                            shape_b, shape_b_frame_abs, surface_world_a, normal_world_a, projection_search, projection_step);
                        if (!projection.valid || std::abs(projection.signed_gap) > projection_search) {
                            continue;
                        }

                        const ChSDFProbeResult probe_a_on_b =
                            shape_a.ProbeLocal(shape_a_frame_abs.TransformPointParentToLocal(projection.projected_world));
                        if (!probe_a_on_b.valid) {
                            continue;
                        }

                        const ChVector3d normal_world_b = shape_b_frame_abs.TransformDirectionLocalToParent(
                                                               SafeNormalized(projection.projected_probe_b.normal_world, -normal_world_a))
                                                               .GetNormalized();
                        const double normal_opposition = -Vdot(normal_world_a, normal_world_b);
                        if (normal_opposition < settings.min_opposed_normal_cosine) {
                            continue;
                        }

                        const double combined_gap = probe_b_at_carrier.distance + probe_a_on_b.distance;
                        if (settings.max_combined_gap >= 0 && combined_gap > settings.max_combined_gap) {
                            continue;
                        }

                        if (item.second.empty()) {
                            continue;
                        }

                        ChSDFBrickPairRegionSample sample;
                        sample.brick_a_index = brick_a_index;
                        sample.brick_b_index = item.second.front()->brick_b_index;
                        sample.carrier_axis = axis;
                        sample.coord = coord;
                        sample.point_world = surface_world_a;
                        sample.point_shape_a = surface_shape_a;
                        sample.point_shape_b = projection.projected_shape_b;
                        sample.surface_world_a = surface_world_a;
                        sample.surface_world_b = projection.projected_world;
                        sample.surface_shape_a = surface_shape_a;
                        sample.surface_shape_b = projection.projected_shape_b;
                        sample.normal_world_a = normal_world_a;
                        sample.normal_world_b = normal_world_b;
                        sample.contact_normal_world = normal_world_a;
                        sample.distance_a = 0.0;
                        sample.distance_b = probe_b_at_carrier.distance;
                        sample.normal_gap = projection.signed_gap;
                        sample.combined_gap = combined_gap;
                        sample.normal_opposition = normal_opposition;
                        sample.area_weight =
                            FaceAreaForAxis(axis, voxel_size_a) / std::max(std::abs(Vdot(normal_shape_a, AxisUnit(axis))), 1.0e-12);
                        sample.probe_a = surface_probe_a;
                        sample.probe_b = projection.projected_probe_b;

                        const CarrierKey key = MakeCarrierKey(axis, surface_probe_a.point_index);
                        const auto it = sample_lookup.find(key);
                        if (it == sample_lookup.end()) {
                            sample_lookup.emplace(key, samples.size());
                            samples.push_back(sample);
                        } else {
                            const std::size_t old_index = it->second;
                            const double old_score =
                                std::abs(samples[old_index].normal_gap) + std::abs(samples[old_index].combined_gap);
                            const double new_score = std::abs(sample.normal_gap) + std::abs(sample.combined_gap);
                            if (new_score < old_score) {
                                samples[old_index] = sample;
                            }
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
    const double sample_spacing =
        info_a.valid && MinComponent(info_a.voxel_size) > 0 ? MinComponent(info_a.voxel_size)
                                                            : ResolveSampleSpacing(info_a, info_b, settings);
    const double neighbor_distance = NeighborDistanceThreshold(sample_spacing, settings.neighbor_mode);

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

            for (std::size_t neighbor_index = 0; neighbor_index < samples.size(); ++neighbor_index) {
                if (visited[neighbor_index] || neighbor_index == current) {
                    continue;
                }

                if ((samples[neighbor_index].point_world - sample.point_world).Length() > neighbor_distance) {
                    continue;
                }

                if (!PairSamplesAreConnected(sample, samples[neighbor_index], settings)) {
                    continue;
                }

                visited[neighbor_index] = 1;
                queue.push(neighbor_index);
            }
        }

        FilterBrickPairRegionByDominantNormal(region, settings.min_region_normal_cosine);
        if (region.samples.size() < settings.min_region_samples) {
            continue;
        }

        region.world_bounds = ChAABB();
        region.patch_bounds = ChAABB();
        region.shape_a_bounds = ChAABB();
        region.shape_b_bounds = ChAABB();
        region.centroid_world = VNULL;
        region.centroid_shape_a = VNULL;
        region.centroid_shape_b = VNULL;
        region.mean_normal_world = VNULL;
        brick_flags_a.clear();
        brick_flags_b.clear();

        for (const auto& sample : region.samples) {
            region.world_bounds += sample.point_world;
            region.shape_a_bounds += sample.surface_shape_a;
            region.shape_b_bounds += sample.surface_shape_b;
            region.centroid_world += sample.point_world;
            region.centroid_shape_a += sample.surface_shape_a;
            region.centroid_shape_b += sample.surface_shape_b;
            region.mean_normal_world += sample.contact_normal_world;
            brick_flags_a[sample.brick_a_index] = 1;
            brick_flags_b[sample.brick_b_index] = 1;
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
