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

#include "chrono/collision/sdf/ChSDFBrickPair.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace chrono {
namespace {

double AxisSeparation(double amin, double amax, double bmin, double bmax) {
    if (amax < bmin) {
        return bmin - amax;
    }
    if (bmax < amin) {
        return amin - bmax;
    }
    return 0.0;
}

double ComputeAABBSeparation(const ChAABB& a, const ChAABB& b) {
    if (a.IsInverted() || b.IsInverted()) {
        return std::numeric_limits<double>::infinity();
    }

    const double dx = AxisSeparation(a.min.x(), a.max.x(), b.min.x(), b.max.x());
    const double dy = AxisSeparation(a.min.y(), a.max.y(), b.min.y(), b.max.y());
    const double dz = AxisSeparation(a.min.z(), a.max.z(), b.min.z(), b.max.z());
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

ChAABB Expanded(const ChAABB& aabb, double margin) {
    if (aabb.IsInverted() || margin <= 0) {
        return aabb;
    }

    return ChAABB(aabb.min - ChVector3d(margin), aabb.max + ChVector3d(margin));
}

bool Overlap1D(double amin, double amax, double bmin, double bmax) {
    return !(amax < bmin || bmax < amin);
}

bool IntersectAABB(const ChAABB& a, const ChAABB& b) {
    return Overlap1D(a.min.x(), a.max.x(), b.min.x(), b.max.x()) &&
           Overlap1D(a.min.y(), a.max.y(), b.min.y(), b.max.y()) &&
           Overlap1D(a.min.z(), a.max.z(), b.min.z(), b.max.z());
}

ChAABB OverlapRegion(const ChAABB& a, const ChAABB& b) {
    if (!IntersectAABB(a, b)) {
        return ChAABB();
    }

    return ChAABB(ChVector3d(std::max(a.min.x(), b.min.x()), std::max(a.min.y(), b.min.y()), std::max(a.min.z(), b.min.z())),
                  ChVector3d(std::min(a.max.x(), b.max.x()), std::min(a.max.y(), b.max.y()), std::min(a.max.z(), b.max.z())));
}

double OverlapVolume(const ChAABB& aabb) {
    if (aabb.IsInverted()) {
        return 0.0;
    }

    const ChVector3d size = aabb.Size();
    return std::max(0.0, size.x()) * std::max(0.0, size.y()) * std::max(0.0, size.z());
}

}  // namespace

std::vector<ChSDFBrickWorldProxy> ChSDFBrickPairBroadphase::BuildWorldBricks(const ChCollisionShapeSDF& shape,
                                                                              const ChFrame<>& shape_frame_abs,
                                                                              const Settings& settings) {
    std::vector<ChSDFBrickWorldProxy> world_bricks;

    if (!shape.IsLoaded()) {
        return world_bricks;
    }

    const auto& bricks = shape.GetLeafBricks();
    world_bricks.reserve(bricks.size());

    for (const auto& brick : bricks) {
        if (!brick.valid || brick.world_bounds.IsInverted()) {
            continue;
        }

        if (settings.max_min_abs_value >= 0 && brick.min_abs_value > settings.max_min_abs_value) {
            continue;
        }

        ChSDFBrickWorldProxy proxy;
        proxy.brick_index = brick.index;
        proxy.brick = brick;
        proxy.world_bounds = brick.world_bounds.Transform(shape_frame_abs);
        proxy.center_world = proxy.world_bounds.Center();
        world_bricks.push_back(proxy);
    }

    return world_bricks;
}

std::vector<ChSDFBrickPairCandidate> ChSDFBrickPairBroadphase::FindBrickPairs(const ChCollisionShapeSDF& shape_a,
                                                                               const ChFrame<>& shape_a_frame_abs,
                                                                               const ChCollisionShapeSDF& shape_b,
                                                                               const ChFrame<>& shape_b_frame_abs,
                                                                               const Settings& settings) {
    std::vector<ChSDFBrickPairCandidate> pairs;

    const auto bricks_a = BuildWorldBricks(shape_a, shape_a_frame_abs, settings);
    const auto bricks_b = BuildWorldBricks(shape_b, shape_b_frame_abs, settings);
    if (bricks_a.empty() || bricks_b.empty()) {
        return pairs;
    }

    pairs.reserve(std::min<std::size_t>(bricks_a.size() * bricks_b.size(), 4096));

    for (const auto& brick_a : bricks_a) {
        const ChAABB query_a = Expanded(brick_a.world_bounds, settings.world_margin);

        for (const auto& brick_b : bricks_b) {
            const ChAABB query_b = Expanded(brick_b.world_bounds, settings.world_margin);

            const double center_distance = (brick_a.center_world - brick_b.center_world).Length();
            if (settings.max_center_distance >= 0 && center_distance > settings.max_center_distance) {
                continue;
            }

            const double separation_distance = ComputeAABBSeparation(query_a, query_b);
            if (settings.max_separation_distance >= 0 && separation_distance > settings.max_separation_distance) {
                continue;
            }

            if (settings.require_world_overlap && !IntersectAABB(query_a, query_b)) {
                continue;
            }

            ChSDFBrickPairCandidate pair;
            pair.brick_a_index = brick_a.brick_index;
            pair.brick_b_index = brick_b.brick_index;
            pair.center_distance = center_distance;
            pair.separation_distance = separation_distance;
            pair.min_abs_value_bound = std::max(brick_a.brick.min_abs_value, brick_b.brick.min_abs_value);
            pair.overlap_world = OverlapRegion(query_a, query_b);
            pair.overlap_volume = OverlapVolume(pair.overlap_world);
            pair.brick_a = brick_a;
            pair.brick_b = brick_b;
            pairs.push_back(pair);
        }
    }

    std::stable_sort(pairs.begin(), pairs.end(), [](const ChSDFBrickPairCandidate& a, const ChSDFBrickPairCandidate& b) {
        if (a.separation_distance != b.separation_distance) {
            return a.separation_distance < b.separation_distance;
        }
        if (a.overlap_volume != b.overlap_volume) {
            return a.overlap_volume > b.overlap_volume;
        }
        if (a.min_abs_value_bound != b.min_abs_value_bound) {
            return a.min_abs_value_bound < b.min_abs_value_bound;
        }
        return a.center_distance < b.center_distance;
    });

    return pairs;
}

}  // end namespace chrono
