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

#include "chrono/collision/sdf/ChSDFPatchCandidate.h"

#include <algorithm>

#include "chrono/core/ChMatrix33.h"

namespace chrono {
namespace {

ChVector3d ProjectOntoTangent(const ChVector3d& direction, const ChVector3d& normal) {
    return direction - normal * Vdot(direction, normal);
}

bool IsDuplicateDirection(const std::vector<ChVector3d>& directions, const ChVector3d& candidate, double cosine_tol) {
    for (const auto& dir : directions) {
        if (std::abs(Vdot(dir, candidate)) >= cosine_tol) {
            return true;
        }
    }

    return false;
}

void TryAddDirection(std::vector<ChVector3d>& directions,
                     const ChVector3d& suggestion,
                     const ChVector3d& normal,
                     const ChSDFPatchCandidateSettings& settings) {
    ChVector3d tangent = ProjectOntoTangent(suggestion, normal);
    const double tangent_length = tangent.Length();
    if (tangent_length <= settings.min_tangent_length) {
        return;
    }

    tangent /= tangent_length;
    if (IsDuplicateDirection(directions, tangent, settings.duplicate_cosine)) {
        return;
    }

    directions.push_back(tangent);
}

double DistanceScore(double abs_distance) {
    return 1.0 / (1.0 + abs_distance);
}

}  // namespace

std::vector<ChSDFPatchCandidate> ChSDFPatchCandidateGenerator::GenerateCandidates(
    const ChCollisionShapeSDF& sdf_shape,
    const ChFrameMoving<>& sdf_frame_abs,
    const ChFrameMoving<>& patch_body_frame_abs,
    const ChSDFPatchCandidateSettings& settings) {
    std::vector<ChSDFPatchCandidate> candidates;

    if (!sdf_shape.IsLoaded()) {
        return candidates;
    }

    const ChVector3d seed_point_abs = patch_body_frame_abs.TransformPointLocalToParent(settings.seed_point_patch);
    const ChVector3d seed_point_shape = sdf_frame_abs.TransformPointParentToLocal(seed_point_abs);
    const ChSDFProbeResult probe = sdf_shape.ProbeLocal(seed_point_shape);

    if (!probe.valid) {
        return candidates;
    }

    const double abs_distance = std::abs(probe.distance);
    if (settings.max_abs_distance > 0 && abs_distance > settings.max_abs_distance) {
        return candidates;
    }

    const double normal_length = probe.normal_world.Length();
    if (normal_length <= settings.min_normal_length) {
        return candidates;
    }

    const ChVector3d normal_shape = probe.normal_world / normal_length;
    const ChVector3d surface_point_shape =
        settings.project_to_surface ? (seed_point_shape - normal_shape * probe.distance) : seed_point_shape;
    const ChVector3d patch_origin_shape = surface_point_shape + normal_shape * settings.surface_offset;

    const ChVector3d patch_seed_speed_abs = patch_body_frame_abs.PointSpeedLocalToParent(settings.seed_point_patch);
    const ChVector3d sdf_surface_speed_abs = sdf_frame_abs.PointSpeedLocalToParent(surface_point_shape);
    const ChVector3d relative_velocity_abs = patch_seed_speed_abs - sdf_surface_speed_abs;
    const ChVector3d relative_velocity_shape = sdf_frame_abs.TransformDirectionParentToLocal(relative_velocity_abs);
    const ChVector3d relative_tangent_shape = ProjectOntoTangent(relative_velocity_shape, normal_shape);
    const double relative_tangent_speed = relative_tangent_shape.Length();
    const ChVector3d relative_tangent_dir =
        relative_tangent_speed > settings.min_tangent_length ? (relative_tangent_shape / relative_tangent_speed) : VNULL;

    std::vector<ChVector3d> tangents;
    tangents.reserve(8);

    if (settings.include_relative_velocity) {
        TryAddDirection(tangents, relative_velocity_shape, normal_shape, settings);
    }

    if (settings.include_patch_axes) {
        TryAddDirection(tangents,
                        sdf_frame_abs.TransformDirectionParentToLocal(patch_body_frame_abs.GetRotMat().GetAxisX()),
                        normal_shape, settings);
        TryAddDirection(tangents,
                        sdf_frame_abs.TransformDirectionParentToLocal(patch_body_frame_abs.GetRotMat().GetAxisY()),
                        normal_shape, settings);
        TryAddDirection(tangents,
                        sdf_frame_abs.TransformDirectionParentToLocal(patch_body_frame_abs.GetRotMat().GetAxisZ()),
                        normal_shape, settings);
    }

    if (settings.include_world_axes) {
        TryAddDirection(tangents, sdf_frame_abs.TransformDirectionParentToLocal(VECT_X), normal_shape, settings);
        TryAddDirection(tangents, sdf_frame_abs.TransformDirectionParentToLocal(VECT_Y), normal_shape, settings);
        TryAddDirection(tangents, sdf_frame_abs.TransformDirectionParentToLocal(VECT_Z), normal_shape, settings);
    }

    if (tangents.empty()) {
        TryAddDirection(tangents, VECT_X, normal_shape, settings);
    }

    for (const auto& tangent_shape : tangents) {
        ChMatrix33<> rot;
        rot.SetFromAxisZ(normal_shape, tangent_shape);

        ChSDFPatchCandidate candidate;
        candidate.valid = true;
        candidate.signed_distance = probe.distance;
        candidate.tangent_speed = relative_tangent_speed;
        candidate.tangent_alignment =
            (relative_tangent_speed > settings.min_tangent_length) ? std::abs(Vdot(tangent_shape, relative_tangent_dir))
                                                                   : 0.0;

        candidate.seed_point_patch = settings.seed_point_patch;
        candidate.seed_point_abs = seed_point_abs;
        candidate.seed_point_shape = seed_point_shape;
        candidate.surface_point_shape = surface_point_shape;
        candidate.normal_shape = normal_shape;
        candidate.tangent_shape = tangent_shape;
        candidate.relative_velocity_shape = relative_velocity_shape;
        candidate.patch_frame_shape = ChFrame<>(patch_origin_shape, rot);
        candidate.patch_frame_abs =
            ChFrame<>(sdf_frame_abs.TransformLocalToParent(ChFrameMoving<>(candidate.patch_frame_shape)).GetCoordsys());
        candidate.patch_frame_patch =
            ChFrame<>(patch_body_frame_abs.TransformParentToLocal(ChFrameMoving<>(candidate.patch_frame_abs)).GetCoordsys());
        candidate.probe = probe;

        candidate.score = DistanceScore(abs_distance) + candidate.tangent_alignment + 0.1 * candidate.tangent_speed;
        if (settings.include_relative_velocity && candidate.tangent_alignment > 0.999) {
            candidate.score += 1.0;
        }

        candidates.push_back(candidate);
    }

    std::stable_sort(candidates.begin(), candidates.end(),
                     [](const ChSDFPatchCandidate& a, const ChSDFPatchCandidate& b) { return a.score > b.score; });

    if (settings.max_candidates > 0 && candidates.size() > settings.max_candidates) {
        candidates.resize(settings.max_candidates);
    }

    return candidates;
}

}  // end namespace chrono
