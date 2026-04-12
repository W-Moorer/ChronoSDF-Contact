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

#ifndef CH_SDF_PATCH_CANDIDATE_H
#define CH_SDF_PATCH_CANDIDATE_H

#include <vector>

#include "chrono/collision/ChCollisionShapeSDF.h"
#include "chrono/core/ChFrameMoving.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// One automatically generated probe patch candidate.
struct ChApi ChSDFPatchCandidate {
    bool valid = false;
    double score = 0;

    double signed_distance = 0;
    double tangent_speed = 0;
    double tangent_alignment = 0;

    ChVector3d seed_point_patch = VNULL;
    ChVector3d seed_point_abs = VNULL;
    ChVector3d seed_point_shape = VNULL;

    ChVector3d surface_point_shape = VNULL;
    ChVector3d normal_shape = VNULL;
    ChVector3d tangent_shape = VNULL;
    ChVector3d relative_velocity_shape = VNULL;

    ChFrame<> patch_frame_shape;
    ChFrame<> patch_frame_patch;
    ChFrame<> patch_frame_abs;

    ChSDFProbeResult probe;
};

/// Settings for automatic patch candidate generation.
/// This first version generates local tangent-plane patches around one probe seed point on the probe body.
struct ChApi ChSDFPatchCandidateSettings {
    ChVector3d seed_point_patch = VNULL;

    double surface_offset = 0;
    double max_abs_distance = -1;
    double min_normal_length = 1e-10;
    double min_tangent_length = 1e-10;
    double duplicate_cosine = 0.995;

    bool project_to_surface = true;
    bool include_relative_velocity = true;
    bool include_patch_axes = true;
    bool include_world_axes = true;

    std::size_t max_candidates = 4;
};

/// Automatic candidate generator for local SDF contact patches.
class ChApi ChSDFPatchCandidateGenerator {
  public:
    /// Generate candidates using explicitly provided moving frames.
    static std::vector<ChSDFPatchCandidate> GenerateCandidates(const ChCollisionShapeSDF& sdf_shape,
                                                               const ChFrameMoving<>& sdf_frame_abs,
                                                               const ChFrameMoving<>& patch_body_frame_abs,
                                                               const ChSDFPatchCandidateSettings& settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
