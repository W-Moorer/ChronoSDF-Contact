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

#ifndef CH_SDF_CONTACT_REGION_H
#define CH_SDF_CONTACT_REGION_H

#include <vector>

#include "chrono/collision/ChCollisionShapeSDF.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// One narrow-band brick retained by the patch-slab broadphase.
struct ChApi ChSDFPatchBrickCandidate {
    std::size_t brick_index = 0;

    double slab_distance = 0;

    ChVector3d center_patch = VNULL;
    ChVector3d size_patch = VNULL;

    ChSDFLeafBrick brick;
};

/// One voxel-center sample used for SDF-only contact-region partitioning.
struct ChApi ChSDFContactRegionSample {
    std::size_t brick_index = 0;

    ChVector3i coord = ChVector3i(0, 0, 0);
    ChVector3d point_shape = VNULL;
    ChVector3d point_patch = VNULL;

    double normal_alignment = 0;

    ChSDFProbeResult probe;
};

/// One connected SDF contact region found inside the patch slab.
struct ChApi ChSDFContactRegion {
    std::size_t region_id = 0;

    ChFrame<> patch_frame_shape;
    ChAABB shape_bounds;
    ChAABB patch_bounds;

    ChVector3d centroid_shape = VNULL;
    ChVector3d centroid_patch = VNULL;
    ChVector3d mean_normal_shape = VNULL;

    std::vector<std::size_t> brick_indices;
    std::vector<ChSDFContactRegionSample> samples;

    ChSDFContactPatchSampler::Settings suggested_patch_settings;

    bool HasSamples() const { return !samples.empty(); }
};

/// SDF-only broadphase and topology partitioning for a rectangular probe patch.
/// The workflow is:
/// - extract sparse NanoVDB leaf bricks
/// - retain the bricks intersecting an oriented patch slab
/// - sample voxel centers inside the slab and near the zero level set
/// - split them into independent connected regions
class ChApi ChSDFContactRegionBuilder {
  public:
    struct Settings {
        double slab_half_thickness = -1;
        double brick_margin = 0;
        double sample_max_abs_distance = -1;
        double max_distance_jump = -1;
        double min_abs_normal_alignment = 0;
        double min_neighbor_normal_cosine = 0.9;

        std::size_t min_region_samples = 1;
        int neighbor_mode = 6;

        bool require_zero_crossing = false;
    };

    /// Retain the sparse NanoVDB bricks intersecting the patch slab.
    static std::vector<ChSDFPatchBrickCandidate> FindPatchBrickCandidates(
        const ChCollisionShapeSDF& shape,
        const ChFrame<>& patch_frame_shape,
        const ChSDFContactPatchSampler::Settings& patch_settings,
        const Settings& settings);

    /// Build connected contact regions from voxel-center samples extracted out of the retained bricks.
    static std::vector<ChSDFContactRegion> BuildPatchRegions(const ChCollisionShapeSDF& shape,
                                                             const ChFrame<>& patch_frame_shape,
                                                             const ChSDFContactPatchSampler::Settings& patch_settings,
                                                             const Settings& settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
