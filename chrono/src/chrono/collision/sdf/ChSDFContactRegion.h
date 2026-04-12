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
#include "chrono/collision/sdf/ChSDFBrickPair.h"

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

/// One dual-SDF sample extracted from a coarse brick pair overlap.
struct ChApi ChSDFBrickPairRegionSample {
    std::size_t brick_a_index = 0;
    std::size_t brick_b_index = 0;

    ChVector3i coord = ChVector3i(0, 0, 0);

    ChVector3d point_world = VNULL;
    ChVector3d point_patch = VNULL;
    ChVector3d point_shape_a = VNULL;
    ChVector3d point_shape_b = VNULL;

    ChVector3d surface_world_a = VNULL;
    ChVector3d surface_world_b = VNULL;
    ChVector3d surface_shape_a = VNULL;
    ChVector3d surface_shape_b = VNULL;

    ChVector3d normal_world_a = VNULL;
    ChVector3d normal_world_b = VNULL;
    ChVector3d contact_normal_world = VNULL;

    double distance_a = 0;
    double distance_b = 0;
    double combined_gap = 0;
    double normal_opposition = 0;

    ChSDFProbeResult probe_a;
    ChSDFProbeResult probe_b;
};

/// One connected region extracted from dual-SDF brick-pair samples.
struct ChApi ChSDFBrickPairRegion {
    std::size_t region_id = 0;

    ChFrame<> contact_frame_world;
    ChFrame<> contact_frame_shape_a;
    ChFrame<> contact_frame_shape_b;

    ChAABB world_bounds;
    ChAABB patch_bounds;
    ChAABB shape_a_bounds;
    ChAABB shape_b_bounds;

    ChVector3d centroid_world = VNULL;
    ChVector3d centroid_shape_a = VNULL;
    ChVector3d centroid_shape_b = VNULL;
    ChVector3d mean_normal_world = VNULL;

    double sample_spacing = 0;

    std::vector<std::size_t> brick_a_indices;
    std::vector<std::size_t> brick_b_indices;
    std::vector<ChSDFBrickPairRegionSample> samples;

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
        double sample_spacing = -1;
        double sample_max_abs_distance = -1;
        double max_combined_gap = -1;
        double max_distance_jump = -1;
        double min_abs_normal_alignment = 0;
        double min_opposed_normal_cosine = 0;
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

    /// Build dual-SDF near-contact samples from the overlap of coarse brick pairs.
    static std::vector<ChSDFBrickPairRegionSample> BuildBrickPairSamples(
        const ChCollisionShapeSDF& shape_a,
        const ChFrame<>& shape_a_frame_abs,
        const ChCollisionShapeSDF& shape_b,
        const ChFrame<>& shape_b_frame_abs,
        const std::vector<ChSDFBrickPairCandidate>& brick_pairs,
        const Settings& settings);

    /// Build connected dual-SDF contact regions from the overlap of coarse brick pairs.
    static std::vector<ChSDFBrickPairRegion> BuildBrickPairRegions(const ChCollisionShapeSDF& shape_a,
                                                                   const ChFrame<>& shape_a_frame_abs,
                                                                   const ChCollisionShapeSDF& shape_b,
                                                                   const ChFrame<>& shape_b_frame_abs,
                                                                   const std::vector<ChSDFBrickPairCandidate>& brick_pairs,
                                                                   const Settings& settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
