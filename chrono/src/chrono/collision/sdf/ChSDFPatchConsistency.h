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

#ifndef CH_SDF_PATCH_CONSISTENCY_H
#define CH_SDF_PATCH_CONSISTENCY_H

#include <cstddef>
#include <vector>

#include "chrono/collision/sdf/ChSDFContactWrench.h"
#include "chrono/collision/sdf/ChSDFSheetRepresentation.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

struct ChApi ChSDFPatchConsistencySettings {
    bool enable = true;
    bool clamp_alpha = true;
    bool use_sample_index_back_projection = true;
    bool add_unassigned_band_residual_patch = true;
    bool use_first_moment_consistent_correction = true;
    bool use_patch_local_redistribution = true;

    double min_alpha = 0.25;
    double max_alpha = 4.0;
    double min_band_area = 1.0e-12;
};

struct ChApi ChSDFPatchBandAggregate {
    std::size_t region_id = 0;
    std::size_t patch_id = 0;

    bool residual_patch = false;

    std::vector<std::size_t> source_sheet_sample_indices;
    std::vector<std::size_t> source_band_sample_indices;

    double band_area = 0;
    double integrated_pressure = 0;

    ChVector3d centroid_world = VNULL;
    ChVector3d pressure_center_world = VNULL;
    ChWrenchd wrench_world_a = {VNULL, VNULL};
    ChWrenchd wrench_world_b = {VNULL, VNULL};
    ChAABB support_bbox_world;
};

struct ChApi ChSDFPatchConsistencyResult {
    std::size_t region_id = 0;
    std::size_t patch_id = 0;

    bool residual_patch = false;

    double band_area = 0;
    double sheet_area = 0;
    double alpha = 1;
    double integrated_pressure_band = 0;
    double integrated_pressure_corrected = 0;

    ChWrenchd wrench_world_a_band = {VNULL, VNULL};
    ChWrenchd wrench_world_b_band = {VNULL, VNULL};
    ChWrenchd wrench_world_a_corrected = {VNULL, VNULL};
    ChWrenchd wrench_world_b_corrected = {VNULL, VNULL};

    ChVector3d centroid_world = VNULL;
    ChVector3d pressure_center_world = VNULL;
    ChVector3d sheet_pressure_center_world = VNULL;
    ChVector3d corrected_pressure_center_world = VNULL;
    ChAABB support_bbox_world;
};

struct ChApi ChSDFPatchConsistencyPairResult {
    std::vector<ChSDFPatchConsistencyResult> patches;

    ChWrenchd wrench_world_a_band = {VNULL, VNULL};
    ChWrenchd wrench_world_b_band = {VNULL, VNULL};
    ChWrenchd wrench_world_a_corrected = {VNULL, VNULL};
    ChWrenchd wrench_world_b_corrected = {VNULL, VNULL};

    double total_band_area = 0;
    double total_sheet_area = 0;
    double total_corrected_area = 0;
    double integrated_pressure_band = 0;
    double integrated_pressure_corrected = 0;
    double mean_alpha = 1;

    std::size_t corrected_patch_count = 0;
};

class ChApi ChSDFPatchConsistencyBridge {
  public:
    static std::vector<ChSDFPatchBandAggregate> BuildPatchBandAggregates(
        const ChSDFBrickPairWrenchResult& band_region,
        const ChSDFSheetRegion& sheet_region,
        const ChFrameMoving<>& shape_a_frame_abs,
        const ChFrameMoving<>& shape_b_frame_abs,
        const ChSDFPatchConsistencySettings& settings);

    static ChSDFPatchConsistencyResult BuildPatchConsistencyResult(
        const ChSDFPatchBandAggregate& band_patch,
        const ChSDFBrickPairWrenchResult& band_region,
        const ChSDFSheetPatch& sheet_patch,
        const ChFrameMoving<>& shape_a_frame_abs,
        const ChFrameMoving<>& shape_b_frame_abs,
        const ChSDFPatchConsistencySettings& settings);

    static ChSDFPatchConsistencyPairResult BuildPairConsistencyResult(
        const ChSDFShapePairContactResult& band_result,
        const ChSDFSheetShapePairResult& sheet_result,
        const ChSDFPatchConsistencySettings& settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
