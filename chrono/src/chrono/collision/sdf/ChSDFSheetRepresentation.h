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

#ifndef CH_SDF_SHEET_REPRESENTATION_H
#define CH_SDF_SHEET_REPRESENTATION_H

#include <cstddef>
#include <vector>

#include "chrono/collision/sdf/ChSDFContactWrench.h"
#include "chrono/core/ChVector2.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// Settings for the zero-thickness sheet layer collapsed out of the band-mechanics cells.
struct ChApi ChSDFSheetCollapseSettings {
    bool enable = true;
    bool use_local_fiber_projection = false;
    double fiber_lateral_tolerance = -1;
    double fiber_normal_cosine = 0.7;
    double min_sheet_sample_area = 0;
    int patch_neighbor_mode = 8;
};

/// One collapsed sheet sample obtained by merging a normal-direction fiber of band cells.
struct ChApi ChSDFSheetFiberSample {
    std::size_t region_id = 0;
    std::size_t fiber_id = 0;

    int dominant_axis = -1;
    ChVector2i lattice_coord = ChVector2i(0, 0);

    double measure_area = 0;
    double footprint_area = 0;

    ChVector3d centroid_world = VNULL;
    ChVector3d pressure_center_world = VNULL;
    ChVector3d normal_world = VNULL;
    ChVector3d force_world = VNULL;

    ChAABB source_bounds_world;

    std::vector<std::size_t> source_sample_indices;

    bool HasSupport() const { return measure_area > 0; }
};

/// One connected component on the collapsed sheet lattice.
struct ChApi ChSDFSheetPatch {
    std::size_t patch_id = 0;
    std::size_t support_columns = 0;

    double measure_area = 0;
    double footprint_area = 0;

    ChVector3d centroid_world = VNULL;
    ChVector3d pressure_center_world = VNULL;
    ChVector3d mean_normal_world = VNULL;

    ChAABB bounds_world;
    ChAABB support_bbox_world;

    std::vector<std::size_t> sample_indices;

    bool HasSamples() const { return !sample_indices.empty(); }
};

/// Sheet-layer counterpart of one band region.
struct ChApi ChSDFSheetRegion {
    std::size_t region_id = 0;
    std::size_t persistent_id = 0;
    std::size_t patch_count = 0;

    int dominant_axis = -1;

    double measure_area = 0;
    double footprint_area = 0;
    double largest_patch_area = 0;

    ChVector3d centroid_world = VNULL;
    ChVector3d pressure_center_world = VNULL;
    ChVector3d mean_normal_world = VNULL;

    ChAABB bounds_world;
    ChAABB support_bbox_world;

    std::vector<ChSDFSheetFiberSample> samples;
    std::vector<ChSDFSheetPatch> patches;

    bool HasSamples() const { return !samples.empty(); }
};

/// Pair-level sheet geometry result built from the band-mechanics output.
struct ChApi ChSDFSheetShapePairResult {
    double occupied_area = 0;
    double band_area = 0;
    double sheet_area = 0;
    double sheet_footprint_area = 0;
    double largest_patch_area = 0;

    ChVector3d sheet_center_world = VNULL;
    ChVector3d pressure_center_world = VNULL;
    ChAABB support_bbox_world;

    std::size_t occupied_cells = 0;
    std::size_t collapsed_samples = 0;
    std::size_t patch_count = 0;

    std::vector<ChSDFSheetRegion> regions;

    bool HasSamples() const { return collapsed_samples > 0; }
};

/// Build a zero-thickness sheet representation from the current band-mechanics result.
class ChApi ChSDFSheetBuilder {
  public:
    /// Dominant world axis used by the minimal column-collapse implementation.
    static int ComputeDominantAxis(const ChVector3d& normal_world);

    /// Project an integer cell coordinate to the two in-plane coordinates orthogonal to the dominant axis.
    static ChVector2i ComputeLatticeCoord(const ChVector3i& coord, int dominant_axis);

    /// Collapse one band region into a minimal dominant-axis sheet representation.
    static ChSDFSheetRegion BuildRegion(const ChSDFBrickPairWrenchResult& band_region,
                                        const ChSDFSheetCollapseSettings& settings);

    /// Collapse multiple band regions independently.
    static std::vector<ChSDFSheetRegion> BuildRegions(const std::vector<ChSDFBrickPairWrenchResult>& band_regions,
                                                      const ChSDFSheetCollapseSettings& settings);

    /// Build the pair-level sheet result from the aggregated band-mechanics output.
    static ChSDFSheetShapePairResult BuildShapePair(const ChSDFShapePairContactResult& band_result,
                                                    const ChSDFSheetCollapseSettings& settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
