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
    /// Default Step 6 v2 path from supplement5.
    bool use_step6_v2 = true;
    /// Legacy v1 path switch. Only used when use_step6_v2 is false.
    bool use_local_fiber_projection = false;
    bool enable_sheet_diagnostics = true;
    bool allow_dominant_axis_fallback = true;

    /// Characteristic sheet spacing h. If <= 0, use the band sample spacing.
    double fiber_lateral_tolerance = -1;
    /// Cosine threshold used in fiber clustering.
    double fiber_normal_cosine = 0.9396926207859084;   // cos(20 deg)
    /// Tangential proximity threshold for the v2 fiber graph. If <= 0, use 1.0 * h.
    double fiber_tangent_tolerance = -1;
    /// Normal-direction coplanarity threshold for the v2 fiber graph. If <= 0, use 0.5 * h.
    double fiber_plane_tolerance = -1;
    /// Connection radius for the v2 patch graph. If <= 0, use 2.0 * h.
    double patch_connection_radius = -1;
    /// Cosine threshold used in the v2 patch graph.
    double patch_normal_cosine = 0.9063077870366499;   // cos(25 deg)
    /// Minimum gradient norm accepted by the local sheet-seed reconstruction.
    double min_gradient_norm = 1.0e-8;
    double min_sheet_sample_area = 0;

    /// Diagnostics/fallback thresholds.
    double min_sheet_area_ratio = 0.5;
    double max_sheet_area_ratio = 1.5;
    std::size_t max_patch_count_before_fallback = 0;
    std::size_t max_fiber_count_before_fallback = 0;

    /// Legacy projected-lattice neighbor mode. Only used by the dominant-axis fallback.
    int patch_neighbor_mode = 8;
};

/// One local sheet seed recovered from a band-mechanics sample.
struct ChApi ChSDFSheetSeed {
    std::size_t region_id = 0;
    std::size_t source_sample_index = 0;

    int source_carrier_axis = -1;
    ChVector3i source_coord = ChVector3i(0, 0, 0);
    ChVector3d seed_world = VNULL;
    ChVector3d seed_normal_world = VNULL;
    ChVector3d band_point_world = VNULL;

    double measure_area = 0;
    double pressure = 0;
    double h_value = 0;

    ChVector3d grad_h_world = VNULL;
    ChVector3d force_world = VNULL;

    ChAABB source_bounds_world;

    bool valid = false;
};

/// Explicit local footprint carried by a recovered sheet sample or patch.
struct ChApi ChSDFSheetLocalFootprint {
    ChVector3d origin_world = VNULL;
    ChVector3d tangent_u_world = VNULL;
    ChVector3d tangent_v_world = VNULL;
    ChVector2d centroid_uv = ChVector2d(0, 0);
    double area = 0;
    std::vector<ChVector2d> polygon_uv;

    bool HasPolygon() const { return polygon_uv.size() >= 3 && area > 1.0e-16; }
};

/// One raw support witness retained before patch-level footprint reconstruction.
struct ChApi ChSDFSheetSupportEvidence {
    int source_carrier_axis = -1;
    ChVector3i source_coord = ChVector3i(0, 0, 0);
    ChVector3d point_world = VNULL;
    ChVector3d seed_world = VNULL;
    ChVector3d normal_world = VNULL;
    double measure_area = 0;
    std::size_t source_sample_index = 0;
};

/// One raw support-layer sample assigned to a unique patch-plane lattice site.
struct ChApi ChSDFPatchPlaneLayerSample {
    int source_carrier_axis = -1;
    ChVector3i source_coord = ChVector3i(0, 0, 0);
    ChVector3d point_world = VNULL;
    ChVector3d seed_world = VNULL;
    ChVector3d normal_world = VNULL;
    double measure_area = 0;
    double normal_depth = 0;
    std::size_t source_sample_index = 0;
};

/// One layered occupancy container before collapsing to a single sheet cell.
struct ChApi ChSDFPatchPlaneLayeredCell {
    int carrier_axis = -1;
    ChVector2i cell_ij = ChVector2i(0, 0);
    ChVector2d center_uv = ChVector2d(0, 0);
    ChVector2d half_extents_uv = ChVector2d(0, 0);
    double measure_area_sum = 0;
    std::vector<ChSDFPatchPlaneLayerSample> samples;
};

/// One discrete support-sheet cell on the local patch plane built from raw support evidence.
struct ChApi ChSDFPatchPlaneSupportCell {
    int carrier_axis = -1;
    ChVector2i cell_ij = ChVector2i(0, 0);
    ChVector2d center_uv = ChVector2d(0, 0);
    ChVector2d half_extents_uv = ChVector2d(0, 0);
    double measure_area = 0;
    std::size_t layer_count = 0;
    bool occupied = false;
    bool first_layer = true;
    bool shell = false;

    ChVector3d representative_point_world = VNULL;
    ChVector3d representative_normal_world = VNULL;
    std::vector<std::size_t> source_sample_indices;
};

/// One collapsed sheet sample obtained by merging a unique fiber cluster of band seeds.
struct ChApi ChSDFSheetFiberSample {
    std::size_t region_id = 0;
    std::size_t fiber_id = 0;

    /// Legacy dominant-axis bookkeeping retained for debugging and fallback inspection.
    int dominant_axis = -1;
    ChVector2i lattice_coord = ChVector2i(0, 0);

    double measure_area = 0;
    /// In v2 this defaults to the measure-preserving area because no explicit support polygon is reconstructed yet.
    double footprint_area = 0;

    std::size_t support_seed_count = 0;

    ChVector3d centroid_world = VNULL;
    ChVector3d pressure_center_world = VNULL;
    ChVector3d normal_world = VNULL;
    ChVector3d force_world = VNULL;

    ChAABB source_bounds_world;
    ChAABB support_bbox_world;
    ChSDFSheetLocalFootprint support_footprint;

    std::vector<ChSDFSheetSupportEvidence> support_evidence;
    std::vector<std::size_t> source_sample_indices;

    bool HasSupport() const { return measure_area > 0; }
};

/// One connected component on the collapsed sheet graph.
struct ChApi ChSDFSheetPatch {
    std::size_t patch_id = 0;
    std::size_t support_columns = 0;
    std::size_t support_seed_count = 0;

    double measure_area = 0;
    double footprint_area = 0;

    ChVector3d centroid_world = VNULL;
    ChVector3d pressure_center_world = VNULL;
    ChVector3d mean_normal_world = VNULL;

    ChAABB bounds_world;
    ChAABB support_bbox_world;
    ChSDFSheetLocalFootprint support_footprint;

    std::size_t layered_cell_count = 0;
    std::size_t max_layer_count_per_cell = 0;
    std::size_t largest_connected_sheet_cells = 0;
    double mean_layer_count_per_cell = 0;
    double sheet_fill_ratio = 0;
    std::vector<ChSDFPatchPlaneSupportCell> support_cells;
    std::vector<std::size_t> sample_indices;

    bool HasSamples() const { return !sample_indices.empty(); }
};

/// Sheet-layer counterpart of one band region.
struct ChApi ChSDFSheetRegion {
    std::size_t region_id = 0;
    std::size_t persistent_id = 0;
    std::size_t patch_count = 0;
    std::size_t fiber_count = 0;

    int dominant_axis = -1;

    double measure_area = 0;
    double footprint_area = 0;
    double largest_patch_area = 0;
    double area_ratio = 0;
    double mean_support_seed_count = 0;
    double normal_spread = 0;

    bool used_fallback = false;

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
    double area_ratio = 0;
    double mean_support_seed_count = 0;
    double normal_spread = 0;

    bool used_fallback = false;

    ChVector3d sheet_center_world = VNULL;
    ChVector3d pressure_center_world = VNULL;
    ChAABB support_bbox_world;

    std::size_t occupied_cells = 0;
    std::size_t collapsed_samples = 0;
    std::size_t patch_count = 0;
    std::size_t fiber_count = 0;
    std::size_t fallback_regions = 0;

    std::vector<ChSDFSheetRegion> regions;

    bool HasSamples() const { return collapsed_samples > 0; }
};

/// Build a zero-thickness sheet representation from the current band-mechanics result.
class ChApi ChSDFSheetBuilder {
  public:
    /// Legacy dominant world axis used by the dominant-axis fallback implementation.
    static int ComputeDominantAxis(const ChVector3d& normal_world);

    /// Legacy projection of an integer cell coordinate to the two in-plane coordinates orthogonal to the dominant axis.
    static ChVector2i ComputeLatticeCoord(const ChVector3i& coord, int dominant_axis);

    /// Collapse one band region into a sheet representation.
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
