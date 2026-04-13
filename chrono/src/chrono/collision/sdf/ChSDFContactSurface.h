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

#ifndef CH_SDF_CONTACT_SURFACE_H
#define CH_SDF_CONTACT_SURFACE_H

#include <array>
#include <cstddef>
#include <vector>

#include "chrono/collision/ChCollisionShapeSDF.h"
#include "chrono/collision/sdf/ChSDFContactRegion.h"
#include "chrono/collision/sdf/ChSDFPotentialField.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

struct ChApi ChSDFLineSupportInterval {
    bool valid = false;

    double s_min = 0;
    double s_max = 0;
    double overlap_thickness = 0;
};

struct ChApi ChSDFRegionChartSettings {
    double cell_size = -1;
    double chart_padding = -1;
    double normal_search_half_length = -1;
    double line_search_step = -1;
    double support_epsilon = 0;
    double equilibrium_tolerance = 1.0e-6;
    std::size_t max_cells_u = 257;
    std::size_t max_cells_v = 257;
};

struct ChApi ChSDFRegionChartNode {
    std::size_t index = 0;

    int iu = 0;
    int iv = 0;

    ChVector3d point_chart = VNULL;
    ChVector3d anchor_world = VNULL;
    ChVector3d point_world_b = VNULL;
    ChVector3d normal_world_b = VNULL;

    ChSDFLineSupportInterval support_interval;

    ChSDFPotentialFieldProbe probe_a;
    ChSDFPotentialFieldProbe probe_b;

    double support_overlap = 0;
    bool supported = false;
};

struct ChApi ChSDFContactSurfacePolygon {
    double area_uv = 0;

    ChVector3d centroid_uv = VNULL;
    ChVector3d centroid_world = VNULL;

    std::vector<ChVector3d> vertices_uv;
    std::vector<ChVector3d> vertices_world;
};

struct ChApi ChSDFRegionChartCell {
    std::size_t index = 0;

    int iu = 0;
    int iv = 0;

    std::array<std::size_t, 4> corner_indices = {{0, 0, 0, 0}};

    bool active = false;

    double clipped_area_uv = 0;

    std::vector<ChSDFContactSurfacePolygon> clipped_polygons;
};

struct ChApi ChSDFContactQuadraturePoint {
    ChVector3d point_world = VNULL;
    ChVector3d point_shape_a = VNULL;
    ChVector3d point_shape_b = VNULL;
    ChVector3d normal_world = VNULL;

    double area_weight = 0;
    double p0_a = 0;
    double p0_b = 0;
    double p_cap = 0;
    double h_value = 0;
    double support_overlap = 0;

    ChSDFPotentialFieldProbe probe_a;
    ChSDFPotentialFieldProbe probe_b;
};

struct ChApi ChSDFContactSurfaceRegion {
    ChSDFBrickPairRegion seed_region;

    ChFrame<> chart_frame_world;
    ChFrame<> chart_frame_shape_a;
    ChFrame<> chart_frame_shape_b;

    ChAABB chart_bounds_uv;

    double cell_size = 0;
    double cell_size_u = 0;
    double cell_size_v = 0;
    double supported_area = 0;

    std::size_t nodes_u = 0;
    std::size_t nodes_v = 0;
    std::size_t cells_u = 0;
    std::size_t cells_v = 0;
    std::size_t supported_cells = 0;

    std::vector<ChSDFRegionChartNode> nodes;
    std::vector<ChSDFRegionChartCell> cells;
    std::vector<ChSDFContactQuadraturePoint> quadrature_points;

    bool HasSupport() const { return supported_cells > 0 && supported_area > 0; }
};

class ChApi ChSDFContactSurfaceBuilder {
  public:
    static ChFrame<> BuildBSideChartFrame(const ChSDFBrickPairRegion& region, const ChFrame<>& shape_b_frame_abs);

    static ChSDFContactSurfaceRegion BuildRegionSurface(const ChCollisionShapeSDF& shape_a,
                                                        const ChFrame<>& shape_a_frame_abs,
                                                        const ChCollisionShapeSDF& shape_b,
                                                        const ChFrame<>& shape_b_frame_abs,
                                                        const ChSDFBrickPairRegion& seed_region,
                                                        const ChSDFRegionChartSettings& chart_settings);

    static std::vector<ChSDFContactSurfaceRegion> BuildRegionSurfaces(
        const ChCollisionShapeSDF& shape_a,
        const ChFrame<>& shape_a_frame_abs,
        const ChCollisionShapeSDF& shape_b,
        const ChFrame<>& shape_b_frame_abs,
        const std::vector<ChSDFBrickPairRegion>& seed_regions,
        const ChSDFRegionChartSettings& chart_settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
