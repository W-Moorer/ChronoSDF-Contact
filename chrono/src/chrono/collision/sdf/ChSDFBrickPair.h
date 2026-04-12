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

#ifndef CH_SDF_BRICK_PAIR_H
#define CH_SDF_BRICK_PAIR_H

#include <vector>

#include "chrono/collision/ChCollisionShapeSDF.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// One retained sparse narrow-band brick expressed in the absolute frame.
struct ChApi ChSDFBrickWorldProxy {
    std::size_t brick_index = 0;

    ChSDFLeafBrick brick;
    ChAABB world_bounds;
    ChVector3d center_world = VNULL;
};

/// One retained coarse brick pair between two SDF shapes.
struct ChApi ChSDFBrickPairCandidate {
    std::size_t brick_a_index = 0;
    std::size_t brick_b_index = 0;

    double center_distance = 0;
    double separation_distance = 0;
    double overlap_volume = 0;
    double min_abs_value_bound = 0;

    ChAABB overlap_world;

    ChSDFBrickWorldProxy brick_a;
    ChSDFBrickWorldProxy brick_b;
};

/// SDF-only coarse broadphase between two sparse narrow-band leaf-brick sets.
/// This layer is the dual-shape counterpart of the single patch-slab filter.
class ChApi ChSDFBrickPairBroadphase {
  public:
    struct Settings {
        double world_margin = 0;
        double max_center_distance = -1;
        double max_separation_distance = 0;
        double max_min_abs_value = -1;

        bool require_world_overlap = true;
    };

    /// Enumerate the sparse bricks of one SDF shape and express them in the absolute frame of the body.
    static std::vector<ChSDFBrickWorldProxy> BuildWorldBricks(const ChCollisionShapeSDF& shape,
                                                              const ChFrame<>& shape_frame_abs,
                                                              const Settings& settings);

    /// Find coarse candidate brick pairs between two SDF shapes.
    static std::vector<ChSDFBrickPairCandidate> FindBrickPairs(const ChCollisionShapeSDF& shape_a,
                                                               const ChFrame<>& shape_a_frame_abs,
                                                               const ChCollisionShapeSDF& shape_b,
                                                               const ChFrame<>& shape_b_frame_abs,
                                                               const Settings& settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
