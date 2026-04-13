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

#ifndef CH_SDF_VOLUME_CONTACT_EVALUATOR_H
#define CH_SDF_VOLUME_CONTACT_EVALUATOR_H

#include <vector>

#include "chrono/collision/sdf/ChSDFBrickPair.h"
#include "chrono/collision/sdf/ChSDFContactRegion.h"
#include "chrono/collision/sdf/ChSDFContactWrench.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// Volume-integral contact evaluator operating on a common voxel grid induced by overlapping brick pairs.
class ChApi ChSDFVolumeContactEvaluator {
  public:
    static ChSDFShapePairContactResult EvaluateBrickPairs(const ChCollisionShapeSDF& shape_a,
                                                          const ChFrameMoving<>& shape_a_frame_abs,
                                                          const ChCollisionShapeSDF& shape_b,
                                                          const ChFrameMoving<>& shape_b_frame_abs,
                                                          const std::vector<ChSDFBrickPairCandidate>& brick_pairs,
                                                          const ChSDFContactRegionBuilder::Settings& region_settings,
                                                          const ChSDFEffectiveMassProperties& body_a,
                                                          const ChSDFEffectiveMassProperties& body_b,
                                                          const ChSDFNormalPressureSettings& pressure_settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
