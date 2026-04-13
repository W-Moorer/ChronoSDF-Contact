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

#ifndef CH_SDF_PFC_EVALUATOR_H
#define CH_SDF_PFC_EVALUATOR_H

#include <vector>

#include "chrono/collision/sdf/ChSDFContactSurface.h"
#include "chrono/collision/sdf/ChSDFContactWrench.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// Pressure-field contact evaluator operating on equal-pressure surface quadrature points.
/// The traction law uses the shape-fixed potential pressure field value p_cap at each quadrature point.
class ChApi ChSDFPFCEvaluator {
  public:
    static ChSDFBrickPairWrenchResult EvaluateRegion(const ChSDFContactSurfaceRegion& surface_region,
                                                     const ChFrameMoving<>& shape_a_frame_abs,
                                                     const ChFrameMoving<>& shape_b_frame_abs,
                                                     const ChSDFEffectiveMassProperties& body_a,
                                                     const ChSDFEffectiveMassProperties& body_b,
                                                     const ChSDFNormalPressureSettings& pressure_settings);

    static ChSDFShapePairContactResult EvaluateRegions(const std::vector<ChSDFContactSurfaceRegion>& surface_regions,
                                                       const ChFrameMoving<>& shape_a_frame_abs,
                                                       const ChFrameMoving<>& shape_b_frame_abs,
                                                       const ChSDFEffectiveMassProperties& body_a,
                                                       const ChSDFEffectiveMassProperties& body_b,
                                                       const ChSDFNormalPressureSettings& pressure_settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
