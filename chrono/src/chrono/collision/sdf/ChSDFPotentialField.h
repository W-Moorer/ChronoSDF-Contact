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

#ifndef CH_SDF_POTENTIAL_FIELD_H
#define CH_SDF_POTENTIAL_FIELD_H

#include "chrono/collision/sdf/ChNanoVDBLevelSet.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// Shape-fixed virtual pressure-field settings derived from the local SDF depth.
struct ChApi ChSDFPotentialFieldSettings {
    /// Pressure-scale coefficient used in the linear depth-to-pressure map p0 = E * depth / depth_scale.
    double modulus = 1.0e7;
    /// Characteristic depth used to normalize the inward distance. If <= 0, use the minimum voxel size.
    double depth_scale = -1;
    /// Optional cap applied to the inward depth before mapping it to p0. If <= 0, no cap is applied.
    double depth_cap = -1;
    /// Offset used to define the support indicator support_value = support_margin - phi.
    double support_margin = 0;
    /// Clamp the potential pressure field to zero outside the zero level set.
    bool clamp_outside_to_zero = true;
};

/// Result of probing the shape-fixed virtual pressure field.
struct ChApi ChSDFPotentialFieldProbe {
    bool valid = false;

    double phi = 0;
    double depth = 0;
    double p0 = 0;
    double support_value = 0;

    ChVector3d point_local = VNULL;
    ChVector3d gradient_local = VNULL;
    ChVector3d normal_local = VNULL;
    ChVector3d grad_p0_local = VNULL;

    ChSDFProbeResult sdf_probe;
};

/// Stateless helper that maps a raw SDF probe to the corresponding virtual pressure-field probe.
class ChApi ChSDFPotentialFieldEvaluator {
  public:
    /// Resolve the depth scale used in the virtual pressure map.
    static double ResolveDepthScale(const ChNanoVDBGridInfo& info, const ChSDFPotentialFieldSettings& settings);

    /// Convert a local SDF probe into a local virtual pressure-field probe.
    static ChSDFPotentialFieldProbe Evaluate(const ChSDFProbeResult& sdf_probe,
                                             const ChNanoVDBGridInfo& info,
                                             const ChSDFPotentialFieldSettings& settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
