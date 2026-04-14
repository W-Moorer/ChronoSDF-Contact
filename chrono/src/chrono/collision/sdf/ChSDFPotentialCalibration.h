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

#ifndef CH_SDF_POTENTIAL_CALIBRATION_H
#define CH_SDF_POTENTIAL_CALIBRATION_H

#include "chrono/collision/sdf/ChSDFPotentialField.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// Helper functions that expose the planar-consistency calibration used by the SDF potential field.
class ChApi ChSDFPotentialCalibration {
  public:
    /// Estimate the per-shape virtual pressure-field modulus that recovers a target linear penalty stiffness
    /// under the symmetric planar-compression assumption.
    static double EstimateSymmetricFieldModulusForLinearPenalty(double reference_stiffness,
                                                                double per_shape_depth_fraction = 0.5,
                                                                double calibration_gain = 1.11);

    /// Return a copy of the input settings with the calibrated modulus applied.
    static ChSDFPotentialFieldSettings MakePlanarConsistentPotentialSettings(
        double reference_stiffness,
        const ChSDFPotentialFieldSettings& base_settings = ChSDFPotentialFieldSettings(),
        double per_shape_depth_fraction = 0.5,
        double calibration_gain = 1.11);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
