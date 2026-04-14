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

#include "chrono/collision/sdf/ChSDFPotentialCalibration.h"

#include <algorithm>

namespace chrono {

double ChSDFPotentialCalibration::EstimateSymmetricFieldModulusForLinearPenalty(double reference_stiffness,
                                                                                double per_shape_depth_fraction,
                                                                                double calibration_gain) {
    const double safe_stiffness = std::max(reference_stiffness, 0.0);
    const double safe_fraction = per_shape_depth_fraction > 1.0e-12 ? per_shape_depth_fraction : 0.5;
    const double safe_gain = calibration_gain > 0 ? calibration_gain : 1.0;
    return safe_stiffness / safe_fraction * safe_gain;
}

ChSDFPotentialFieldSettings ChSDFPotentialCalibration::MakePlanarConsistentPotentialSettings(
    double reference_stiffness,
    const ChSDFPotentialFieldSettings& base_settings,
    double per_shape_depth_fraction,
    double calibration_gain) {
    auto settings = base_settings;
    settings.modulus = EstimateSymmetricFieldModulusForLinearPenalty(reference_stiffness, per_shape_depth_fraction,
                                                                     calibration_gain);
    return settings;
}

}  // end namespace chrono
