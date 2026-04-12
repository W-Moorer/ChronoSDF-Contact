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

#ifndef CH_SDF_CONTACT_WRENCH_H
#define CH_SDF_CONTACT_WRENCH_H

#include <cstddef>
#include <vector>

#include "chrono/collision/ChCollisionShapeSDF.h"
#include "chrono/core/ChFrame.h"
#include "chrono/core/ChVector3.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// Relative patch kinematics expressed in the local frame of an SDF shape.
/// The linear velocity is the relative velocity of the patch origin with respect to the SDF shape.
/// The angular velocity describes the relative rigid motion of the probe patch about the same origin.
struct ChApi ChSDFPatchKinematics {
    ChVector3d linear_velocity = VNULL;
    ChVector3d angular_velocity = VNULL;

    /// Evaluate the relative velocity at a local point on the patch.
    ChVector3d PointVelocity(const ChVector3d& point_shape, const ChFrame<>& patch_frame_shape) const;
};

/// Normal-only pressure law used to aggregate a distributed SDF contact patch.
/// Coefficients are interpreted as pressure-density coefficients:
///   pressure = stiffness * penetration + damping * closing_speed + adhesion_pressure
struct ChApi ChSDFNormalPressureSettings {
    double stiffness = 1.0e7;
    double damping = 0;
    double distance_offset = 0;
    double adhesion_pressure = 0;
    double max_pressure = -1;

    bool use_only_closing_speed = true;
    bool clamp_negative_pressure = true;
};

/// One sample contribution in the distributed patch-to-wrench aggregation.
struct ChApi ChSDFContactWrenchSample {
    ChSDFContactPatchSample patch_sample;

    bool active = false;

    double quadrature_area = 0;
    double signed_gap = 0;
    double penetration = 0;
    double normal_speed = 0;
    double pressure = 0;

    ChVector3d point_shape = VNULL;
    ChVector3d point_patch = VNULL;
    ChVector3d force_shape = VNULL;
    ChVector3d torque_shape = VNULL;
    ChVector3d force_patch = VNULL;
    ChVector3d torque_patch = VNULL;
};

/// Aggregated distributed contact result for one sampled patch.
/// All vectors are expressed in the local frame of the SDF shape unless noted otherwise.
struct ChApi ChSDFContactWrenchResult {
    ChFrame<> patch_frame;

    /// Net wrench acting on the probe patch, expressed in the SDF shape frame and reduced to the SDF shape origin.
    ChWrenchd wrench_shape = {VNULL, VNULL};
    /// Same physical wrench acting on the probe patch, but reduced to the patch origin and expressed in the patch frame.
    ChWrenchd wrench_patch = {VNULL, VNULL};

    std::size_t total_samples = 0;
    std::size_t accepted_samples = 0;
    std::size_t active_samples = 0;

    double patch_area = 0;
    double accepted_area = 0;
    double active_area = 0;
    double integrated_pressure = 0;
    double mean_pressure = 0;
    double max_penetration = 0;
    double max_pressure = 0;

    /// Pressure-weighted application point of the active patch force, expressed in the SDF shape frame.
    ChVector3d pressure_center_shape = VNULL;

    std::vector<ChSDFContactWrenchSample> samples;

    bool HasActiveContact() const { return active_samples > 0; }
};

/// Distributed contact evaluator that turns a sampled SDF patch into a net wrench.
class ChApi ChSDFContactWrenchEvaluator {
  public:
    /// Evaluate a pre-sampled patch against a normal-only penalty law.
    static ChSDFContactWrenchResult EvaluatePatch(const ChSDFContactPatch& patch,
                                                  const ChSDFContactPatchSampler::Settings& patch_settings,
                                                  const ChSDFPatchKinematics& relative_kinematics,
                                                  const ChSDFNormalPressureSettings& pressure_settings);

    /// Sample and evaluate a patch directly from an SDF collision shape.
    static ChSDFContactWrenchResult EvaluatePatchLocal(const ChCollisionShapeSDF& shape,
                                                       const ChFrame<>& patch_frame_local,
                                                       const ChSDFContactPatchSampler::Settings& patch_settings,
                                                       const ChSDFPatchKinematics& relative_kinematics,
                                                       const ChSDFNormalPressureSettings& pressure_settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
