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
#include "chrono/collision/sdf/ChSDFContactRegion.h"
#include "chrono/core/ChFrame.h"
#include "chrono/core/ChFrameMoving.h"
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

/// Normal pressure law with optional regularized tangential traction.
/// Coefficients are interpreted as pressure-density coefficients:
///   pressure = stiffness * penetration + damping * closing_speed + adhesion_pressure
struct ChApi ChSDFNormalPressureSettings {
    double stiffness = 1.0e7;
    double damping = 0;
    double distance_offset = 0;
    double adhesion_pressure = 0;
    double max_pressure = -1;
    /// Exponential attenuation applied to the local gradient-quality indicator q = ||grad phi|| - 1|.
    double gradient_quality_gain = 0;
    /// Rational attenuation applied to a local curvature proxy extracted from neighboring normal variation.
    double curvature_gain = 0;
    /// Rational attenuation applied to the local resolution length h.
    double resolution_scale_gain = 0;
    /// Nominal resolution length used to normalize h. If <= 0, the raw local length is used.
    double reference_resolution_length = -1;
    double min_stiffness_scale = 0;
    double max_stiffness_scale = -1;
    /// Regularized Coulomb friction coefficient used for the tangential traction density.
    double friction_coefficient = 0;
    /// Velocity scale used in the smooth traction law t = -mu * p * vt / sqrt(|vt|^2 + eps^2).
    double tangential_velocity_regularization = 1.0e-3;

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
    double tangential_speed = 0;
    double gradient_quality = 0;
    double curvature_proxy = 0;
    double resolution_length = 0;
    double stiffness_scale = 1;
    double local_stiffness = 0;
    double local_damping = 0;
    double pressure = 0;

    ChVector3d point_shape = VNULL;
    ChVector3d point_patch = VNULL;
    ChVector3d tangential_velocity_shape = VNULL;
    /// Tangential traction density expressed in the SDF shape frame.
    ChVector3d traction_shape = VNULL;
    /// Tangential traction density expressed in the patch frame.
    ChVector3d traction_patch = VNULL;
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
    double mean_local_stiffness = 0;
    double max_local_stiffness = 0;
    double max_penetration = 0;
    double max_pressure = 0;

    /// Pressure-weighted application point of the active patch force, expressed in the SDF shape frame.
    ChVector3d pressure_center_shape = VNULL;

    std::vector<ChSDFContactWrenchSample> samples;

    bool HasActiveContact() const { return active_samples > 0; }
};

/// One sample contribution in the distributed dual-SDF region-to-wrench aggregation.
struct ChApi ChSDFBrickPairWrenchSample {
    ChSDFBrickPairRegionSample region_sample;

    bool active = false;

    double quadrature_area = 0;
    double signed_gap = 0;
    double penetration = 0;
    double normal_speed = 0;
    double tangential_speed = 0;
    double gradient_quality = 0;
    double curvature_proxy = 0;
    double resolution_length = 0;
    double stiffness_scale = 1;
    double local_stiffness = 0;
    double local_damping = 0;
    double pressure = 0;

    ChVector3d point_patch = VNULL;
    ChVector3d tangential_velocity_world = VNULL;
    /// Tangential traction density expressed in the absolute frame.
    ChVector3d traction_world = VNULL;
    ChVector3d force_world = VNULL;
    ChVector3d force_shape_a = VNULL;
    ChVector3d torque_shape_a = VNULL;
    ChVector3d force_shape_b = VNULL;
    ChVector3d torque_shape_b = VNULL;
};

/// Aggregated distributed contact result for one dual-SDF connected region.
struct ChApi ChSDFBrickPairWrenchResult {
    ChSDFBrickPairRegion region;

    /// Net wrench acting on shape A, reduced to the origin of shape A and expressed in the local frame of shape A.
    ChWrenchd wrench_shape_a = {VNULL, VNULL};
    /// Net wrench acting on shape B, reduced to the origin of shape B and expressed in the local frame of shape B.
    ChWrenchd wrench_shape_b = {VNULL, VNULL};

    /// Same physical wrenches expressed in the absolute frame and reduced to the corresponding shape origins.
    ChWrenchd wrench_world_a = {VNULL, VNULL};
    ChWrenchd wrench_world_b = {VNULL, VNULL};

    std::size_t total_samples = 0;
    std::size_t active_samples = 0;

    double patch_area = 0;
    double active_area = 0;
    double integrated_pressure = 0;
    double mean_pressure = 0;
    double mean_local_stiffness = 0;
    double max_local_stiffness = 0;
    double max_penetration = 0;
    double max_pressure = 0;

    ChVector3d pressure_center_world = VNULL;

    std::vector<ChSDFBrickPairWrenchSample> samples;

    bool HasActiveContact() const { return active_samples > 0; }
};

/// Aggregated distributed contact result across all connected regions of one SDF-shape pair.
struct ChApi ChSDFShapePairContactResult {
    bool valid = false;

    ChFrameMoving<> shape_a_frame_abs;
    ChFrameMoving<> shape_b_frame_abs;

    ChWrenchd wrench_shape_a = {VNULL, VNULL};
    ChWrenchd wrench_shape_b = {VNULL, VNULL};
    ChWrenchd wrench_world_a = {VNULL, VNULL};
    ChWrenchd wrench_world_b = {VNULL, VNULL};

    std::size_t total_regions = 0;
    std::size_t active_regions = 0;

    double active_area = 0;
    double integrated_pressure = 0;
    double mean_pressure = 0;
    double mean_local_stiffness = 0;
    double max_local_stiffness = 0;
    double max_penetration = 0;
    double max_pressure = 0;

    std::vector<ChSDFBrickPairWrenchResult> regions;

    bool HasActiveContact() const { return active_regions > 0; }
};

/// Distributed contact evaluator that turns a sampled SDF patch into a net wrench.
class ChApi ChSDFContactWrenchEvaluator {
  public:
    /// Evaluate a pre-sampled patch against a penalty law with optional tangential traction.
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

    /// Evaluate one parameterized dual-SDF connected region against the same pressure and traction law.
    static ChSDFBrickPairWrenchResult EvaluateBrickPairRegion(const ChSDFBrickPairRegion& region,
                                                              const ChFrameMoving<>& shape_a_frame_abs,
                                                              const ChFrameMoving<>& shape_b_frame_abs,
                                                              const ChSDFNormalPressureSettings& pressure_settings);

    /// Evaluate multiple parameterized dual-SDF connected regions and aggregate the resulting wrenches.
    static ChSDFShapePairContactResult EvaluateBrickPairRegions(
        const std::vector<ChSDFBrickPairRegion>& regions,
        const ChFrameMoving<>& shape_a_frame_abs,
        const ChFrameMoving<>& shape_b_frame_abs,
        const ChSDFNormalPressureSettings& pressure_settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
