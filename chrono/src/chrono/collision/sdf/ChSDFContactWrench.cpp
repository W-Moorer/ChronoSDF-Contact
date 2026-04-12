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

#include "chrono/collision/sdf/ChSDFContactWrench.h"

#include <algorithm>

namespace chrono {
namespace {

double GridStep(std::size_t count, double half_length) {
    return (count > 1) ? (2.0 * half_length / static_cast<double>(count - 1)) : (2.0 * half_length);
}

double TrapezoidWeight(std::size_t index, std::size_t count) {
    if (count <= 1) {
        return 1.0;
    }

    return (index == 0 || index + 1 == count) ? 0.5 : 1.0;
}

double ClampPressure(double pressure, const ChSDFNormalPressureSettings& settings) {
    if (settings.clamp_negative_pressure) {
        pressure = std::max(0.0, pressure);
    }

    if (settings.max_pressure >= 0) {
        pressure = std::min(pressure, settings.max_pressure);
    }

    return pressure;
}

}  // namespace

ChVector3d ChSDFPatchKinematics::PointVelocity(const ChVector3d& point_shape, const ChFrame<>& patch_frame_shape) const {
    return linear_velocity + Vcross(angular_velocity, point_shape - patch_frame_shape.GetPos());
}

ChSDFContactWrenchResult ChSDFContactWrenchEvaluator::EvaluatePatch(
    const ChSDFContactPatch& patch,
    const ChSDFContactPatchSampler::Settings& patch_settings,
    const ChSDFPatchKinematics& relative_kinematics,
    const ChSDFNormalPressureSettings& pressure_settings) {
    ChSDFContactWrenchResult result;
    result.patch_frame = patch.patch_frame;
    result.total_samples = patch.total_samples;
    result.accepted_samples = patch.accepted_samples;
    result.patch_area = 4.0 * patch_settings.half_length_u * patch_settings.half_length_v;

    if (patch.samples.empty()) {
        return result;
    }

    const double du = GridStep(patch_settings.samples_u, patch_settings.half_length_u);
    const double dv = GridStep(patch_settings.samples_v, patch_settings.half_length_v);

    double pressure_weight_sum = 0;
    ChVector3d pressure_center_sum = VNULL;

    result.samples.reserve(patch.samples.size());

    for (const auto& patch_sample : patch.samples) {
        ChSDFContactWrenchSample sample;
        sample.patch_sample = patch_sample;

        const double weight_u = TrapezoidWeight(patch_sample.iu, patch_settings.samples_u);
        const double weight_v = TrapezoidWeight(patch_sample.iv, patch_settings.samples_v);
        sample.quadrature_area = du * dv * weight_u * weight_v;

        sample.point_shape = patch_sample.probe.point_world;
        sample.point_patch = patch.patch_frame.TransformPointParentToLocal(sample.point_shape);

        sample.signed_gap = patch_sample.probe.distance - pressure_settings.distance_offset;
        sample.penetration = std::max(-sample.signed_gap, 0.0);

        const ChVector3d point_velocity =
            relative_kinematics.PointVelocity(sample.point_shape, patch.patch_frame);
        sample.normal_speed = -Vdot(point_velocity, patch_sample.probe.normal_world);

        const double damping_speed =
            pressure_settings.use_only_closing_speed ? std::max(sample.normal_speed, 0.0) : sample.normal_speed;

        sample.pressure = pressure_settings.stiffness * sample.penetration +
                          pressure_settings.damping * damping_speed + pressure_settings.adhesion_pressure;
        sample.pressure = ClampPressure(sample.pressure, pressure_settings);

        if (sample.pressure > 0) {
            sample.active = true;
            sample.force_shape = patch_sample.probe.normal_world * (sample.pressure * sample.quadrature_area);

            ChWrenchd point_wrench_shape = {sample.force_shape, Vcross(sample.point_shape, sample.force_shape)};
            ChWrenchd point_wrench_patch = patch.patch_frame.TransformWrenchParentToLocal(point_wrench_shape);

            sample.torque_shape = point_wrench_shape.torque;
            sample.force_patch = point_wrench_patch.force;
            sample.torque_patch = point_wrench_patch.torque;

            result.wrench_shape.force += point_wrench_shape.force;
            result.wrench_shape.torque += point_wrench_shape.torque;
            result.wrench_patch.force += point_wrench_patch.force;
            result.wrench_patch.torque += point_wrench_patch.torque;

            result.active_samples++;
            result.active_area += sample.quadrature_area;
            result.integrated_pressure += sample.pressure * sample.quadrature_area;
            result.max_penetration = std::max(result.max_penetration, sample.penetration);
            result.max_pressure = std::max(result.max_pressure, sample.pressure);

            const double sample_weight = sample.pressure * sample.quadrature_area;
            pressure_weight_sum += sample_weight;
            pressure_center_sum += sample.point_shape * sample_weight;
        }

        result.accepted_area += sample.quadrature_area;
        result.samples.push_back(sample);
    }

    if (result.active_area > 0) {
        result.mean_pressure = result.integrated_pressure / result.active_area;
    }

    if (pressure_weight_sum > 0) {
        result.pressure_center_shape = pressure_center_sum / pressure_weight_sum;
    }

    return result;
}

ChSDFContactWrenchResult ChSDFContactWrenchEvaluator::EvaluatePatchLocal(
    const ChCollisionShapeSDF& shape,
    const ChFrame<>& patch_frame_local,
    const ChSDFContactPatchSampler::Settings& patch_settings,
    const ChSDFPatchKinematics& relative_kinematics,
    const ChSDFNormalPressureSettings& pressure_settings) {
    const auto patch = shape.SamplePatchLocal(patch_frame_local, patch_settings);
    return EvaluatePatch(patch, patch_settings, relative_kinematics, pressure_settings);
}

}  // end namespace chrono
