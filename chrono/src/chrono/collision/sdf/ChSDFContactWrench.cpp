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

ChVector3d ComputeTangentialTraction(const ChVector3d& tangential_velocity,
                                     double pressure,
                                     const ChSDFNormalPressureSettings& settings) {
    if (pressure <= 0 || settings.friction_coefficient <= 0) {
        return VNULL;
    }

    const double speed2 = tangential_velocity.Length2();
    if (speed2 <= 1.0e-24) {
        return VNULL;
    }

    const double regularization =
        std::max(settings.tangential_velocity_regularization, 0.0);

    if (regularization > 0) {
        const double denom = std::sqrt(speed2 + regularization * regularization);
        return tangential_velocity * (-settings.friction_coefficient * pressure / denom);
    }

    const double speed = std::sqrt(speed2);
    if (speed <= 1.0e-12) {
        return VNULL;
    }

    return tangential_velocity * (-settings.friction_coefficient * pressure / speed);
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
        const double normal_velocity_component = Vdot(point_velocity, patch_sample.probe.normal_world);
        sample.normal_speed = -normal_velocity_component;
        sample.tangential_velocity_shape =
            point_velocity - patch_sample.probe.normal_world * normal_velocity_component;
        sample.tangential_speed = sample.tangential_velocity_shape.Length();

        const double damping_speed =
            pressure_settings.use_only_closing_speed ? std::max(sample.normal_speed, 0.0) : sample.normal_speed;

        sample.pressure = pressure_settings.stiffness * sample.penetration +
                          pressure_settings.damping * damping_speed + pressure_settings.adhesion_pressure;
        sample.pressure = ClampPressure(sample.pressure, pressure_settings);

        if (sample.pressure > 0) {
            sample.active = true;
            sample.traction_shape =
                ComputeTangentialTraction(sample.tangential_velocity_shape, sample.pressure, pressure_settings);

            const ChVector3d normal_force_shape =
                patch_sample.probe.normal_world * (sample.pressure * sample.quadrature_area);
            const ChVector3d tangential_force_shape = sample.traction_shape * sample.quadrature_area;
            sample.force_shape = normal_force_shape + tangential_force_shape;

            ChWrenchd point_wrench_shape = {sample.force_shape, Vcross(sample.point_shape, sample.force_shape)};
            ChWrenchd point_wrench_patch = patch.patch_frame.TransformWrenchParentToLocal(point_wrench_shape);

            sample.torque_shape = point_wrench_shape.torque;
            sample.traction_patch = patch.patch_frame.TransformDirectionParentToLocal(sample.traction_shape);
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

ChSDFBrickPairWrenchResult ChSDFContactWrenchEvaluator::EvaluateBrickPairRegion(
    const ChSDFBrickPairRegion& region,
    const ChFrameMoving<>& shape_a_frame_abs,
    const ChFrameMoving<>& shape_b_frame_abs,
    const ChSDFNormalPressureSettings& pressure_settings) {
    ChSDFBrickPairWrenchResult result;
    result.region = region;
    result.total_samples = region.samples.size();
    result.patch_area = 4.0 * region.suggested_patch_settings.half_length_u * region.suggested_patch_settings.half_length_v;

    if (region.samples.empty()) {
        return result;
    }

    const double du = region.sample_spacing > 0 ? region.sample_spacing : GridStep(region.suggested_patch_settings.samples_u,
                                                                                   region.suggested_patch_settings.half_length_u);
    const double dv = region.sample_spacing > 0 ? region.sample_spacing : GridStep(region.suggested_patch_settings.samples_v,
                                                                                   region.suggested_patch_settings.half_length_v);
    const double quadrature_area = du * dv;

    double pressure_weight_sum = 0;
    ChVector3d pressure_center_sum = VNULL;

    result.samples.reserve(region.samples.size());

    for (const auto& region_sample : region.samples) {
        ChSDFBrickPairWrenchSample sample;
        sample.region_sample = region_sample;
        sample.point_patch = region_sample.point_patch;
        sample.quadrature_area = quadrature_area;

        sample.signed_gap = region_sample.combined_gap - pressure_settings.distance_offset;
        sample.penetration = std::max(-sample.signed_gap, 0.0);

        const ChVector3d speed_world_a = shape_a_frame_abs.PointSpeedLocalToParent(region_sample.surface_shape_a);
        const ChVector3d speed_world_b = shape_b_frame_abs.PointSpeedLocalToParent(region_sample.surface_shape_b);
        const ChVector3d relative_speed_world = speed_world_b - speed_world_a;
        const double normal_velocity_component = Vdot(relative_speed_world, region_sample.contact_normal_world);
        sample.normal_speed = -normal_velocity_component;
        sample.tangential_velocity_world =
            relative_speed_world - region_sample.contact_normal_world * normal_velocity_component;
        sample.tangential_speed = sample.tangential_velocity_world.Length();

        const double damping_speed =
            pressure_settings.use_only_closing_speed ? std::max(sample.normal_speed, 0.0) : sample.normal_speed;
        sample.pressure = pressure_settings.stiffness * sample.penetration +
                          pressure_settings.damping * damping_speed + pressure_settings.adhesion_pressure;
        sample.pressure = ClampPressure(sample.pressure, pressure_settings);

        if (sample.pressure > 0) {
            sample.active = true;
            sample.traction_world =
                ComputeTangentialTraction(sample.tangential_velocity_world, sample.pressure, pressure_settings);

            const ChVector3d normal_force_world_b =
                region_sample.contact_normal_world * (sample.pressure * sample.quadrature_area);
            const ChVector3d tangential_force_world_b = sample.traction_world * sample.quadrature_area;
            const ChVector3d force_world_b = normal_force_world_b + tangential_force_world_b;
            const ChVector3d force_world_a = -force_world_b;

            const ChVector3d torque_world_a =
                Vcross(region_sample.point_world - shape_a_frame_abs.GetPos(), force_world_a);
            const ChVector3d torque_world_b =
                Vcross(region_sample.point_world - shape_b_frame_abs.GetPos(), force_world_b);

            sample.force_world = force_world_b;
            sample.force_shape_a = shape_a_frame_abs.TransformDirectionParentToLocal(force_world_a);
            sample.torque_shape_a = shape_a_frame_abs.TransformDirectionParentToLocal(torque_world_a);
            sample.force_shape_b = shape_b_frame_abs.TransformDirectionParentToLocal(force_world_b);
            sample.torque_shape_b = shape_b_frame_abs.TransformDirectionParentToLocal(torque_world_b);

            result.wrench_world_a.force += force_world_a;
            result.wrench_world_a.torque += torque_world_a;
            result.wrench_world_b.force += force_world_b;
            result.wrench_world_b.torque += torque_world_b;
            result.wrench_shape_a.force += sample.force_shape_a;
            result.wrench_shape_a.torque += sample.torque_shape_a;
            result.wrench_shape_b.force += sample.force_shape_b;
            result.wrench_shape_b.torque += sample.torque_shape_b;

            result.active_samples++;
            result.active_area += sample.quadrature_area;
            result.integrated_pressure += sample.pressure * sample.quadrature_area;
            result.max_penetration = std::max(result.max_penetration, sample.penetration);
            result.max_pressure = std::max(result.max_pressure, sample.pressure);

            const double sample_weight = sample.pressure * sample.quadrature_area;
            pressure_weight_sum += sample_weight;
            pressure_center_sum += region_sample.point_world * sample_weight;
        }

        result.samples.push_back(sample);
    }

    if (result.active_area > 0) {
        result.mean_pressure = result.integrated_pressure / result.active_area;
    }

    if (pressure_weight_sum > 0) {
        result.pressure_center_world = pressure_center_sum / pressure_weight_sum;
    }

    return result;
}

ChSDFShapePairContactResult ChSDFContactWrenchEvaluator::EvaluateBrickPairRegions(
    const std::vector<ChSDFBrickPairRegion>& regions,
    const ChFrameMoving<>& shape_a_frame_abs,
    const ChFrameMoving<>& shape_b_frame_abs,
    const ChSDFNormalPressureSettings& pressure_settings) {
    ChSDFShapePairContactResult result;
    result.valid = true;
    result.shape_a_frame_abs = shape_a_frame_abs;
    result.shape_b_frame_abs = shape_b_frame_abs;
    result.total_regions = regions.size();
    result.regions.reserve(regions.size());

    for (const auto& region : regions) {
        auto region_result = EvaluateBrickPairRegion(region, shape_a_frame_abs, shape_b_frame_abs, pressure_settings);

        if (region_result.HasActiveContact()) {
            result.active_regions++;
            result.active_area += region_result.active_area;
            result.integrated_pressure += region_result.integrated_pressure;
            result.max_penetration = std::max(result.max_penetration, region_result.max_penetration);
            result.max_pressure = std::max(result.max_pressure, region_result.max_pressure);
        }

        result.wrench_world_a.force += region_result.wrench_world_a.force;
        result.wrench_world_a.torque += region_result.wrench_world_a.torque;
        result.wrench_world_b.force += region_result.wrench_world_b.force;
        result.wrench_world_b.torque += region_result.wrench_world_b.torque;
        result.wrench_shape_a.force += region_result.wrench_shape_a.force;
        result.wrench_shape_a.torque += region_result.wrench_shape_a.torque;
        result.wrench_shape_b.force += region_result.wrench_shape_b.force;
        result.wrench_shape_b.torque += region_result.wrench_shape_b.torque;

        result.regions.push_back(std::move(region_result));
    }

    if (result.active_area > 0) {
        result.mean_pressure = result.integrated_pressure / result.active_area;
    }

    return result;
}

}  // end namespace chrono
