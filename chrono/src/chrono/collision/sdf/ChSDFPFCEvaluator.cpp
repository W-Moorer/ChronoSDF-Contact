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

#include "chrono/collision/sdf/ChSDFPFCEvaluator.h"

#include <algorithm>
#include <cmath>

namespace chrono {
namespace {

double ClampPressure(double pressure, const ChSDFNormalPressureSettings& settings) {
    if (settings.clamp_negative_pressure) {
        pressure = std::max(0.0, pressure);
    }

    if (settings.max_pressure >= 0) {
        pressure = std::min(pressure, settings.max_pressure);
    }

    return pressure;
}

double GradientQuality(const ChSDFPotentialFieldProbe& probe) {
    return probe.valid ? std::abs(probe.gradient_local.Length() - 1.0) : 0.0;
}

double ResolutionMetric(double resolution_length, const ChSDFNormalPressureSettings& settings) {
    if (resolution_length <= 0 || settings.resolution_scale_gain <= 0) {
        return 0.0;
    }

    if (settings.reference_resolution_length > 0) {
        return std::max(resolution_length / settings.reference_resolution_length - 1.0, 0.0);
    }

    return resolution_length;
}

double LocalStiffnessScale(double gradient_quality,
                           double curvature_proxy,
                           double resolution_length,
                           const ChSDFNormalPressureSettings& settings) {
    double scale = 1.0;

    if (settings.gradient_quality_gain > 0) {
        scale *= std::exp(-settings.gradient_quality_gain * std::max(gradient_quality, 0.0));
    }

    if (settings.curvature_gain > 0) {
        scale *= 1.0 / (1.0 + settings.curvature_gain * std::max(curvature_proxy, 0.0));
    }

    if (settings.resolution_scale_gain > 0) {
        scale *= 1.0 / (1.0 + settings.resolution_scale_gain * ResolutionMetric(resolution_length, settings));
    }

    if (settings.min_stiffness_scale > 0) {
        scale = std::max(scale, settings.min_stiffness_scale);
    }

    if (settings.max_stiffness_scale >= 0) {
        scale = std::min(scale, settings.max_stiffness_scale);
    }

    return scale;
}

double ClampEffectiveMass(double effective_mass, const ChSDFNormalPressureSettings& settings) {
    if (settings.min_effective_mass > 0) {
        effective_mass = std::max(effective_mass, settings.min_effective_mass);
    }

    if (settings.max_effective_mass >= 0) {
        effective_mass = std::min(effective_mass, settings.max_effective_mass);
    }

    return effective_mass;
}

double ResolveLocalDamping(double local_stiffness,
                           double effective_mass,
                           double stiffness_scale,
                           const ChSDFNormalPressureSettings& settings) {
    if (settings.damping_ratio >= 0) {
        if (local_stiffness <= 0 || effective_mass <= 0) {
            return 0.0;
        }

        return 2.0 * settings.damping_ratio * std::sqrt(local_stiffness * effective_mass);
    }

    return settings.damping * std::sqrt(std::max(stiffness_scale, 0.0));
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

    const double regularization = std::max(settings.tangential_velocity_regularization, 0.0);
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

ChVector3d TransformPotentialGradientWorld(const ChFrameMoving<>& shape_frame_abs,
                                           const ChSDFPotentialFieldProbe& probe) {
    return shape_frame_abs.TransformDirectionLocalToParent(probe.grad_p0_local);
}

double ResolveContactNormalStiffness(const ChSDFContactQuadraturePoint& quadrature_point,
                                     const ChFrameMoving<>& shape_a_frame_abs,
                                     const ChFrameMoving<>& shape_b_frame_abs) {
    const ChVector3d normal_world = quadrature_point.normal_world;
    const ChVector3d grad_a_world = TransformPotentialGradientWorld(shape_a_frame_abs, quadrature_point.probe_a);
    const ChVector3d grad_b_world = TransformPotentialGradientWorld(shape_b_frame_abs, quadrature_point.probe_b);

    const double stiffness_a = std::abs(Vdot(grad_a_world, normal_world));
    const double stiffness_b = std::abs(Vdot(grad_b_world, normal_world));
    const double stiffness = 0.5 * (stiffness_a + stiffness_b);
    if (stiffness > 1.0e-12) {
        return stiffness;
    }

    if (quadrature_point.support_overlap > 1.0e-12) {
        return std::max(quadrature_point.p_cap / quadrature_point.support_overlap, 0.0);
    }

    return 0.0;
}

}  // namespace

ChSDFBrickPairWrenchResult ChSDFPFCEvaluator::EvaluateRegion(const ChSDFContactSurfaceRegion& surface_region,
                                                             const ChFrameMoving<>& shape_a_frame_abs,
                                                             const ChFrameMoving<>& shape_b_frame_abs,
                                                             const ChSDFEffectiveMassProperties& body_a,
                                                             const ChSDFEffectiveMassProperties& body_b,
                                                             const ChSDFNormalPressureSettings& pressure_settings) {
    ChSDFBrickPairWrenchResult result;
    result.region = surface_region.seed_region;
    result.total_samples = surface_region.quadrature_points.size();

    if (surface_region.quadrature_points.empty()) {
        return result;
    }

    const double resolution_length = surface_region.cell_size > 0 ? surface_region.cell_size : 1.0e-3;

    double pressure_weight_sum = 0;
    ChVector3d pressure_center_sum = VNULL;
    double local_stiffness_area_sum = 0;
    double effective_mass_area_sum = 0;

    result.samples.reserve(surface_region.quadrature_points.size());

    for (const auto& quadrature_point : surface_region.quadrature_points) {
        ChSDFBrickPairWrenchSample sample;
        sample.quadrature_point = quadrature_point;
        sample.point_patch = surface_region.chart_frame_world.TransformPointParentToLocal(quadrature_point.point_world);
        sample.quadrature_area = std::max(quadrature_point.area_weight, 0.0);
        result.patch_area += sample.quadrature_area;

        const ChVector3d speed_world_a = shape_a_frame_abs.PointSpeedLocalToParent(quadrature_point.point_shape_a);
        const ChVector3d speed_world_b = shape_b_frame_abs.PointSpeedLocalToParent(quadrature_point.point_shape_b);
        const ChVector3d relative_speed_world = speed_world_b - speed_world_a;
        const double normal_velocity_component = Vdot(relative_speed_world, quadrature_point.normal_world);
        sample.normal_speed = -normal_velocity_component;
        sample.tangential_velocity_world =
            relative_speed_world - quadrature_point.normal_world * normal_velocity_component;
        sample.tangential_speed = sample.tangential_velocity_world.Length();

        sample.gradient_quality =
            0.5 * (GradientQuality(quadrature_point.probe_a) + GradientQuality(quadrature_point.probe_b));
        sample.curvature_proxy = 0.0;
        sample.resolution_length = resolution_length;
        sample.stiffness_scale = LocalStiffnessScale(sample.gradient_quality, sample.curvature_proxy,
                                                     sample.resolution_length, pressure_settings);

        sample.effective_mass = ClampEffectiveMass(
            ChSDFContactWrenchEvaluator::EstimateEffectiveMass(body_a, body_b, quadrature_point.point_world,
                                                               quadrature_point.normal_world,
                                                               pressure_settings.fallback_effective_mass),
            pressure_settings);

        const double base_stiffness =
            ResolveContactNormalStiffness(quadrature_point, shape_a_frame_abs, shape_b_frame_abs);
        sample.local_stiffness = base_stiffness * sample.stiffness_scale;
        sample.local_damping =
            ResolveLocalDamping(sample.local_stiffness, sample.effective_mass, sample.stiffness_scale, pressure_settings);

        const double damping_speed =
            pressure_settings.use_only_closing_speed ? std::max(sample.normal_speed, 0.0) : sample.normal_speed;
        const double base_pressure = quadrature_point.p_cap * sample.stiffness_scale;
        sample.pressure = ClampPressure(base_pressure + sample.local_damping * damping_speed +
                                            pressure_settings.adhesion_pressure,
                                        pressure_settings);

        sample.penetration =
            sample.local_stiffness > 1.0e-12 ? std::max(base_pressure / sample.local_stiffness, 0.0) : 0.0;
        sample.signed_gap = -sample.penetration;

        if (sample.pressure > 0) {
            sample.active = true;
            sample.traction_world =
                ComputeTangentialTraction(sample.tangential_velocity_world, sample.pressure, pressure_settings);

            const ChVector3d normal_force_world_b =
                quadrature_point.normal_world * (sample.pressure * sample.quadrature_area);
            const ChVector3d tangential_force_world_b = sample.traction_world * sample.quadrature_area;
            const ChVector3d force_world_b = normal_force_world_b + tangential_force_world_b;
            const ChVector3d force_world_a = -force_world_b;

            const ChVector3d torque_world_a =
                Vcross(quadrature_point.point_world - shape_a_frame_abs.GetPos(), force_world_a);
            const ChVector3d torque_world_b =
                Vcross(quadrature_point.point_world - shape_b_frame_abs.GetPos(), force_world_b);

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
            local_stiffness_area_sum += sample.local_stiffness * sample.quadrature_area;
            effective_mass_area_sum += sample.effective_mass * sample.quadrature_area;
            result.max_local_stiffness = std::max(result.max_local_stiffness, sample.local_stiffness);
            result.max_effective_mass = std::max(result.max_effective_mass, sample.effective_mass);
            result.max_penetration = std::max(result.max_penetration, sample.penetration);
            result.max_pressure = std::max(result.max_pressure, sample.pressure);

            const double sample_weight = sample.pressure * sample.quadrature_area;
            pressure_weight_sum += sample_weight;
            pressure_center_sum += quadrature_point.point_world * sample_weight;
        }

        result.samples.push_back(std::move(sample));
    }

    if (result.active_area > 0) {
        result.mean_pressure = result.integrated_pressure / result.active_area;
        result.mean_local_stiffness = local_stiffness_area_sum / result.active_area;
        result.mean_effective_mass = effective_mass_area_sum / result.active_area;
    }

    if (pressure_weight_sum > 0) {
        result.pressure_center_world = pressure_center_sum / pressure_weight_sum;
    }

    return result;
}

ChSDFShapePairContactResult ChSDFPFCEvaluator::EvaluateRegions(
    const std::vector<ChSDFContactSurfaceRegion>& surface_regions,
    const ChFrameMoving<>& shape_a_frame_abs,
    const ChFrameMoving<>& shape_b_frame_abs,
    const ChSDFEffectiveMassProperties& body_a,
    const ChSDFEffectiveMassProperties& body_b,
    const ChSDFNormalPressureSettings& pressure_settings) {
    ChSDFShapePairContactResult result;
    result.valid = true;
    result.shape_a_frame_abs = shape_a_frame_abs;
    result.shape_b_frame_abs = shape_b_frame_abs;
    result.total_regions = surface_regions.size();
    result.regions.reserve(surface_regions.size());

    double local_stiffness_area_sum = 0;
    double effective_mass_area_sum = 0;

    for (const auto& surface_region : surface_regions) {
        auto region_result =
            EvaluateRegion(surface_region, shape_a_frame_abs, shape_b_frame_abs, body_a, body_b, pressure_settings);

        if (region_result.HasActiveContact()) {
            result.active_regions++;
            result.active_area += region_result.active_area;
            result.integrated_pressure += region_result.integrated_pressure;
            local_stiffness_area_sum += region_result.mean_local_stiffness * region_result.active_area;
            effective_mass_area_sum += region_result.mean_effective_mass * region_result.active_area;
            result.max_local_stiffness = std::max(result.max_local_stiffness, region_result.max_local_stiffness);
            result.max_effective_mass = std::max(result.max_effective_mass, region_result.max_effective_mass);
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
        result.mean_local_stiffness = local_stiffness_area_sum / result.active_area;
        result.mean_effective_mass = effective_mass_area_sum / result.active_area;
    }

    return result;
}

}  // end namespace chrono
