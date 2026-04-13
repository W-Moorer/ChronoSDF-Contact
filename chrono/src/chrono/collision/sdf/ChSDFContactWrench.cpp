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
#include <cmath>
#include <unordered_map>

#include "chrono/physics/ChBody.h"

namespace chrono {
namespace {

struct CoordKey3 {
    int x = 0;
    int y = 0;
    int z = 0;

    bool operator==(const CoordKey3& other) const { return x == other.x && y == other.y && z == other.z; }
};

struct CoordKey3Hasher {
    std::size_t operator()(const CoordKey3& key) const {
        std::size_t seed = 0;
        seed ^= static_cast<std::size_t>(key.x) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
        seed ^= static_cast<std::size_t>(key.y) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
        seed ^= static_cast<std::size_t>(key.z) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
        return seed;
    }
};

struct CoordKey2 {
    int u = 0;
    int v = 0;

    bool operator==(const CoordKey2& other) const { return u == other.u && v == other.v; }
};

struct CoordKey2Hasher {
    std::size_t operator()(const CoordKey2& key) const {
        std::size_t seed = 0;
        seed ^= static_cast<std::size_t>(key.u) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
        seed ^= static_cast<std::size_t>(key.v) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
        return seed;
    }
};

struct SurfacePatchSample {
    const ChSDFBrickPairRegionSample* sample = nullptr;
    CoordKey2 key;
};

struct SurfacePatchSelection {
    std::vector<SurfacePatchSample> samples;
    int min_u = 0;
    int max_u = -1;
    int min_v = 0;
    int max_v = -1;
};

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

double GradientQuality(const ChSDFProbeResult& probe) {
    return probe.valid ? std::abs(probe.gradient_world.Length() - 1.0) : 0.0;
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

double InverseEffectiveMassContribution(const ChSDFEffectiveMassProperties& body,
                                        const ChVector3d& point_world,
                                        const ChVector3d& normal_world) {
    if (!body.valid || body.fixed || body.mass <= 1.0e-20) {
        return 0.0;
    }

    const ChVector3d r = point_world - body.com_world;
    const ChVector3d rxn = Vcross(r, normal_world);
    const ChVector3d ang = body.inv_inertia_world * rxn;

    return (1.0 / body.mass) + Vdot(normal_world, Vcross(ang, r));
}

double EstimateEffectiveMassImpl(const ChSDFEffectiveMassProperties& body_a,
                                 const ChSDFEffectiveMassProperties& body_b,
                                 const ChVector3d& point_world,
                                 const ChVector3d& normal_world,
                                 double fallback_effective_mass) {
    const double inverse_mass = InverseEffectiveMassContribution(body_a, point_world, normal_world) +
                                InverseEffectiveMassContribution(body_b, point_world, normal_world);

    if (inverse_mass > 1.0e-20) {
        return 1.0 / inverse_mass;
    }

    return fallback_effective_mass > 0 ? fallback_effective_mass : 0.0;
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

std::vector<double> EstimatePatchCurvatureProxies(const ChSDFContactPatch& patch,
                                                  const ChSDFContactPatchSampler::Settings& patch_settings) {
    std::vector<double> curvature(patch.samples.size(), 0.0);
    if (patch.samples.empty() || patch_settings.samples_u == 0 || patch_settings.samples_v == 0) {
        return curvature;
    }

    std::vector<int> lookup(patch_settings.samples_u * patch_settings.samples_v, -1);
    for (std::size_t i = 0; i < patch.samples.size(); ++i) {
        const auto& sample = patch.samples[i];
        if (sample.iu < patch_settings.samples_u && sample.iv < patch_settings.samples_v) {
            lookup[sample.iv * patch_settings.samples_u + sample.iu] = static_cast<int>(i);
        }
    }

    const double du = std::abs(GridStep(patch_settings.samples_u, patch_settings.half_length_u));
    const double dv = std::abs(GridStep(patch_settings.samples_v, patch_settings.half_length_v));

    for (std::size_t i = 0; i < patch.samples.size(); ++i) {
        const auto& sample = patch.samples[i];
        double sum = 0;
        int count = 0;

        const struct Neighbor2D {
            int du;
            int dv;
            double distance;
        } neighbors[] = {
            {-1, 0, du}, {1, 0, du}, {0, -1, dv}, {0, 1, dv},
        };

        for (const auto& neighbor : neighbors) {
            const int nu = static_cast<int>(sample.iu) + neighbor.du;
            const int nv = static_cast<int>(sample.iv) + neighbor.dv;
            if (nu < 0 || nv < 0 || nu >= static_cast<int>(patch_settings.samples_u) ||
                nv >= static_cast<int>(patch_settings.samples_v) || neighbor.distance <= 1.0e-12) {
                continue;
            }

            const int neighbor_index = lookup[nv * static_cast<int>(patch_settings.samples_u) + nu];
            if (neighbor_index < 0) {
                continue;
            }

            const ChVector3d delta_normal =
                patch.samples[static_cast<std::size_t>(neighbor_index)].probe.normal_world - sample.probe.normal_world;
            sum += delta_normal.Length() / neighbor.distance;
            ++count;
        }

        if (count > 0) {
            curvature[i] = sum / static_cast<double>(count);
        }
    }

    return curvature;
}

std::vector<double> EstimateRegionCurvatureProxies(const ChSDFBrickPairRegion& region) {
    std::vector<double> curvature(region.samples.size(), 0.0);
    if (region.samples.empty() || region.sample_spacing <= 1.0e-12) {
        return curvature;
    }

    std::unordered_map<CoordKey3, std::size_t, CoordKey3Hasher> lookup;
    lookup.reserve(region.samples.size());
    for (std::size_t i = 0; i < region.samples.size(); ++i) {
        const auto& coord = region.samples[i].coord;
        lookup.emplace(CoordKey3{coord.x(), coord.y(), coord.z()}, i);
    }

    const CoordKey3 neighbors[] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};

    for (std::size_t i = 0; i < region.samples.size(); ++i) {
        const auto& sample = region.samples[i];
        double sum = 0;
        int count = 0;

        for (const auto& offset : neighbors) {
            const CoordKey3 key{sample.coord.x() + offset.x, sample.coord.y() + offset.y, sample.coord.z() + offset.z};
            const auto it = lookup.find(key);
            if (it == lookup.end()) {
                continue;
            }

            const ChVector3d delta_normal =
                region.samples[it->second].contact_normal_world - sample.contact_normal_world;
            sum += delta_normal.Length() / region.sample_spacing;
            ++count;
        }

        if (count > 0) {
            curvature[i] = sum / static_cast<double>(count);
        }
    }

    return curvature;
}

SurfacePatchSelection SelectSurfacePatchSamples(const ChSDFBrickPairRegion& region, double du, double dv) {
    SurfacePatchSelection selection;
    if (region.samples.empty()) {
        return selection;
    }

    const double safe_du = std::abs(du) > 1.0e-12 ? std::abs(du) : 1.0;
    const double safe_dv = std::abs(dv) > 1.0e-12 ? std::abs(dv) : 1.0;

    struct CandidateChoice {
        std::size_t sample_index = 0;
        double plane_distance = std::numeric_limits<double>::infinity();
        double gap_score = std::numeric_limits<double>::infinity();
    };

    std::unordered_map<CoordKey2, CandidateChoice, CoordKey2Hasher> best_by_cell;
    best_by_cell.reserve(region.samples.size());

    for (std::size_t i = 0; i < region.samples.size(); ++i) {
        const auto& sample = region.samples[i];
        const CoordKey2 key{static_cast<int>(std::llround(sample.point_patch.x() / safe_du)),
                            static_cast<int>(std::llround(sample.point_patch.y() / safe_dv))};
        const double plane_distance = std::abs(sample.point_patch.z());
        const double gap_score = std::abs(sample.combined_gap);

        const auto it = best_by_cell.find(key);
        if (it == best_by_cell.end() || plane_distance < it->second.plane_distance - 1.0e-12 ||
            (std::abs(plane_distance - it->second.plane_distance) <= 1.0e-12 && gap_score < it->second.gap_score)) {
            best_by_cell[key] = CandidateChoice{i, plane_distance, gap_score};
        }
    }

    if (best_by_cell.empty()) {
        return selection;
    }

    selection.min_u = std::numeric_limits<int>::max();
    selection.max_u = std::numeric_limits<int>::min();
    selection.min_v = std::numeric_limits<int>::max();
    selection.max_v = std::numeric_limits<int>::min();
    selection.samples.reserve(best_by_cell.size());

    for (const auto& item : best_by_cell) {
        selection.min_u = std::min(selection.min_u, item.first.u);
        selection.max_u = std::max(selection.max_u, item.first.u);
        selection.min_v = std::min(selection.min_v, item.first.v);
        selection.max_v = std::max(selection.max_v, item.first.v);
        selection.samples.push_back(SurfacePatchSample{&region.samples[item.second.sample_index], item.first});
    }

    std::stable_sort(selection.samples.begin(), selection.samples.end(),
                     [](const SurfacePatchSample& a, const SurfacePatchSample& b) {
                         if (a.key.v != b.key.v) {
                             return a.key.v < b.key.v;
                         }
                         return a.key.u < b.key.u;
                     });

    return selection;
}

std::vector<double> EstimateSurfaceRegionCurvatureProxies(const SurfacePatchSelection& selection,
                                                          double du,
                                                          double dv) {
    std::vector<double> curvature(selection.samples.size(), 0.0);
    if (selection.samples.empty()) {
        return curvature;
    }

    const double safe_du = std::abs(du) > 1.0e-12 ? std::abs(du) : 1.0;
    const double safe_dv = std::abs(dv) > 1.0e-12 ? std::abs(dv) : 1.0;

    std::unordered_map<CoordKey2, std::size_t, CoordKey2Hasher> lookup;
    lookup.reserve(selection.samples.size());
    for (std::size_t i = 0; i < selection.samples.size(); ++i) {
        lookup.emplace(selection.samples[i].key, i);
    }

    const struct Neighbor2D {
        int du;
        int dv;
        double distance;
    } neighbors[] = {
        {-1, 0, safe_du}, {1, 0, safe_du}, {0, -1, safe_dv}, {0, 1, safe_dv},
    };

    for (std::size_t i = 0; i < selection.samples.size(); ++i) {
        const auto& sample = selection.samples[i];
        double sum = 0;
        int count = 0;

        for (const auto& neighbor : neighbors) {
            const CoordKey2 key{sample.key.u + neighbor.du, sample.key.v + neighbor.dv};
            const auto it = lookup.find(key);
            if (it == lookup.end() || neighbor.distance <= 1.0e-12) {
                continue;
            }

            const ChVector3d delta_normal =
                selection.samples[it->second].sample->contact_normal_world - sample.sample->contact_normal_world;
            sum += delta_normal.Length() / neighbor.distance;
            ++count;
        }

        if (count > 0) {
            curvature[i] = sum / static_cast<double>(count);
        }
    }

    return curvature;
}

}  // namespace

ChSDFEffectiveMassProperties ChSDFContactWrenchEvaluator::MakeEffectiveMassProperties(ChBody* body) {
    ChSDFEffectiveMassProperties properties;
    if (!body) {
        return properties;
    }

    properties.valid = true;
    properties.fixed = body->IsFixed();
    properties.mass = body->GetMass();
    properties.com_world = body->GetFrameCOMToAbs().GetPos();

    const ChMatrix33<>& rot = body->GetFrameCOMToAbs().GetRotMat();
    properties.inv_inertia_world = rot * body->GetInvInertia() * rot.transpose();

    return properties;
}

double ChSDFContactWrenchEvaluator::EstimateEffectiveMass(const ChSDFEffectiveMassProperties& body_a,
                                                          const ChSDFEffectiveMassProperties& body_b,
                                                          const ChVector3d& point_world,
                                                          const ChVector3d& normal_world,
                                                          double fallback_effective_mass) {
    return EstimateEffectiveMassImpl(body_a, body_b, point_world, normal_world, fallback_effective_mass);
}

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
    const double resolution_length = 0.5 * (std::abs(du) + std::abs(dv));
    const auto curvature_proxies = EstimatePatchCurvatureProxies(patch, patch_settings);

    double pressure_weight_sum = 0;
    ChVector3d pressure_center_sum = VNULL;
    double local_stiffness_area_sum = 0;
    double effective_mass_area_sum = 0;

    result.samples.reserve(patch.samples.size());

    for (std::size_t sample_index = 0; sample_index < patch.samples.size(); ++sample_index) {
        const auto& patch_sample = patch.samples[sample_index];
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

        sample.gradient_quality = GradientQuality(patch_sample.probe);
        sample.curvature_proxy = curvature_proxies[sample_index];
        sample.resolution_length = resolution_length;
        sample.stiffness_scale = LocalStiffnessScale(sample.gradient_quality, sample.curvature_proxy,
                                                     sample.resolution_length, pressure_settings);
        sample.effective_mass = ClampEffectiveMass(
            EstimateEffectiveMassImpl(relative_kinematics.body_a, relative_kinematics.body_b, sample.point_shape,
                                      patch_sample.probe.normal_world,
                                      relative_kinematics.fallback_effective_mass > 0 ? relative_kinematics.fallback_effective_mass
                                                                                       : pressure_settings.fallback_effective_mass),
            pressure_settings);
        sample.local_stiffness = pressure_settings.stiffness * sample.stiffness_scale;
        sample.local_damping =
            ResolveLocalDamping(sample.local_stiffness, sample.effective_mass, sample.stiffness_scale, pressure_settings);

        sample.pressure = sample.local_stiffness * sample.penetration +
                          sample.local_damping * damping_speed + pressure_settings.adhesion_pressure;
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
            local_stiffness_area_sum += sample.local_stiffness * sample.quadrature_area;
            effective_mass_area_sum += sample.effective_mass * sample.quadrature_area;
            result.max_local_stiffness = std::max(result.max_local_stiffness, sample.local_stiffness);
            result.max_effective_mass = std::max(result.max_effective_mass, sample.effective_mass);
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
        result.mean_local_stiffness = local_stiffness_area_sum / result.active_area;
        result.mean_effective_mass = effective_mass_area_sum / result.active_area;
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
    const ChSDFEffectiveMassProperties& body_a,
    const ChSDFEffectiveMassProperties& body_b,
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
    const double quadrature_area = std::abs(du * dv);
    const double resolution_length = region.sample_spacing > 0 ? region.sample_spacing : 0.5 * (std::abs(du) + std::abs(dv));
    const auto surface_selection = SelectSurfacePatchSamples(region, du, dv);
    const auto curvature_proxies = EstimateSurfaceRegionCurvatureProxies(surface_selection, du, dv);

    double pressure_weight_sum = 0;
    ChVector3d pressure_center_sum = VNULL;
    double local_stiffness_area_sum = 0;
    double effective_mass_area_sum = 0;

    result.samples.reserve(surface_selection.samples.size());

    for (std::size_t sample_index = 0; sample_index < surface_selection.samples.size(); ++sample_index) {
        const auto& surface_sample = surface_selection.samples[sample_index];
        const auto& region_sample = *surface_sample.sample;
        ChSDFBrickPairWrenchSample sample;
        sample.region_sample = region_sample;
        sample.point_patch = region_sample.point_patch;
        const std::size_t count_u = surface_selection.max_u >= surface_selection.min_u
                                        ? static_cast<std::size_t>(surface_selection.max_u - surface_selection.min_u + 1)
                                        : 1;
        const std::size_t count_v = surface_selection.max_v >= surface_selection.min_v
                                        ? static_cast<std::size_t>(surface_selection.max_v - surface_selection.min_v + 1)
                                        : 1;
        const std::size_t index_u = static_cast<std::size_t>(surface_sample.key.u - surface_selection.min_u);
        const std::size_t index_v = static_cast<std::size_t>(surface_sample.key.v - surface_selection.min_v);
        sample.quadrature_area =
            quadrature_area * TrapezoidWeight(index_u, count_u) * TrapezoidWeight(index_v, count_v);

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
        sample.gradient_quality = 0.5 * (GradientQuality(region_sample.probe_a) + GradientQuality(region_sample.probe_b));
        sample.curvature_proxy = curvature_proxies[sample_index];
        sample.resolution_length = resolution_length;
        sample.stiffness_scale = LocalStiffnessScale(sample.gradient_quality, sample.curvature_proxy,
                                                     sample.resolution_length, pressure_settings);
        sample.effective_mass = ClampEffectiveMass(
            EstimateEffectiveMassImpl(body_a, body_b, region_sample.point_world, region_sample.contact_normal_world,
                                      pressure_settings.fallback_effective_mass),
            pressure_settings);
        sample.local_stiffness = pressure_settings.stiffness * sample.stiffness_scale;
        sample.local_damping =
            ResolveLocalDamping(sample.local_stiffness, sample.effective_mass, sample.stiffness_scale, pressure_settings);
        sample.pressure = sample.local_stiffness * sample.penetration +
                          sample.local_damping * damping_speed + pressure_settings.adhesion_pressure;
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
            local_stiffness_area_sum += sample.local_stiffness * sample.quadrature_area;
            effective_mass_area_sum += sample.effective_mass * sample.quadrature_area;
            result.max_local_stiffness = std::max(result.max_local_stiffness, sample.local_stiffness);
            result.max_effective_mass = std::max(result.max_effective_mass, sample.effective_mass);
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
        result.mean_local_stiffness = local_stiffness_area_sum / result.active_area;
        result.mean_effective_mass = effective_mass_area_sum / result.active_area;
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
    const ChSDFEffectiveMassProperties& body_a,
    const ChSDFEffectiveMassProperties& body_b,
    const ChSDFNormalPressureSettings& pressure_settings) {
    ChSDFShapePairContactResult result;
    result.valid = true;
    result.shape_a_frame_abs = shape_a_frame_abs;
    result.shape_b_frame_abs = shape_b_frame_abs;
    result.total_regions = regions.size();
    result.regions.reserve(regions.size());
    double local_stiffness_area_sum = 0;
    double effective_mass_area_sum = 0;

    for (const auto& region : regions) {
        auto region_result =
            EvaluateBrickPairRegion(region, shape_a_frame_abs, shape_b_frame_abs, body_a, body_b, pressure_settings);

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
