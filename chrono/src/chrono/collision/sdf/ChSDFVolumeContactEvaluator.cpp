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

#include "chrono/collision/sdf/ChSDFVolumeContactEvaluator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace chrono {
namespace {

struct GridCoordKey {
    int x = 0;
    int y = 0;
    int z = 0;

    bool operator==(const GridCoordKey& other) const { return x == other.x && y == other.y && z == other.z; }
};

struct GridCoordHash {
    std::size_t operator()(const GridCoordKey& key) const {
        const std::size_t hx = static_cast<std::size_t>(static_cast<uint32_t>(key.x));
        const std::size_t hy = static_cast<std::size_t>(static_cast<uint32_t>(key.y));
        const std::size_t hz = static_cast<std::size_t>(static_cast<uint32_t>(key.z));
        return hx * 73856093u ^ hy * 19349663u ^ hz * 83492791u;
    }
};

struct IndexRange {
    int min = 0;
    int max = -1;

    bool IsValid() const { return min <= max; }
};

struct GridCellCandidate {
    GridCoordKey coord;
    ChVector3d center_world = VNULL;
    std::size_t brick_a_index = 0;
    std::size_t brick_b_index = 0;
};

struct ActiveVoxelContribution {
    GridCoordKey coord;
    ChSDFBrickPairWrenchSample sample;
};

double MinPositiveComponent(const ChVector3d& v) {
    double result = std::numeric_limits<double>::infinity();
    if (v.x() > 0) {
        result = std::min(result, v.x());
    }
    if (v.y() > 0) {
        result = std::min(result, v.y());
    }
    if (v.z() > 0) {
        result = std::min(result, v.z());
    }
    return std::isfinite(result) ? result : 0.0;
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
    return speed > 1.0e-12 ? tangential_velocity * (-settings.friction_coefficient * pressure / speed) : VNULL;
}

ChVector3d TransformPotentialGradientWorld(const ChFrameMoving<>& shape_frame_abs,
                                           const ChSDFPotentialFieldProbe& probe) {
    return shape_frame_abs.TransformDirectionLocalToParent(probe.grad_p0_local);
}

double ResolveSpacing(const ChCollisionShapeSDF& shape_a,
                      const ChCollisionShapeSDF& shape_b,
                      const ChSDFContactRegionBuilder::Settings& region_settings,
                      const std::vector<ChSDFBrickPairCandidate>& brick_pairs) {
    if (region_settings.sample_spacing > 0) {
        return region_settings.sample_spacing;
    }

    double spacing = std::numeric_limits<double>::infinity();
    const auto info_a = shape_a.GetGridInfo();
    const auto info_b = shape_b.GetGridInfo();
    if (info_a.valid) {
        spacing = std::min(spacing, MinPositiveComponent(info_a.voxel_size));
    }
    if (info_b.valid) {
        spacing = std::min(spacing, MinPositiveComponent(info_b.voxel_size));
    }

    for (const auto& pair : brick_pairs) {
        spacing = std::min(spacing, MinPositiveComponent(pair.brick_a.brick.voxel_size));
        spacing = std::min(spacing, MinPositiveComponent(pair.brick_b.brick.voxel_size));
    }

    return std::isfinite(spacing) && spacing > 0 ? spacing : 1.0e-3;
}

IndexRange ComputeCenterIndexRange(double min_value, double max_value, double center_origin, double spacing) {
    IndexRange range;
    if (max_value < min_value || spacing <= 0) {
        return range;
    }

    range.min = static_cast<int>(std::ceil((min_value - center_origin) / spacing));
    range.max = static_cast<int>(std::floor((max_value - center_origin) / spacing));
    return range;
}

double CenterOrigin(double min_value, double spacing) {
    return std::floor(min_value / spacing) * spacing + 0.5 * spacing;
}

ChVector3d GridCenter(const GridCoordKey& coord, const ChVector3d& origin, double spacing) {
    return origin + ChVector3d(coord.x * spacing, coord.y * spacing, coord.z * spacing);
}

double RegularizedDelta(double h_value, double eta_h) {
    if (eta_h <= 1.0e-12 || std::abs(h_value) > eta_h) {
        return 0.0;
    }

    return 0.5 / eta_h * (1.0 + std::cos(CH_PI * h_value / eta_h));
}

double ResolvePressureBandWidth(double grad_h_norm, double spacing) {
    const double geometric_band = 2.0 * std::max(spacing, 1.0e-8);
    return std::max(grad_h_norm * geometric_band, 1.0e-12);
}

std::vector<ChVector3i> BuildNeighborOffsets(int neighbor_mode) {
    std::vector<ChVector3i> offsets;
    offsets.reserve(26);

    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) {
                    continue;
                }

                const int l1 = std::abs(dx) + std::abs(dy) + std::abs(dz);
                const bool keep = (neighbor_mode == 6) ? (l1 == 1) : ((neighbor_mode == 18) ? (l1 <= 2) : true);
                if (keep) {
                    offsets.emplace_back(dx, dy, dz);
                }
            }
        }
    }

    return offsets;
}

void FinalizeAggregateResult(ChSDFShapePairContactResult& result) {
    result.total_regions = result.regions.size();
    result.active_regions = 0;
    result.active_area = 0;
    result.integrated_pressure = 0;
    result.mean_pressure = 0;
    result.mean_local_stiffness = 0;
    result.mean_effective_mass = 0;
    result.max_local_stiffness = 0;
    result.max_effective_mass = 0;
    result.max_penetration = 0;
    result.max_pressure = 0;
    result.wrench_shape_a = {VNULL, VNULL};
    result.wrench_shape_b = {VNULL, VNULL};
    result.wrench_world_a = {VNULL, VNULL};
    result.wrench_world_b = {VNULL, VNULL};

    double local_stiffness_area_sum = 0;
    double effective_mass_area_sum = 0;

    for (const auto& region_result : result.regions) {
        result.wrench_shape_a.force += region_result.wrench_shape_a.force;
        result.wrench_shape_a.torque += region_result.wrench_shape_a.torque;
        result.wrench_shape_b.force += region_result.wrench_shape_b.force;
        result.wrench_shape_b.torque += region_result.wrench_shape_b.torque;
        result.wrench_world_a.force += region_result.wrench_world_a.force;
        result.wrench_world_a.torque += region_result.wrench_world_a.torque;
        result.wrench_world_b.force += region_result.wrench_world_b.force;
        result.wrench_world_b.torque += region_result.wrench_world_b.torque;

        if (!region_result.HasActiveContact()) {
            continue;
        }

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

    if (result.active_area > 0) {
        result.mean_pressure = result.integrated_pressure / result.active_area;
        result.mean_local_stiffness = local_stiffness_area_sum / result.active_area;
        result.mean_effective_mass = effective_mass_area_sum / result.active_area;
    }
}

}  // namespace

ChSDFShapePairContactResult ChSDFVolumeContactEvaluator::EvaluateBrickPairs(
    const ChCollisionShapeSDF& shape_a,
    const ChFrameMoving<>& shape_a_frame_abs,
    const ChCollisionShapeSDF& shape_b,
    const ChFrameMoving<>& shape_b_frame_abs,
    const std::vector<ChSDFBrickPairCandidate>& brick_pairs,
    const ChSDFContactRegionBuilder::Settings& region_settings,
    const ChSDFEffectiveMassProperties& body_a,
    const ChSDFEffectiveMassProperties& body_b,
    const ChSDFNormalPressureSettings& pressure_settings) {
    ChSDFShapePairContactResult result;
    result.valid = true;
    result.shape_a_frame_abs = shape_a_frame_abs;
    result.shape_b_frame_abs = shape_b_frame_abs;

    if (brick_pairs.empty()) {
        return result;
    }

    const double spacing = ResolveSpacing(shape_a, shape_b, region_settings, brick_pairs);
    if (spacing <= 0) {
        return result;
    }

    ChAABB union_bounds;
    bool has_union = false;
    for (const auto& pair : brick_pairs) {
        if (pair.overlap_world.IsInverted()) {
            continue;
        }
        union_bounds += pair.overlap_world;
        has_union = true;
    }

    if (!has_union || union_bounds.IsInverted()) {
        return result;
    }

    const ChVector3d grid_origin(CenterOrigin(union_bounds.min.x(), spacing), CenterOrigin(union_bounds.min.y(), spacing),
                                 CenterOrigin(union_bounds.min.z(), spacing));

    std::unordered_map<GridCoordKey, GridCellCandidate, GridCoordHash> cells;
    for (const auto& pair : brick_pairs) {
        if (pair.overlap_world.IsInverted()) {
            continue;
        }

        const IndexRange ix =
            ComputeCenterIndexRange(pair.overlap_world.min.x(), pair.overlap_world.max.x(), grid_origin.x(), spacing);
        const IndexRange iy =
            ComputeCenterIndexRange(pair.overlap_world.min.y(), pair.overlap_world.max.y(), grid_origin.y(), spacing);
        const IndexRange iz =
            ComputeCenterIndexRange(pair.overlap_world.min.z(), pair.overlap_world.max.z(), grid_origin.z(), spacing);
        if (!ix.IsValid() || !iy.IsValid() || !iz.IsValid()) {
            continue;
        }

        for (int kz = iz.min; kz <= iz.max; ++kz) {
            for (int ky = iy.min; ky <= iy.max; ++ky) {
                for (int kx = ix.min; kx <= ix.max; ++kx) {
                    GridCoordKey coord{kx, ky, kz};
                    auto [it, inserted] = cells.emplace(coord, GridCellCandidate());
                    if (inserted) {
                        it->second.coord = coord;
                        it->second.center_world = GridCenter(coord, grid_origin, spacing);
                        it->second.brick_a_index = pair.brick_a_index;
                        it->second.brick_b_index = pair.brick_b_index;
                    }
                }
            }
        }
    }

    if (cells.empty()) {
        return result;
    }

    const double cell_volume = spacing * spacing * spacing;
    std::vector<ActiveVoxelContribution> active_samples;
    active_samples.reserve(cells.size());

    for (const auto& item : cells) {
        const GridCellCandidate& cell = item.second;
        const ChVector3d point_shape_a = shape_a_frame_abs.TransformPointParentToLocal(cell.center_world);
        const ChVector3d point_shape_b = shape_b_frame_abs.TransformPointParentToLocal(cell.center_world);
        const ChSDFPotentialFieldProbe probe_a = shape_a.ProbePotentialLocal(point_shape_a);
        const ChSDFPotentialFieldProbe probe_b = shape_b.ProbePotentialLocal(point_shape_b);

        if (!probe_a.valid || !probe_b.valid || probe_a.phi > 0 || probe_b.phi > 0) {
            continue;
        }

        const double p_a = probe_a.p0;
        const double p_b = probe_b.p0;
        const double h_value = p_a - p_b;

        const ChVector3d grad_p_a_world = TransformPotentialGradientWorld(shape_a_frame_abs, probe_a);
        const ChVector3d grad_p_b_world = TransformPotentialGradientWorld(shape_b_frame_abs, probe_b);
        const ChVector3d grad_h_world = grad_p_a_world - grad_p_b_world;
        const double grad_h_norm = grad_h_world.Length();
        if (grad_h_norm <= 1.0e-12) {
            continue;
        }

        const double eta_h = ResolvePressureBandWidth(grad_h_norm, spacing);
        if (std::abs(h_value) > eta_h) {
            continue;
        }

        const double area_weight = RegularizedDelta(h_value, eta_h) * grad_h_norm * cell_volume;
        if (area_weight <= 1.0e-16) {
            continue;
        }

        const ChVector3d normal_world = -(grad_h_world / grad_h_norm);
        const ChVector3d point_world = cell.center_world - grad_h_world * (h_value / (grad_h_norm * grad_h_norm));
        const ChVector3d point_surface_a = shape_a_frame_abs.TransformPointParentToLocal(point_world);
        const ChVector3d point_surface_b = shape_b_frame_abs.TransformPointParentToLocal(point_world);

        ChSDFBrickPairWrenchSample sample;
        sample.active = true;
        sample.quadrature_area = area_weight;
        sample.point_patch = VNULL;
        sample.region_sample.brick_a_index = cell.brick_a_index;
        sample.region_sample.brick_b_index = cell.brick_b_index;
        sample.region_sample.coord = ChVector3i(cell.coord.x, cell.coord.y, cell.coord.z);
        sample.region_sample.point_world = point_world;
        sample.region_sample.point_shape_a = point_shape_a;
        sample.region_sample.point_shape_b = point_shape_b;
        sample.region_sample.surface_world_a = point_world;
        sample.region_sample.surface_world_b = point_world;
        sample.region_sample.surface_shape_a = point_surface_a;
        sample.region_sample.surface_shape_b = point_surface_b;
        sample.region_sample.normal_world_a = shape_a_frame_abs.TransformDirectionLocalToParent(probe_a.normal_local);
        sample.region_sample.normal_world_b = shape_b_frame_abs.TransformDirectionLocalToParent(probe_b.normal_local);
        sample.region_sample.contact_normal_world = normal_world;
        sample.region_sample.distance_a = probe_a.phi;
        sample.region_sample.distance_b = probe_b.phi;
        sample.region_sample.normal_gap = -h_value / grad_h_norm;
        sample.region_sample.combined_gap = probe_a.phi + probe_b.phi;
        sample.region_sample.normal_opposition =
            -Vdot(sample.region_sample.normal_world_a, sample.region_sample.normal_world_b);
        sample.region_sample.area_weight = area_weight;

        sample.quadrature_point.point_world = point_world;
        sample.quadrature_point.point_shape_a = point_surface_a;
        sample.quadrature_point.point_shape_b = point_surface_b;
        sample.quadrature_point.normal_world = normal_world;
        sample.quadrature_point.area_weight = area_weight;
        sample.quadrature_point.p0_a = p_a;
        sample.quadrature_point.p0_b = p_b;
        sample.quadrature_point.p_cap = 0.5 * (p_a + p_b);
        sample.quadrature_point.h_value = h_value;
        sample.quadrature_point.support_overlap = 2.0 * eta_h / grad_h_norm;
        sample.quadrature_point.probe_a = probe_a;
        sample.quadrature_point.probe_b = probe_b;

        const ChVector3d speed_world_a = shape_a_frame_abs.PointSpeedLocalToParent(point_surface_a);
        const ChVector3d speed_world_b = shape_b_frame_abs.PointSpeedLocalToParent(point_surface_b);
        const ChVector3d relative_speed_world = speed_world_b - speed_world_a;
        const double normal_velocity_component = Vdot(relative_speed_world, normal_world);
        sample.normal_speed = -normal_velocity_component;
        sample.tangential_velocity_world = relative_speed_world - normal_world * normal_velocity_component;
        sample.tangential_speed = sample.tangential_velocity_world.Length();
        sample.gradient_quality = 0.5 * (GradientQuality(probe_a) + GradientQuality(probe_b));
        sample.curvature_proxy = 0.0;
        sample.resolution_length = spacing;
        sample.stiffness_scale =
            LocalStiffnessScale(sample.gradient_quality, sample.curvature_proxy, sample.resolution_length, pressure_settings);
        sample.effective_mass = ClampEffectiveMass(
            ChSDFContactWrenchEvaluator::EstimateEffectiveMass(body_a, body_b, point_world, normal_world,
                                                               pressure_settings.fallback_effective_mass),
            pressure_settings);
        sample.local_stiffness = grad_h_norm * sample.stiffness_scale;
        sample.local_damping = pressure_settings.damping * sample.local_stiffness;

        const double base_pressure = sample.quadrature_point.p_cap * sample.stiffness_scale;
        const double damping_speed =
            pressure_settings.use_only_closing_speed ? std::max(sample.normal_speed, 0.0) : sample.normal_speed;
        double pressure = base_pressure;
        if (pressure_settings.damping_ratio >= 0) {
            const double local_damping =
                2.0 * pressure_settings.damping_ratio *
                std::sqrt(std::max(sample.local_stiffness, 0.0) * std::max(sample.effective_mass, 0.0));
            sample.local_damping = local_damping;
            pressure += local_damping * damping_speed;
        } else {
            pressure *= (1.0 + pressure_settings.damping * damping_speed);
        }
        sample.pressure = ClampPressure(pressure + pressure_settings.adhesion_pressure, pressure_settings);
        sample.penetration =
            sample.local_stiffness > 1.0e-12 ? std::max(base_pressure / sample.local_stiffness, 0.0) : 0.0;
        sample.signed_gap = -sample.penetration;

        if (sample.pressure <= 0) {
            continue;
        }

        sample.traction_world = ComputeTangentialTraction(sample.tangential_velocity_world, sample.pressure, pressure_settings);
        const ChVector3d normal_force_world_b = normal_world * (sample.pressure * sample.quadrature_area);
        const ChVector3d tangential_force_world_b = sample.traction_world * sample.quadrature_area;
        const ChVector3d force_world_b = normal_force_world_b + tangential_force_world_b;
        const ChVector3d force_world_a = -force_world_b;
        const ChVector3d torque_world_a = Vcross(point_world - shape_a_frame_abs.GetPos(), force_world_a);
        const ChVector3d torque_world_b = Vcross(point_world - shape_b_frame_abs.GetPos(), force_world_b);
        sample.force_world = force_world_b;
        sample.force_shape_a = shape_a_frame_abs.TransformDirectionParentToLocal(force_world_a);
        sample.torque_shape_a = shape_a_frame_abs.TransformDirectionParentToLocal(torque_world_a);
        sample.force_shape_b = shape_b_frame_abs.TransformDirectionParentToLocal(force_world_b);
        sample.torque_shape_b = shape_b_frame_abs.TransformDirectionParentToLocal(torque_world_b);

        active_samples.push_back({cell.coord, std::move(sample)});
    }

    if (active_samples.empty()) {
        return result;
    }

    std::unordered_map<GridCoordKey, std::size_t, GridCoordHash> index_by_coord;
    index_by_coord.reserve(active_samples.size());
    for (std::size_t i = 0; i < active_samples.size(); ++i) {
        index_by_coord.emplace(active_samples[i].coord, i);
    }

    const auto neighbor_offsets = BuildNeighborOffsets(region_settings.neighbor_mode);
    std::vector<char> visited(active_samples.size(), 0);
    std::vector<ChSDFBrickPairWrenchResult> regions;
    regions.reserve(active_samples.size());

    for (std::size_t seed = 0; seed < active_samples.size(); ++seed) {
        if (visited[seed]) {
            continue;
        }

        std::queue<std::size_t> frontier;
        frontier.push(seed);
        visited[seed] = 1;

        std::vector<std::size_t> component;
        component.reserve(64);

        while (!frontier.empty()) {
            const std::size_t current = frontier.front();
            frontier.pop();
            component.push_back(current);

            const GridCoordKey coord = active_samples[current].coord;
            for (const auto& offset : neighbor_offsets) {
                const GridCoordKey neighbor{coord.x + offset.x(), coord.y + offset.y(), coord.z + offset.z()};
                auto it = index_by_coord.find(neighbor);
                if (it == index_by_coord.end() || visited[it->second]) {
                    continue;
                }

                visited[it->second] = 1;
                frontier.push(it->second);
            }
        }

        if (component.size() < region_settings.min_region_samples) {
            continue;
        }

        ChSDFBrickPairWrenchResult region_result;
        region_result.total_samples = component.size();

        std::unordered_set<std::size_t> brick_a_indices;
        std::unordered_set<std::size_t> brick_b_indices;
        double pressure_weight_sum = 0;
        ChVector3d pressure_center_sum = VNULL;
        double local_stiffness_area_sum = 0;
        double effective_mass_area_sum = 0;

        for (const std::size_t index : component) {
            const auto& voxel = active_samples[index];
            const auto& sample = voxel.sample;
            region_result.samples.push_back(sample);
            region_result.region.samples.push_back(sample.region_sample);
            brick_a_indices.insert(sample.region_sample.brick_a_index);
            brick_b_indices.insert(sample.region_sample.brick_b_index);

            region_result.patch_area += sample.quadrature_area;
            region_result.active_samples++;
            region_result.active_area += sample.quadrature_area;
            region_result.integrated_pressure += sample.pressure * sample.quadrature_area;
            local_stiffness_area_sum += sample.local_stiffness * sample.quadrature_area;
            effective_mass_area_sum += sample.effective_mass * sample.quadrature_area;
            region_result.max_local_stiffness = std::max(region_result.max_local_stiffness, sample.local_stiffness);
            region_result.max_effective_mass = std::max(region_result.max_effective_mass, sample.effective_mass);
            region_result.max_penetration = std::max(region_result.max_penetration, sample.penetration);
            region_result.max_pressure = std::max(region_result.max_pressure, sample.pressure);

            region_result.wrench_world_a.force -= sample.force_world;
            region_result.wrench_world_a.torque +=
                shape_a_frame_abs.TransformDirectionLocalToParent(sample.torque_shape_a);
            region_result.wrench_world_b.force += sample.force_world;
            region_result.wrench_world_b.torque +=
                shape_b_frame_abs.TransformDirectionLocalToParent(sample.torque_shape_b);
            region_result.wrench_shape_a.force += sample.force_shape_a;
            region_result.wrench_shape_a.torque += sample.torque_shape_a;
            region_result.wrench_shape_b.force += sample.force_shape_b;
            region_result.wrench_shape_b.torque += sample.torque_shape_b;

            const double sample_weight = sample.pressure * sample.quadrature_area;
            pressure_weight_sum += sample_weight;
            pressure_center_sum += sample.region_sample.point_world * sample_weight;
        }

        region_result.region.sample_spacing = spacing;
        region_result.region.brick_a_indices.assign(brick_a_indices.begin(), brick_a_indices.end());
        region_result.region.brick_b_indices.assign(brick_b_indices.begin(), brick_b_indices.end());
        std::sort(region_result.region.brick_a_indices.begin(), region_result.region.brick_a_indices.end());
        std::sort(region_result.region.brick_b_indices.begin(), region_result.region.brick_b_indices.end());
        ChSDFContactRegionBuilder::ReparameterizeBrickPairRegion(
            region_result.region, ChFrame<>(shape_a_frame_abs.GetCoordsys()), ChFrame<>(shape_b_frame_abs.GetCoordsys()));

        if (region_result.active_area > 0) {
            region_result.mean_pressure = region_result.integrated_pressure / region_result.active_area;
            region_result.mean_local_stiffness = local_stiffness_area_sum / region_result.active_area;
            region_result.mean_effective_mass = effective_mass_area_sum / region_result.active_area;
        }
        if (pressure_weight_sum > 0) {
            region_result.pressure_center_world = pressure_center_sum / pressure_weight_sum;
        }

        regions.push_back(std::move(region_result));
    }

    std::stable_sort(regions.begin(), regions.end(), [](const ChSDFBrickPairWrenchResult& a, const ChSDFBrickPairWrenchResult& b) {
        return a.samples.size() > b.samples.size();
    });

    for (std::size_t i = 0; i < regions.size(); ++i) {
        regions[i].region.region_id = i;
    }

    result.regions = std::move(regions);
    FinalizeAggregateResult(result);
    return result;
}

}  // end namespace chrono
