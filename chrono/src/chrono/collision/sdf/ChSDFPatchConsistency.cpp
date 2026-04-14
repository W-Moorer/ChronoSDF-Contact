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

#include "chrono/collision/sdf/ChSDFPatchConsistency.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>

namespace chrono {
namespace {

bool IsFiniteVector(const ChVector3d& v) {
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
}

ChVector3d SafeNormalized(const ChVector3d& v, const ChVector3d& fallback = VECT_Y) {
    const double length = v.Length();
    return length > 1.0e-12 ? v / length : fallback;
}

void BuildOrthonormalBasis(const ChVector3d& normal, ChVector3d& tangent_u, ChVector3d& tangent_v) {
    const ChVector3d n = SafeNormalized(normal, VECT_Y);
    const ChVector3d reference = (std::abs(n.y()) < 0.9) ? VECT_Y : VECT_X;
    tangent_u = SafeNormalized(Vcross(reference, n), VECT_X);
    tangent_v = SafeNormalized(Vcross(n, tangent_u), VECT_Z);
}

double ClampValue(double value, double lower, double upper) {
    return std::max(lower, std::min(value, upper));
}

ChVector3d ResolveSamplePointWorld(const ChSDFBrickPairWrenchSample& sample) {
    if (IsFiniteVector(sample.quadrature_point.point_world)) {
        return sample.quadrature_point.point_world;
    }
    if (IsFiniteVector(sample.region_sample.point_world)) {
        return sample.region_sample.point_world;
    }
    if (IsFiniteVector(sample.region_sample.surface_world_b)) {
        return sample.region_sample.surface_world_b;
    }
    return VNULL;
}

void AccumulateSampleIntoPatch(ChSDFPatchBandAggregate& aggregate,
                               const ChSDFBrickPairWrenchSample& sample,
                               const ChFrameMoving<>& shape_a_frame_abs,
                               const ChFrameMoving<>& shape_b_frame_abs) {
    if (!sample.active || sample.quadrature_area <= 1.0e-16) {
        return;
    }

    const ChVector3d point_world = ResolveSamplePointWorld(sample);
    const double pressure_weight = std::max(sample.pressure * sample.quadrature_area, 0.0);

    aggregate.band_area += sample.quadrature_area;
    aggregate.integrated_pressure += pressure_weight;
    aggregate.centroid_world += point_world * sample.quadrature_area;
    aggregate.pressure_center_world += point_world * pressure_weight;
    aggregate.wrench_world_b.force += sample.force_world;
    aggregate.wrench_world_b.torque += Vcross(point_world - shape_b_frame_abs.GetPos(), sample.force_world);
    aggregate.wrench_world_a.force -= sample.force_world;
    aggregate.wrench_world_a.torque += Vcross(point_world - shape_a_frame_abs.GetPos(), -sample.force_world);
    aggregate.support_bbox_world += point_world;
}

ChSDFPatchBandAggregate FinalizeAggregate(ChSDFPatchBandAggregate aggregate) {
    if (aggregate.band_area > 0) {
        aggregate.centroid_world /= aggregate.band_area;
    }
    if (aggregate.integrated_pressure > 0) {
        aggregate.pressure_center_world /= aggregate.integrated_pressure;
    } else {
        aggregate.pressure_center_world = aggregate.centroid_world;
    }
    return aggregate;
}

ChSDFPatchConsistencyResult MakePassthroughPatchResult(const ChSDFPatchBandAggregate& aggregate) {
    ChSDFPatchConsistencyResult result;
    result.region_id = aggregate.region_id;
    result.patch_id = aggregate.patch_id;
    result.residual_patch = aggregate.residual_patch;
    result.band_area = aggregate.band_area;
    result.sheet_area = aggregate.band_area;
    result.alpha = 1.0;
    result.integrated_pressure_band = aggregate.integrated_pressure;
    result.integrated_pressure_corrected = aggregate.integrated_pressure;
    result.wrench_world_a_band = aggregate.wrench_world_a;
    result.wrench_world_b_band = aggregate.wrench_world_b;
    result.wrench_world_a_corrected = aggregate.wrench_world_a;
    result.wrench_world_b_corrected = aggregate.wrench_world_b;
    result.centroid_world = aggregate.centroid_world;
    result.pressure_center_world = aggregate.pressure_center_world;
    result.sheet_pressure_center_world = aggregate.pressure_center_world;
    result.corrected_pressure_center_world = aggregate.pressure_center_world;
    result.support_bbox_world = aggregate.support_bbox_world;
    return result;
}

ChWrenchd BuildFirstMomentCorrectedWrench(const ChWrenchd& wrench_band,
                                          const ChVector3d& origin_world,
                                          const ChVector3d& band_pressure_center_world,
                                          const ChVector3d& corrected_pressure_center_world,
                                          double alpha) {
    ChWrenchd corrected = {alpha * wrench_band.force, VNULL};

    const ChVector3d band_arm = band_pressure_center_world - origin_world;
    const ChVector3d corrected_arm = corrected_pressure_center_world - origin_world;
    const ChVector3d intrinsic_band = wrench_band.torque - Vcross(band_arm, wrench_band.force);
    corrected.torque = alpha * intrinsic_band + Vcross(corrected_arm, corrected.force);
    return corrected;
}

struct ProjectedRect {
    double min_u = 0;
    double max_u = 0;
    double min_v = 0;
    double max_v = 0;
    bool valid = false;
};

ProjectedRect MakeCenteredRect(double center_u, double center_v, double area) {
    ProjectedRect rect;
    const double half_extent = 0.5 * std::sqrt(std::max(area, 0.0));
    rect.min_u = center_u - half_extent;
    rect.max_u = center_u + half_extent;
    rect.min_v = center_v - half_extent;
    rect.max_v = center_v + half_extent;
    rect.valid = std::isfinite(rect.min_u) && std::isfinite(rect.max_u) && std::isfinite(rect.min_v) &&
                 std::isfinite(rect.max_v);
    return rect;
}

ProjectedRect ProjectSupportRect(const ChAABB& bbox,
                                 const ChVector3d& origin_world,
                                 const ChVector3d& tangent_u,
                                 const ChVector3d& tangent_v) {
    ProjectedRect rect;
    if (bbox.IsInverted()) {
        return rect;
    }

    rect.min_u = std::numeric_limits<double>::infinity();
    rect.max_u = -std::numeric_limits<double>::infinity();
    rect.min_v = std::numeric_limits<double>::infinity();
    rect.max_v = -std::numeric_limits<double>::infinity();

    const std::array<ChVector3d, 8> corners = {
        ChVector3d(bbox.min.x(), bbox.min.y(), bbox.min.z()), ChVector3d(bbox.min.x(), bbox.min.y(), bbox.max.z()),
        ChVector3d(bbox.min.x(), bbox.max.y(), bbox.min.z()), ChVector3d(bbox.min.x(), bbox.max.y(), bbox.max.z()),
        ChVector3d(bbox.max.x(), bbox.min.y(), bbox.min.z()), ChVector3d(bbox.max.x(), bbox.min.y(), bbox.max.z()),
        ChVector3d(bbox.max.x(), bbox.max.y(), bbox.min.z()), ChVector3d(bbox.max.x(), bbox.max.y(), bbox.max.z())};

    for (const auto& corner : corners) {
        const ChVector3d rel = corner - origin_world;
        const double u = Vdot(rel, tangent_u);
        const double v = Vdot(rel, tangent_v);
        rect.min_u = std::min(rect.min_u, u);
        rect.max_u = std::max(rect.max_u, u);
        rect.min_v = std::min(rect.min_v, v);
        rect.max_v = std::max(rect.max_v, v);
    }

    rect.valid = std::isfinite(rect.min_u) && std::isfinite(rect.max_u) && std::isfinite(rect.min_v) &&
                 std::isfinite(rect.max_v);
    return rect;
}

ProjectedRect RescaleRectToArea(const ProjectedRect& rect, double target_area) {
    if (!rect.valid) {
        return rect;
    }

    ProjectedRect scaled = rect;
    const double span_u = std::max(rect.max_u - rect.min_u, 0.0);
    const double span_v = std::max(rect.max_v - rect.min_v, 0.0);
    const double rect_area = span_u * span_v;
    if (rect_area <= 1.0e-16 || target_area <= 1.0e-16) {
        return scaled;
    }

    const double center_u = 0.5 * (rect.min_u + rect.max_u);
    const double center_v = 0.5 * (rect.min_v + rect.max_v);
    const double scale = std::sqrt(target_area / rect_area);
    const double half_u = 0.5 * span_u * scale;
    const double half_v = 0.5 * span_v * scale;
    scaled.min_u = center_u - half_u;
    scaled.max_u = center_u + half_u;
    scaled.min_v = center_v - half_v;
    scaled.max_v = center_v + half_v;
    return scaled;
}

ProjectedRect ProjectSheetSampleSupport(const ChSDFSheetFiberSample& sample,
                                        const ChVector3d& patch_origin_world,
                                        const ChVector3d& tangent_u,
                                        const ChVector3d& tangent_v) {
    auto rect = ProjectSupportRect(sample.support_bbox_world, patch_origin_world, tangent_u, tangent_v);
    const double target_area = std::max(sample.footprint_area, sample.measure_area);
    if (rect.valid) {
        rect = RescaleRectToArea(rect, target_area);
        return rect;
    }

    const ChVector3d rel = sample.centroid_world - patch_origin_world;
    return MakeCenteredRect(Vdot(rel, tangent_u), Vdot(rel, tangent_v), target_area);
}

double OverlapInterval(double a0, double a1, double b0, double b1, double& min_overlap, double& max_overlap) {
    min_overlap = std::max(a0, b0);
    max_overlap = std::min(a1, b1);
    return std::max(max_overlap - min_overlap, 0.0);
}

struct RedistributionContribution {
    double corrected_area = 0;
    ChVector3d corrected_point_world = VNULL;
};

RedistributionContribution ComputeRedistributedContribution(const ChSDFBrickPairWrenchSample& sample,
                                                            const ChVector3d& patch_origin_world,
                                                            const std::vector<ProjectedRect>& support_rects,
                                                            const ChVector3d& tangent_u,
                                                            const ChVector3d& tangent_v,
                                                            double sample_half_extent) {
    RedistributionContribution contribution;
    if (!sample.active || sample.quadrature_area <= 1.0e-16 || sample_half_extent <= 1.0e-12 || support_rects.empty()) {
        return contribution;
    }

    const ChVector3d point_world = ResolveSamplePointWorld(sample);
    const ChVector3d rel = point_world - patch_origin_world;
    const double center_u = Vdot(rel, tangent_u);
    const double center_v = Vdot(rel, tangent_v);
    const double half = sample_half_extent;

    double overlap_area_sum = 0;
    double overlap_center_u_sum = 0;
    double overlap_center_v_sum = 0;
    for (const auto& support_rect : support_rects) {
        if (!support_rect.valid) {
            continue;
        }

        double min_u = 0;
        double max_u = 0;
        double min_v = 0;
        double max_v = 0;
        const double overlap_u =
            OverlapInterval(center_u - half, center_u + half, support_rect.min_u, support_rect.max_u, min_u, max_u);
        const double overlap_v =
            OverlapInterval(center_v - half, center_v + half, support_rect.min_v, support_rect.max_v, min_v, max_v);
        const double overlap_area = overlap_u * overlap_v;
        if (overlap_area <= 1.0e-16) {
            continue;
        }

        overlap_area_sum += overlap_area;
        overlap_center_u_sum += overlap_area * 0.5 * (min_u + max_u);
        overlap_center_v_sum += overlap_area * 0.5 * (min_v + max_v);
    }

    const double corrected_area = std::min(overlap_area_sum, sample.quadrature_area);
    if (corrected_area <= 1.0e-16) {
        return contribution;
    }

    contribution.corrected_area = corrected_area;
    const double inv_overlap = 1.0 / overlap_area_sum;
    const double overlap_center_u = overlap_center_u_sum * inv_overlap;
    const double overlap_center_v = overlap_center_v_sum * inv_overlap;
    contribution.corrected_point_world =
        patch_origin_world + tangent_u * overlap_center_u + tangent_v * overlap_center_v;
    return contribution;
}

}  // namespace

std::vector<ChSDFPatchBandAggregate> ChSDFPatchConsistencyBridge::BuildPatchBandAggregates(
    const ChSDFBrickPairWrenchResult& band_region,
    const ChSDFSheetRegion& sheet_region,
    const ChFrameMoving<>& shape_a_frame_abs,
    const ChFrameMoving<>& shape_b_frame_abs,
    const ChSDFPatchConsistencySettings& settings) {
    std::vector<ChSDFPatchBandAggregate> aggregates;
    if (!band_region.HasActiveContact()) {
        return aggregates;
    }

    std::vector<char> assigned(band_region.samples.size(), 0);
    aggregates.reserve(sheet_region.patches.size() + 1);

    for (const auto& patch : sheet_region.patches) {
        ChSDFPatchBandAggregate aggregate;
        aggregate.region_id = sheet_region.region_id;
        aggregate.patch_id = patch.patch_id;
        aggregate.source_sheet_sample_indices = patch.sample_indices;

        std::vector<char> patch_seen(band_region.samples.size(), 0);
        for (const auto sheet_sample_index : patch.sample_indices) {
            if (sheet_sample_index >= sheet_region.samples.size()) {
                continue;
            }

            const auto& sheet_sample = sheet_region.samples[sheet_sample_index];
            for (const auto band_sample_index : sheet_sample.source_sample_indices) {
                if (band_sample_index >= band_region.samples.size() || patch_seen[band_sample_index]) {
                    continue;
                }
                patch_seen[band_sample_index] = 1;
                assigned[band_sample_index] = 1;
                aggregate.source_band_sample_indices.push_back(band_sample_index);
                AccumulateSampleIntoPatch(aggregate, band_region.samples[band_sample_index], shape_a_frame_abs,
                                          shape_b_frame_abs);
            }
        }

        if (aggregate.band_area > settings.min_band_area) {
            aggregates.push_back(FinalizeAggregate(std::move(aggregate)));
        }
    }

    if (settings.add_unassigned_band_residual_patch) {
        ChSDFPatchBandAggregate residual;
        residual.region_id = band_region.region.region_id;
        residual.patch_id = aggregates.size();
        residual.residual_patch = true;

        for (std::size_t band_sample_index = 0; band_sample_index < band_region.samples.size(); ++band_sample_index) {
            if (assigned[band_sample_index]) {
                continue;
            }
            const auto& sample = band_region.samples[band_sample_index];
            if (!sample.active || sample.quadrature_area <= settings.min_band_area) {
                continue;
            }
            residual.source_band_sample_indices.push_back(band_sample_index);
            AccumulateSampleIntoPatch(residual, sample, shape_a_frame_abs, shape_b_frame_abs);
        }

        if (residual.band_area > settings.min_band_area) {
            aggregates.push_back(FinalizeAggregate(std::move(residual)));
        }
    }

    return aggregates;
}

ChSDFPatchConsistencyResult ChSDFPatchConsistencyBridge::BuildPatchConsistencyResult(
    const ChSDFPatchBandAggregate& band_patch,
    const ChSDFBrickPairWrenchResult& band_region,
    const ChSDFSheetRegion& sheet_region,
    const ChSDFSheetPatch& sheet_patch,
    const ChFrameMoving<>& shape_a_frame_abs,
    const ChFrameMoving<>& shape_b_frame_abs,
    const ChSDFPatchConsistencySettings& settings) {
    if (band_patch.band_area <= settings.min_band_area) {
        return MakePassthroughPatchResult(band_patch);
    }

    ChSDFPatchConsistencyResult result;
    result.region_id = band_patch.region_id;
    result.patch_id = band_patch.patch_id;
    result.residual_patch = band_patch.residual_patch;
    result.band_area = band_patch.band_area;
    result.sheet_area = sheet_patch.footprint_area > 0 ? sheet_patch.footprint_area : sheet_patch.measure_area;
    result.integrated_pressure_band = band_patch.integrated_pressure;
    result.wrench_world_a_band = band_patch.wrench_world_a;
    result.wrench_world_b_band = band_patch.wrench_world_b;
    result.centroid_world = band_patch.centroid_world;
    result.pressure_center_world = band_patch.pressure_center_world;
    result.sheet_pressure_center_world =
        IsFiniteVector(sheet_patch.pressure_center_world) ? sheet_patch.pressure_center_world : sheet_patch.centroid_world;
    result.corrected_pressure_center_world =
        IsFiniteVector(sheet_patch.centroid_world)
            ? sheet_patch.centroid_world
            : (IsFiniteVector(result.sheet_pressure_center_world) ? result.sheet_pressure_center_world
                                                                  : result.pressure_center_world);
    result.support_bbox_world = sheet_patch.support_bbox_world.IsInverted() ? band_patch.support_bbox_world
                                                                            : sheet_patch.support_bbox_world;

    double alpha = result.sheet_area > 0 ? (result.sheet_area / result.band_area) : 1.0;
    if (settings.clamp_alpha) {
        alpha = ClampValue(alpha, settings.min_alpha, settings.max_alpha);
    }

    result.alpha = alpha;
    result.integrated_pressure_corrected = alpha * result.integrated_pressure_band;

    if (settings.use_patch_local_redistribution) {
        const ChVector3d patch_normal = SafeNormalized(sheet_patch.mean_normal_world, VECT_Y);
        ChVector3d tangent_u;
        ChVector3d tangent_v;
        BuildOrthonormalBasis(patch_normal, tangent_u, tangent_v);

        std::vector<ProjectedRect> support_rects;
        support_rects.reserve(sheet_patch.sample_indices.size());
        for (const auto sample_index : sheet_patch.sample_indices) {
            if (sample_index >= sheet_region.samples.size()) {
                continue;
            }
            const auto rect =
                ProjectSheetSampleSupport(sheet_region.samples[sample_index], sheet_patch.centroid_world, tangent_u, tangent_v);
            if (rect.valid) {
                support_rects.push_back(rect);
            }
        }

        ChWrenchd redistributed_a = {VNULL, VNULL};
        ChWrenchd redistributed_b = {VNULL, VNULL};
        ChVector3d corrected_pressure_center_sum = VNULL;
        double corrected_pressure_weight_sum = 0;
        double raw_redistributed_area = 0;

        struct LocalContribution {
            std::size_t sample_index = 0;
            RedistributionContribution geom;
        };
        std::vector<LocalContribution> local;
        local.reserve(band_patch.source_band_sample_indices.size());

        for (const auto sample_index : band_patch.source_band_sample_indices) {
            if (sample_index >= band_region.samples.size()) {
                continue;
            }
            const auto& band_sample = band_region.samples[sample_index];
            double sample_half_extent = 0.5 * std::sqrt(std::max(band_sample.quadrature_area, 0.0));
            if (sample_half_extent <= 1.0e-12 && band_region.region.sample_spacing > 1.0e-12) {
                sample_half_extent = 0.5 * band_region.region.sample_spacing;
            }
            const auto geom =
                ComputeRedistributedContribution(band_sample, sheet_patch.centroid_world, support_rects, tangent_u,
                                                 tangent_v, sample_half_extent);
            if (geom.corrected_area <= 1.0e-16 || !IsFiniteVector(geom.corrected_point_world)) {
                continue;
            }
            raw_redistributed_area += geom.corrected_area;
            local.push_back(LocalContribution{sample_index, geom});
        }

        if (!local.empty() && raw_redistributed_area > 1.0e-16) {
            const double area_scale = result.sheet_area / raw_redistributed_area;
            for (const auto& item : local) {
                const auto& sample = band_region.samples[item.sample_index];
                const double corrected_area = item.geom.corrected_area * area_scale;
                const double scale = corrected_area / sample.quadrature_area;
                const ChVector3d corrected_force = scale * sample.force_world;
                redistributed_b.force += corrected_force;
                redistributed_b.torque +=
                    Vcross(item.geom.corrected_point_world - shape_b_frame_abs.GetPos(), corrected_force);
                redistributed_a.force -= corrected_force;
                redistributed_a.torque +=
                    Vcross(item.geom.corrected_point_world - shape_a_frame_abs.GetPos(), -corrected_force);

                const double pressure_weight = std::max(sample.pressure * corrected_area, 0.0);
                corrected_pressure_weight_sum += pressure_weight;
                corrected_pressure_center_sum += item.geom.corrected_point_world * pressure_weight;
            }

            result.wrench_world_a_corrected = redistributed_a;
            result.wrench_world_b_corrected = redistributed_b;
            result.integrated_pressure_corrected = alpha * result.integrated_pressure_band;
            result.corrected_pressure_center_world =
                corrected_pressure_weight_sum > 0 ? corrected_pressure_center_sum / corrected_pressure_weight_sum
                                                  : result.corrected_pressure_center_world;
            return result;
        }
    }

    if (settings.use_first_moment_consistent_correction) {
        result.wrench_world_a_corrected =
            BuildFirstMomentCorrectedWrench(result.wrench_world_a_band, shape_a_frame_abs.GetPos(),
                                            result.pressure_center_world, result.corrected_pressure_center_world, alpha);
        result.wrench_world_b_corrected =
            BuildFirstMomentCorrectedWrench(result.wrench_world_b_band, shape_b_frame_abs.GetPos(),
                                            result.pressure_center_world, result.corrected_pressure_center_world, alpha);
    } else {
        result.wrench_world_a_corrected.force = alpha * result.wrench_world_a_band.force;
        result.wrench_world_a_corrected.torque = alpha * result.wrench_world_a_band.torque;
        result.wrench_world_b_corrected.force = alpha * result.wrench_world_b_band.force;
        result.wrench_world_b_corrected.torque = alpha * result.wrench_world_b_band.torque;
    }
    return result;
}

ChSDFPatchConsistencyPairResult ChSDFPatchConsistencyBridge::BuildPairConsistencyResult(
    const ChSDFShapePairContactResult& band_result,
    const ChSDFSheetShapePairResult& sheet_result,
    const ChSDFPatchConsistencySettings& settings) {
    ChSDFPatchConsistencyPairResult result;
    result.wrench_world_a_band = band_result.wrench_world_a;
    result.wrench_world_b_band = band_result.wrench_world_b;

    if (!settings.enable || !band_result.valid) {
        result.wrench_world_a_corrected = result.wrench_world_a_band;
        result.wrench_world_b_corrected = result.wrench_world_b_band;
        result.total_band_area = band_result.active_area;
        result.total_sheet_area = band_result.active_area;
        result.total_corrected_area = band_result.active_area;
        result.integrated_pressure_band = band_result.integrated_pressure;
        result.integrated_pressure_corrected = band_result.integrated_pressure;
        return result;
    }

    std::unordered_map<std::size_t, const ChSDFSheetRegion*> sheet_region_by_id;
    sheet_region_by_id.reserve(sheet_result.regions.size());
    for (const auto& region : sheet_result.regions) {
        sheet_region_by_id.emplace(region.region_id, &region);
    }

    double alpha_weight_sum = 0;
    for (const auto& band_region : band_result.regions) {
        if (!band_region.HasActiveContact()) {
            continue;
        }

        result.total_band_area += band_region.active_area;
        result.integrated_pressure_band += band_region.integrated_pressure;

        const ChSDFSheetRegion* sheet_region = nullptr;
        const auto sheet_it = sheet_region_by_id.find(band_region.region.region_id);
        if (sheet_it != sheet_region_by_id.end()) {
            sheet_region = sheet_it->second;
        }

        if (!sheet_region) {
            ChSDFPatchBandAggregate passthrough;
            passthrough.region_id = band_region.region.region_id;
            passthrough.patch_id = 0;
            for (std::size_t i = 0; i < band_region.samples.size(); ++i) {
                if (!band_region.samples[i].active || band_region.samples[i].quadrature_area <= settings.min_band_area) {
                    continue;
                }
                passthrough.source_band_sample_indices.push_back(i);
                AccumulateSampleIntoPatch(passthrough, band_region.samples[i], band_result.shape_a_frame_abs,
                                          band_result.shape_b_frame_abs);
            }
            if (passthrough.band_area > settings.min_band_area) {
                auto patch_result = MakePassthroughPatchResult(FinalizeAggregate(std::move(passthrough)));
                result.total_sheet_area += patch_result.sheet_area;
                result.total_corrected_area += patch_result.alpha * patch_result.band_area;
                result.integrated_pressure_corrected += patch_result.integrated_pressure_corrected;
                alpha_weight_sum += patch_result.alpha * patch_result.band_area;
                result.wrench_world_a_corrected.force += patch_result.wrench_world_a_corrected.force;
                result.wrench_world_a_corrected.torque += patch_result.wrench_world_a_corrected.torque;
                result.wrench_world_b_corrected.force += patch_result.wrench_world_b_corrected.force;
                result.wrench_world_b_corrected.torque += patch_result.wrench_world_b_corrected.torque;
                result.patches.push_back(std::move(patch_result));
            }
            continue;
        }

        auto aggregates = BuildPatchBandAggregates(band_region, *sheet_region, band_result.shape_a_frame_abs,
                                                   band_result.shape_b_frame_abs, settings);
        std::unordered_map<std::size_t, const ChSDFSheetPatch*> patch_by_id;
        patch_by_id.reserve(sheet_region->patches.size());
        for (const auto& patch : sheet_region->patches) {
            patch_by_id.emplace(patch.patch_id, &patch);
        }

        for (const auto& aggregate : aggregates) {
            ChSDFPatchConsistencyResult patch_result;
            if (aggregate.residual_patch) {
                patch_result = MakePassthroughPatchResult(aggregate);
            } else {
                const auto patch_it = patch_by_id.find(aggregate.patch_id);
                if (patch_it != patch_by_id.end()) {
                    patch_result = BuildPatchConsistencyResult(aggregate, band_region, *sheet_region, *patch_it->second,
                                                               band_result.shape_a_frame_abs,
                                                               band_result.shape_b_frame_abs, settings);
                } else {
                    patch_result = MakePassthroughPatchResult(aggregate);
                }
            }

            result.total_sheet_area += patch_result.sheet_area;
            result.total_corrected_area += patch_result.alpha * patch_result.band_area;
            result.integrated_pressure_corrected += patch_result.integrated_pressure_corrected;
            alpha_weight_sum += patch_result.alpha * patch_result.band_area;
            result.wrench_world_a_corrected.force += patch_result.wrench_world_a_corrected.force;
            result.wrench_world_a_corrected.torque += patch_result.wrench_world_a_corrected.torque;
            result.wrench_world_b_corrected.force += patch_result.wrench_world_b_corrected.force;
            result.wrench_world_b_corrected.torque += patch_result.wrench_world_b_corrected.torque;
            result.patches.push_back(std::move(patch_result));
        }
    }

    result.corrected_patch_count = result.patches.size();
    if (result.total_band_area > 0) {
        result.mean_alpha = alpha_weight_sum / result.total_band_area;
    }
    if (result.corrected_patch_count == 0) {
        result.wrench_world_a_corrected = result.wrench_world_a_band;
        result.wrench_world_b_corrected = result.wrench_world_b_band;
        result.total_sheet_area = result.total_band_area;
        result.total_corrected_area = result.total_band_area;
        result.integrated_pressure_corrected = result.integrated_pressure_band;
        result.mean_alpha = 1.0;
    }

    return result;
}

}  // end namespace chrono
