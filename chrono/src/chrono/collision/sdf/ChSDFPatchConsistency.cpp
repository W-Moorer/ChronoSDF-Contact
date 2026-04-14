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

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace chrono {
namespace {

bool IsFiniteVector(const ChVector3d& v) {
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
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

double SignedPolygonArea(const std::vector<ChVector2d>& polygon) {
    if (polygon.size() < 3) {
        return 0;
    }
    double twice_area = 0;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        twice_area += a.x() * b.y() - a.y() * b.x();
    }
    return 0.5 * twice_area;
}

double PolygonArea(const std::vector<ChVector2d>& polygon) {
    return std::abs(SignedPolygonArea(polygon));
}

ChVector2d PolygonCentroid(const std::vector<ChVector2d>& polygon) {
    if (polygon.empty()) {
        return ChVector2d(0, 0);
    }
    if (polygon.size() < 3) {
        ChVector2d centroid(0, 0);
        for (const auto& point : polygon) {
            centroid += point;
        }
        return centroid / static_cast<double>(polygon.size());
    }

    const double signed_area = SignedPolygonArea(polygon);
    if (std::abs(signed_area) <= 1.0e-16) {
        ChVector2d centroid(0, 0);
        for (const auto& point : polygon) {
            centroid += point;
        }
        return centroid / static_cast<double>(polygon.size());
    }

    ChVector2d centroid(0, 0);
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        const double cross = a.x() * b.y() - a.y() * b.x();
        centroid += (a + b) * cross;
    }
    return centroid / (6.0 * signed_area);
}

ChVector2d IntersectSegmentWithVertical(const ChVector2d& a, const ChVector2d& b, double x_clip) {
    const double dx = b.x() - a.x();
    if (std::abs(dx) <= 1.0e-16) {
        return ChVector2d(x_clip, a.y());
    }
    const double t = (x_clip - a.x()) / dx;
    return ChVector2d(x_clip, a.y() + t * (b.y() - a.y()));
}

ChVector2d IntersectSegmentWithHorizontal(const ChVector2d& a, const ChVector2d& b, double y_clip) {
    const double dy = b.y() - a.y();
    if (std::abs(dy) <= 1.0e-16) {
        return ChVector2d(a.x(), y_clip);
    }
    const double t = (y_clip - a.y()) / dy;
    return ChVector2d(a.x() + t * (b.x() - a.x()), y_clip);
}

std::vector<ChVector2d> ClipPolygonAgainstLeft(const std::vector<ChVector2d>& polygon, double min_u) {
    std::vector<ChVector2d> clipped;
    if (polygon.empty()) {
        return clipped;
    }
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto& current = polygon[i];
        const auto& prev = polygon[(i + polygon.size() - 1) % polygon.size()];
        const bool inside_current = current.x() >= min_u - 1.0e-12;
        const bool inside_prev = prev.x() >= min_u - 1.0e-12;
        if (inside_current) {
            if (!inside_prev) {
                clipped.push_back(IntersectSegmentWithVertical(prev, current, min_u));
            }
            clipped.push_back(current);
        } else if (inside_prev) {
            clipped.push_back(IntersectSegmentWithVertical(prev, current, min_u));
        }
    }
    return clipped;
}

std::vector<ChVector2d> ClipPolygonAgainstRight(const std::vector<ChVector2d>& polygon, double max_u) {
    std::vector<ChVector2d> clipped;
    if (polygon.empty()) {
        return clipped;
    }
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto& current = polygon[i];
        const auto& prev = polygon[(i + polygon.size() - 1) % polygon.size()];
        const bool inside_current = current.x() <= max_u + 1.0e-12;
        const bool inside_prev = prev.x() <= max_u + 1.0e-12;
        if (inside_current) {
            if (!inside_prev) {
                clipped.push_back(IntersectSegmentWithVertical(prev, current, max_u));
            }
            clipped.push_back(current);
        } else if (inside_prev) {
            clipped.push_back(IntersectSegmentWithVertical(prev, current, max_u));
        }
    }
    return clipped;
}

std::vector<ChVector2d> ClipPolygonAgainstBottom(const std::vector<ChVector2d>& polygon, double min_v) {
    std::vector<ChVector2d> clipped;
    if (polygon.empty()) {
        return clipped;
    }
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto& current = polygon[i];
        const auto& prev = polygon[(i + polygon.size() - 1) % polygon.size()];
        const bool inside_current = current.y() >= min_v - 1.0e-12;
        const bool inside_prev = prev.y() >= min_v - 1.0e-12;
        if (inside_current) {
            if (!inside_prev) {
                clipped.push_back(IntersectSegmentWithHorizontal(prev, current, min_v));
            }
            clipped.push_back(current);
        } else if (inside_prev) {
            clipped.push_back(IntersectSegmentWithHorizontal(prev, current, min_v));
        }
    }
    return clipped;
}

std::vector<ChVector2d> ClipPolygonAgainstTop(const std::vector<ChVector2d>& polygon, double max_v) {
    std::vector<ChVector2d> clipped;
    if (polygon.empty()) {
        return clipped;
    }
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto& current = polygon[i];
        const auto& prev = polygon[(i + polygon.size() - 1) % polygon.size()];
        const bool inside_current = current.y() <= max_v + 1.0e-12;
        const bool inside_prev = prev.y() <= max_v + 1.0e-12;
        if (inside_current) {
            if (!inside_prev) {
                clipped.push_back(IntersectSegmentWithHorizontal(prev, current, max_v));
            }
            clipped.push_back(current);
        } else if (inside_prev) {
            clipped.push_back(IntersectSegmentWithHorizontal(prev, current, max_v));
        }
    }
    return clipped;
}

std::vector<ChVector2d> ClipPolygonToRect(const std::vector<ChVector2d>& polygon,
                                          double min_u,
                                          double max_u,
                                          double min_v,
                                          double max_v) {
    auto clipped = ClipPolygonAgainstLeft(polygon, min_u);
    clipped = ClipPolygonAgainstRight(clipped, max_u);
    clipped = ClipPolygonAgainstBottom(clipped, min_v);
    clipped = ClipPolygonAgainstTop(clipped, max_v);
    return clipped;
}

struct RedistributionContribution {
    double corrected_area = 0;
    ChVector3d corrected_point_world = VNULL;
};

RedistributionContribution ComputeRedistributedContribution(const ChSDFBrickPairWrenchSample& sample,
                                                            const ChSDFSheetLocalFootprint& footprint,
                                                            double sample_half_extent) {
    RedistributionContribution contribution;
    if (!sample.active || sample.quadrature_area <= 1.0e-16 || sample_half_extent <= 1.0e-12 || !footprint.HasPolygon()) {
        return contribution;
    }

    const ChVector3d point_world = ResolveSamplePointWorld(sample);
    const ChVector3d rel = point_world - footprint.origin_world;
    const double center_u = Vdot(rel, footprint.tangent_u_world);
    const double center_v = Vdot(rel, footprint.tangent_v_world);
    const double half = sample_half_extent;

    const auto clipped = ClipPolygonToRect(footprint.polygon_uv, center_u - half, center_u + half, center_v - half,
                                           center_v + half);
    const double corrected_area = std::min(PolygonArea(clipped), sample.quadrature_area);
    if (corrected_area <= 1.0e-16) {
        return contribution;
    }

    contribution.corrected_area = corrected_area;
    const ChVector2d overlap_center = PolygonCentroid(clipped);
    contribution.corrected_point_world =
        footprint.origin_world + footprint.tangent_u_world * overlap_center.x() +
        footprint.tangent_v_world * overlap_center.y();
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

    if (settings.use_patch_local_redistribution && sheet_patch.support_footprint.HasPolygon()) {
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
                ComputeRedistributedContribution(band_sample, sheet_patch.support_footprint, sample_half_extent);
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
                    patch_result = BuildPatchConsistencyResult(aggregate, band_region, *patch_it->second,
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
