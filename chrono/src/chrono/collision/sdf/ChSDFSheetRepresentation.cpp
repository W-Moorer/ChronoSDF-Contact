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

#include "chrono/collision/sdf/ChSDFSheetRepresentation.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <unordered_map>

namespace chrono {
namespace {

struct FiberKey {
    int a = 0;
    int b = 0;

    bool operator==(const FiberKey& other) const { return a == other.a && b == other.b; }
};

struct FiberKeyHash {
    std::size_t operator()(const FiberKey& key) const {
        const std::size_t ha = static_cast<std::size_t>(static_cast<uint32_t>(key.a));
        const std::size_t hb = static_cast<std::size_t>(static_cast<uint32_t>(key.b));
        return ha * 73856093u ^ hb * 19349663u;
    }
};

struct FiberAccumulator {
    int dominant_axis = -1;
    ChVector2i lattice_coord = ChVector2i(0, 0);

    double measure_area = 0;
    double pressure_weight = 0;

    ChVector3d area_centroid_sum = VNULL;
    ChVector3d pressure_center_sum = VNULL;
    ChVector3d normal_sum = VNULL;
    ChVector3d force_sum = VNULL;

    ChAABB bounds_world;
    std::vector<std::size_t> source_sample_indices;
};

ChVector3d SafeNormalized(const ChVector3d& v, const ChVector3d& fallback = VECT_Y) {
    const double length = v.Length();
    return length > 1.0e-12 ? v / length : fallback;
}

ChVector3d ResolveRegionNormal(const ChSDFBrickPairWrenchResult& band_region) {
    ChVector3d normal_sum = band_region.region.mean_normal_world;
    if (normal_sum.Length2() > 1.0e-24) {
        return SafeNormalized(normal_sum);
    }

    normal_sum = VNULL;
    for (const auto& sample : band_region.samples) {
        if (!sample.active || sample.quadrature_area <= 1.0e-16) {
            continue;
        }
        normal_sum += sample.region_sample.contact_normal_world * sample.quadrature_area;
    }

    return SafeNormalized(normal_sum);
}

double ResolveRegionSpacing(const ChSDFBrickPairWrenchResult& band_region) {
    if (band_region.region.sample_spacing > 0) {
        return band_region.region.sample_spacing;
    }

    for (const auto& sample : band_region.samples) {
        if (sample.resolution_length > 0) {
            return sample.resolution_length;
        }
    }

    return 0.0;
}

std::vector<ChVector2i> BuildPatchNeighborOffsets(int neighbor_mode) {
    std::vector<ChVector2i> offsets;
    offsets.reserve(8);

    for (int dv = -1; dv <= 1; ++dv) {
        for (int du = -1; du <= 1; ++du) {
            if (du == 0 && dv == 0) {
                continue;
            }

            const int l1 = std::abs(du) + std::abs(dv);
            const bool keep = (neighbor_mode == 4) ? (l1 == 1) : true;
            if (keep) {
                offsets.emplace_back(du, dv);
            }
        }
    }

    return offsets;
}

void FinalizePatch(ChSDFSheetPatch& patch,
                   const std::vector<ChSDFSheetFiberSample>& samples,
                   const ChVector3d& fallback_normal) {
    patch.measure_area = 0;
    patch.footprint_area = 0;
    patch.centroid_world = VNULL;
    patch.mean_normal_world = VNULL;
    patch.bounds_world = ChAABB();

    if (patch.sample_indices.empty()) {
        return;
    }

    ChVector3d centroid_sum = VNULL;
    ChVector3d normal_sum = VNULL;

    for (const std::size_t sample_index : patch.sample_indices) {
        const auto& sample = samples[sample_index];
        patch.measure_area += sample.measure_area;
        patch.footprint_area += sample.footprint_area;
        centroid_sum += sample.centroid_world * sample.measure_area;
        normal_sum += sample.normal_world * sample.measure_area;
        patch.bounds_world += sample.source_bounds_world;
    }

    if (patch.measure_area > 0) {
        patch.centroid_world = centroid_sum / patch.measure_area;
    }
    patch.mean_normal_world = SafeNormalized(normal_sum, fallback_normal);
}

std::vector<ChSDFSheetPatch> BuildPatches(const std::vector<ChSDFSheetFiberSample>& samples,
                                          int neighbor_mode,
                                          const ChVector3d& fallback_normal) {
    std::vector<ChSDFSheetPatch> patches;
    if (samples.empty()) {
        return patches;
    }

    std::unordered_map<FiberKey, std::size_t, FiberKeyHash> sample_by_coord;
    sample_by_coord.reserve(samples.size());
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto& sample = samples[i];
        sample_by_coord.emplace(FiberKey{sample.lattice_coord.x(), sample.lattice_coord.y()}, i);
    }

    const auto neighbor_offsets = BuildPatchNeighborOffsets(neighbor_mode);
    std::vector<char> visited(samples.size(), 0);

    for (std::size_t seed = 0; seed < samples.size(); ++seed) {
        if (visited[seed]) {
            continue;
        }

        ChSDFSheetPatch patch;
        patch.patch_id = patches.size();

        std::queue<std::size_t> frontier;
        frontier.push(seed);
        visited[seed] = 1;

        while (!frontier.empty()) {
            const std::size_t current = frontier.front();
            frontier.pop();
            patch.sample_indices.push_back(current);

            const auto& coord = samples[current].lattice_coord;
            for (const auto& offset : neighbor_offsets) {
                FiberKey neighbor{coord.x() + offset.x(), coord.y() + offset.y()};
                auto it = sample_by_coord.find(neighbor);
                if (it == sample_by_coord.end() || visited[it->second]) {
                    continue;
                }

                visited[it->second] = 1;
                frontier.push(it->second);
            }
        }

        FinalizePatch(patch, samples, fallback_normal);
        patches.push_back(std::move(patch));
    }

    return patches;
}

}  // namespace

int ChSDFSheetBuilder::ComputeDominantAxis(const ChVector3d& normal_world) {
    const ChVector3d abs_normal(std::abs(normal_world.x()), std::abs(normal_world.y()), std::abs(normal_world.z()));
    if (abs_normal.x() >= abs_normal.y() && abs_normal.x() >= abs_normal.z()) {
        return 0;
    }
    if (abs_normal.y() >= abs_normal.x() && abs_normal.y() >= abs_normal.z()) {
        return 1;
    }
    return 2;
}

ChVector2i ChSDFSheetBuilder::ComputeLatticeCoord(const ChVector3i& coord, int dominant_axis) {
    switch (dominant_axis) {
        case 0:
            return ChVector2i(coord.y(), coord.z());
        case 1:
            return ChVector2i(coord.x(), coord.z());
        case 2:
        default:
            return ChVector2i(coord.x(), coord.y());
    }
}

ChSDFSheetRegion ChSDFSheetBuilder::BuildRegion(const ChSDFBrickPairWrenchResult& band_region,
                                                const ChSDFSheetCollapseSettings& settings) {
    ChSDFSheetRegion region;
    region.region_id = band_region.region.region_id;
    region.persistent_id = band_region.region.persistent_id;

    if (!settings.enable || !band_region.HasActiveContact() || band_region.samples.empty()) {
        return region;
    }

    const ChVector3d fallback_normal = ResolveRegionNormal(band_region);
    region.dominant_axis = ComputeDominantAxis(fallback_normal);

    std::unordered_map<FiberKey, FiberAccumulator, FiberKeyHash> fibers;
    fibers.reserve(band_region.samples.size());

    for (std::size_t sample_index = 0; sample_index < band_region.samples.size(); ++sample_index) {
        const auto& sample = band_region.samples[sample_index];
        if (!sample.active || sample.quadrature_area <= 1.0e-16) {
            continue;
        }

        const ChVector2i lattice_coord = ComputeLatticeCoord(sample.region_sample.coord, region.dominant_axis);
        const FiberKey key{lattice_coord.x(), lattice_coord.y()};
        auto& fiber = fibers[key];
        fiber.dominant_axis = region.dominant_axis;
        fiber.lattice_coord = lattice_coord;
        fiber.measure_area += sample.quadrature_area;
        fiber.area_centroid_sum += sample.region_sample.point_world * sample.quadrature_area;
        fiber.normal_sum += sample.region_sample.contact_normal_world * sample.quadrature_area;
        fiber.force_sum += sample.force_world;
        fiber.bounds_world += sample.region_sample.point_world;
        fiber.source_sample_indices.push_back(sample_index);

        const double pressure_weight = std::max(sample.pressure * sample.quadrature_area, 0.0);
        fiber.pressure_weight += pressure_weight;
        fiber.pressure_center_sum += sample.region_sample.point_world * pressure_weight;
    }

    if (fibers.empty()) {
        return region;
    }

    const double spacing = ResolveRegionSpacing(band_region);
    const double footprint_area = spacing > 0 ? spacing * spacing : 0.0;

    region.samples.reserve(fibers.size());

    ChVector3d centroid_sum = VNULL;
    ChVector3d pressure_center_sum = VNULL;
    ChVector3d normal_sum = VNULL;
    double pressure_weight_sum = 0;

    std::vector<std::pair<FiberKey, FiberAccumulator>> ordered_fibers;
    ordered_fibers.reserve(fibers.size());
    for (auto& item : fibers) {
        ordered_fibers.push_back({item.first, std::move(item.second)});
    }
    std::stable_sort(ordered_fibers.begin(), ordered_fibers.end(),
                     [](const auto& a, const auto& b) { return a.first.a != b.first.a ? a.first.a < b.first.a : a.first.b < b.first.b; });

    for (const auto& item : ordered_fibers) {
        const auto& fiber_data = item.second;
        if (fiber_data.measure_area <= std::max(settings.min_sheet_sample_area, 0.0)) {
            continue;
        }

        ChSDFSheetFiberSample sample;
        sample.region_id = region.region_id;
        sample.fiber_id = region.samples.size();
        sample.dominant_axis = region.dominant_axis;
        sample.lattice_coord = fiber_data.lattice_coord;
        sample.measure_area = fiber_data.measure_area;
        sample.footprint_area = footprint_area;
        sample.centroid_world = fiber_data.area_centroid_sum / fiber_data.measure_area;
        sample.pressure_center_world =
            fiber_data.pressure_weight > 0 ? fiber_data.pressure_center_sum / fiber_data.pressure_weight : sample.centroid_world;
        sample.normal_world = SafeNormalized(fiber_data.normal_sum, fallback_normal);
        sample.force_world = fiber_data.force_sum;
        sample.source_bounds_world = fiber_data.bounds_world;
        sample.source_sample_indices = fiber_data.source_sample_indices;

        region.measure_area += sample.measure_area;
        region.footprint_area += sample.footprint_area;
        centroid_sum += sample.centroid_world * sample.measure_area;
        normal_sum += sample.normal_world * sample.measure_area;
        region.bounds_world += sample.source_bounds_world;

        const double region_pressure_weight = sample.force_world.Length();
        pressure_weight_sum += region_pressure_weight;
        pressure_center_sum += sample.pressure_center_world * region_pressure_weight;

        region.samples.push_back(std::move(sample));
    }

    if (region.samples.empty()) {
        return region;
    }

    region.centroid_world = centroid_sum / region.measure_area;
    region.mean_normal_world = SafeNormalized(normal_sum, fallback_normal);
    region.pressure_center_world =
        pressure_weight_sum > 0 ? pressure_center_sum / pressure_weight_sum : region.centroid_world;
    region.patches = BuildPatches(region.samples, settings.patch_neighbor_mode, region.mean_normal_world);

    return region;
}

std::vector<ChSDFSheetRegion> ChSDFSheetBuilder::BuildRegions(
    const std::vector<ChSDFBrickPairWrenchResult>& band_regions,
    const ChSDFSheetCollapseSettings& settings) {
    std::vector<ChSDFSheetRegion> regions;
    regions.reserve(band_regions.size());

    for (const auto& band_region : band_regions) {
        auto region = BuildRegion(band_region, settings);
        if (region.HasSamples()) {
            regions.push_back(std::move(region));
        }
    }

    return regions;
}

ChSDFSheetShapePairResult ChSDFSheetBuilder::BuildShapePair(const ChSDFShapePairContactResult& band_result,
                                                            const ChSDFSheetCollapseSettings& settings) {
    ChSDFSheetShapePairResult result;
    result.band_area = band_result.active_area;

    if (!settings.enable || !band_result.HasActiveContact()) {
        return result;
    }

    ChVector3d sheet_center_sum = VNULL;
    ChVector3d pressure_center_sum = VNULL;
    double pressure_weight_sum = 0;

    result.regions.reserve(band_result.regions.size());

    for (const auto& band_region : band_result.regions) {
        if (!band_region.HasActiveContact()) {
            continue;
        }

        const double spacing = ResolveRegionSpacing(band_region);
        if (spacing > 0) {
            result.occupied_area += static_cast<double>(band_region.active_samples) * spacing * spacing;
        }
        result.occupied_cells += band_region.active_samples;

        auto sheet_region = BuildRegion(band_region, settings);
        if (!sheet_region.HasSamples()) {
            continue;
        }

        result.sheet_area += sheet_region.measure_area;
        result.sheet_footprint_area += sheet_region.footprint_area;
        result.collapsed_samples += sheet_region.samples.size();
        sheet_center_sum += sheet_region.centroid_world * sheet_region.measure_area;

        const double region_pressure_weight = band_region.integrated_pressure;
        pressure_weight_sum += region_pressure_weight;
        pressure_center_sum += sheet_region.pressure_center_world * region_pressure_weight;

        result.regions.push_back(std::move(sheet_region));
    }

    if (result.sheet_area > 0) {
        result.sheet_center_world = sheet_center_sum / result.sheet_area;
    }
    result.pressure_center_world =
        pressure_weight_sum > 0 ? pressure_center_sum / pressure_weight_sum : result.sheet_center_world;

    return result;
}

}  // end namespace chrono
