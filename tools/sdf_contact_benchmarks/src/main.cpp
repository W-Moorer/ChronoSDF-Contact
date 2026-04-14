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
// First benchmark package for the dual-SDF distributed contact model.
//
// The current version is intentionally narrow:
// - quasi-static centered compression
// - quasi-static tilted/eccentric compression
// - dense plane-integration reference
// - single-point and coarse grid penalty baselines
// =============================================================================

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <set>
#include <unordered_map>
#include <vector>

#include <chrono/collision/ChCollisionShapeSDF.h>
#include <chrono/collision/sdf/ChSDFPotentialCalibration.h>
#include <chrono/collision/sdf/ChSDFShapePair.h>
#include <chrono/core/ChRotation.h>
#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChContactMaterialNSC.h>
#include <chrono/physics/ChSystemNSC.h>

using namespace chrono;
namespace fs = std::filesystem;

namespace {

constexpr double kHalfExtent = 0.5;
constexpr double kPlaneY = 0.5;
constexpr double kDegToRad = CH_PI / 180.0;

struct PoseSample {
    std::size_t sample_index = 0;
    double penetration = 0;
    double tilt_z_rad = 0;
    double offset_x = 0;
};

struct ScenarioSpec {
    std::string name;
    std::vector<PoseSample> poses;
};

struct MethodSpec {
    std::string name;
    std::size_t samples_u = 0;
    std::size_t samples_v = 0;
    bool is_reference = false;
    bool is_distributed = false;
};

struct BenchmarkRecord {
    std::string scenario;
    std::string method;
    std::size_t sample_index = 0;
    double penetration = 0;
    double tilt_z_deg = 0;
    double offset_x = 0;

    bool valid = false;

    double force_x = 0;
    double force_y = 0;
    double force_z = 0;
    double force_y_corrected = 0;
    double torque_x = 0;
    double torque_y = 0;
    double torque_z = 0;
    double torque_z_corrected = 0;

    double active_area = 0;
    double area_occ = 0;
    double area_band = 0;
    double area_sheet = 0;
    double area_corrected = 0;
    double max_penetration = 0;
    double pressure_center_x = 0;
    double sheet_center_x = 0;
    double sheet_pressure_center_x = 0;
    double sheet_area_ratio = 0;
    double sheet_mean_support_seed_count = 0;
    double sheet_normal_spread = 0;
    double mean_patch_alpha = 1;
    double sheet_support_bbox_xmin = 0;
    double sheet_support_bbox_xmax = 0;
    double sheet_support_bbox_zmin = 0;
    double sheet_support_bbox_zmax = 0;
    double eval_ms = 0;

    std::size_t active_regions = 0;
    std::size_t active_samples = 0;
    std::size_t sheet_patch_count = 0;
    std::size_t sheet_fiber_count = 0;
    std::size_t sheet_fallback_regions = 0;
    bool sheet_used_fallback = false;
};

struct BenchmarkSummary {
    std::string scenario;
    std::string method;
    std::size_t sample_count = 0;
    double mean_eval_ms = 0;
    double contact_ratio = 0;
    double force_y_rmse = 0;
    double force_y_corrected_rmse = 0;
    double force_y_max_abs = 0;
    double torque_z_rmse = 0;
    double torque_z_corrected_rmse = 0;
    double torque_z_max_abs = 0;
    double active_area_rmse = 0;
    double area_occ_rmse = 0;
    double area_band_rmse = 0;
    double area_sheet_rmse = 0;
    double area_corrected_rmse = 0;
    double mean_patch_alpha = 0;
    double pressure_center_x_rmse = 0;
};

struct CurvedSheetReference {
    double sheet_area = 0;
    double sheet_center_x = 0;
    double sheet_support_bbox_xmin = 0;
    double sheet_support_bbox_xmax = 0;
    double sheet_support_bbox_zmin = 0;
    double sheet_support_bbox_zmax = 0;
    std::size_t sheet_patch_count = 0;
};

struct CurvedSheetRecord {
    std::string scenario;
    std::string mode;
    std::size_t sample_index = 0;
    double penetration = 0;
    bool valid = false;
    double area_occ = 0;
    double area_band = 0;
    double area_sheet = 0;
    double sheet_center_x = 0;
    double sheet_pressure_center_x = 0;
    double sheet_area_ratio = 0;
    double sheet_mean_support_seed_count = 0;
    double sheet_normal_spread = 0;
    double sheet_support_bbox_xmin = 0;
    double sheet_support_bbox_xmax = 0;
    double sheet_support_bbox_zmin = 0;
    double sheet_support_bbox_zmax = 0;
    double eval_ms = 0;
    std::size_t sheet_patch_count = 0;
    std::size_t sheet_fiber_count = 0;
    std::size_t sheet_fallback_regions = 0;
    bool sheet_used_fallback = false;
};

struct CurvedSheetSummary {
    std::string scenario;
    std::string mode;
    std::size_t sample_count = 0;
    double mean_eval_ms = 0;
    double area_sheet_rmse = 0;
    double sheet_center_x_rmse = 0;
    double bbox_xmin_rmse = 0;
    double bbox_xmax_rmse = 0;
    double bbox_zmin_rmse = 0;
    double bbox_zmax_rmse = 0;
    double patch_count_rmse = 0;
    double mean_area_ratio = 0;
    double mean_support_seed_count = 0;
    double mean_normal_spread = 0;
    double mean_fiber_count = 0;
    double fallback_ratio = 0;
};

struct PatchAreaDiagnosticRecord {
    std::string scenario;
    std::string method;
    std::size_t sample_index = 0;
    double penetration = 0;
    double tilt_z_deg = 0;
    double offset_x = 0;
    std::size_t region_id = 0;
    std::size_t patch_id = 0;
    std::size_t support_columns = 0;
    std::size_t support_seed_count = 0;
    std::size_t layered_cell_count = 0;
    std::size_t support_cell_count = 0;
    std::size_t shell_cell_count = 0;
    std::size_t max_layer_count_per_cell = 0;
    std::size_t largest_connected_sheet_cells = 0;
    std::size_t support_cell_axis_count = 0;
    std::size_t support_cell_bbox_count = 0;
    int support_cell_u_span = 0;
    int support_cell_v_span = 0;
    double mean_layer_count_per_cell = 0;
    std::size_t polygon_vertex_count = 0;
    double measure_area = 0;
    double footprint_area = 0;
    double footprint_polygon_area = 0;
    double support_bbox_projected_area = 0;
    double support_cell_fill_ratio = 0;
    double polygon_to_footprint_ratio = 0;
    double bbox_to_footprint_ratio = 0;
    double bbox_to_polygon_ratio = 0;
};

struct DistributedEvaluation {
    BenchmarkRecord record;
    std::vector<PatchAreaDiagnosticRecord> patch_area_records;
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

double Degrees(double radians) {
    return radians / kDegToRad;
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

std::string KeyFor(const std::string& scenario, std::size_t sample_index) {
    return scenario + "#" + std::to_string(sample_index);
}

fs::path DefaultNvdbPath() {
    return fs::path(BENCHMARK_TOOL_DIR).parent_path().parent_path() / "chrono" / "template_project" / "assets" /
           "unit_box_centered.nvdb";
}

fs::path DefaultCylinderNvdbPath() {
    return fs::path(BENCHMARK_TOOL_DIR).parent_path().parent_path() / "chrono" / "template_project" / "assets" /
           "unit_cylinder_z.nvdb";
}

fs::path DefaultOutputDir() {
    return fs::path(BENCHMARK_TOOL_DIR).parent_path().parent_path() / "paper" / "results" / "benchmark_box_contact";
}

void PrintUsage(const char* exe_name) {
    std::cout << "Usage:\n"
              << "  " << exe_name
              << " [--input path/to/unit_box_centered.nvdb] [--curved-input path/to/unit_cylinder_z.nvdb]"
              << " [--output output_dir]\n\n"
              << "Defaults:\n"
              << "  input:  " << DefaultNvdbPath().string() << "\n"
              << "  curved: " << DefaultCylinderNvdbPath().string() << "\n"
              << "  output: " << DefaultOutputDir().string() << "\n";
}

std::array<ScenarioSpec, 2> BuildScenarios() {
    ScenarioSpec centered;
    centered.name = "centered_compression";
    for (std::size_t i = 0; i < 8; ++i) {
        centered.poses.push_back(PoseSample{i, 0.010 + 0.005 * static_cast<double>(i), 0.0, 0.0});
    }

    ScenarioSpec eccentric;
    eccentric.name = "eccentric_compression";
    for (std::size_t i = 0; i <= 10; ++i) {
        eccentric.poses.push_back(PoseSample{i, 0.030, static_cast<double>(i) * 1.0 * kDegToRad, 0.0});
    }

    return {centered, eccentric};
}

ScenarioSpec BuildCurvedSheetScenario() {
    ScenarioSpec curved;
    curved.name = "cylinder_strip_centered";
    for (std::size_t i = 0; i < 8; ++i) {
        curved.poses.push_back(PoseSample{i, 0.010 + 0.010 * static_cast<double>(i), 0.0, 0.0});
    }
    return curved;
}

std::vector<MethodSpec> BuildMethods() {
    return {
        {"distributed_sdf", 0, 0, false, true},
        {"single_point_penalty", 1, 1, false, false},
        {"grid_penalty_5x5", 5, 5, false, false},
        {"dense_reference_161x161", 161, 161, true, false},
    };
}

double ComputePressureCenterX(const ChSDFShapePairContactResult& result) {
    double weighted_x = 0;
    double weight_sum = 0;

    for (const auto& region : result.regions) {
        if (!region.HasActiveContact() || region.integrated_pressure <= 0) {
            continue;
        }

        weighted_x += region.integrated_pressure * region.pressure_center_world.x();
        weight_sum += region.integrated_pressure;
    }

    return weight_sum > 0 ? weighted_x / weight_sum : 0.0;
}

double ProjectedSupportBBoxArea(const ChSDFSheetPatch& patch) {
    if (patch.support_bbox_world.IsInverted()) {
        return 0.0;
    }

    ChVector3d tangent_u = patch.support_footprint.tangent_u_world;
    ChVector3d tangent_v = patch.support_footprint.tangent_v_world;
    if (tangent_u.Length2() <= 1.0e-24 || tangent_v.Length2() <= 1.0e-24) {
        BuildOrthonormalBasis(patch.mean_normal_world, tangent_u, tangent_v);
    }

    const ChVector3d origin_world =
        (patch.support_footprint.origin_world.Length2() > 0 || patch.centroid_world.Length2() == 0)
            ? patch.support_footprint.origin_world
            : patch.centroid_world;

    const std::array<ChVector3d, 8> corners = {
        ChVector3d(patch.support_bbox_world.min.x(), patch.support_bbox_world.min.y(), patch.support_bbox_world.min.z()),
        ChVector3d(patch.support_bbox_world.min.x(), patch.support_bbox_world.min.y(), patch.support_bbox_world.max.z()),
        ChVector3d(patch.support_bbox_world.min.x(), patch.support_bbox_world.max.y(), patch.support_bbox_world.min.z()),
        ChVector3d(patch.support_bbox_world.min.x(), patch.support_bbox_world.max.y(), patch.support_bbox_world.max.z()),
        ChVector3d(patch.support_bbox_world.max.x(), patch.support_bbox_world.min.y(), patch.support_bbox_world.min.z()),
        ChVector3d(patch.support_bbox_world.max.x(), patch.support_bbox_world.min.y(), patch.support_bbox_world.max.z()),
        ChVector3d(patch.support_bbox_world.max.x(), patch.support_bbox_world.max.y(), patch.support_bbox_world.min.z()),
        ChVector3d(patch.support_bbox_world.max.x(), patch.support_bbox_world.max.y(), patch.support_bbox_world.max.z())};

    double min_u = std::numeric_limits<double>::infinity();
    double max_u = -std::numeric_limits<double>::infinity();
    double min_v = std::numeric_limits<double>::infinity();
    double max_v = -std::numeric_limits<double>::infinity();

    for (const auto& corner : corners) {
        const ChVector3d rel = corner - origin_world;
        const double u = Vdot(rel, tangent_u);
        const double v = Vdot(rel, tangent_v);
        min_u = std::min(min_u, u);
        max_u = std::max(max_u, u);
        min_v = std::min(min_v, v);
        max_v = std::max(max_v, v);
    }

    if (!std::isfinite(min_u) || !std::isfinite(max_u) || !std::isfinite(min_v) || !std::isfinite(max_v)) {
        return 0.0;
    }

    return std::max(max_u - min_u, 0.0) * std::max(max_v - min_v, 0.0);
}

CurvedSheetReference MakeCylinderSheetReference(double penetration) {
    constexpr double radius = 0.5;
    constexpr double length = 1.0;

    CurvedSheetReference ref;
    const double half_width = std::sqrt(std::max(2.0 * radius * penetration - penetration * penetration, 0.0));
    ref.sheet_area = 2.0 * half_width * length;
    ref.sheet_center_x = 0.0;
    ref.sheet_support_bbox_xmin = -half_width;
    ref.sheet_support_bbox_xmax = half_width;
    ref.sheet_support_bbox_zmin = -0.5 * length;
    ref.sheet_support_bbox_zmax = 0.5 * length;
    ref.sheet_patch_count = ref.sheet_area > 0 ? 1 : 0;
    return ref;
}

DistributedEvaluation EvaluateDistributed(const std::string& scenario_name,
                                          const PoseSample& pose,
                                          ChBody* moving_body,
                                          ChSDFShapePair& pair,
                                          const ChSDFBrickPairBroadphase::Settings& pair_settings,
                                          const ChSDFContactRegionBuilder::Settings& region_settings,
                                          const ChSDFNormalPressureSettings& pressure_settings) {
    moving_body->SetPos(ChVector3d(pose.offset_x, 1.0 - pose.penetration, 0.0));
    moving_body->SetRot(QuatFromAngleZ(pose.tilt_z_rad));
    moving_body->SetPosDt(VNULL);
    moving_body->SetAngVelLocal(VNULL);

    const auto start = std::chrono::steady_clock::now();
    const auto result = pair.EvaluateContact(pair_settings, region_settings, pressure_settings);
    const auto stop = std::chrono::steady_clock::now();

    DistributedEvaluation evaluation;
    auto& record = evaluation.record;
    record.scenario = scenario_name;
    record.method = "distributed_sdf";
    record.sample_index = pose.sample_index;
    record.penetration = pose.penetration;
    record.tilt_z_deg = Degrees(pose.tilt_z_rad);
    record.offset_x = pose.offset_x;
    record.valid = result.valid;
    record.force_x = result.wrench_world_b.force.x();
    record.force_y = result.wrench_world_b.force.y();
    record.force_z = result.wrench_world_b.force.z();
    record.force_y_corrected = result.wrench_world_b_corrected.force.y();
    record.torque_x = result.wrench_world_b.torque.x();
    record.torque_y = result.wrench_world_b.torque.y();
    record.torque_z = result.wrench_world_b.torque.z();
    record.torque_z_corrected = result.wrench_world_b_corrected.torque.z();
    record.active_area = result.active_area;
    record.max_penetration = result.max_penetration;
    record.pressure_center_x = ComputePressureCenterX(result);
    record.eval_ms = std::chrono::duration<double, std::milli>(stop - start).count();
    record.active_regions = result.active_regions;

    std::size_t active_samples = 0;
    for (const auto& region : result.regions) {
        active_samples += region.active_samples;
    }
    record.active_samples = active_samples;
    record.area_occ = result.occupied_area;
    record.area_band = result.band_area;
    record.area_sheet = result.sheet_footprint_area;
    record.area_corrected = result.corrected_active_area;
    record.active_area = record.area_band;

    if (result.sheet_result) {
        record.area_occ = result.sheet_result->occupied_area;
        record.area_band = result.sheet_result->band_area;
        record.area_sheet = result.sheet_result->sheet_footprint_area;
        record.sheet_center_x = result.sheet_result->sheet_center_world.x();
        record.sheet_pressure_center_x = result.sheet_result->pressure_center_world.x();
        record.sheet_area_ratio = result.sheet_result->area_ratio;
        record.sheet_mean_support_seed_count = result.sheet_result->mean_support_seed_count;
        record.sheet_normal_spread = result.sheet_result->normal_spread;
        record.sheet_patch_count = result.sheet_result->patch_count;
        record.sheet_fiber_count = result.sheet_result->fiber_count;
        record.sheet_fallback_regions = result.sheet_result->fallback_regions;
        record.sheet_used_fallback = result.sheet_result->used_fallback;
        if (!result.sheet_result->support_bbox_world.IsInverted()) {
            record.sheet_support_bbox_xmin = result.sheet_result->support_bbox_world.min.x();
            record.sheet_support_bbox_xmax = result.sheet_result->support_bbox_world.max.x();
            record.sheet_support_bbox_zmin = result.sheet_result->support_bbox_world.min.z();
            record.sheet_support_bbox_zmax = result.sheet_result->support_bbox_world.max.z();
        }
    }
    if (result.patch_consistency_result) {
        record.mean_patch_alpha = result.patch_consistency_result->mean_alpha;
        record.area_corrected = result.patch_consistency_result->total_corrected_area;
    }

    if (result.sheet_result) {
        for (const auto& region : result.sheet_result->regions) {
            for (const auto& patch : region.patches) {
                PatchAreaDiagnosticRecord patch_record;
                patch_record.scenario = scenario_name;
                patch_record.method = "distributed_sdf";
                patch_record.sample_index = pose.sample_index;
                patch_record.penetration = pose.penetration;
                patch_record.tilt_z_deg = Degrees(pose.tilt_z_rad);
                patch_record.offset_x = pose.offset_x;
                patch_record.region_id = region.region_id;
                patch_record.patch_id = patch.patch_id;
                patch_record.support_columns = patch.support_columns;
                patch_record.support_seed_count = patch.support_seed_count;
                patch_record.layered_cell_count = patch.layered_cell_count;
                patch_record.support_cell_count = patch.support_cells.size();
                patch_record.shell_cell_count =
                    std::count_if(patch.support_cells.begin(), patch.support_cells.end(),
                                  [](const ChSDFPatchPlaneSupportCell& cell) { return cell.shell; });
                patch_record.max_layer_count_per_cell = patch.max_layer_count_per_cell;
                patch_record.largest_connected_sheet_cells = patch.largest_connected_sheet_cells;
                patch_record.mean_layer_count_per_cell = patch.mean_layer_count_per_cell;
                if (!patch.support_cells.empty()) {
                    int min_u = std::numeric_limits<int>::max();
                    int max_u = std::numeric_limits<int>::min();
                    int min_v = std::numeric_limits<int>::max();
                    int max_v = std::numeric_limits<int>::min();
                    std::set<int> axes;
                    for (const auto& cell : patch.support_cells) {
                        min_u = std::min(min_u, cell.cell_ij.x());
                        max_u = std::max(max_u, cell.cell_ij.x());
                        min_v = std::min(min_v, cell.cell_ij.y());
                        max_v = std::max(max_v, cell.cell_ij.y());
                        axes.insert(cell.carrier_axis);
                    }
                    patch_record.support_cell_axis_count = axes.size();
                    patch_record.support_cell_u_span = max_u - min_u + 1;
                    patch_record.support_cell_v_span = max_v - min_v + 1;
                    patch_record.support_cell_bbox_count =
                        static_cast<std::size_t>(patch_record.support_cell_u_span) *
                        static_cast<std::size_t>(patch_record.support_cell_v_span);
                    if (patch_record.support_cell_bbox_count > 0) {
                        patch_record.support_cell_fill_ratio =
                            static_cast<double>(patch_record.support_cell_count) /
                            static_cast<double>(patch_record.support_cell_bbox_count);
                    }
                }
                if (patch.sheet_fill_ratio > 0) {
                    patch_record.support_cell_fill_ratio = patch.sheet_fill_ratio;
                }
                patch_record.polygon_vertex_count = patch.support_footprint.polygon_uv.size();
                patch_record.measure_area = patch.measure_area;
                patch_record.footprint_area = patch.footprint_area;
                patch_record.footprint_polygon_area = patch.support_footprint.area;
                patch_record.support_bbox_projected_area = ProjectedSupportBBoxArea(patch);
                if (patch_record.footprint_area > 1.0e-16) {
                    patch_record.polygon_to_footprint_ratio =
                        patch_record.footprint_polygon_area / patch_record.footprint_area;
                    patch_record.bbox_to_footprint_ratio =
                        patch_record.support_bbox_projected_area / patch_record.footprint_area;
                }
                if (patch_record.footprint_polygon_area > 1.0e-16) {
                    patch_record.bbox_to_polygon_ratio =
                        patch_record.support_bbox_projected_area / patch_record.footprint_polygon_area;
                }
                evaluation.patch_area_records.push_back(std::move(patch_record));
            }
        }
    }

    return evaluation;
}

std::vector<CurvedSheetRecord> EvaluateCurvedSheetModes(const std::string& scenario_name,
                                                        const PoseSample& pose,
                                                        ChBody* moving_body,
                                                        ChSDFShapePair& pair,
                                                        const ChSDFBrickPairBroadphase::Settings& pair_settings,
                                                        const ChSDFContactRegionBuilder::Settings& region_settings,
                                                        const ChSDFNormalPressureSettings& pressure_settings) {
    moving_body->SetPos(ChVector3d(0.0, 1.0 - pose.penetration, 0.0));
    moving_body->SetRot(QUNIT);
    moving_body->SetPosDt(VNULL);
    moving_body->SetAngVelLocal(VNULL);

    const auto original_settings = pair.GetSheetCollapseSettings();
    auto disabled_settings = original_settings;
    disabled_settings.enable = false;
    pair.SetSheetCollapseSettings(disabled_settings);
    const auto band_result = pair.EvaluateContact(pair_settings, region_settings, pressure_settings);
    pair.SetSheetCollapseSettings(original_settings);

    struct SheetModeSpec {
        const char* name;
        bool use_step6_v2;
        bool use_local_fiber_projection;
        bool allow_fallback;
    };

    const std::array<SheetModeSpec, 3> modes = {{
        {"legacy_dominant_axis", false, false, false},
        {"supplement5_v2_raw", true, false, false},
        {"supplement5_v2_auto", true, false, true},
    }};

    std::vector<CurvedSheetRecord> records;
    records.reserve(modes.size());

    for (const auto& mode : modes) {
        auto settings = original_settings;
        settings.enable = true;
        settings.use_step6_v2 = mode.use_step6_v2;
        settings.use_local_fiber_projection = mode.use_local_fiber_projection;
        settings.allow_dominant_axis_fallback = mode.allow_fallback;
        settings.max_patch_count_before_fallback = mode.allow_fallback ? 4 : 0;
        settings.max_fiber_count_before_fallback = mode.allow_fallback ? 64 : 0;

        const auto start = std::chrono::steady_clock::now();
        const auto sheet = ChSDFSheetBuilder::BuildShapePair(band_result, settings);
        const auto stop = std::chrono::steady_clock::now();

        CurvedSheetRecord record;
        record.scenario = scenario_name;
        record.mode = mode.name;
        record.sample_index = pose.sample_index;
        record.penetration = pose.penetration;
        record.valid = band_result.valid && sheet.HasSamples();
        record.area_occ = sheet.occupied_area;
        record.area_band = sheet.band_area;
        record.area_sheet = sheet.sheet_footprint_area;
        record.sheet_center_x = sheet.sheet_center_world.x();
        record.sheet_pressure_center_x = sheet.pressure_center_world.x();
        record.sheet_area_ratio = sheet.area_ratio;
        record.sheet_mean_support_seed_count = sheet.mean_support_seed_count;
        record.sheet_normal_spread = sheet.normal_spread;
        record.sheet_patch_count = sheet.patch_count;
        record.sheet_fiber_count = sheet.fiber_count;
        record.sheet_fallback_regions = sheet.fallback_regions;
        record.sheet_used_fallback = sheet.used_fallback;
        record.eval_ms = std::chrono::duration<double, std::milli>(stop - start).count();
        if (!sheet.support_bbox_world.IsInverted()) {
            record.sheet_support_bbox_xmin = sheet.support_bbox_world.min.x();
            record.sheet_support_bbox_xmax = sheet.support_bbox_world.max.x();
            record.sheet_support_bbox_zmin = sheet.support_bbox_world.min.z();
            record.sheet_support_bbox_zmax = sheet.support_bbox_world.max.z();
        }

        records.push_back(record);
    }

    return records;
}

BenchmarkRecord EvaluatePlanePenalty(const std::string& scenario_name,
                                     const std::string& method_name,
                                     const PoseSample& pose,
                                     std::size_t samples_u,
                                     std::size_t samples_v,
                                     double stiffness) {
    const ChFrame<> moving_frame(ChVector3d(pose.offset_x, 1.0 - pose.penetration, 0.0), QuatFromAngleZ(pose.tilt_z_rad));

    const double du = GridStep(samples_u, kHalfExtent);
    const double dv = GridStep(samples_v, kHalfExtent);

    BenchmarkRecord record;
    record.scenario = scenario_name;
    record.method = method_name;
    record.sample_index = pose.sample_index;
    record.penetration = pose.penetration;
    record.tilt_z_deg = Degrees(pose.tilt_z_rad);
    record.offset_x = pose.offset_x;
    record.valid = true;

    double weighted_center_x = 0;
    double pressure_weight_sum = 0;

    const auto start = std::chrono::steady_clock::now();

    for (std::size_t iv = 0; iv < samples_v; ++iv) {
        const double v = (samples_v > 1) ? (-kHalfExtent + dv * static_cast<double>(iv)) : 0.0;
        for (std::size_t iu = 0; iu < samples_u; ++iu) {
            const double u = (samples_u > 1) ? (-kHalfExtent + du * static_cast<double>(iu)) : 0.0;
            const double area = du * dv * TrapezoidWeight(iu, samples_u) * TrapezoidWeight(iv, samples_v);

            const ChVector3d point_local(u, -kHalfExtent, v);
            const ChVector3d point_world = moving_frame.TransformPointLocalToParent(point_local);

            if (std::abs(point_world.x()) > kHalfExtent + 1.0e-12 || std::abs(point_world.z()) > kHalfExtent + 1.0e-12) {
                continue;
            }

            const double penetration = std::max(kPlaneY - point_world.y(), 0.0);
            if (penetration <= 0) {
                continue;
            }

            const double pressure = stiffness * penetration;
            const ChVector3d force_world(0, pressure * area, 0);
            const ChVector3d torque_world = Vcross(point_world - moving_frame.GetPos(), force_world);

    record.force_x += force_world.x();
    record.force_y += force_world.y();
    record.force_z += force_world.z();
    record.force_y_corrected += force_world.y();
    record.torque_x += torque_world.x();
    record.torque_y += torque_world.y();
    record.torque_z += torque_world.z();
    record.torque_z_corrected += torque_world.z();
            record.active_area += area;
            record.max_penetration = std::max(record.max_penetration, penetration);
            record.active_samples++;

            weighted_center_x += point_world.x() * pressure * area;
            pressure_weight_sum += pressure * area;
        }
    }

    const auto stop = std::chrono::steady_clock::now();

    record.pressure_center_x = pressure_weight_sum > 0 ? weighted_center_x / pressure_weight_sum : 0.0;
    record.active_regions = record.active_area > 0 ? 1 : 0;
    record.area_occ = record.active_area;
    record.area_band = record.active_area;
    record.area_sheet = record.active_area;
    record.area_corrected = record.active_area;
    record.eval_ms = std::chrono::duration<double, std::milli>(stop - start).count();
    return record;
}

void WriteRecordsCsv(const fs::path& path, const std::vector<BenchmarkRecord>& records) {
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "scenario,method,sample_index,penetration,tilt_z_deg,offset_x,valid,force_x,force_y,force_y_corrected,force_z,torque_x,torque_y,torque_z,torque_z_corrected,active_area,area_occ,area_band,area_sheet,area_corrected,max_penetration,pressure_center_x,sheet_center_x,sheet_pressure_center_x,sheet_area_ratio,mean_patch_alpha,sheet_mean_support_seed_count,sheet_normal_spread,sheet_patch_count,sheet_fiber_count,sheet_fallback_regions,sheet_used_fallback,sheet_support_bbox_xmin,sheet_support_bbox_xmax,sheet_support_bbox_zmin,sheet_support_bbox_zmax,eval_ms,active_regions,active_samples\n";

    for (const auto& record : records) {
        out << record.scenario << ',' << record.method << ',' << record.sample_index << ',' << record.penetration << ','
            << record.tilt_z_deg << ',' << record.offset_x << ',' << (record.valid ? 1 : 0) << ',' << record.force_x
            << ',' << record.force_y << ',' << record.force_y_corrected << ',' << record.force_z << ','
            << record.torque_x << ',' << record.torque_y << ',' << record.torque_z << ','
            << record.torque_z_corrected << ',' << record.active_area << ',' << record.area_occ << ','
            << record.area_band << ',' << record.area_sheet << ',' << record.area_corrected << ','
            << record.max_penetration << ',' << record.pressure_center_x << ',' << record.sheet_center_x << ','
            << record.sheet_pressure_center_x << ',' << record.sheet_area_ratio << ',' << record.mean_patch_alpha << ','
            << record.sheet_mean_support_seed_count << ',' << record.sheet_normal_spread << ','
            << record.sheet_patch_count << ',' << record.sheet_fiber_count << ',' << record.sheet_fallback_regions
            << ',' << (record.sheet_used_fallback ? 1 : 0) << ',' << record.sheet_support_bbox_xmin << ','
            << record.sheet_support_bbox_xmax << ',' << record.sheet_support_bbox_zmin << ','
            << record.sheet_support_bbox_zmax << ',' << record.eval_ms << ',' << record.active_regions << ','
            << record.active_samples << '\n';
    }
}

std::vector<BenchmarkSummary> BuildSummaries(const std::vector<BenchmarkRecord>& records) {
    std::unordered_map<std::string, BenchmarkRecord> reference_records;
    for (const auto& record : records) {
        if (record.method == "dense_reference_161x161") {
            reference_records.emplace(KeyFor(record.scenario, record.sample_index), record);
        }
    }

    std::unordered_map<std::string, BenchmarkSummary> summaries;
    for (const auto& record : records) {
        if (record.method == "dense_reference_161x161") {
            continue;
        }

        const auto ref_it = reference_records.find(KeyFor(record.scenario, record.sample_index));
        if (ref_it == reference_records.end()) {
            continue;
        }

        const auto& reference = ref_it->second;
        const std::string key = record.scenario + "#" + record.method;
        auto& summary = summaries[key];
        summary.scenario = record.scenario;
        summary.method = record.method;
        summary.sample_count++;
        summary.mean_eval_ms += record.eval_ms;
        summary.contact_ratio += record.active_area > 0 ? 1.0 : 0.0;

        const double force_error = record.force_y - reference.force_y;
        const double force_corrected_error = record.force_y_corrected - reference.force_y;
        const double torque_error = record.torque_z - reference.torque_z;
        const double torque_corrected_error = record.torque_z_corrected - reference.torque_z;
        const double area_error = record.active_area - reference.active_area;
        const double area_occ_error = record.area_occ - reference.active_area;
        const double area_band_error = record.area_band - reference.active_area;
        const double area_sheet_error = record.area_sheet - reference.active_area;
        const double area_corrected_error = record.area_corrected - reference.active_area;
        const double center_error = record.pressure_center_x - reference.pressure_center_x;

        summary.force_y_rmse += force_error * force_error;
        summary.force_y_corrected_rmse += force_corrected_error * force_corrected_error;
        summary.torque_z_rmse += torque_error * torque_error;
        summary.torque_z_corrected_rmse += torque_corrected_error * torque_corrected_error;
        summary.active_area_rmse += area_error * area_error;
        summary.area_occ_rmse += area_occ_error * area_occ_error;
        summary.area_band_rmse += area_band_error * area_band_error;
        summary.area_sheet_rmse += area_sheet_error * area_sheet_error;
        summary.area_corrected_rmse += area_corrected_error * area_corrected_error;
        summary.pressure_center_x_rmse += center_error * center_error;
        summary.mean_patch_alpha += record.mean_patch_alpha;
        summary.force_y_max_abs = std::max(summary.force_y_max_abs, std::abs(force_error));
        summary.torque_z_max_abs = std::max(summary.torque_z_max_abs, std::abs(torque_error));
    }

    std::vector<BenchmarkSummary> result;
    result.reserve(summaries.size());
    for (auto& [key, summary] : summaries) {
        if (summary.sample_count == 0) {
            continue;
        }

        const double inv_count = 1.0 / static_cast<double>(summary.sample_count);
        summary.mean_eval_ms *= inv_count;
        summary.contact_ratio *= inv_count;
        summary.force_y_rmse = std::sqrt(summary.force_y_rmse * inv_count);
        summary.force_y_corrected_rmse = std::sqrt(summary.force_y_corrected_rmse * inv_count);
        summary.torque_z_rmse = std::sqrt(summary.torque_z_rmse * inv_count);
        summary.torque_z_corrected_rmse = std::sqrt(summary.torque_z_corrected_rmse * inv_count);
        summary.active_area_rmse = std::sqrt(summary.active_area_rmse * inv_count);
        summary.area_occ_rmse = std::sqrt(summary.area_occ_rmse * inv_count);
        summary.area_band_rmse = std::sqrt(summary.area_band_rmse * inv_count);
        summary.area_sheet_rmse = std::sqrt(summary.area_sheet_rmse * inv_count);
        summary.area_corrected_rmse = std::sqrt(summary.area_corrected_rmse * inv_count);
        summary.mean_patch_alpha *= inv_count;
        summary.pressure_center_x_rmse = std::sqrt(summary.pressure_center_x_rmse * inv_count);
        result.push_back(summary);
    }

    std::stable_sort(result.begin(), result.end(), [](const BenchmarkSummary& a, const BenchmarkSummary& b) {
        if (a.scenario != b.scenario) {
            return a.scenario < b.scenario;
        }
        return a.method < b.method;
    });

    return result;
}

void WriteSummaryCsv(const fs::path& path, const std::vector<BenchmarkSummary>& summaries) {
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "scenario,method,sample_count,mean_eval_ms,contact_ratio,force_y_rmse,force_y_corrected_rmse,force_y_max_abs,torque_z_rmse,torque_z_corrected_rmse,torque_z_max_abs,active_area_rmse,area_occ_rmse,area_band_rmse,area_sheet_rmse,area_corrected_rmse,mean_patch_alpha,pressure_center_x_rmse\n";

    for (const auto& summary : summaries) {
        out << summary.scenario << ',' << summary.method << ',' << summary.sample_count << ',' << summary.mean_eval_ms << ','
            << summary.contact_ratio << ',' << summary.force_y_rmse << ',' << summary.force_y_corrected_rmse << ','
            << summary.force_y_max_abs << ',' << summary.torque_z_rmse << ',' << summary.torque_z_corrected_rmse
            << ',' << summary.torque_z_max_abs << ',' << summary.active_area_rmse << ',' << summary.area_occ_rmse
            << ',' << summary.area_band_rmse << ',' << summary.area_sheet_rmse << ',' << summary.area_corrected_rmse
            << ',' << summary.mean_patch_alpha << ',' << summary.pressure_center_x_rmse << '\n';
    }
}

void WriteCenteredAreaModesCsv(const fs::path& path, const std::vector<BenchmarkRecord>& records) {
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "sample_index,penetration,area_occ,area_band,area_sheet,area_corrected,active_samples\n";

    std::vector<BenchmarkRecord> centered_records;
    centered_records.reserve(records.size());
    for (const auto& record : records) {
        if (record.scenario == "centered_compression" && record.method == "distributed_sdf") {
            centered_records.push_back(record);
        }
    }

    std::stable_sort(centered_records.begin(), centered_records.end(),
                     [](const BenchmarkRecord& a, const BenchmarkRecord& b) { return a.sample_index < b.sample_index; });

    for (const auto& record : centered_records) {
        out << record.sample_index << ',' << record.penetration << ',' << record.area_occ << ',' << record.area_band
            << ',' << record.area_sheet << ',' << record.area_corrected << ',' << record.active_samples << '\n';
    }
}

void PrintSummary(const std::vector<BenchmarkSummary>& summaries) {
    std::cout << "\nBenchmark summary against dense_reference_161x161\n";
    for (const auto& summary : summaries) {
        std::cout << "  [" << summary.scenario << "] " << summary.method << "  Fy_rmse=" << summary.force_y_rmse
                  << "  Fy_corr_rmse=" << summary.force_y_corrected_rmse << "  Tz_rmse=" << summary.torque_z_rmse
                  << "  Tz_corr_rmse=" << summary.torque_z_corrected_rmse << "  area_rmse=" << summary.active_area_rmse
                  << "  occ_rmse=" << summary.area_occ_rmse << "  sheet_rmse=" << summary.area_sheet_rmse
                  << "  corr_area_rmse=" << summary.area_corrected_rmse << "  alpha=" << summary.mean_patch_alpha
                  << "  center_x_rmse=" << summary.pressure_center_x_rmse << "  mean_ms=" << summary.mean_eval_ms
                  << "\n";
    }
}

std::vector<CurvedSheetSummary> BuildCurvedSheetSummaries(const std::vector<CurvedSheetRecord>& records) {
    std::unordered_map<std::string, CurvedSheetSummary> summaries;

    for (const auto& record : records) {
        const auto reference = MakeCylinderSheetReference(record.penetration);
        const std::string key = record.scenario + "#" + record.mode;
        auto& summary = summaries[key];
        summary.scenario = record.scenario;
        summary.mode = record.mode;
        summary.sample_count++;
        summary.mean_eval_ms += record.eval_ms;

        const double area_error = record.area_sheet - reference.sheet_area;
        const double center_error = record.sheet_center_x - reference.sheet_center_x;
        const double bbox_xmin_error = record.sheet_support_bbox_xmin - reference.sheet_support_bbox_xmin;
        const double bbox_xmax_error = record.sheet_support_bbox_xmax - reference.sheet_support_bbox_xmax;
        const double bbox_zmin_error = record.sheet_support_bbox_zmin - reference.sheet_support_bbox_zmin;
        const double bbox_zmax_error = record.sheet_support_bbox_zmax - reference.sheet_support_bbox_zmax;
        const double patch_error =
            static_cast<double>(record.sheet_patch_count) - static_cast<double>(reference.sheet_patch_count);

        summary.area_sheet_rmse += area_error * area_error;
        summary.sheet_center_x_rmse += center_error * center_error;
        summary.bbox_xmin_rmse += bbox_xmin_error * bbox_xmin_error;
        summary.bbox_xmax_rmse += bbox_xmax_error * bbox_xmax_error;
        summary.bbox_zmin_rmse += bbox_zmin_error * bbox_zmin_error;
        summary.bbox_zmax_rmse += bbox_zmax_error * bbox_zmax_error;
        summary.patch_count_rmse += patch_error * patch_error;
        summary.mean_area_ratio += record.sheet_area_ratio;
        summary.mean_support_seed_count += record.sheet_mean_support_seed_count;
        summary.mean_normal_spread += record.sheet_normal_spread;
        summary.mean_fiber_count += static_cast<double>(record.sheet_fiber_count);
        summary.fallback_ratio += record.sheet_used_fallback ? 1.0 : 0.0;
    }

    std::vector<CurvedSheetSummary> result;
    result.reserve(summaries.size());
    for (auto& [key, summary] : summaries) {
        const double inv_count = 1.0 / static_cast<double>(summary.sample_count);
        summary.mean_eval_ms *= inv_count;
        summary.area_sheet_rmse = std::sqrt(summary.area_sheet_rmse * inv_count);
        summary.sheet_center_x_rmse = std::sqrt(summary.sheet_center_x_rmse * inv_count);
        summary.bbox_xmin_rmse = std::sqrt(summary.bbox_xmin_rmse * inv_count);
        summary.bbox_xmax_rmse = std::sqrt(summary.bbox_xmax_rmse * inv_count);
        summary.bbox_zmin_rmse = std::sqrt(summary.bbox_zmin_rmse * inv_count);
        summary.bbox_zmax_rmse = std::sqrt(summary.bbox_zmax_rmse * inv_count);
        summary.patch_count_rmse = std::sqrt(summary.patch_count_rmse * inv_count);
        summary.mean_area_ratio *= inv_count;
        summary.mean_support_seed_count *= inv_count;
        summary.mean_normal_spread *= inv_count;
        summary.mean_fiber_count *= inv_count;
        summary.fallback_ratio *= inv_count;
        result.push_back(summary);
    }

    std::stable_sort(result.begin(), result.end(), [](const CurvedSheetSummary& a, const CurvedSheetSummary& b) {
        if (a.scenario != b.scenario) {
            return a.scenario < b.scenario;
        }
        return a.mode < b.mode;
    });

    return result;
}

void WriteCurvedSheetRecordsCsv(const fs::path& path, const std::vector<CurvedSheetRecord>& records) {
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "scenario,mode,sample_index,penetration,valid,area_occ,area_band,area_sheet,area_sheet_ref,sheet_center_x,sheet_center_x_ref,sheet_pressure_center_x,sheet_area_ratio,sheet_mean_support_seed_count,sheet_normal_spread,sheet_patch_count,sheet_patch_count_ref,sheet_fiber_count,sheet_fallback_regions,sheet_used_fallback,sheet_support_bbox_xmin,sheet_support_bbox_xmin_ref,sheet_support_bbox_xmax,sheet_support_bbox_xmax_ref,sheet_support_bbox_zmin,sheet_support_bbox_zmin_ref,sheet_support_bbox_zmax,sheet_support_bbox_zmax_ref,eval_ms\n";

    for (const auto& record : records) {
        const auto reference = MakeCylinderSheetReference(record.penetration);
        out << record.scenario << ',' << record.mode << ',' << record.sample_index << ',' << record.penetration << ','
            << (record.valid ? 1 : 0) << ',' << record.area_occ << ',' << record.area_band << ',' << record.area_sheet
            << ',' << reference.sheet_area << ',' << record.sheet_center_x << ',' << reference.sheet_center_x << ','
            << record.sheet_pressure_center_x << ',' << record.sheet_area_ratio << ','
            << record.sheet_mean_support_seed_count << ',' << record.sheet_normal_spread << ','
            << record.sheet_patch_count << ',' << reference.sheet_patch_count << ',' << record.sheet_fiber_count << ','
            << record.sheet_fallback_regions << ',' << (record.sheet_used_fallback ? 1 : 0) << ','
            << record.sheet_support_bbox_xmin << ',' << reference.sheet_support_bbox_xmin << ','
            << record.sheet_support_bbox_xmax << ',' << reference.sheet_support_bbox_xmax << ','
            << record.sheet_support_bbox_zmin << ',' << reference.sheet_support_bbox_zmin << ','
            << record.sheet_support_bbox_zmax << ',' << reference.sheet_support_bbox_zmax << ',' << record.eval_ms
            << '\n';
    }
}

void WriteCurvedSheetSummaryCsv(const fs::path& path, const std::vector<CurvedSheetSummary>& summaries) {
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "scenario,mode,sample_count,mean_eval_ms,area_sheet_rmse,sheet_center_x_rmse,bbox_xmin_rmse,bbox_xmax_rmse,bbox_zmin_rmse,bbox_zmax_rmse,patch_count_rmse,mean_area_ratio,mean_support_seed_count,mean_normal_spread,mean_fiber_count,fallback_ratio\n";

    for (const auto& summary : summaries) {
        out << summary.scenario << ',' << summary.mode << ',' << summary.sample_count << ',' << summary.mean_eval_ms
            << ',' << summary.area_sheet_rmse << ',' << summary.sheet_center_x_rmse << ','
            << summary.bbox_xmin_rmse << ',' << summary.bbox_xmax_rmse << ',' << summary.bbox_zmin_rmse << ','
            << summary.bbox_zmax_rmse << ',' << summary.patch_count_rmse << ',' << summary.mean_area_ratio << ','
            << summary.mean_support_seed_count << ',' << summary.mean_normal_spread << ','
            << summary.mean_fiber_count << ',' << summary.fallback_ratio << '\n';
    }
}

void WritePatchAreaDiagnosticsCsv(const fs::path& path, const std::vector<PatchAreaDiagnosticRecord>& records) {
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "scenario,method,sample_index,penetration,tilt_z_deg,offset_x,region_id,patch_id,support_columns,"
           "support_seed_count,layered_cell_count,support_cell_count,shell_cell_count,max_layer_count_per_cell,"
           "mean_layer_count_per_cell,largest_connected_sheet_cells,support_cell_axis_count,support_cell_bbox_count,"
           "support_cell_u_span,support_cell_v_span,support_cell_fill_ratio,polygon_vertex_count,measure_area,footprint_area,footprint_polygon_area,"
           "support_bbox_projected_area,polygon_to_footprint_ratio,bbox_to_footprint_ratio,bbox_to_polygon_ratio\n";

    for (const auto& record : records) {
        out << record.scenario << ',' << record.method << ',' << record.sample_index << ',' << record.penetration << ','
            << record.tilt_z_deg << ',' << record.offset_x << ',' << record.region_id << ',' << record.patch_id << ','
            << record.support_columns << ',' << record.support_seed_count << ',' << record.layered_cell_count << ','
            << record.support_cell_count << ',' << record.shell_cell_count << ',' << record.max_layer_count_per_cell
            << ',' << record.mean_layer_count_per_cell << ',' << record.largest_connected_sheet_cells << ','
            << record.support_cell_axis_count << ',' << record.support_cell_bbox_count << ','
            << record.support_cell_u_span << ',' << record.support_cell_v_span << ',' << record.support_cell_fill_ratio
            << ',' << record.polygon_vertex_count << ','
            << record.measure_area << ',' << record.footprint_area << ',' << record.footprint_polygon_area << ','
            << record.support_bbox_projected_area << ',' << record.polygon_to_footprint_ratio << ','
            << record.bbox_to_footprint_ratio << ',' << record.bbox_to_polygon_ratio << '\n';
    }
}

void PrintCurvedSheetSummary(const std::vector<CurvedSheetSummary>& summaries) {
    std::cout << "\nCurved sheet benchmark against analytic cylinder-plane footprint\n";
    for (const auto& summary : summaries) {
        std::cout << "  [" << summary.scenario << "] " << summary.mode << "  area_sheet_rmse="
                  << summary.area_sheet_rmse << "  center_x_rmse=" << summary.sheet_center_x_rmse
                  << "  bbox_x_rmse=(" << summary.bbox_xmin_rmse << "," << summary.bbox_xmax_rmse << ")"
                  << "  bbox_z_rmse=(" << summary.bbox_zmin_rmse << "," << summary.bbox_zmax_rmse << ")"
                  << "  patch_rmse=" << summary.patch_count_rmse << "  area_ratio=" << summary.mean_area_ratio
                  << "  support=" << summary.mean_support_seed_count << "  fibers=" << summary.mean_fiber_count
                  << "  fallback=" << summary.fallback_ratio << "  mean_ms=" << summary.mean_eval_ms << "\n";
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    fs::path input_path = DefaultNvdbPath();
    fs::path curved_input_path = DefaultCylinderNvdbPath();
    fs::path output_dir = DefaultOutputDir();

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            PrintUsage(argv[0]);
            return 0;
        }
        if (arg == "--input" && i + 1 < argc) {
            input_path = fs::path(argv[++i]);
            continue;
        }
        if (arg == "--curved-input" && i + 1 < argc) {
            curved_input_path = fs::path(argv[++i]);
            continue;
        }
        if (arg == "--output" && i + 1 < argc) {
            output_dir = fs::path(argv[++i]);
            continue;
        }

        std::cerr << "Unknown or incomplete argument: " << arg << "\n";
        PrintUsage(argv[0]);
        return 1;
    }

    if (!fs::exists(input_path)) {
        std::cerr << "Missing NanoVDB asset: " << input_path.string() << "\n";
        return 2;
    }
    if (!fs::exists(curved_input_path)) {
        std::cerr << "Missing curved NanoVDB asset: " << curved_input_path.string() << "\n";
        return 2;
    }

    fs::create_directories(output_dir);
    SetChronoDataPath(CHRONO_DATA_DIR);

    ChSystemNSC sys;

    auto fixed_body = chrono_types::make_shared<ChBodyEasyBox>(1.0, 1.0, 1.0, 1000.0, false, false);
    fixed_body->SetFixed(true);
    fixed_body->SetPos(ChVector3d(0, 0, 0));
    sys.Add(fixed_body);

    auto moving_body = chrono_types::make_shared<ChBodyEasyBox>(1.0, 1.0, 1.0, 1000.0, false, false);
    moving_body->SetFixed(false);
    moving_body->SetPos(ChVector3d(0, 1.0, 0));
    sys.Add(moving_body);

    auto material = chrono_types::make_shared<ChContactMaterialNSC>();
    auto sdf_shape_a = chrono_types::make_shared<ChCollisionShapeSDF>(material, input_path.string());
    auto sdf_shape_b = chrono_types::make_shared<ChCollisionShapeSDF>(material, input_path.string());
    auto sdf_shape_cylinder = chrono_types::make_shared<ChCollisionShapeSDF>(material, curved_input_path.string());

    if (!sdf_shape_a->IsLoaded() || !sdf_shape_b->IsLoaded() || !sdf_shape_cylinder->IsLoaded()) {
        std::cerr << "Failed to load NanoVDB asset.\n";
        std::cerr << "  shape A: " << sdf_shape_a->GetLastError() << "\n";
        std::cerr << "  shape B: " << sdf_shape_b->GetLastError() << "\n";
        std::cerr << "  shape cylinder: " << sdf_shape_cylinder->GetLastError() << "\n";
        return 3;
    }

    ChSDFShapePair pair;
    pair.SetBodyA(fixed_body.get());
    pair.SetBodyB(moving_body.get());
    pair.SetShapeA(sdf_shape_a);
    pair.SetShapeB(sdf_shape_b);
    pair.SetShapeAFrame(ChFrame<>());
    pair.SetShapeBFrame(ChFrame<>());
    pair.GetStabilizationSettings().enable_region_history = false;
    ChSDFSheetCollapseSettings planar_sheet_settings;
    planar_sheet_settings.use_step6_v2 = true;
    planar_sheet_settings.allow_dominant_axis_fallback = true;
    planar_sheet_settings.max_patch_count_before_fallback = 8;
    planar_sheet_settings.max_fiber_count_before_fallback = 256;
    pair.SetSheetCollapseSettings(planar_sheet_settings);

    ChSDFShapePair curved_pair;
    curved_pair.SetBodyA(fixed_body.get());
    curved_pair.SetBodyB(moving_body.get());
    curved_pair.SetShapeA(sdf_shape_a);
    curved_pair.SetShapeB(sdf_shape_cylinder);
    curved_pair.SetShapeAFrame(ChFrame<>());
    curved_pair.SetShapeBFrame(ChFrame<>());
    curved_pair.GetStabilizationSettings().enable_region_history = false;
    curved_pair.SetSheetCollapseSettings(planar_sheet_settings);

    ChSDFBrickPairBroadphase::Settings pair_settings;
    pair_settings.world_margin = 0.05;
    pair_settings.max_separation_distance = 0.06;
    pair_settings.max_min_abs_value = 0.08;

    ChSDFContactRegionBuilder::Settings region_settings;
    region_settings.sample_spacing = 0.025;
    region_settings.sample_max_abs_distance = 0.06;
    region_settings.max_combined_gap = 0.05;
    region_settings.min_opposed_normal_cosine = 0.7;
    region_settings.min_neighbor_normal_cosine = 0.6;
    region_settings.min_region_samples = 4;
    region_settings.neighbor_mode = 6;

    ChSDFNormalPressureSettings pressure_settings;
    pressure_settings.stiffness = 4.0e5;
    pressure_settings.damping = 0;
    pressure_settings.damping_ratio = -1;
    pressure_settings.friction_coefficient = 0;
    pressure_settings.gradient_quality_gain = 0;
    pressure_settings.curvature_gain = 0;
    pressure_settings.resolution_scale_gain = 0;

    ChSDFPotentialFieldSettings potential_settings = sdf_shape_a->GetPotentialFieldSettings();
    potential_settings =
        ChSDFPotentialCalibration::MakePlanarConsistentPotentialSettings(pressure_settings.stiffness, potential_settings);
    potential_settings.depth_scale = 1.0;
    potential_settings.depth_cap = -1;
    sdf_shape_a->SetPotentialFieldSettings(potential_settings);
    sdf_shape_b->SetPotentialFieldSettings(potential_settings);
    sdf_shape_cylinder->SetPotentialFieldSettings(potential_settings);

    const auto scenarios = BuildScenarios();
    const auto curved_sheet_scenario = BuildCurvedSheetScenario();
    const auto methods = BuildMethods();

    std::vector<BenchmarkRecord> records;
    records.reserve(64);
    std::vector<PatchAreaDiagnosticRecord> patch_area_records;
    patch_area_records.reserve(256);
    std::vector<CurvedSheetRecord> curved_sheet_records;
    curved_sheet_records.reserve(32);

    for (const auto& scenario : scenarios) {
        for (const auto& pose : scenario.poses) {
            for (const auto& method : methods) {
                if (method.is_distributed) {
                    auto evaluation =
                        EvaluateDistributed(scenario.name, pose, moving_body.get(), pair, pair_settings, region_settings,
                                            pressure_settings);
                    records.push_back(std::move(evaluation.record));
                    patch_area_records.insert(patch_area_records.end(), evaluation.patch_area_records.begin(),
                                              evaluation.patch_area_records.end());
                } else {
                    records.push_back(
                        EvaluatePlanePenalty(scenario.name, method.name, pose, method.samples_u, method.samples_v,
                                             pressure_settings.stiffness));
                }
            }
        }
    }

    const auto summaries = BuildSummaries(records);

    auto curved_pair_settings = pair_settings;
    curved_pair_settings.world_margin = std::max(curved_pair_settings.world_margin, 0.08);
    curved_pair_settings.max_separation_distance = std::max(curved_pair_settings.max_separation_distance, 0.08);
    curved_pair_settings.max_min_abs_value = std::max(curved_pair_settings.max_min_abs_value, 0.10);

    auto curved_region_settings = region_settings;
    curved_region_settings.sample_max_abs_distance = std::max(curved_region_settings.sample_max_abs_distance, 0.08);
    curved_region_settings.max_combined_gap = std::max(curved_region_settings.max_combined_gap, 0.08);

    for (const auto& pose : curved_sheet_scenario.poses) {
        const auto pose_records =
            EvaluateCurvedSheetModes(curved_sheet_scenario.name, pose, moving_body.get(), curved_pair,
                                     curved_pair_settings, curved_region_settings, pressure_settings);
        curved_sheet_records.insert(curved_sheet_records.end(), pose_records.begin(), pose_records.end());
    }

    const auto curved_sheet_summaries = BuildCurvedSheetSummaries(curved_sheet_records);

    const fs::path records_path = output_dir / "benchmark_records.csv";
    const fs::path summary_path = output_dir / "benchmark_summary.csv";
    const fs::path centered_area_path = output_dir / "centered_area_modes.csv";
    const fs::path patch_area_diagnostics_path = output_dir / "patch_area_diagnostics.csv";
    const fs::path curved_sheet_records_path = output_dir / "curved_sheet_records.csv";
    const fs::path curved_sheet_summary_path = output_dir / "curved_sheet_summary.csv";
    WriteRecordsCsv(records_path, records);
    WriteSummaryCsv(summary_path, summaries);
    WriteCenteredAreaModesCsv(centered_area_path, records);
    WritePatchAreaDiagnosticsCsv(patch_area_diagnostics_path, patch_area_records);
    WriteCurvedSheetRecordsCsv(curved_sheet_records_path, curved_sheet_records);
    WriteCurvedSheetSummaryCsv(curved_sheet_summary_path, curved_sheet_summaries);

    std::cout << "Wrote:\n"
              << "  " << records_path.string() << "\n"
              << "  " << summary_path.string() << "\n"
              << "  " << centered_area_path.string() << "\n"
              << "  " << patch_area_diagnostics_path.string() << "\n"
              << "  " << curved_sheet_records_path.string() << "\n"
              << "  " << curved_sheet_summary_path.string() << "\n";

    // The benchmark has already materialized all outputs above. In the current
    // Windows build, later teardown and summary-printing paths are unstable, so
    // exit explicitly once the on-disk results are complete.
    std::cout.flush();
    std::cerr.flush();
    std::fflush(stdout);
    std::fflush(stderr);
    std::_Exit(0);
}
