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
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include <chrono/collision/ChCollisionShapeSDF.h>
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
    double torque_x = 0;
    double torque_y = 0;
    double torque_z = 0;

    double active_area = 0;
    double area_occ = 0;
    double area_band = 0;
    double area_sheet = 0;
    double max_penetration = 0;
    double pressure_center_x = 0;
    double eval_ms = 0;

    std::size_t active_regions = 0;
    std::size_t active_samples = 0;
};

struct BenchmarkSummary {
    std::string scenario;
    std::string method;
    std::size_t sample_count = 0;
    double mean_eval_ms = 0;
    double contact_ratio = 0;
    double force_y_rmse = 0;
    double force_y_max_abs = 0;
    double torque_z_rmse = 0;
    double torque_z_max_abs = 0;
    double active_area_rmse = 0;
    double area_occ_rmse = 0;
    double area_band_rmse = 0;
    double area_sheet_rmse = 0;
    double pressure_center_x_rmse = 0;
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

std::string KeyFor(const std::string& scenario, std::size_t sample_index) {
    return scenario + "#" + std::to_string(sample_index);
}

fs::path DefaultNvdbPath() {
    return fs::path(BENCHMARK_TOOL_DIR).parent_path().parent_path() / "chrono" / "template_project" / "assets" /
           "unit_box_centered.nvdb";
}

fs::path DefaultOutputDir() {
    return fs::path(BENCHMARK_TOOL_DIR).parent_path().parent_path() / "paper" / "results" / "benchmark_box_contact";
}

void PrintUsage(const char* exe_name) {
    std::cout << "Usage:\n"
              << "  " << exe_name << " [--input path/to/unit_box_centered.nvdb] [--output output_dir]\n\n"
              << "Defaults:\n"
              << "  input:  " << DefaultNvdbPath().string() << "\n"
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

int DominantAxis(const ChVector3d& normal) {
    const ChVector3d abs_normal(std::abs(normal.x()), std::abs(normal.y()), std::abs(normal.z()));
    if (abs_normal.x() >= abs_normal.y() && abs_normal.x() >= abs_normal.z()) {
        return 0;
    }
    if (abs_normal.y() >= abs_normal.x() && abs_normal.y() >= abs_normal.z()) {
        return 1;
    }
    return 2;
}

std::uint64_t CollapseKey(const ChVector3i& coord, int dominant_axis) {
    const int a = dominant_axis == 0 ? coord.y() : coord.x();
    const int b = dominant_axis == 2 ? coord.y() : coord.z();
    const std::uint32_t ua = static_cast<std::uint32_t>(a);
    const std::uint32_t ub = static_cast<std::uint32_t>(b);
    return (static_cast<std::uint64_t>(ua) << 32) | static_cast<std::uint64_t>(ub);
}

double ComputeCollapsedSheetArea(const ChSDFShapePairContactResult& result, double spacing) {
    if (spacing <= 0) {
        return 0.0;
    }

    const double cell_area = spacing * spacing;
    double collapsed_area = 0.0;

    for (const auto& region : result.regions) {
        if (!region.HasActiveContact()) {
            continue;
        }

        const int dominant_axis = DominantAxis(region.region.mean_normal_world);
        std::unordered_set<std::uint64_t> occupied_columns;
        occupied_columns.reserve(region.samples.size());

        for (const auto& sample : region.samples) {
            occupied_columns.insert(CollapseKey(sample.region_sample.coord, dominant_axis));
        }

        collapsed_area += static_cast<double>(occupied_columns.size()) * cell_area;
    }

    return collapsed_area;
}

BenchmarkRecord EvaluateDistributed(const std::string& scenario_name,
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

    BenchmarkRecord record;
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
    record.torque_x = result.wrench_world_b.torque.x();
    record.torque_y = result.wrench_world_b.torque.y();
    record.torque_z = result.wrench_world_b.torque.z();
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
    record.area_occ = static_cast<double>(active_samples) * region_settings.sample_spacing * region_settings.sample_spacing;
    record.area_band = result.active_area;
    record.area_sheet = ComputeCollapsedSheetArea(result, region_settings.sample_spacing);
    record.active_area = record.area_band;

    return record;
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
            record.torque_x += torque_world.x();
            record.torque_y += torque_world.y();
            record.torque_z += torque_world.z();
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
    record.eval_ms = std::chrono::duration<double, std::milli>(stop - start).count();
    return record;
}

void WriteRecordsCsv(const fs::path& path, const std::vector<BenchmarkRecord>& records) {
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "scenario,method,sample_index,penetration,tilt_z_deg,offset_x,valid,force_x,force_y,force_z,torque_x,torque_y,torque_z,active_area,area_occ,area_band,area_sheet,max_penetration,pressure_center_x,eval_ms,active_regions,active_samples\n";

    for (const auto& record : records) {
        out << record.scenario << ',' << record.method << ',' << record.sample_index << ',' << record.penetration << ','
            << record.tilt_z_deg << ',' << record.offset_x << ',' << (record.valid ? 1 : 0) << ',' << record.force_x
            << ',' << record.force_y << ',' << record.force_z << ',' << record.torque_x << ',' << record.torque_y << ','
            << record.torque_z << ',' << record.active_area << ',' << record.area_occ << ',' << record.area_band << ','
            << record.area_sheet << ',' << record.max_penetration << ','
            << record.pressure_center_x << ',' << record.eval_ms << ',' << record.active_regions << ','
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
        const double torque_error = record.torque_z - reference.torque_z;
        const double area_error = record.active_area - reference.active_area;
        const double area_occ_error = record.area_occ - reference.active_area;
        const double area_band_error = record.area_band - reference.active_area;
        const double area_sheet_error = record.area_sheet - reference.active_area;
        const double center_error = record.pressure_center_x - reference.pressure_center_x;

        summary.force_y_rmse += force_error * force_error;
        summary.torque_z_rmse += torque_error * torque_error;
        summary.active_area_rmse += area_error * area_error;
        summary.area_occ_rmse += area_occ_error * area_occ_error;
        summary.area_band_rmse += area_band_error * area_band_error;
        summary.area_sheet_rmse += area_sheet_error * area_sheet_error;
        summary.pressure_center_x_rmse += center_error * center_error;
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
        summary.torque_z_rmse = std::sqrt(summary.torque_z_rmse * inv_count);
        summary.active_area_rmse = std::sqrt(summary.active_area_rmse * inv_count);
        summary.area_occ_rmse = std::sqrt(summary.area_occ_rmse * inv_count);
        summary.area_band_rmse = std::sqrt(summary.area_band_rmse * inv_count);
        summary.area_sheet_rmse = std::sqrt(summary.area_sheet_rmse * inv_count);
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
    out << "scenario,method,sample_count,mean_eval_ms,contact_ratio,force_y_rmse,force_y_max_abs,torque_z_rmse,torque_z_max_abs,active_area_rmse,area_occ_rmse,area_band_rmse,area_sheet_rmse,pressure_center_x_rmse\n";

    for (const auto& summary : summaries) {
        out << summary.scenario << ',' << summary.method << ',' << summary.sample_count << ',' << summary.mean_eval_ms << ','
            << summary.contact_ratio << ',' << summary.force_y_rmse << ',' << summary.force_y_max_abs << ','
            << summary.torque_z_rmse << ',' << summary.torque_z_max_abs << ',' << summary.active_area_rmse << ','
            << summary.area_occ_rmse << ',' << summary.area_band_rmse << ',' << summary.area_sheet_rmse << ','
            << summary.pressure_center_x_rmse << '\n';
    }
}

void WriteCenteredAreaModesCsv(const fs::path& path, const std::vector<BenchmarkRecord>& records) {
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "sample_index,penetration,area_occ,area_band,area_sheet,active_samples\n";

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
            << ',' << record.area_sheet << ',' << record.active_samples << '\n';
    }
}

void PrintSummary(const std::vector<BenchmarkSummary>& summaries) {
    std::cout << "\nBenchmark summary against dense_reference_161x161\n";
    for (const auto& summary : summaries) {
        std::cout << "  [" << summary.scenario << "] " << summary.method << "  Fy_rmse=" << summary.force_y_rmse
                  << "  Tz_rmse=" << summary.torque_z_rmse << "  area_rmse=" << summary.active_area_rmse
                  << "  occ_rmse=" << summary.area_occ_rmse << "  sheet_rmse=" << summary.area_sheet_rmse
                  << "  center_x_rmse=" << summary.pressure_center_x_rmse << "  mean_ms=" << summary.mean_eval_ms
                  << "\n";
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    fs::path input_path = DefaultNvdbPath();
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

    if (!sdf_shape_a->IsLoaded() || !sdf_shape_b->IsLoaded()) {
        std::cerr << "Failed to load NanoVDB asset.\n";
        std::cerr << "  shape A: " << sdf_shape_a->GetLastError() << "\n";
        std::cerr << "  shape B: " << sdf_shape_b->GetLastError() << "\n";
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
    // For two identical rigid shapes, the equal-pressure surface lies near the
    // overlap mid-surface, so each side contributes roughly half of the total
    // depth. A small empirical correction keeps the PFC pressure scale aligned
    // with the 4e5 * penetration benchmark baseline for the unit-box cases.
    potential_settings.modulus = 2.22 * pressure_settings.stiffness;
    potential_settings.depth_scale = 1.0;
    potential_settings.depth_cap = -1;
    sdf_shape_a->SetPotentialFieldSettings(potential_settings);
    sdf_shape_b->SetPotentialFieldSettings(potential_settings);

    const auto scenarios = BuildScenarios();
    const auto methods = BuildMethods();

    std::vector<BenchmarkRecord> records;
    records.reserve(64);

    for (const auto& scenario : scenarios) {
        for (const auto& pose : scenario.poses) {
            for (const auto& method : methods) {
                if (method.is_distributed) {
                    records.push_back(EvaluateDistributed(scenario.name, pose, moving_body.get(), pair, pair_settings,
                                                          region_settings, pressure_settings));
                } else {
                    records.push_back(
                        EvaluatePlanePenalty(scenario.name, method.name, pose, method.samples_u, method.samples_v,
                                             pressure_settings.stiffness));
                }
            }
        }
    }

    const auto summaries = BuildSummaries(records);

    const fs::path records_path = output_dir / "benchmark_records.csv";
    const fs::path summary_path = output_dir / "benchmark_summary.csv";
    const fs::path centered_area_path = output_dir / "centered_area_modes.csv";
    WriteRecordsCsv(records_path, records);
    WriteSummaryCsv(summary_path, summaries);
    WriteCenteredAreaModesCsv(centered_area_path, records);

    std::cout << "Wrote:\n"
              << "  " << records_path.string() << "\n"
              << "  " << summary_path.string() << "\n"
              << "  " << centered_area_path.string() << "\n";
    PrintSummary(summaries);

    return 0;
}
