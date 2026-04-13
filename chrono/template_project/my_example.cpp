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
// Minimal dual-SDF contact demo.
//
// Usage:
//   my_demo [path/to/unit_box_centered.nvdb]
//
// If no argument is provided, this example looks for:
//   TEMPLATE_PROJECT_DIR/assets/unit_box_centered.nvdb
//
// You can generate that file with the standalone preprocessing tool:
//   tools/mesh_to_vdb_nvdb/build/Release/mesh_to_vdb_nvdb.exe
//       --input chrono/template_project/assets/unit_box_centered.obj
//       --voxel-size 0.05
//       --band-width 4
// =============================================================================

#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include <chrono/collision/ChCollisionShapeSDF.h>
#include <chrono/collision/sdf/ChSDFShapePair.h>
#include <chrono/physics/ChBodyEasy.h>
#include <chrono/physics/ChContactMaterialNSC.h>
#include <chrono/physics/ChSystemNSC.h>

using namespace chrono;
namespace fs = std::filesystem;

namespace {

fs::path DefaultNvdbPath() {
    return fs::path(TEMPLATE_PROJECT_DIR) / "assets" / "unit_box_centered.nvdb";
}

void PrintUsage(const char* exe_name) {
    std::cout << "Usage:\n"
              << "  " << exe_name << " [path/to/unit_box_centered.nvdb]\n\n"
              << "Default path:\n"
              << "  " << DefaultNvdbPath().string() << "\n";
}

bool IsFiniteWrench(const ChWrenchd& wrench) {
    const double values[] = {wrench.force.x(),  wrench.force.y(),  wrench.force.z(),
                             wrench.torque.x(), wrench.torque.y(), wrench.torque.z()};
    for (double value : values) {
        if (!std::isfinite(value)) {
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    const fs::path nvdb_path = (argc > 1) ? fs::path(argv[1]) : DefaultNvdbPath();
    if (argc > 1 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        PrintUsage(argv[0]);
        return 0;
    }

    if (!fs::exists(nvdb_path)) {
        std::cerr << "Missing NanoVDB file: " << nvdb_path.string() << "\n";
        PrintUsage(argv[0]);
        return 1;
    }

    SetChronoDataPath(CHRONO_DATA_DIR);

    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, -9.81, 0));

    auto fixed_body = chrono_types::make_shared<ChBodyEasyBox>(1.0, 1.0, 1.0, 1000.0, false, false);
    fixed_body->SetFixed(true);
    fixed_body->SetPos(ChVector3d(0, 0, 0));
    sys.Add(fixed_body);

    auto moving_body = chrono_types::make_shared<ChBodyEasyBox>(1.0, 1.0, 1.0, 1000.0, false, false);
    moving_body->SetPos(ChVector3d(0, 0.98, 0));
    moving_body->SetLinVel(ChVector3d(0.20, -0.05, 0));
    sys.Add(moving_body);

    auto material = chrono_types::make_shared<ChContactMaterialNSC>();
    auto sdf_shape_a = chrono_types::make_shared<ChCollisionShapeSDF>(material, nvdb_path.string());
    auto sdf_shape_b = chrono_types::make_shared<ChCollisionShapeSDF>(material, nvdb_path.string());

    if (!sdf_shape_a->IsLoaded() || !sdf_shape_b->IsLoaded()) {
        std::cerr << "Failed to load NanoVDB SDF.\n";
        std::cerr << "  shape A: " << sdf_shape_a->GetLastError() << "\n";
        std::cerr << "  shape B: " << sdf_shape_b->GetLastError() << "\n";
        return 2;
    }

    ChSDFShapePair pair;
    pair.SetBodyA(fixed_body.get());
    pair.SetBodyB(moving_body.get());
    pair.SetShapeA(sdf_shape_a);
    pair.SetShapeB(sdf_shape_b);
    pair.SetShapeAFrame(ChFrame<>());
    pair.SetShapeBFrame(ChFrame<>());
    pair.CreateAccumulators();

    ChSDFBrickPairBroadphase::Settings pair_settings;
    pair_settings.world_margin = 0.10;
    pair_settings.max_separation_distance = 0.10;
    pair_settings.max_min_abs_value = 0.12;

    ChSDFContactRegionBuilder::Settings region_settings;
    region_settings.sample_spacing = 0.05;
    region_settings.sample_max_abs_distance = 0.08;
    region_settings.max_combined_gap = 0.08;
    region_settings.min_opposed_normal_cosine = 0.8;
    region_settings.min_neighbor_normal_cosine = 0.7;
    region_settings.min_region_samples = 9;
    region_settings.neighbor_mode = 6;

    ChSDFNormalPressureSettings pressure_settings;
    pressure_settings.stiffness = 4.0e5;
    pressure_settings.damping = 4.0e3;
    pressure_settings.max_pressure = 2.0e6;
    pressure_settings.gradient_quality_gain = 2.0;
    pressure_settings.resolution_scale_gain = 1.0;
    pressure_settings.reference_resolution_length = 0.025;
    pressure_settings.min_stiffness_scale = 0.25;
    pressure_settings.friction_coefficient = 0.20;
    pressure_settings.tangential_velocity_regularization = 0.02;

    const double step_size = 1.0e-3;
    const double end_time = 0.25;

    std::size_t contact_steps = 0;
    double max_upward_force = 0;
    double max_abs_tangential_force = 0;
    double min_height = std::numeric_limits<double>::infinity();
    bool finite = true;
    ChSDFShapePairContactResult last_result;

    while (sys.GetChTime() < end_time) {
        pair.EmptyAccumulators();
        const auto result = pair.EvaluateAndApply(pair_settings, region_settings, pressure_settings);
        last_result = result;

        if (result.HasActiveContact()) {
            ++contact_steps;
            max_upward_force = std::max(max_upward_force, result.wrench_world_b.force.y());
            max_abs_tangential_force = std::max(max_abs_tangential_force, std::abs(result.wrench_world_b.force.x()));
        }

        finite = finite && IsFiniteWrench(result.wrench_world_a) && IsFiniteWrench(result.wrench_world_b) &&
                 std::isfinite(moving_body->GetPos().x()) && std::isfinite(moving_body->GetPos().y()) &&
                 std::isfinite(moving_body->GetLinVel().x()) && std::isfinite(moving_body->GetLinVel().y());
        min_height = std::min(min_height, moving_body->GetPos().y());

        if (static_cast<int>(std::round(sys.GetChTime() / step_size)) % 50 == 0) {
            std::cout << "t=" << sys.GetChTime() << "  x=" << moving_body->GetPos().x()
                      << "  y=" << moving_body->GetPos().y() << "  vx=" << moving_body->GetLinVel().x()
                      << "  vy=" << moving_body->GetLinVel().y() << "  active_regions=" << result.active_regions
                      << "  k_mean=" << result.mean_local_stiffness << "  k_max=" << result.max_local_stiffness
                      << "  Fx=" << result.wrench_world_b.force.x() << "  Fy=" << result.wrench_world_b.force.y()
                      << "\n";
        }

        sys.DoStepDynamics(step_size);
    }

    std::cout << "\nFinal state\n";
    std::cout << "  time:          " << sys.GetChTime() << "\n";
    std::cout << "  moving x:      " << moving_body->GetPos().x() << "\n";
    std::cout << "  moving y:      " << moving_body->GetPos().y() << "\n";
    std::cout << "  moving vx:     " << moving_body->GetLinVel().x() << "\n";
    std::cout << "  moving vy:     " << moving_body->GetLinVel().y() << "\n";
    std::cout << "  min y:         " << min_height << "\n";
    std::cout << "  contact steps: " << contact_steps << "\n";
    std::cout << "  mean k:        " << last_result.mean_local_stiffness << "\n";
    std::cout << "  max k:         " << last_result.max_local_stiffness << "\n";
    std::cout << "  max |Fx|:      " << max_abs_tangential_force << "\n";
    std::cout << "  max Fy:        " << max_upward_force << "\n";

    if (!finite) {
        std::cerr << "Non-finite state detected during integration.\n";
        return 3;
    }

    if (contact_steps == 0) {
        std::cerr << "The demo did not detect any active contact regions.\n";
        return 4;
    }

    if (max_abs_tangential_force <= 1.0) {
        std::cerr << "The demo did not detect any meaningful tangential traction.\n";
        return 5;
    }

    return 0;
}
