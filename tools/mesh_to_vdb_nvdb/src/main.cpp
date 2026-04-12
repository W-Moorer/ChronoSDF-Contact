#include "MeshIO.h"

#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <openvdb/io/File.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>

#include <nanovdb/NanoVDB.h>
#include <nanovdb/io/IO.h>
#include <nanovdb/tools/CreateNanoGrid.h>

namespace fs = std::filesystem;

namespace {

struct Arguments {
    fs::path input_path;
    fs::path vdb_output;
    fs::path nvdb_output;
    std::string grid_name;
    double voxel_size = 0.005;
    double exterior_band_width = 3.0;
    double interior_band_width = 3.0;
    double scale = 1.0;
    bool allow_open = false;
};

void PrintUsage() {
    std::cout << "Usage:\n"
                 "  mesh_to_vdb_nvdb --input <mesh.obj|mesh.stl> --voxel-size <meters>\n"
                 "                   [--vdb-out <path>] [--nvdb-out <path>]\n"
                 "                   [--grid-name <name>] [--band-width <voxels>]\n"
                 "                   [--ex-band <voxels>] [--in-band <voxels>]\n"
                 "                   [--scale <factor>] [--allow-open]\n\n"
                 "Notes:\n"
                 "  - OBJ and binary STL are supported.\n"
                 "  - Band widths are specified in voxel units.\n"
                 "  - Without --allow-open, the tool rejects non-watertight meshes.\n";
}

std::string RequireValue(const std::string& option, int argc, char* argv[], int& index) {
    if (index + 1 >= argc) {
        throw std::runtime_error("Missing value for " + option);
    }
    ++index;
    return argv[index];
}

double ParseDouble(const std::string& text, const std::string& option) {
    try {
        std::size_t consumed = 0;
        const double value = std::stod(text, &consumed);
        if (consumed != text.size()) {
            throw std::runtime_error("");
        }
        return value;
    } catch (...) {
        throw std::runtime_error("Invalid numeric value for " + option + ": " + text);
    }
}

Arguments ParseArguments(int argc, char* argv[]) {
    Arguments args;

    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--help" || option == "-h") {
            PrintUsage();
            std::exit(0);
        } else if (option == "--input" || option == "-i") {
            args.input_path = RequireValue(option, argc, argv, i);
        } else if (option == "--vdb-out") {
            args.vdb_output = RequireValue(option, argc, argv, i);
        } else if (option == "--nvdb-out") {
            args.nvdb_output = RequireValue(option, argc, argv, i);
        } else if (option == "--grid-name") {
            args.grid_name = RequireValue(option, argc, argv, i);
        } else if (option == "--voxel-size") {
            args.voxel_size = ParseDouble(RequireValue(option, argc, argv, i), option);
        } else if (option == "--band-width") {
            const double width = ParseDouble(RequireValue(option, argc, argv, i), option);
            args.exterior_band_width = width;
            args.interior_band_width = width;
        } else if (option == "--ex-band") {
            args.exterior_band_width = ParseDouble(RequireValue(option, argc, argv, i), option);
        } else if (option == "--in-band") {
            args.interior_band_width = ParseDouble(RequireValue(option, argc, argv, i), option);
        } else if (option == "--scale") {
            args.scale = ParseDouble(RequireValue(option, argc, argv, i), option);
        } else if (option == "--allow-open") {
            args.allow_open = true;
        } else {
            throw std::runtime_error("Unknown option: " + option);
        }
    }

    if (args.input_path.empty()) {
        throw std::runtime_error("Missing required option --input");
    }
    if (args.voxel_size <= 0.0) {
        throw std::runtime_error("--voxel-size must be positive");
    }
    if (args.exterior_band_width <= 0.0 || args.interior_band_width <= 0.0) {
        throw std::runtime_error("Band widths must be positive");
    }
    if (args.scale <= 0.0) {
        throw std::runtime_error("--scale must be positive");
    }

    if (args.vdb_output.empty()) {
        args.vdb_output = args.input_path;
        args.vdb_output.replace_extension(".vdb");
    }
    if (args.nvdb_output.empty()) {
        args.nvdb_output = args.input_path;
        args.nvdb_output.replace_extension(".nvdb");
    }
    if (args.grid_name.empty()) {
        args.grid_name = args.input_path.stem().string();
    }

    return args;
}

std::vector<openvdb::Vec3s> ConvertPoints(const sdfprep::TriangleMeshData& mesh) {
    std::vector<openvdb::Vec3s> points;
    points.reserve(mesh.points.size());
    for (const auto& p : mesh.points) {
        points.emplace_back(p.x, p.y, p.z);
    }
    return points;
}

std::vector<openvdb::Vec3I> ConvertTriangles(const sdfprep::TriangleMeshData& mesh) {
    std::vector<openvdb::Vec3I> triangles;
    triangles.reserve(mesh.triangles.size());
    for (const auto& tri : mesh.triangles) {
        triangles.emplace_back(tri.v0, tri.v1, tri.v2);
    }
    return triangles;
}

void WriteVdb(const openvdb::GridBase::Ptr& grid, const fs::path& output_path) {
    if (!output_path.parent_path().empty()) {
        fs::create_directories(output_path.parent_path());
    }
    openvdb::GridPtrVec grids;
    grids.push_back(grid);

    openvdb::io::File file(output_path.string());
    file.write(grids);
    file.close();
}

auto CreateNanoHandle(const openvdb::GridBase::Ptr& grid) {
    return nanovdb::tools::openToNanoVDB(grid);
}

template <class HandleT>
void ValidateNanoFile(const fs::path& nvdb_path, const HandleT&) {
    auto loaded = nanovdb::io::readGrid(nvdb_path.string());
    if (loaded.template grid<float>() == nullptr) {
        throw std::runtime_error("NanoVDB validation failed for '" + nvdb_path.string() + "'");
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Arguments args = ParseArguments(argc, argv);

        openvdb::initialize();

        const auto mesh = sdfprep::LoadTriangleMesh(args.input_path, args.scale);
        if (mesh.points.empty() || mesh.triangles.empty()) {
            throw std::runtime_error("Input mesh is empty after loading.");
        }

        const auto bbox = sdfprep::ComputeBoundingBox(mesh);
        const auto topology = sdfprep::AnalyzeTopology(mesh);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Input mesh: " << args.input_path.string() << "\n";
        std::cout << "Vertices:   " << mesh.points.size() << "\n";
        std::cout << "Triangles:  " << mesh.triangles.size() << "\n";
        if (bbox.valid) {
            std::cout << "BBox min:   (" << bbox.min.x << ", " << bbox.min.y << ", " << bbox.min.z << ")\n";
            std::cout << "BBox max:   (" << bbox.max.x << ", " << bbox.max.y << ", " << bbox.max.z << ")\n";
        }
        std::cout << "Boundary edges:     " << topology.boundary_edges << "\n";
        std::cout << "Non-manifold edges: " << topology.non_manifold_edges << "\n";
        std::cout << "Degenerate faces:   " << topology.degenerate_faces << "\n";

        if (!args.allow_open && topology.boundary_edges > 0) {
            throw std::runtime_error(
                "Mesh has open boundary edges. Fix the mesh or rerun with --allow-open to force signed conversion.");
        }

        const auto points = ConvertPoints(mesh);
        const auto triangles = ConvertTriangles(mesh);
        const std::vector<openvdb::Vec4I> quads;

        auto transform = openvdb::math::Transform::createLinearTransform(args.voxel_size);
        auto grid = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(
            *transform, points, triangles, quads, static_cast<float>(args.exterior_band_width),
            static_cast<float>(args.interior_band_width));
        grid->setGridClass(openvdb::GRID_LEVEL_SET);
        grid->setName(args.grid_name);

        std::cout << "Active voxels: " << grid->activeVoxelCount() << "\n";

        WriteVdb(grid, args.vdb_output);
        std::cout << "Wrote VDB:   " << args.vdb_output.string() << "\n";

        auto nano_handle = CreateNanoHandle(grid);
        nanovdb::io::writeGrid(args.nvdb_output.string(), nano_handle);
        ValidateNanoFile(args.nvdb_output, nano_handle);
        std::cout << "Wrote NVDB:  " << args.nvdb_output.string() << "\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "mesh_to_vdb_nvdb failed: " << e.what() << "\n";
        return 1;
    }
}
