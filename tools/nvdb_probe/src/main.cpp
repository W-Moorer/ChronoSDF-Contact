#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <nanovdb/NanoVDB.h>
#include <nanovdb/io/IO.h>
#include <nanovdb/math/SampleFromVoxels.h>

namespace fs = std::filesystem;

namespace {

struct ProbePoint {
    nanovdb::Vec3d world;
};

struct Arguments {
    fs::path input_path;
    std::string grid_name;
    std::vector<ProbePoint> points;
    bool metadata_only = false;
};

void PrintUsage() {
    std::cout
        << "Usage:\n"
           "  nvdb_probe --input <grid.nvdb> [--grid-name <name>] --point <x> <y> <z> [--point <x> <y> <z> ...]\n"
           "  nvdb_probe --input <grid.nvdb> [--grid-name <name>] --metadata-only\n\n"
           "Notes:\n"
           "  - Coordinates passed to --point are interpreted in world space.\n"
           "  - The tool performs trilinear sampling and reports value plus gradient.\n";
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

ProbePoint ParsePoint(int argc, char* argv[], int& index) {
    const double x = ParseDouble(RequireValue("--point", argc, argv, index), "--point.x");
    const double y = ParseDouble(RequireValue("--point", argc, argv, index), "--point.y");
    const double z = ParseDouble(RequireValue("--point", argc, argv, index), "--point.z");
    return ProbePoint{nanovdb::Vec3d(x, y, z)};
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
        } else if (option == "--grid-name") {
            args.grid_name = RequireValue(option, argc, argv, i);
        } else if (option == "--point") {
            args.points.push_back(ParsePoint(argc, argv, i));
        } else if (option == "--metadata-only") {
            args.metadata_only = true;
        } else {
            throw std::runtime_error("Unknown option: " + option);
        }
    }

    if (args.input_path.empty()) {
        throw std::runtime_error("Missing required option --input");
    }
    if (!fs::exists(args.input_path)) {
        throw std::runtime_error("Input file does not exist: " + args.input_path.string());
    }
    if (!args.metadata_only && args.points.empty()) {
        throw std::runtime_error("At least one --point is required unless --metadata-only is used.");
    }

    return args;
}

template <typename Vec3T>
double Norm(const Vec3T& v) {
    const double x = static_cast<double>(v[0]);
    const double y = static_cast<double>(v[1]);
    const double z = static_cast<double>(v[2]);
    return std::sqrt(x * x + y * y + z * z);
}

template <typename Vec3T>
Vec3T NormalizeOrZero(const Vec3T& v) {
    using Scalar = std::remove_cv_t<std::remove_reference_t<decltype(v[0])>>;
    const double norm = Norm(v);
    if (norm <= 0.0) {
        return Vec3T(Scalar(0));
    }
    return Vec3T(Scalar(v[0] / norm), Scalar(v[1] / norm), Scalar(v[2] / norm));
}

template <typename Vec3T>
std::string FormatVec3(const Vec3T& v, int precision = 6) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << "(" << static_cast<double>(v[0]) << ", "
        << static_cast<double>(v[1]) << ", " << static_cast<double>(v[2]) << ")";
    return out.str();
}

std::string GridClassName(const nanovdb::GridClass grid_class) {
    switch (grid_class) {
        case nanovdb::GridClass::LevelSet:
            return "LevelSet";
        case nanovdb::GridClass::FogVolume:
            return "FogVolume";
        case nanovdb::GridClass::VoxelVolume:
            return "VoxelVolume";
        case nanovdb::GridClass::PointIndex:
            return "PointIndex";
        case nanovdb::GridClass::PointData:
            return "PointData";
        case nanovdb::GridClass::Topology:
            return "Topology";
        case nanovdb::GridClass::IndexGrid:
            return "IndexGrid";
        case nanovdb::GridClass::Staggered:
            return "Staggered";
        default:
            return "Unknown";
    }
}

void PrintMetadata(const nanovdb::GridMetaData& meta) {
    std::cout << "Grid name:           " << meta.shortGridName() << "\n";
    std::cout << "Grid class:          " << GridClassName(meta.gridClass()) << "\n";
    std::cout << "Active voxels:       " << meta.activeVoxelCount() << "\n";
    std::cout << "Index bbox min/max:  " << FormatVec3(meta.indexBBox().min()) << " / "
              << FormatVec3(meta.indexBBox().max()) << "\n";
    std::cout << "World bbox min/max:  " << FormatVec3(meta.worldBBox().min()) << " / "
              << FormatVec3(meta.worldBBox().max()) << "\n";
    std::cout << "Voxel size:          " << FormatVec3(meta.voxelSize()) << "\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Arguments args = ParseArguments(argc, argv);

        auto handle = args.grid_name.empty() ? nanovdb::io::readGrid(args.input_path.string())
                                             : nanovdb::io::readGrid(args.input_path.string(), args.grid_name);

        const auto* grid = handle.grid<float>();
        if (grid == nullptr) {
            throw std::runtime_error("The requested grid is not a float NanoVDB grid.");
        }

        const nanovdb::GridMetaData meta(*grid);
        PrintMetadata(meta);
        if (!grid->isLevelSet()) {
            std::cout << "Warning: grid class is not LevelSet. Distances and normals may not be meaningful.\n";
        }

        if (args.metadata_only) {
            return 0;
        }

        const auto accessor = grid->getAccessor();
        auto sampler = nanovdb::math::createSampler<1>(accessor);
        const auto index_bbox_real = meta.indexBBox().template asReal<double>();

        for (std::size_t i = 0; i < args.points.size(); ++i) {
            const auto& world = args.points[i].world;
            const auto index = grid->worldToIndexF(world);
            const float value = sampler(index);
            const auto gradient_index = sampler.gradient(index);
            const auto gradient_world = grid->indexToWorldGradF(gradient_index);
            const auto normal_world = NormalizeOrZero(gradient_world);
            const bool zero_crossing = sampler.zeroCrossing(index);
            const bool inside_world_bbox = meta.worldBBox().isInside(world);
            const bool inside_index_bbox = index_bbox_real.isInside(index);

            std::cout << "\nProbe " << (i + 1) << "\n";
            std::cout << "  world xyz:         " << FormatVec3(world) << "\n";
            std::cout << "  index xyz:         " << FormatVec3(index) << "\n";
            std::cout << "  sdf value:         " << std::fixed << std::setprecision(6) << value << "\n";
            std::cout << "  grad (index):      " << FormatVec3(gradient_index) << "\n";
            std::cout << "  grad (world):      " << FormatVec3(gradient_world) << "\n";
            std::cout << "  |grad world|:      " << std::fixed << std::setprecision(6) << Norm(gradient_world) << "\n";
            std::cout << "  normal (world):    " << FormatVec3(normal_world) << "\n";
            std::cout << "  zero crossing:     " << (zero_crossing ? "true" : "false") << "\n";
            std::cout << "  in world bbox:     " << (inside_world_bbox ? "true" : "false") << "\n";
            std::cout << "  in index bbox:     " << (inside_index_bbox ? "true" : "false") << "\n";
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "nvdb_probe failed: " << e.what() << "\n";
        return 1;
    }
}
