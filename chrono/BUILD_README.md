# Chrono Minimal MBS Build Guide

This fork is trimmed for the core Chrono multibody dynamics framework on Windows.

## Environment

- CMake: 3.31 or newer
- Visual Studio: 2022
- Generator: `Visual Studio 17 2022`
- Architecture: `x64`
- Eigen: 3.4.0 or newer

The local workspace already includes Eigen at:

```text
E:/workspace/ChronoSDF-Contact/eigen-3.4.0
```

## Default Build Behavior

`chrono/src/CMakeLists.txt` now enables `CHRONO_BUILD_MINIMAL_MBS=ON` by default.

That means the default configuration:

- builds `Chrono_core`
- skips optional modules
- skips demos, tests, and benchmarks
- skips MPI, CUDA, and Thrust discovery

## Configure

```powershell
cmake -S E:/workspace/ChronoSDF-Contact/chrono `
      -B E:/workspace/ChronoSDF-Contact/chrono/build_mbs_clean `
      -G "Visual Studio 17 2022" `
      -A x64 `
      -DEIGEN3_INCLUDE_DIR="E:/workspace/ChronoSDF-Contact/eigen-3.4.0"
```

## Build

```powershell
cmake --build E:/workspace/ChronoSDF-Contact/chrono/build_mbs_clean `
      --config Release `
      --target Chrono_core `
      -j 8
```

## Optional Install

```powershell
cmake --install E:/workspace/ChronoSDF-Contact/chrono/build_mbs_clean `
      --config Release `
      --prefix E:/workspace/ChronoSDF-Contact/chrono/install_mbs_clean
```

## Optional NanoVDB SDF Utilities

The trimmed core now includes optional NanoVDB-based sparse signed-distance utilities under `chrono/collision/sdf`.

To enable them on Windows with `vcpkg`:

```powershell
C:\vcpkg\vcpkg.exe install openvdb[nanovdb]:x64-windows

cmake -S E:/workspace/ChronoSDF-Contact/chrono `
      -B E:/workspace/ChronoSDF-Contact/chrono/build_mbs_nanovdb `
      -G "Visual Studio 17 2022" `
      -A x64 `
      -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake `
      -DVCPKG_TARGET_TRIPLET=x64-windows `
      -DCH_ENABLE_NANOVDB_SDF=ON
```

Then build:

```powershell
cmake --build E:/workspace/ChronoSDF-Contact/chrono/build_mbs_nanovdb `
      --config Release `
      --target Chrono_core
```

## Scope

This is now a hard-trimmed fork. `CHRONO_BUILD_MINIMAL_MBS=OFF` is no longer supported in this workspace.
