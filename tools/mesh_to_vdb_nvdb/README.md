# mesh_to_vdb_nvdb

Independent preprocessing tool for:

`triangle mesh -> narrow-band OpenVDB level set (.vdb) -> NanoVDB runtime grid (.nvdb)`

This tool is intentionally decoupled from the trimmed Chrono core build. It reuses the vendored mesh readers already present in this repository and links directly against OpenVDB/NanoVDB.

## Inputs

- `.obj`
- `.stl` (binary STL only)

## Outputs

- `.vdb`: authoring/debug version of the sparse narrow-band signed distance field
- `.nvdb`: lightweight runtime version intended for later collision/query integration

## Dependencies

One practical Windows path is `vcpkg`:

```powershell
C:\vcpkg\vcpkg.exe install openvdb[nanovdb]:x64-windows
```

## Configure

```powershell
cmake -S tools/mesh_to_vdb_nvdb -B tools/mesh_to_vdb_nvdb/build `
  -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake `
  -DVCPKG_TARGET_TRIPLET=x64-windows
```

## Build

```powershell
cmake --build tools/mesh_to_vdb_nvdb/build --config Release
```

## Usage

```powershell
tools/mesh_to_vdb_nvdb/build/Release/mesh_to_vdb_nvdb.exe `
  --input path/to/model.obj `
  --voxel-size 0.002 `
  --band-width 4
```

Optional arguments:

- `--vdb-out <path>`
- `--nvdb-out <path>`
- `--grid-name <name>`
- `--ex-band <voxels>`
- `--in-band <voxels>`
- `--scale <factor>`
- `--allow-open`

## Notes

- Band widths are expressed in voxel units because that matches the OpenVDB API.
- Without `--allow-open`, the tool rejects meshes with open boundary edges.
- Non-manifold edges and degenerate faces are still reported because they are useful diagnostics, but OpenVDB can often process them.
- For reliable signed distance fields, use watertight triangle meshes whenever possible.
