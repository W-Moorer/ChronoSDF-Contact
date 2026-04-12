# nvdb_probe

Independent NanoVDB query tool for probing a sparse signed distance field at world-space coordinates.

## Configure

```powershell
cmake -S tools/nvdb_probe -B tools/nvdb_probe/build `
  -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake `
  -DVCPKG_TARGET_TRIPLET=x64-windows
```

## Build

```powershell
cmake --build tools/nvdb_probe/build --config Release
```

## Usage

```powershell
tools/nvdb_probe/build/Release/nvdb_probe.exe `
  --input path/to/model.nvdb `
  --point 0.10 0.10 0.10 `
  --point 0.50 0.50 0.50
```

Optional arguments:

- `--grid-name <name>`
- `--metadata-only`

## Output

For each query point the tool prints:

- world-space coordinate
- mapped index-space coordinate
- trilinearly sampled SDF value
- gradient in index space
- gradient in world space
- normalized world-space normal
- local zero-crossing flag
- bbox inclusion checks

## Notes

- `--point` coordinates are interpreted in world space.
- This tool expects a float NanoVDB grid, which matches the output of `mesh_to_vdb_nvdb`.
- For level-set grids, negative values are typically inside and positive values outside.
