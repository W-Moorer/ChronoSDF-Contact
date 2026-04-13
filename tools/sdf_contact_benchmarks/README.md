# SDF Contact Benchmarks

This tool runs a first benchmark package for the trimmed Chrono SDF contact pipeline.

Current scenarios:

- `centered_compression`: quasi-static penetration sweep for two aligned unit boxes
- `eccentric_compression`: quasi-static tilt sweep for a unit box compressed against a fixed unit box

Current methods:

- `distributed_sdf`
- `single_point_penalty`
- `grid_penalty_5x5`
- `dense_reference_161x161`

Outputs:

- `benchmark_records.csv`: one row per scenario / sample / method
- `benchmark_summary.csv`: RMSE and max-absolute-error summaries against the dense reference

The default NanoVDB asset is:

- `chrono/template_project/assets/unit_box_centered.nvdb`

Build example:

```powershell
cmake -S tools/sdf_contact_benchmarks -B tools/sdf_contact_benchmarks/build `
  -DChrono_DIR=E:/workspace/ChronoSDF-Contact/chrono/build_mbs_nanovdb/cmake `
  -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake `
  -DVCPKG_TARGET_TRIPLET=x64-windows
cmake --build tools/sdf_contact_benchmarks/build --config Release
```

Run example:

```powershell
tools/sdf_contact_benchmarks/build/Release/sdf_contact_benchmarks.exe `
  --input chrono/template_project/assets/unit_box_centered.nvdb `
  --output paper/results/benchmark_box_contact
```
