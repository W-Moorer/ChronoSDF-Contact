# Supplement 6 Mapping

This note maps `reference/supplement6.tex` to the current `ChronoSDF-Contact` codebase.
It focuses on the missing `D` and `E` stages:

- `D`: patch-wise band-sheet consistency bridge
- `E`: benchmark gate / acceptance workflow

The goal is to turn the current pipeline

`band mechanics -> sheet representation -> reporting`

into

`band mechanics -> sheet representation -> patch-wise consistency bridge -> corrected wrench -> benchmark gate`

## 1. Current State

The current code already has:

- `A`: preprocessing
  - `tools/mesh_to_vdb_nvdb`
  - `tools/nvdb_probe`
- `B`: runtime broadphase + band mechanics
  - `ChSDFBrickPair*`
  - `ChSDFVolumeContactEvaluator`
- `C`: sheet representation
  - `ChSDFSheetRepresentation`
- `E` partial: benchmark runners
  - `tools/sdf_contact_benchmarks`

The main missing piece is:

- `D`: patch-wise bridge from `A_band^m` / `F_band^m` / `tau_band^m`
  to `A_sheet^m` / `alpha_m` / corrected patch wrench

## 2. Target Theory

For each recovered sheet patch `m`:

- `A_band^m = sum_i DeltaA_i` over band samples assigned to patch `m`
- `A_sheet^m = recovered patch footprint area`
- `alpha_m = A_sheet^m / A_band^m`

Then correct patch wrench by:

- `F_m_corr = alpha_m * F_m_band`
- `tau_m_corr = alpha_m * tau_m_band`

Pair-level corrected wrench:

- `F_corr = sum_m F_m_corr`
- `tau_corr = sum_m tau_m_corr`

This is the minimum implementation required by Supplement 6.

## 3. New Data Structures

Add a new module:

- `chrono/src/chrono/collision/sdf/ChSDFPatchConsistency.h`
- `chrono/src/chrono/collision/sdf/ChSDFPatchConsistency.cpp`

Recommended structures:

### `ChSDFPatchBandAggregate`

Purpose:
- band-side aggregate for one recovered sheet patch

Fields:
- `std::size_t region_id`
- `std::size_t patch_id`
- `std::vector<std::size_t> source_sheet_sample_indices`
- `std::vector<std::size_t> source_band_sample_indices`
- `double band_area`
- `double integrated_pressure`
- `ChVector3d centroid_world`
- `ChVector3d pressure_center_world`
- `ChWrenchd wrench_world_b`
- `ChAABB support_bbox_world`

### `ChSDFPatchConsistencyResult`

Purpose:
- final bridge result for one patch

Fields:
- `std::size_t region_id`
- `std::size_t patch_id`
- `double band_area`
- `double sheet_area`
- `double alpha`
- `ChWrenchd wrench_world_b_band`
- `ChWrenchd wrench_world_b_corrected`
- `ChVector3d pressure_center_world`
- `ChAABB support_bbox_world`

### `ChSDFPatchConsistencySettings`

Purpose:
- control bridge behavior

Fields:
- `bool enable = true`
- `bool clamp_alpha = true`
- `double min_alpha = 0.25`
- `double max_alpha = 4.0`
- `double min_band_area = 1.0e-12`
- `bool use_patch_bbox_for_assignment = true`
- `bool use_sample_index_back_projection = true`

### `ChSDFPatchConsistencyPairResult`

Purpose:
- pair-level bridge result

Fields:
- `std::vector<ChSDFPatchConsistencyResult> patches`
- `ChWrenchd wrench_world_b_band`
- `ChWrenchd wrench_world_b_corrected`
- `double total_band_area`
- `double total_sheet_area`
- `double mean_alpha`
- `std::size_t corrected_patch_count`

## 4. Required Result Extensions

Extend `ChSDFShapePairContactResult` in:

- `chrono/src/chrono/collision/sdf/ChSDFContactWrench.h`

Add:

- `ChWrenchd wrench_world_b_band`
- `ChWrenchd wrench_world_b_corrected`
- `double corrected_active_area`
- `double corrected_integrated_pressure`
- `std::shared_ptr<ChSDFPatchConsistencyPairResult> patch_consistency_result`

Keep existing fields unchanged so current callers do not break.

Recommended policy:

- preserve current `wrench_world_b` as the existing band result until the bridge is validated
- add corrected fields separately first
- only after validation consider promoting corrected wrench to default

## 5. Band-to-Patch Assignment

This is the critical algorithmic step.

Input:
- one `ChSDFBrickPairWrenchResult` band region
- one `ChSDFSheetRegion`

Need:
- assign each active band sample to one sheet patch

Recommended staged implementation:

### Stage 1: sample-index back projection

Use the already preserved lineage:

- `ChSDFSheetFiberSample::source_sample_indices`
- `ChSDFSheetPatch::sample_indices`

Algorithm:
- for each patch
- gather all sheet sample indices in that patch
- union all `source_sample_indices`
- these become the band sample membership of the patch

This is the safest first implementation because it does not depend on extra geometry heuristics.

### Stage 2: fallback geometric assignment

If lineage is incomplete:

- assign a band sample to the nearest patch whose:
  - support bbox contains or nearly contains the sample point
  - normal compatibility is acceptable
  - tangent distance is minimal

This is only fallback logic, not the primary path.

## 6. Core Functions To Add

In `ChSDFPatchConsistency.cpp` add:

### `BuildPatchBandAggregates`

Signature sketch:

```cpp
std::vector<ChSDFPatchBandAggregate> BuildPatchBandAggregates(
    const ChSDFBrickPairWrenchResult& band_region,
    const ChSDFSheetRegion& sheet_region,
    const ChSDFPatchConsistencySettings& settings);
```

Responsibilities:
- iterate patches
- recover member band samples
- accumulate:
  - area
  - wrench
  - pressure center
  - bbox

### `BuildPatchConsistencyResult`

Signature sketch:

```cpp
ChSDFPatchConsistencyResult BuildPatchConsistencyResult(
    const ChSDFPatchBandAggregate& band_patch,
    const ChSDFSheetPatch& sheet_patch,
    const ChSDFPatchConsistencySettings& settings);
```

Responsibilities:
- compute `alpha`
- clamp if enabled
- scale patch wrench

### `BuildPairConsistencyResult`

Signature sketch:

```cpp
ChSDFPatchConsistencyPairResult BuildPairConsistencyResult(
    const ChSDFShapePairContactResult& band_result,
    const ChSDFSheetShapePairResult& sheet_result,
    const ChSDFPatchConsistencySettings& settings);
```

Responsibilities:
- match regions
- build patch aggregates
- compute corrected patch results
- sum pair-level corrected wrench

## 7. Region Matching Rule

Current `sheet_result.regions` are built from current band regions, so region matching should be direct:

- match by `region_id`
- secondarily verify `persistent_id` if needed

No global nearest-neighbor region matching should be introduced unless this direct mapping fails.

## 8. Integration Into `ChSDFShapePair`

Modify:

- `chrono/src/chrono/collision/sdf/ChSDFShapePair.h`
- `chrono/src/chrono/collision/sdf/ChSDFShapePair.cpp`

Add:

- `SetPatchConsistencySettings(...)`
- `GetPatchConsistencySettings()`

Recommended `EvaluateContact()` order:

1. `FindBrickPairs`
2. `EvaluateBrickPairs` -> band result
3. `ApplyActivationHysteresis`
4. `UpdateSheetResult`
5. `UpdatePatchConsistencyResult`
6. `UpdateHistory`

`UpdatePatchConsistencyResult(...)` should:

- build pair-level patch bridge result
- write:
  - `wrench_world_b_band`
  - `wrench_world_b_corrected`
  - `corrected_active_area`
  - `patch_consistency_result`

Do not change `Apply(...)` yet.
First validate by benchmark only.

## 9. Benchmark Changes

Modify:

- `tools/sdf_contact_benchmarks/src/main.cpp`

### Add new output columns

For main benchmark:
- `force_y_corrected`
- `torque_z_corrected`
- `area_corrected`
- `mean_patch_alpha`
- `patch_count_corrected`

For patch CSV if added later:
- `region_id`
- `patch_id`
- `band_area`
- `sheet_area`
- `alpha`
- `force_y_band`
- `force_y_corrected`
- `torque_z_band`
- `torque_z_corrected`

### Add corrected summary metrics

For centered/eccentric summary:
- `force_y_corr_rmse`
- `torque_z_corr_rmse`
- `area_corr_rmse`
- `alpha_mean`
- `alpha_rmse`

### Gate logic matching Supplement 6

Add explicit benchmark gate checks:

#### Centered gate
- `sheet_area_rmse` below threshold
- `pressure_center_x_rmse` below threshold
- `torque_z_corr_rmse` below threshold
- `force_y_corr_rmse` below threshold

#### Eccentric gate
- `area_sheet_rmse < area_band_rmse`
- `force_y_corr_rmse < force_y_band_rmse`
- `torque_z_corr_rmse < torque_z_band_rmse`

#### Cylinder gate
- `patch_count_rmse == 0`
- `area_sheet_rmse` below threshold
- `bbox_x/bbox_z rmse` below threshold

The thresholds should be written in code as explicit constants, not implied by discussion.

## 10. Planar Consistency Calibration Layer

Supplement 6 explicitly asks for a more formal planar consistency calibration.

Current benchmark uses:
- `potential.modulus = 2.22 * stiffness`

To satisfy the document better, add a helper:

- `ChSDFPotentialCalibration.h/.cpp`

Provide:

```cpp
double EstimateSymmetricFieldModulusForLinearPenalty(double reference_stiffness);
```

First implementation can still return the currently validated benchmark value,
but the calibration should become an explicit, named layer rather than an inline constant.

This keeps benchmark setup and theory mapping clearer.

## 11. Recommended Implementation Order

### Pass 1
- add patch consistency data structures
- add patch-wise lineage-based band aggregation
- compute `alpha_m`
- export corrected patch/pair wrench
- extend benchmark CSV

Goal:
- observe whether corrected `Fy/Tz` move in the expected direction

### Pass 2
- add benchmark gate logic
- add explicit planar calibration helper
- compare:
  - band
  - corrected
  - sheet geometry

Goal:
- make Supplement 6 executable as a pass/fail workflow

### Pass 3
- if corrected wrench is validated, optionally promote corrected wrench to default `Apply()`
- otherwise keep dual-output mode and continue improving bridge details

## 12. Files To Change

### New
- `chrono/src/chrono/collision/sdf/ChSDFPatchConsistency.h`
- `chrono/src/chrono/collision/sdf/ChSDFPatchConsistency.cpp`
- optionally `chrono/src/chrono/collision/sdf/ChSDFPotentialCalibration.h`
- optionally `chrono/src/chrono/collision/sdf/ChSDFPotentialCalibration.cpp`

### Modify
- `chrono/src/chrono/collision/sdf/ChSDFContactWrench.h`
- `chrono/src/chrono/collision/sdf/ChSDFShapePair.h`
- `chrono/src/chrono/collision/sdf/ChSDFShapePair.cpp`
- `chrono/src/chrono/CMakeLists.txt`
- `tools/sdf_contact_benchmarks/src/main.cpp`

## 13. Immediate Acceptance Criteria

After implementing Pass 1, the minimum expected observations are:

- centered:
  - corrected `Fy_rmse` should improve over band `Fy_rmse`
- eccentric:
  - corrected `Fy_rmse` and/or `Tz_rmse` should improve over band
- cylinder:
  - corrected bridge should not break `patch_count_rmse == 0`

If these are not true, then the bridge implementation is incomplete or patch assignment is wrong.

## 14. Bottom Line

To satisfy `supplement6`, the current project must stop treating `sheet` as output-only geometry.

The required change is:

- recover patch geometry from sheet
- compute patch-wise `alpha_m = A_sheet^m / A_band^m`
- feed that back into patch wrench
- benchmark the corrected wrench explicitly

That is the missing closure layer in the current codebase.
