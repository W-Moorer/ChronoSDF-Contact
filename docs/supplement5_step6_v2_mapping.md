# Supplement5 Step 6 v2 Mapping for Current Chrono Codebase

## 1. Goal

This document maps `reference/supplement5.tex` onto the current
`ChronoSDF-Contact` implementation.

Supplement5 does not change `Band mechanics`.
It rewrites only the `Sheet representation` layer.

The central upgrade is:

- v1:
  region-level projection plane -> quantized `(u,v)` bins -> per-bin accumulators -> lattice BFS
- v2:
  band sample -> local sheet seed -> 3D geometric fiber clustering -> measure-preserving collapse -> world-space patch graph

In the current codebase, this means:

- keep `ChSDFVolumeContactEvaluator` unchanged as the force path
- replace the core of `ChSDFSheetBuilder`
- keep `ChSDFShapePair::EvaluateContact()` orchestration unchanged except for consuming the new sheet builder output

## 2. Current State and Why It Fails

Current Step 6 implementation lives in:

- `chrono/src/chrono/collision/sdf/ChSDFSheetRepresentation.h`
- `chrono/src/chrono/collision/sdf/ChSDFSheetRepresentation.cpp`

Current v1 behavior:

- one region-wide `fallback_normal`
- one region-wide tangent basis
- one region-wide `projection_origin`
- projected `(u,v)` quantization via integer lattice bins
- multiple accumulators allowed in the same projected bin
- `footprint_area = lateral_tolerance^2`
- patch BFS on the projected integer lattice

This is exactly the representation that Supplement5 replaces.

Observed failure modes already confirmed by the current benchmark:

- curved strips fragment into many `sheet patches`
- `area_sheet` depends on the number of collapsed fibers
- `support_bbox_world` does not really reflect sheet geometry
- `local_fiber` can perform worse than `dominant_axis`

## 3. Supplement5 Objects Mapped to Current Data

Supplement5 uses the following logical objects:

- `band sample`
- `local sheet seed`
- `fiber cluster`
- `collapsed sheet point`
- `patch`

Mapped to current code:

- `band sample`
  - current type: `ChSDFBrickPairWrenchSample`
  - source fields:
    - `region_sample.point_world`
    - `region_sample.contact_normal_world`
    - `quadrature_area`
    - `force_world`
    - `pressure`
    - `region_sample.h_value`
    - `region_sample.grad_h_world`
- `local sheet seed`
  - new type to add: `ChSDFSheetSeed`
- `fiber cluster`
  - internal clustering object, may remain implementation-private
- `collapsed sheet point`
  - current public type to reuse and redefine: `ChSDFSheetFiberSample`
- `patch`
  - current public type to keep: `ChSDFSheetPatch`

The important semantic change is:

- `ChSDFSheetFiberSample` should no longer mean "one projected lattice bin"
- it should mean "one measure-preserving collapsed point from one unique 3D fiber cluster"

## 4. Required Public Type Changes

File:

- `chrono/src/chrono/collision/sdf/ChSDFSheetRepresentation.h`

### 4.1 Extend `ChSDFSheetCollapseSettings`

Keep:

- `enable`
- `min_sheet_sample_area`

Deprecate or mark legacy:

- `use_local_fiber_projection`
- `fiber_lateral_tolerance`
- `patch_neighbor_mode`

Add v2 settings:

- `double fiber_normal_cosine = cos(20 deg)`
- `double fiber_tangent_tolerance = -1`
- `double fiber_plane_tolerance = -1`
- `double patch_connection_radius = -1`
- `double patch_normal_cosine = cos(25 deg)`
- `double min_gradient_norm = 1.0e-8`
- `bool enable_sheet_diagnostics = true`
- `bool allow_dominant_axis_fallback = true`
- `double min_sheet_area_ratio = 0.5`
- `double max_sheet_area_ratio = 1.5`
- `std::size_t max_patch_count_before_fallback = 0`

Resolution rule:

- if `fiber_tangent_tolerance <= 0`, use `1.0 * spacing`
- if `fiber_plane_tolerance <= 0`, use `0.5 * spacing`
- if `patch_connection_radius <= 0`, use `2.0 * spacing`

### 4.2 Add `ChSDFSheetSeed`

Add a new public or file-local struct:

- `std::size_t region_id`
- `std::size_t source_sample_index`
- `ChVector3d seed_world`
- `ChVector3d seed_normal_world`
- `double measure_area`
- `ChVector3d force_world`
- `double pressure`
- `double h_value`
- `ChVector3d grad_h_world`
- `ChVector3d band_point_world`
- `ChAABB source_bounds_world`

Semantics:

- one band sample projected onto its own local first-order `h=0` sheet point

### 4.3 Redefine `ChSDFSheetFiberSample`

Keep the type name to minimize API churn.

Change semantics:

- one collapsed point produced from one unique fiber cluster

Keep fields:

- `measure_area`
- `centroid_world`
- `pressure_center_world`
- `normal_world`
- `force_world`
- `source_sample_indices`

Change the meaning of:

- `footprint_area`
  - no longer fixed `h^2`
  - in v2 it should default to `measure_area`
  - later it can be separated again when support polygon recovery exists

Add:

- `std::size_t support_seed_count = 0`
- `ChAABB support_bbox_world`

### 4.4 Extend `ChSDFSheetPatch`

Keep:

- `patch_id`
- `measure_area`
- `centroid_world`
- `pressure_center_world`
- `mean_normal_world`
- `bounds_world`

Change:

- `support_bbox_world`
  - compute from collapsed sheet points, not source band bounds

Add:

- `std::size_t support_seed_count = 0`

### 4.5 Extend `ChSDFSheetRegion` and `ChSDFSheetShapePairResult`

Add diagnostics:

- `double area_ratio = 0`
- `std::size_t fiber_count = 0`
- `double mean_support_seed_count = 0`
- `double normal_spread = 0`
- `bool used_fallback = false`

Reason:

- Supplement5 explicitly asks for diagnostics and fallback policy

## 5. Functions to Retire or Rewrite

File:

- `chrono/src/chrono/collision/sdf/ChSDFSheetRepresentation.cpp`

### 5.1 Retire v1 projection/bin logic

These functions should be removed from the v2 main path:

- `ResolveLateralTolerance(...)`
- `QuantizeProjectedCoord(...)`
- `ComputeDominantAxis(...)`
- `ComputeLatticeCoord(...)`
- the current lattice-based `BuildPatches(...)`

They may remain as:

- legacy fallback implementation
- benchmark comparison path

But they should not remain the default Step 6 path.

### 5.2 Rewrite `BuildRegion(...)`

Current:

- `BuildRegion(...)` directly projects samples into lattice bins and collapses per bin

Target:

- `BuildRegion(...)` should become a staged v2 pipeline:

1. `BuildSeeds(...)`
2. `BuildFiberClusters(...)`
3. `CollapseFiberClusters(...)`
4. `BuildPatchGraph(...)`
5. `FinalizeRegion(...)`

## 6. New Internal Algorithm Stages

### 6.1 `BuildSeeds(const ChSDFBrickPairWrenchResult&)`

Input:

- active `band_region.samples`

For each active sample:

- read
  - `x_i = region_sample.point_world`
  - `h_i = region_sample.h_value`
  - `grad_h_i = region_sample.grad_h_world`
  - `A_i = quadrature_area`
  - `f_i = force_world`
- compute local seed:
  - `x_i_sigma = x_i - h_i / ||grad_h_i||^2 * grad_h_i`
- if `||grad_h_i|| < min_gradient_norm`
  - fallback to `x_i_sigma = x_i`
  - mark low-quality seed if needed

Output:

- `std::vector<ChSDFSheetSeed>`

This directly implements Supplement5 Eq. `x_i^Sigma`.

### 6.2 `BuildFiberClusters(const std::vector<ChSDFSheetSeed>&, settings)`

This is the core replacement for the old binning logic.

Do not build projected bins.

Instead:

- insert seeds into a spatial hash grid or kd-tree using `seed_world`
- for each seed, search neighbors within radius `fiber_tangent_tolerance + fiber_plane_tolerance`
- connect two seeds when all hold:
  - `dot(n_i, n_j) >= fiber_normal_cosine`
  - tangential distance `<= fiber_tangent_tolerance`
  - normal-direction distance `<= fiber_plane_tolerance`

Distance decomposition:

- `n_bar = normalize(n_i + n_j)`
- `d = x_j_sigma - x_i_sigma`
- `d_n = abs(dot(n_bar, d))`
- `d_t = ||(I - n_bar n_bar^T) d||`

Cluster construction:

- union-find is preferred
- BFS/DFS on the implicit graph is also acceptable

Output:

- `std::vector<std::vector<std::size_t>>` or an equivalent cluster index map

### 6.3 `CollapseFiberClusters(...)`

For each cluster `F_alpha`:

- `A_alpha = sum A_i`
- `x_alpha = sum(A_i * x_i_sigma) / A_alpha`
- `n_alpha = normalize(sum(A_i * n_i_sigma))`
- `f_alpha = sum f_i`
- `p_alpha = sum(A_i * p_i) / A_alpha`

This stage must be measure-preserving.

Required semantics:

- `sheet_area = sum A_alpha = sum A_i`
- collapse must not depend on cluster count

Public mapping:

- one collapsed cluster becomes one `ChSDFSheetFiberSample`

### 6.4 `BuildPatchGraph(...)`

This replaces current lattice BFS.

Patch adjacency must be world-space.

Connect two collapsed sheet points `alpha`, `beta` when:

- `||x_beta - x_alpha|| <= patch_connection_radius`
- `dot(n_alpha, n_beta) >= patch_normal_cosine`

Optional extra filter:

- tangential separation limit in the average patch tangent plane

Cluster connected components to form `ChSDFSheetPatch`.

### 6.5 `FinalizeRegion(...)`

Compute:

- `measure_area`
- `footprint_area`
- `centroid_world`
- `pressure_center_world`
- `mean_normal_world`
- `support_bbox_world`
- `patch_count`
- `largest_patch_area`
- diagnostics

Crucial bbox rule:

- `support_bbox_world` must be computed from collapsed point positions
- do not inherit `region.bounds_world` from band samples

## 7. Pair-Level Aggregation Changes

File:

- `chrono/src/chrono/collision/sdf/ChSDFSheetRepresentation.cpp`

Function:

- `ChSDFSheetBuilder::BuildShapePair(...)`

Required changes:

- `sheet_area = sum(region.measure_area)`
- `sheet_footprint_area = sum(region.footprint_area)`
- `sheet_center_world = area-weighted centroid of collapsed points`
- `pressure_center_world = force-weighted average of collapsed pressure centers`
- `support_bbox_world = AABB(collapsed point positions)`
- `patch_count = sum(region patch_count)`

Important semantic correction:

- pair-level bbox must no longer be inherited from band region bounds

## 8. `ChSDFShapePair` Integration

Files:

- `chrono/src/chrono/collision/sdf/ChSDFShapePair.h`
- `chrono/src/chrono/collision/sdf/ChSDFShapePair.cpp`

No dataflow rewrite is needed.

Keep:

- `EvaluateContact()`
  - `FindBrickPairs -> EvaluateBrickPairs -> ApplyActivationHysteresis -> UpdateSheetResult -> UpdateHistory`
- `Apply()`
  - still uses only band wrench

Only the internal behavior of:

- `BuildSheetRepresentation(...)`

needs to change because it calls the rewritten `ChSDFSheetBuilder`.

## 9. Benchmark Output Changes

File:

- `tools/sdf_contact_benchmarks/src/main.cpp`

No schema rewrite is required for the existing CSVs.

But Supplement5 diagnostics should be added for curved-sheet debugging:

- `sheet_fiber_count`
- `sheet_mean_support_seed_count`
- `sheet_area_ratio = area_sheet / area_band`
- `sheet_used_fallback`
- `sheet_normal_spread`

This should be added first to:

- `curved_sheet_records.csv`
- `curved_sheet_summary.csv`

Reason:

- the curved benchmark is currently the fastest way to validate Step 6 v2

## 10. Suggested Rollout Plan

### Step A

Add new data types and diagnostics fields without changing the active algorithm.

### Step B

Implement `BuildSeeds(...)` and export a debug path:

- `band sample count`
- `seed count`
- `mean projection distance`

### Step C

Implement `BuildFiberClusters(...)` with spatial hash + union-find.

Do not remove v1 yet.

### Step D

Implement `CollapseFiberClusters(...)` and switch:

- `footprint_area = measure_area`

for the v2 path.

### Step E

Implement world-space `BuildPatchGraph(...)`.

### Step F

Add diagnostics and automatic fallback:

- if patch count explodes
- if `sheet_area / band_area` is out of range
- if mean support count is too low

then fall back to the current dominant-axis path.

### Step G

Make v2 the default path and keep v1 only as explicit legacy mode.

## 11. Minimal Concrete API Direction

To minimize churn, the recommended public API strategy is:

- keep `ChSDFSheetBuilder::BuildRegion(...)`
- keep `ChSDFSheetBuilder::BuildShapePair(...)`
- keep existing public result types
- replace only their internal semantics and add fields

Recommended private helper additions in `ChSDFSheetRepresentation.cpp`:

- `BuildSeeds(...)`
- `BuildSeedSpatialHash(...)`
- `BuildFiberClusters(...)`
- `CollapseFiberClusters(...)`
- `BuildPatchAdjacency(...)`
- `FinalizePatch(...)`
- `FinalizeRegionDiagnostics(...)`
- `ShouldFallbackToDominantAxis(...)`
- `BuildRegionLegacyV1(...)`
- `BuildRegionV2(...)`

This preserves downstream code while allowing Step 6 v2 to become the new default.

## 12. Final Mapping Summary

Supplement5's direct message to the current codebase is:

- stop treating `fiber` as a projected bin
- stop treating `sheet area` as `number of bins * h^2`
- stop building `patch` on a quantized plane lattice
- build local sheet seeds first
- cluster fibers in 3D
- collapse with measure conservation
- build patch connectivity in world space
- compute bbox from collapsed sheet points

That is the smallest change set that is actually faithful to Supplement5.
