# Supplement4 Dual-Representation Mapping for Current Chrono Codebase

## 1. Goal

This document maps `reference/supplement4.tex` onto the current
`ChronoSDF-Contact` codebase.

The supplement's core claim is:

- `Band mechanics` remains the force-and-wrench computation path
- `Sheet representation` is added as a separate geometric/topological layer
- the two layers must not be conflated

In the current codebase, the force path already exists. The missing part is a
first-class `Sheet` layer and a clean split between:

- occupied-cell area
- band measure area
- collapsed sheet area

The intended outcome is:

- keep the current `Band` force calculation stable
- add a `Sheet` collapse layer without changing the force law
- expose `A_occ`, `A_band`, `A_sheet`, footprint, patch, and center outputs
  with correct semantics

## 2. Current Code Mapping

### 2.1 Broadphase and common grid

Current files:

- `chrono/src/chrono/collision/sdf/ChSDFBrickPair.h`
- `chrono/src/chrono/collision/sdf/ChSDFBrickPair.cpp`

Current responsibility:

- sparse leaf-brick broadphase
- world-space overlap AABB construction
- coarse candidate brick pairs

This already corresponds to the supplement's pre-band candidate search layer.

### 2.2 Band mechanics

Current files:

- `chrono/src/chrono/collision/sdf/ChSDFVolumeContactEvaluator.h`
- `chrono/src/chrono/collision/sdf/ChSDFVolumeContactEvaluator.cpp`
- `chrono/src/chrono/collision/sdf/ChSDFContactWrench.h`

Current responsibility:

- build a common integration grid from overlapping brick pairs
- evaluate active cells with `2x2x2` quadrature
- compute
  - `quadrature_area`
  - `force_world`
  - `point_world`
  - `contact_normal_world`
  - `pressure`
- aggregate to region and pair wrenches

This is already the supplement's `Band mechanics` layer.

Current struct mapping:

- `ChSDFBrickPairWrenchSample` = one band cell contribution
- `ChSDFBrickPairWrenchResult` = one band region
- `ChSDFShapePairContactResult` = pair-level band mechanics aggregate

### 2.3 Shape-pair orchestration

Current files:

- `chrono/src/chrono/collision/sdf/ChSDFShapePair.h`
- `chrono/src/chrono/collision/sdf/ChSDFShapePair.cpp`

Current responsibility:

- `FindBrickPairs()`
- `EvaluateContact()`
- `Apply()`
- activation hysteresis and region history

This is the right top-level location for supplement4's dual-layer orchestration.

### 2.4 Benchmark consumer

Current file:

- `tools/sdf_contact_benchmarks/src/main.cpp`

Current responsibility:

- run centered/eccentric box benchmarks
- compare distributed SDF against penalty/grid/reference baselines
- export CSV summaries

This is the right place to expose supplement4's multi-area semantics.

## 3. What Supplement4 Adds

Supplement4 does not replace the current band evaluator.

It adds a second layer:

- `Band` answers:
  - force
  - torque
  - pressure-weighted wrench quantities
- `Sheet` answers:
  - zero-thickness area
  - footprint
  - patch topology
  - support region
  - center-like geometric quantities

Therefore the current code should be split as:

- keep `ChSDFVolumeContactEvaluator` as the source of band mechanics
- add a new `Sheet` builder that consumes band-region samples
- expose both layers in the top-level result

## 4. Recommended Type-Level Design

## 4.1 Reuse existing band structs

Do not create a brand-new band representation file yet.

The current structs are already close to the supplement's notation:

- `ChSDFBrickPairWrenchSample`
  - `quadrature_area` maps to `ΔA_c`
  - `force_world` maps to `Δf_c`
  - `region_sample.point_world` maps to `x_c^band`
  - `region_sample.contact_normal_world` maps to `n_c`
  - `pressure` and `local_stiffness` carry local mechanics metadata
- `ChSDFBrickPairWrenchResult`
  - region-level band aggregate
- `ChSDFShapePairContactResult`
  - pair-level band aggregate

Required change:

- document these structs as the canonical `Band layer` outputs
- stop treating `active_area` as generic "contact area"
- treat it explicitly as `A_band`

## 4.2 Add a first-class sheet representation file

New files:

- `chrono/src/chrono/collision/sdf/ChSDFSheetRepresentation.h`
- `chrono/src/chrono/collision/sdf/ChSDFSheetRepresentation.cpp`

Recommended types:

- `struct ChSDFSheetCollapseSettings`
- `struct ChSDFSheetFiberSample`
- `struct ChSDFSheetPatch`
- `struct ChSDFSheetRegion`
- `struct ChSDFSheetShapePairResult`
- `class ChSDFSheetBuilder`

Recommended field set:

### `ChSDFSheetCollapseSettings`

- `bool enable = true`
- `bool use_local_fiber_projection = false`
- `double fiber_lateral_tolerance = -1`
- `double fiber_normal_cosine = 0.7`
- `double min_sheet_sample_area = 0`
- `int patch_neighbor_mode = 8`

Minimal behavior:

- Stage 1:
  region-level dominant-axis column collapse using integer grid coordinates
- Stage 2:
  local-fiber collapse using projected tangent-plane coordinates

### `ChSDFSheetFiberSample`

- `std::size_t region_id`
- `std::vector<std::size_t> source_sample_indices`
- `double area`
- `ChVector3d centroid_world`
- `ChVector3d normal_world`
- `ChVector3d force_world`
- `ChVector3d pressure_center_world`
- `ChAABB source_bounds_world`
- `int dominant_axis`
- `ChVector2i lattice_coord`

Semantics:

- one collapsed sheet sample
- measure-preserving merge of one fiber's band cells

### `ChSDFSheetPatch`

- `std::size_t patch_id`
- `double area`
- `ChVector3d centroid_world`
- `ChVector3d mean_normal_world`
- `ChAABB bounds_world`
- `std::vector<std::size_t> sample_indices`

Semantics:

- one connected component on the collapsed sheet

### `ChSDFSheetRegion`

- `std::size_t region_id`
- `std::size_t persistent_id`
- `double area`
- `ChVector3d centroid_world`
- `ChVector3d mean_normal_world`
- `ChAABB bounds_world`
- `std::vector<ChSDFSheetFiberSample> samples`
- `std::vector<ChSDFSheetPatch> patches`

Semantics:

- sheet counterpart of one band region

### `ChSDFSheetShapePairResult`

- `double occupied_area`
- `double band_area`
- `double sheet_area`
- `ChVector3d sheet_center_world`
- `std::size_t occupied_cells`
- `std::size_t collapsed_samples`
- `std::vector<ChSDFSheetRegion> regions`

Semantics:

- pair-level geometry layer output

## 4.3 Extend the top-level contact result

File:

- `chrono/src/chrono/collision/sdf/ChSDFContactWrench.h`

Recommended additions to `ChSDFShapePairContactResult`:

- `double occupied_area = 0`
- `double band_area = 0`
- `double sheet_area = 0`
- `ChVector3d sheet_center_world = VNULL`
- `ChSDFSheetShapePairResult sheet`

Compatibility rule:

- keep `active_area` as a backward-compatible alias of `band_area`

This preserves existing downstream users while making the semantics explicit.

## 5. Recommended Data Flow

## 5.1 Current mainline

Current:

`FindBrickPairs -> EvaluateBrickPairs -> ApplyActivationHysteresis -> UpdateHistory`

## 5.2 Target supplement4 mainline

Target:

`FindBrickPairs -> EvaluateBand -> BuildSheet -> Copy band/sheet summaries -> ApplyActivationHysteresis -> UpdateHistory`

Mapped to current files:

1. `ChSDFShapePair::FindBrickPairs()`
2. `ChSDFVolumeContactEvaluator::EvaluateBrickPairs()`
   returns band mechanics aggregate
3. `ChSDFSheetBuilder::BuildFromBandResult(...)`
   consumes the band result and constructs the sheet layer
4. `ChSDFShapePair::EvaluateContact()`
   attaches sheet outputs to `ChSDFShapePairContactResult`
5. `ChSDFShapePair::Apply()`
   still uses only top-level band wrench fields

## 5.3 Why `Apply()` should stay unchanged

Supplement4 is explicit:

- `Band` remains the mechanical truth for force integration
- `Sheet` is not a replacement force path

Therefore:

- `Apply()` should continue using `wrench_world_a` and `wrench_world_b`
- sheet outputs are for geometry, benchmarking, reduction, and patch logic

## 6. Sheet Collapse Algorithm Mapped to Current Data

## 6.1 Minimal implementation

Input:

- `ChSDFBrickPairWrenchResult region`
- its `region.samples`

Fiber key in the minimal version:

- choose `dominant_axis = argmax(|region.region.mean_normal_world|)`
- collapse all samples with the same two orthogonal integer coordinates

Examples:

- dominant axis `Y`
  - fiber key = `(coord.x, coord.z)`
- dominant axis `X`
  - fiber key = `(coord.y, coord.z)`
- dominant axis `Z`
  - fiber key = `(coord.x, coord.y)`

Per-fiber collapse:

- `A_tilde = sum(sample.quadrature_area)`
- `x_tilde = sum(A_i * point_i) / A_tilde`
- `n_tilde = normalize(sum(A_i * normal_i))`
- optional `f_tilde = sum(force_world_i)`

This directly implements supplement4's measure-preserving collapse while
keeping the current force computation untouched.

## 6.2 Patch formation on sheet samples

After collapse:

- build adjacency in the 2D lattice induced by the dominant-axis plane
- cluster neighboring fiber samples
- produce `ChSDFSheetPatch`

This is the right place to derive:

- footprint
- patch count
- support region
- sheet center

## 6.3 Future upgrade path

The minimal dominant-axis column collapse is enough for:

- centered compression
- flat or nearly flat contact
- immediate benchmark diagnosis

Later upgrades:

- replace dominant-axis columns with local tangent-plane fibers
- use per-sample normals and projected tangent coordinates
- allow curved-sheet clustering

## 7. File-Level Change List

## 7.1 `ChSDFVolumeContactEvaluator.*`

Keep:

- current band mechanics algorithm
- current force and torque accumulation

Add:

- comments clarifying that this class produces `Band layer` outputs
- optional helper accessors if needed for occupied-cell counts

Do not:

- move collapse into this file
- mix geometry collapse logic into the force evaluator

## 7.2 `ChSDFContactWrench.h`

Extend:

- `ChSDFShapePairContactResult`
  with `occupied_area`, `band_area`, `sheet_area`, `sheet_center_world`,
  and embedded `sheet` result

Keep:

- existing wrench fields
- existing `regions` vector as the canonical band result

## 7.3 `ChSDFShapePair.*`

Add member:

- `ChSDFSheetCollapseSettings m_sheet_settings;`

Add API:

- `void SetSheetCollapseSettings(const ChSDFSheetCollapseSettings&)`
- `const ChSDFSheetCollapseSettings& GetSheetCollapseSettings() const`

Change `EvaluateContact()`:

- call `ChSDFVolumeContactEvaluator::EvaluateBrickPairs(...)`
- call `ChSDFSheetBuilder::BuildFromBandResult(...)`
- copy:
  - `occupied_area`
  - `band_area`
  - `sheet_area`
  - `sheet_center_world`
  - `sheet`
- keep hysteresis and `Apply()` tied to band mechanics fields

## 7.4 `tools/sdf_contact_benchmarks/src/main.cpp`

This file already became the first supplement4 consumer.

Target output semantics:

- `active_area` = compatibility alias to `A_band`
- `area_occ` = `A_occ`
- `area_band` = `A_band`
- `area_sheet` = `A_sheet`

Recommended next additions:

- `sheet_center_x`
- `sheet_patch_count`
- `sheet_support_bbox_xmin/xmax/zmin/zmax`

## 7.5 `ChSDFContactSurface.*`

Current status:

- this file belongs to the older chart/clipped-cell PFC path

Supplement4 mapping:

- do not make `ChSDFContactSurfaceBuilder` the main sheet implementation
- keep it as legacy/experimental geometry code

If retained:

- mark it explicitly as legacy PFC surface reconstruction
- do not route supplement4 sheet collapse through it

Reason:

- supplement4's sheet is a collapse of band cells
- it is not the same object as the old chart-based equal-pressure surface

## 8. Benchmark Semantics After the Refactor

For every pair result:

- `occupied_area`
  - reflects active-cell occupancy and layer switching
- `band_area`
  - reflects band measure from `ΔA_c`
- `sheet_area`
  - reflects collapsed zero-thickness geometry

Interpretation:

- if `occupied_area` jumps but `sheet_area` stays stable:
  the problem is representation-layer occupancy, not true contact geometry
- if `band_area` jumps with force while `sheet_area` stays stable:
  force mechanics and geometric representation are decoupled, as intended
- if `sheet_area` also jumps:
  collapse/fiber construction still needs work

This is exactly the diagnostic logic supplement4 wants.

## 9. Implementation Order

### Step 1

Add `ChSDFSheetRepresentation.*` and implement minimal dominant-axis
measure-preserving collapse from `ChSDFBrickPairWrenchResult`.

### Step 2

Extend `ChSDFShapePairContactResult` with `occupied_area`, `band_area`,
`sheet_area`, and embedded sheet output.

### Step 3

Update `ChSDFShapePair::EvaluateContact()` to execute:

`FindBrickPairs -> Band evaluation -> Sheet collapse -> result packaging`

### Step 4

Update benchmark output to consume the top-level sheet fields instead of
recomputing them externally.

### Step 5

Add sheet-patch clustering and sheet-based pressure center / footprint outputs.

### Step 6

Upgrade from dominant-axis column collapse to local-fiber collapse.

## 10. Final Mapping Summary

In the current codebase, supplement4 should be implemented as:

- existing `ChSDFVolumeContactEvaluator` = `Band mechanics`
- new `ChSDFSheetBuilder` = `Sheet representation`
- existing `ChSDFShapePair` = dual-layer orchestrator
- existing benchmark tool = first consumer of `A_occ / A_band / A_sheet`

The key rule is:

- do not rewrite the force path first
- do not read geometric semantics directly from band cells
- collapse band cells into a sheet, then read area/patch/footprint from the
  sheet layer

That is the cleanest supplement4-aligned upgrade with the least disruption to
the current Chrono implementation.
