# SDF-Equivalent PFC Design for Current Chrono Codebase

## 1. Scope

This document defines how to turn the current dual-SDF contact path into an
SDF-equivalent implementation of the paper's pressure field contact (PFC)
model, while preserving the existing brick-pair broadphase and rigid-body
wrench application path.

The current code already has:

- sparse brick-pair broadphase
- connected-region extraction
- shape-pair wrench aggregation and application

The missing middle layer is the paper's real core:

- object-fixed potential pressure fields `p0_A(x)` and `p0_B(x)`
- an explicit contact-surface object for each connected region
- separation between support domain and scalar field
- clipped-cell / marching-squares integration on a 2D local chart

## 2. Current vs target pipeline

### Current pipeline

`FindBrickPairs -> BuildBrickPairSamples -> BuildBrickPairRegions -> EvaluateBrickPairRegion -> Apply`

Current primary state:

- carrier sample on A surface
- projection to B
- `normal_gap`
- `area_weight`

Current normal traction:

`pressure = stiffness * penetration(normal_gap) + damping * closing_speed`

### Target PFC pipeline

`FindBrickPairs -> BuildBrickPairRegions(seed only) -> BuildRegionCharts -> BuildContactSurfaceRegions -> EvaluatePFCRegions -> Apply`

Target primary state:

- object-fixed potential fields `p0_A(x)`, `p0_B(x)`
- field residual `h(x) = p0_A(x) - p0_B(x)`
- support overlap scalar `omega(u,v)` on a local chart
- clipped support polygon per chart cell
- quadrature points lifted to the equal-pressure contact surface

Target normal traction:

`pressure = p_cap + damping_term`

where `p_cap = p0_A(x_cap) = p0_B(x_cap)` and `x_cap` is the equal-pressure
point on the local contact surface.

## 3. File-level design

### 3.1 Add a shape-fixed potential-field layer

New file:

- `chrono/src/chrono/collision/sdf/ChSDFPotentialField.h`
- `chrono/src/chrono/collision/sdf/ChSDFPotentialField.cpp`

New types:

- `struct ChSDFPotentialFieldSettings`
- `struct ChSDFPotentialFieldProbe`
- `class ChSDFPotentialFieldEvaluator`

Proposed responsibilities:

- map raw signed distance `phi(x)` to inward depth `d(x) = max(-phi(x), 0)`
- map depth to object-fixed virtual pressure field `p0(x)`
- provide `grad_p0(x)` from the SDF gradient

Recommended minimal field definition:

- `depth = clamp(max(-phi, 0), 0, depth_cap)`
- `strain = depth / depth_scale`
- `p0 = modulus * strain`

Recommended settings:

- `double modulus`
- `double depth_scale`
- `double depth_cap`
- `double support_margin`
- `bool clamp_outside_to_zero`

Recommended probe output:

- `bool valid`
- `double phi`
- `double depth`
- `double p0`
- `double support_value`
- `ChVector3d normal_world`
- `ChVector3d grad_p0_world`
- `ChSDFProbeResult sdf_probe`

Changes to existing shape API:

File:

- `chrono/src/chrono/collision/ChCollisionShapeSDF.h`

Add members:

- `void SetPotentialFieldSettings(const ChSDFPotentialFieldSettings& settings);`
- `const ChSDFPotentialFieldSettings& GetPotentialFieldSettings() const;`
- `ChSDFPotentialFieldProbe ProbePotentialLocal(const ChVector3d& point_local) const;`

Rationale:

- `p0_A` and `p0_B` are object-fixed fields, so they belong with the shape, not
  with the contact evaluator.

## 3.2 Keep the current region builder as seed generation only

File:

- `chrono/src/chrono/collision/sdf/ChSDFContactRegion.h`
- `chrono/src/chrono/collision/sdf/ChSDFContactRegion.cpp`

Keep:

- `ChSDFBrickPairRegionSample`
- `ChSDFBrickPairRegion`
- `BuildBrickPairSamples(...)`
- `BuildBrickPairRegions(...)`

But change their role:

- `BuildBrickPairSamples(...)` is no longer the final contact quadrature source.
- `BuildBrickPairRegions(...)` becomes a seed-region extractor only.
- `normal_gap` becomes an auxiliary broadphase/debug quantity only.

No longer treat these fields as constitutive variables:

- `normal_gap`
- `combined_gap`
- `area_weight`

Still useful as seed/debug fields:

- `surface_world_a`
- `surface_world_b`
- `normal_world_a`
- `normal_world_b`
- `contact_normal_world`
- `coord`
- `brick_a_index`
- `brick_b_index`

## 3.3 Add an explicit region contact-surface layer

New file:

- `chrono/src/chrono/collision/sdf/ChSDFContactSurface.h`
- `chrono/src/chrono/collision/sdf/ChSDFContactSurface.cpp`

New types:

- `struct ChSDFLineSupportInterval`
- `struct ChSDFRegionChartSettings`
- `struct ChSDFRegionChartNode`
- `struct ChSDFRegionChartCell`
- `struct ChSDFContactQuadraturePoint`
- `struct ChSDFContactSurfaceRegion`
- `class ChSDFContactSurfaceBuilder`

### `ChSDFLineSupportInterval`

Purpose:

- represent the overlap of the two shapes along one chart-normal line

Fields:

- `bool valid`
- `double s_min`
- `double s_max`
- `double overlap_thickness`

Interpretation:

- `overlap_thickness` is the support scalar `omega`
- support indicator is `omega > 0`
- this replaces the old use of `normal_gap` as the main activity test

### `ChSDFRegionChartSettings`

Purpose:

- control the B-side local chart resolution and probing window

Fields:

- `double cell_size`
- `double chart_padding`
- `double normal_search_half_length`
- `double support_epsilon`
- `double equilibrium_tolerance`
- `std::size_t max_cells_u`
- `std::size_t max_cells_v`

### `ChSDFRegionChartNode`

Purpose:

- per-node data on the 2D chart before cell clipping

Fields:

- `int iu`
- `int iv`
- `ChVector3d anchor_world`
- `ChVector3d point_world_b`
- `ChVector3d normal_world_b`
- `ChSDFLineSupportInterval support_interval`
- `ChSDFPotentialFieldProbe probe_a`
- `ChSDFPotentialFieldProbe probe_b`
- `double support_overlap`
- `bool supported`

Important note:

- node probes are evaluated along the chart-normal line, not directly on the B
  surface only

### `ChSDFRegionChartCell`

Purpose:

- one UV cell to be clipped by support

Fields:

- `int iu`
- `int iv`
- `std::array<std::size_t, 4> corner_indices`
- `bool active`
- `double clipped_area_uv`
- `std::vector<ChVector3d> clipped_polygon_world`
- `std::vector<ChVector3d> clipped_polygon_uv`

Support clipping rule:

- use marching squares on the scalar `support_overlap`
- clip the UV cell by the zero-contour of `support_overlap`

### `ChSDFContactQuadraturePoint`

Purpose:

- one final integration sample on the equal-pressure surface

Fields:

- `ChVector3d point_world`
- `ChVector3d point_shape_a`
- `ChVector3d point_shape_b`
- `ChVector3d normal_world`
- `double area_weight`
- `double p0_a`
- `double p0_b`
- `double p_cap`
- `double h_value`
- `double support_overlap`
- `ChSDFPotentialFieldProbe probe_a`
- `ChSDFPotentialFieldProbe probe_b`

Interpretation:

- `point_world` is `x_cap`
- `p_cap` is the paper's contact pressure field value on the contact surface
- `h_value` should be near zero after the 1D equilibrium solve

### `ChSDFContactSurfaceRegion`

Purpose:

- region-level equivalent of Drake's `ContactSurface`

Fields:

- `ChSDFBrickPairRegion seed_region`
- `ChFrame<> chart_frame_world`
- `ChFrame<> chart_frame_shape_a`
- `ChFrame<> chart_frame_shape_b`
- `ChAABB chart_bounds_uv`
- `double cell_size`
- `std::vector<ChSDFRegionChartNode> nodes`
- `std::vector<ChSDFRegionChartCell> cells`
- `std::vector<ChSDFContactQuadraturePoint> quadrature_points`

### `ChSDFContactSurfaceBuilder`

Core methods:

- `static ChFrame<> BuildBSideChartFrame(const ChSDFBrickPairRegion& region, const ChFrame<>& shape_b_frame_abs);`
- `static ChSDFContactSurfaceRegion BuildRegionSurface(const ChCollisionShapeSDF& shape_a, const ChFrame<>& shape_a_frame_abs, const ChCollisionShapeSDF& shape_b, const ChFrame<>& shape_b_frame_abs, const ChSDFBrickPairRegion& seed_region, const ChSDFRegionChartSettings& chart_settings);`
- `static std::vector<ChSDFContactSurfaceRegion> BuildRegionSurfaces(...);`

Detailed runtime for one region:

1. Build a B-side chart frame from B-side seed surfels using PCA in the tangent
   plane.
2. Rasterize a UV chart around the seed-region footprint.
3. For each chart node, trace one line along the chart normal.
4. Along that line, recover the local inside intervals for shape A and shape B.
5. Compute `omega = measure(I_A intersect I_B)`.
6. Mark the node as supported if `omega > support_epsilon`.
7. For each active UV cell, clip the support polygon with marching squares.
8. For each quadrature point inside the clipped polygon, solve
   `h(s) = p0_A(x(s)) - p0_B(x(s)) = 0` along the chart-normal line segment.
9. Store the lifted equilibrium point as a `ChSDFContactQuadraturePoint`.

## 3.4 Add a dedicated PFC evaluator

New file:

- `chrono/src/chrono/collision/sdf/ChSDFPFCEvaluator.h`
- `chrono/src/chrono/collision/sdf/ChSDFPFCEvaluator.cpp`

New types:

- `struct ChSDFPFCSettings`
- `struct ChSDFPFCQuadratureSample`
- `struct ChSDFPFCRegionResult`
- `class ChSDFPFCEvaluator`

### `ChSDFPFCSettings`

Purpose:

- hold traction-law settings for the PFC path

Fields:

- `double dissipation`
- `double damping_ratio`
- `double adhesion_pressure`
- `double friction_coefficient`
- `double tangential_velocity_regularization`
- `double min_effective_mass`
- `double max_effective_mass`
- `double fallback_effective_mass`
- `bool clamp_negative_pressure`

Important change:

- `stiffness` and `distance_offset` are no longer the main normal constitutive
  parameters for the PFC path
- pressure comes from `p_cap`, not from `penetration(normal_gap)`

### `ChSDFPFCQuadratureSample`

Fields:

- `ChSDFContactQuadraturePoint quadrature_point`
- `bool active`
- `double normal_speed`
- `double tangential_speed`
- `double effective_mass`
- `double damping_pressure`
- `double total_pressure`
- `ChVector3d tangential_velocity_world`
- `ChVector3d traction_world`
- `ChVector3d force_world`
- `ChVector3d torque_world_a`
- `ChVector3d torque_world_b`

### `ChSDFPFCRegionResult`

Fields:

- `ChSDFContactSurfaceRegion surface_region`
- `ChWrenchd wrench_shape_a`
- `ChWrenchd wrench_shape_b`
- `ChWrenchd wrench_world_a`
- `ChWrenchd wrench_world_b`
- `double active_area`
- `double integrated_pressure`
- `double max_pressure`
- `std::vector<ChSDFPFCQuadratureSample> samples`

### `ChSDFPFCEvaluator`

Core methods:

- `static ChSDFPFCRegionResult EvaluateRegion(const ChSDFContactSurfaceRegion& surface_region, const ChFrameMoving<>& shape_a_frame_abs, const ChFrameMoving<>& shape_b_frame_abs, const ChSDFEffectiveMassProperties& body_a, const ChSDFEffectiveMassProperties& body_b, const ChSDFPFCSettings& settings);`
- `static ChSDFShapePairContactResult EvaluateRegions(const std::vector<ChSDFContactSurfaceRegion>& surfaces, const ChFrameMoving<>& shape_a_frame_abs, const ChFrameMoving<>& shape_b_frame_abs, const ChSDFEffectiveMassProperties& body_a, const ChSDFEffectiveMassProperties& body_b, const ChSDFPFCSettings& settings);`

Normal traction for one quadrature point:

- `p_cap = 0.5 * (p0_a + p0_b)`
- `damping_pressure = c_eff * max(closing_speed, 0)`
- `total_pressure = clamp(p_cap + damping_pressure + adhesion_pressure)`

Normal direction:

- preferred: `normalize(grad_p0_A - grad_p0_B)`
- fallback: normalized average of opposing surface normals

## 3.5 Keep the top-level shape-pair object, but reroute its middle stage

File:

- `chrono/src/chrono/collision/sdf/ChSDFShapePair.h`
- `chrono/src/chrono/collision/sdf/ChSDFShapePair.cpp`

Keep unchanged:

- `FindBrickPairs(...)`
- `Apply(...)`
- `EvaluateAndApply(...)`

Change:

- `EvaluateContact(...)`

Recommended new internal flow:

1. `FindBrickPairs(...)`
2. `BuildContactRegions(...)`
3. `StabilizeRegions(...)`
4. `BuildRegionSurfaces(...)`
5. `EvaluateRegions(...)`
6. `ApplyActivationHysteresis(...)`
7. `UpdateHistory(...)`

Recommended API additions:

- `ChSDFShapePairContactResult EvaluatePFCContact(const ChSDFBrickPairBroadphase::Settings& pair_settings, const ChSDFContactRegionBuilder::Settings& region_settings, const ChSDFRegionChartSettings& chart_settings, const ChSDFPFCSettings& pfc_settings);`
- `ChSDFShapePairContactResult EvaluatePFCAndApply(...);`

Compatibility recommendation:

- keep the current `EvaluateContact(...)` as the legacy penalty path
- add the PFC path in parallel first
- only replace the default call site after the PFC path is benchmarked

## 4. Function replacement map

### Keep as-is

- `ChSDFShapePair::FindBrickPairs(...)`
- `ChSDFShapePair::Apply(...)`
- `ChSDFShapePair::EvaluateAndApply(...)` for the legacy path
- `ChSDFContactWrenchEvaluator::MakeEffectiveMassProperties(...)`
- `ChSDFContactWrenchEvaluator::EstimateEffectiveMass(...)`

### Keep, but downgrade to seed/debug role

- `ChSDFContactRegionBuilder::BuildBrickPairSamples(...)`
- `ChSDFContactRegionBuilder::BuildBrickPairRegions(...)`
- `ChSDFContactRegionBuilder::ReparameterizeBrickPairRegion(...)`

### Replace in the new PFC path

- Replace `ChSDFContactWrenchEvaluator::EvaluateBrickPairRegion(...)`
  with `ChSDFPFCEvaluator::EvaluateRegion(...)`
- Replace `ChSDFContactWrenchEvaluator::EvaluateBrickPairRegions(...)`
  with `ChSDFPFCEvaluator::EvaluateRegions(...)`
- Insert `ChSDFContactSurfaceBuilder::BuildRegionSurface(...)`
  between region extraction and wrench evaluation

### Preserve only for legacy benchmark comparison

- `ChSDFContactWrenchEvaluator::EvaluatePatch(...)`
- `ChSDFContactWrenchEvaluator::EvaluatePatchLocal(...)`

## 5. Variable migration

### Existing seed variables that remain

- `surface_world_a`
- `surface_world_b`
- `surface_shape_a`
- `surface_shape_b`
- `normal_world_a`
- `normal_world_b`
- `contact_normal_world`
- `coord`
- `point_patch`

### Variables that lose primary constitutive meaning

- `normal_gap`
- `combined_gap`
- `area_weight`
- `distance_b`

### New primary variables

- `phi_a(x)`
- `phi_b(x)`
- `p0_a(x)`
- `p0_b(x)`
- `h(x) = p0_a(x) - p0_b(x)`
- `omega(u,v)` = support overlap thickness along the chart-normal line
- `x_cap(u,v)` = equal-pressure point on the contact surface
- `p_cap(u,v)` = contact pressure field value on the surface
- `grad_p0_a(x_cap)`
- `grad_p0_b(x_cap)`

### Concrete replacement table

| Current variable | Current meaning | New role in PFC path |
| --- | --- | --- |
| `normal_gap` | penetration proxy along A normal | debug seed only; optional root bracket initializer |
| `combined_gap` | heuristic near-contact score | seed/filter score only |
| `area_weight` | sample-based face area | replaced by clipped-cell quadrature area |
| `point_world` | carrier sample point | seed point only; final integration uses `x_cap` |
| `contact_normal_world` | carrier normal | seed chart normal only; final traction uses `n_cap` |
| `pressure` | `k * penetration + c * vn` | `p_cap + damping_pressure` |

## 6. Why this is closer to the paper and Drake

The paper and Drake separate:

- object-fixed field definition
- contact-surface geometry
- traction integration

The current Chrono path collapses all three into a surfel sample with
`normal_gap`.

The proposed SDF-equivalent PFC path restores the same separation without
copying Drake's tetrahedral data structure:

- `ProbePotentialLocal(...)` replaces Drake's volume-field query
- `ChSDFContactSurfaceRegion` replaces Drake's `ContactSurface`
- clipped UV cells replace tetrahedron-plane polygon clipping
- `h(x) = p0_A(x) - p0_B(x)` replaces the equal-pressure plane residual

## 7. Implementation order

### Step 1

Add `ChSDFPotentialFieldSettings`, `ChSDFPotentialFieldProbe`, and shape-side
potential-field probing.

### Step 2

Add `ChSDFContactSurfaceBuilder` and make it output region charts and clipped
support cells, without changing the force model yet.

### Step 3

Add 1D equilibrium solve on each quadrature point and store
`p0_a`, `p0_b`, `h_value`, `p_cap`, and `x_cap`.

### Step 4

Add `ChSDFPFCEvaluator` and compute normal traction from `p_cap`.

### Step 5

Add friction, effective mass, and damping on the new quadrature points.

### Step 6

Switch `ChSDFShapePair` to the PFC path behind a new API or mode flag, while
keeping the current penalty path for baseline comparison.

## 8. Minimal acceptance criteria

The rewrite is only complete when all of the following are true:

- no region-level force integration reads `region_sample.normal_gap` as the main
  normal constitutive input
- all final integration points carry `p0_a`, `p0_b`, `h_value`, and `p_cap`
- active area comes from clipped-cell polygons, not raw surfel area sums
- final normal direction is derived from the equal-pressure surface, not just the
  seed carrier normal
- `ChSDFShapePair` can evaluate and apply the new PFC path end-to-end
