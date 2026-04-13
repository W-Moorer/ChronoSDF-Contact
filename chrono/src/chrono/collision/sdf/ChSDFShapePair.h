// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================

#ifndef CH_SDF_SHAPE_PAIR_H
#define CH_SDF_SHAPE_PAIR_H

#include <limits>
#include <memory>
#include <vector>

#include "chrono/collision/sdf/ChSDFBrickPair.h"
#include "chrono/collision/sdf/ChSDFContactSurface.h"
#include "chrono/collision/sdf/ChSDFContactWrench.h"
#include "chrono/collision/sdf/ChSDFContactRegion.h"
#include "chrono/physics/ChBody.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// Minimal rigid-body wrapper for dual-SDF brick-pair broadphase.
/// This object does not generate contacts by itself. It only exposes the coarse candidate block pairs
/// between two SDF shapes attached to two rigid bodies.
class ChApi ChSDFShapePair {
  public:
    static constexpr unsigned int kInvalidAccumulator = std::numeric_limits<unsigned int>::max();

    struct StabilizationSettings {
        bool enable_region_history = true;
        bool enable_normal_filter = true;
        double normal_filter_weight = 0.35;
        double max_match_distance = 0.10;
        double min_match_normal_cosine = 0.5;
        std::size_t max_missed_steps = 2;
        double activation_area = 0;
        double deactivation_area = 0;
        double activation_penetration = 0;
        double deactivation_penetration = 0;
    };

    void SetBodyA(ChBody* body) { m_body_a = body; }
    void SetBodyB(ChBody* body) { m_body_b = body; }

    void SetShapeA(std::shared_ptr<ChCollisionShapeSDF> shape) { m_shape_a = shape; }
    void SetShapeB(std::shared_ptr<ChCollisionShapeSDF> shape) { m_shape_b = shape; }

    /// Set the SDF shape frames in the local reference frame used by the owning body collision models.
    void SetShapeAFrame(const ChFrame<>& frame) { m_shape_a_frame = frame; }
    void SetShapeBFrame(const ChFrame<>& frame) { m_shape_b_frame = frame; }
    void SetChartSettings(const ChSDFRegionChartSettings& settings) { m_chart_settings = settings; }

    ChBody* GetBodyA() const { return m_body_a; }
    ChBody* GetBodyB() const { return m_body_b; }

    std::shared_ptr<ChCollisionShapeSDF> GetShapeA() const { return m_shape_a; }
    std::shared_ptr<ChCollisionShapeSDF> GetShapeB() const { return m_shape_b; }

    const ChFrame<>& GetShapeAFrame() const { return m_shape_a_frame; }
    const ChFrame<>& GetShapeBFrame() const { return m_shape_b_frame; }
    const ChSDFRegionChartSettings& GetChartSettings() const { return m_chart_settings; }
    StabilizationSettings& GetStabilizationSettings() { return m_stabilization_settings; }
    const StabilizationSettings& GetStabilizationSettings() const { return m_stabilization_settings; }

    /// Return true if all required references are available and both SDF grids are loaded.
    bool IsReady() const;

    /// Create one accumulator on each body if not already available.
    void CreateAccumulators();

    /// Clear the bound accumulators, if any.
    void EmptyAccumulators();

    void SetAccumulatorAIndex(unsigned int idx) { m_accumulator_a = idx; }
    void SetAccumulatorBIndex(unsigned int idx) { m_accumulator_b = idx; }

    unsigned int GetAccumulatorAIndex() const { return m_accumulator_a; }
    unsigned int GetAccumulatorBIndex() const { return m_accumulator_b; }

    /// Return the absolute frames of the two SDF shapes.
    ChFrame<> GetShapeAFrameAbs() const;
    ChFrame<> GetShapeBFrameAbs() const;
    ChFrameMoving<> GetShapeAFrameAbsMoving() const;
    ChFrameMoving<> GetShapeBFrameAbsMoving() const;

    /// Enumerate coarse sparse brick pairs in the absolute frame.
    std::vector<ChSDFBrickPairCandidate> FindBrickPairs(const ChSDFBrickPairBroadphase::Settings& settings) const;

    /// Build dual-SDF region samples from the retained brick pairs.
    std::vector<ChSDFBrickPairRegionSample> BuildRegionSamples(const ChSDFBrickPairBroadphase::Settings& pair_settings,
                                                               const ChSDFContactRegionBuilder::Settings& region_settings) const;

    /// Build connected dual-SDF contact regions from the retained brick pairs.
    std::vector<ChSDFBrickPairRegion> BuildContactRegions(const ChSDFBrickPairBroadphase::Settings& pair_settings,
                                                          const ChSDFContactRegionBuilder::Settings& region_settings) const;

    /// Build B-side local charts and clipped support cells from the retained brick pairs.
    std::vector<ChSDFContactSurfaceRegion> BuildContactSurfaces(
        const ChSDFBrickPairBroadphase::Settings& pair_settings,
        const ChSDFContactRegionBuilder::Settings& region_settings,
        const ChSDFRegionChartSettings& chart_settings) const;

    /// Build connected dual-SDF contact regions and evaluate the resulting distributed contact wrenches.
    ChSDFShapePairContactResult EvaluateContact(const ChSDFBrickPairBroadphase::Settings& pair_settings,
                                                const ChSDFContactRegionBuilder::Settings& region_settings,
                                                const ChSDFNormalPressureSettings& pressure_settings);

    /// Apply a previously evaluated shape-pair contact result to the currently bound body accumulators.
    bool Apply(const ChSDFShapePairContactResult& result);

    /// Evaluate and immediately apply the resulting distributed contact wrenches to the two bodies.
    ChSDFShapePairContactResult EvaluateAndApply(const ChSDFBrickPairBroadphase::Settings& pair_settings,
                                                 const ChSDFContactRegionBuilder::Settings& region_settings,
                                                 const ChSDFNormalPressureSettings& pressure_settings);

    /// Drop any stored region history state.
    void ClearHistory();

  private:
    struct RegionHistoryState {
        std::size_t persistent_id = 0;
        ChVector3d centroid_world = VNULL;
        ChVector3d filtered_normal_world = VECT_Z;
        double last_active_area = 0;
        double last_max_penetration = 0;
        std::size_t age = 0;
        std::size_t missed_steps = 0;
        bool active = false;
    };

    std::vector<ChSDFBrickPairRegion> StabilizeRegions(const std::vector<ChSDFBrickPairRegion>& regions,
                                                       const ChFrame<>& shape_a_frame_abs,
                                                       const ChFrame<>& shape_b_frame_abs);
    void ApplyActivationHysteresis(ChSDFShapePairContactResult& result);
    void RebuildAggregateResult(ChSDFShapePairContactResult& result) const;
    void UpdateHistory(const ChSDFShapePairContactResult& result);

    ChBody* m_body_a = nullptr;
    ChBody* m_body_b = nullptr;

    std::shared_ptr<ChCollisionShapeSDF> m_shape_a;
    std::shared_ptr<ChCollisionShapeSDF> m_shape_b;

    ChFrame<> m_shape_a_frame;
    ChFrame<> m_shape_b_frame;
    ChSDFRegionChartSettings m_chart_settings;

    StabilizationSettings m_stabilization_settings;
    std::vector<RegionHistoryState> m_region_history;
    std::size_t m_next_persistent_region_id = 1;

    unsigned int m_accumulator_a = kInvalidAccumulator;
    unsigned int m_accumulator_b = kInvalidAccumulator;
};

/// @} chrono_collision

}  // end namespace chrono

#endif
