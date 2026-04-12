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

#include <memory>
#include <vector>

#include "chrono/collision/sdf/ChSDFBrickPair.h"
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
    void SetBodyA(ChBody* body) { m_body_a = body; }
    void SetBodyB(ChBody* body) { m_body_b = body; }

    void SetShapeA(std::shared_ptr<ChCollisionShapeSDF> shape) { m_shape_a = shape; }
    void SetShapeB(std::shared_ptr<ChCollisionShapeSDF> shape) { m_shape_b = shape; }

    /// Set the SDF shape frames in the local reference frame used by the owning body collision models.
    void SetShapeAFrame(const ChFrame<>& frame) { m_shape_a_frame = frame; }
    void SetShapeBFrame(const ChFrame<>& frame) { m_shape_b_frame = frame; }

    ChBody* GetBodyA() const { return m_body_a; }
    ChBody* GetBodyB() const { return m_body_b; }

    std::shared_ptr<ChCollisionShapeSDF> GetShapeA() const { return m_shape_a; }
    std::shared_ptr<ChCollisionShapeSDF> GetShapeB() const { return m_shape_b; }

    const ChFrame<>& GetShapeAFrame() const { return m_shape_a_frame; }
    const ChFrame<>& GetShapeBFrame() const { return m_shape_b_frame; }

    /// Return true if all required references are available and both SDF grids are loaded.
    bool IsReady() const;

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

    /// Build connected dual-SDF contact regions and evaluate the resulting distributed contact wrenches.
    ChSDFShapePairContactResult EvaluateContact(const ChSDFBrickPairBroadphase::Settings& pair_settings,
                                                const ChSDFContactRegionBuilder::Settings& region_settings,
                                                const ChSDFNormalPressureSettings& pressure_settings) const;

  private:
    ChBody* m_body_a = nullptr;
    ChBody* m_body_b = nullptr;

    std::shared_ptr<ChCollisionShapeSDF> m_shape_a;
    std::shared_ptr<ChCollisionShapeSDF> m_shape_b;

    ChFrame<> m_shape_a_frame;
    ChFrame<> m_shape_b_frame;
};

/// @} chrono_collision

}  // end namespace chrono

#endif
