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

#include "chrono/collision/sdf/ChSDFShapePair.h"

#include "chrono/core/ChFrameMoving.h"

namespace chrono {

bool ChSDFShapePair::IsReady() const {
    return m_body_a && m_body_b && m_shape_a && m_shape_b && m_shape_a->IsLoaded() && m_shape_b->IsLoaded();
}

ChFrame<> ChSDFShapePair::GetShapeAFrameAbs() const {
    if (!m_body_a) {
        return m_shape_a_frame;
    }

    return ChFrame<>(m_body_a->GetFrameRefToAbs().TransformLocalToParent(ChFrameMoving<>(m_shape_a_frame)).GetCoordsys());
}

ChFrame<> ChSDFShapePair::GetShapeBFrameAbs() const {
    if (!m_body_b) {
        return m_shape_b_frame;
    }

    return ChFrame<>(m_body_b->GetFrameRefToAbs().TransformLocalToParent(ChFrameMoving<>(m_shape_b_frame)).GetCoordsys());
}

std::vector<ChSDFBrickPairCandidate> ChSDFShapePair::FindBrickPairs(
    const ChSDFBrickPairBroadphase::Settings& settings) const {
    if (!IsReady()) {
        return {};
    }

    return ChSDFBrickPairBroadphase::FindBrickPairs(*m_shape_a, GetShapeAFrameAbs(), *m_shape_b, GetShapeBFrameAbs(),
                                                    settings);
}

std::vector<ChSDFBrickPairRegionSample> ChSDFShapePair::BuildRegionSamples(
    const ChSDFBrickPairBroadphase::Settings& pair_settings,
    const ChSDFContactRegionBuilder::Settings& region_settings) const {
    if (!IsReady()) {
        return {};
    }

    const auto brick_pairs = FindBrickPairs(pair_settings);
    return ChSDFContactRegionBuilder::BuildBrickPairSamples(*m_shape_a, GetShapeAFrameAbs(), *m_shape_b,
                                                            GetShapeBFrameAbs(), brick_pairs, region_settings);
}

std::vector<ChSDFBrickPairRegion> ChSDFShapePair::BuildContactRegions(
    const ChSDFBrickPairBroadphase::Settings& pair_settings,
    const ChSDFContactRegionBuilder::Settings& region_settings) const {
    if (!IsReady()) {
        return {};
    }

    const auto brick_pairs = FindBrickPairs(pair_settings);
    return ChSDFContactRegionBuilder::BuildBrickPairRegions(*m_shape_a, GetShapeAFrameAbs(), *m_shape_b,
                                                            GetShapeBFrameAbs(), brick_pairs, region_settings);
}

}  // end namespace chrono
