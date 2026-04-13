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

void ChSDFShapePair::CreateAccumulators() {
    if (m_body_a && m_accumulator_a == kInvalidAccumulator) {
        m_accumulator_a = m_body_a->AddAccumulator();
    }

    if (m_body_b && m_accumulator_b == kInvalidAccumulator) {
        m_accumulator_b = m_body_b->AddAccumulator();
    }
}

void ChSDFShapePair::EmptyAccumulators() {
    if (m_body_a && m_accumulator_a != kInvalidAccumulator) {
        m_body_a->EmptyAccumulator(m_accumulator_a);
    }

    if (m_body_b && m_accumulator_b != kInvalidAccumulator) {
        m_body_b->EmptyAccumulator(m_accumulator_b);
    }
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

ChFrameMoving<> ChSDFShapePair::GetShapeAFrameAbsMoving() const {
    if (!m_body_a) {
        return ChFrameMoving<>(m_shape_a_frame);
    }

    return m_body_a->GetFrameRefToAbs().TransformLocalToParent(ChFrameMoving<>(m_shape_a_frame));
}

ChFrameMoving<> ChSDFShapePair::GetShapeBFrameAbsMoving() const {
    if (!m_body_b) {
        return ChFrameMoving<>(m_shape_b_frame);
    }

    return m_body_b->GetFrameRefToAbs().TransformLocalToParent(ChFrameMoving<>(m_shape_b_frame));
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

ChSDFShapePairContactResult ChSDFShapePair::EvaluateContact(
    const ChSDFBrickPairBroadphase::Settings& pair_settings,
    const ChSDFContactRegionBuilder::Settings& region_settings,
    const ChSDFNormalPressureSettings& pressure_settings) const {
    if (!IsReady()) {
        return {};
    }

    const auto regions = BuildContactRegions(pair_settings, region_settings);
    return ChSDFContactWrenchEvaluator::EvaluateBrickPairRegions(
        regions, GetShapeAFrameAbsMoving(), GetShapeBFrameAbsMoving(),
        ChSDFContactWrenchEvaluator::MakeEffectiveMassProperties(m_body_a),
        ChSDFContactWrenchEvaluator::MakeEffectiveMassProperties(m_body_b), pressure_settings);
}

bool ChSDFShapePair::Apply(const ChSDFShapePairContactResult& result) {
    if (!result.valid || !IsReady()) {
        return false;
    }

    if (m_accumulator_a == kInvalidAccumulator || m_accumulator_b == kInvalidAccumulator) {
        return false;
    }

    m_body_a->AccumulateForce(m_accumulator_a, result.wrench_world_a.force, result.shape_a_frame_abs.GetPos(), false);
    m_body_a->AccumulateTorque(m_accumulator_a, result.wrench_world_a.torque, false);

    m_body_b->AccumulateForce(m_accumulator_b, result.wrench_world_b.force, result.shape_b_frame_abs.GetPos(), false);
    m_body_b->AccumulateTorque(m_accumulator_b, result.wrench_world_b.torque, false);

    return true;
}

ChSDFShapePairContactResult ChSDFShapePair::EvaluateAndApply(
    const ChSDFBrickPairBroadphase::Settings& pair_settings,
    const ChSDFContactRegionBuilder::Settings& region_settings,
    const ChSDFNormalPressureSettings& pressure_settings) {
    auto result = EvaluateContact(pair_settings, region_settings, pressure_settings);
    Apply(result);
    return result;
}

}  // end namespace chrono
