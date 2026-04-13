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

#include "chrono/collision/sdf/ChSDFContactPair.h"

namespace chrono {
namespace {

ChWrenchd NegateWrench(const ChWrenchd& wrench) {
    return {-wrench.force, -wrench.torque};
}

}  // namespace

bool ChSDFContactPair::IsReady() const {
    return m_sdf_body && m_patch_body && m_sdf_shape && m_sdf_shape->IsLoaded();
}

void ChSDFContactPair::CreateAccumulators() {
    if (m_sdf_body && m_sdf_accumulator == kInvalidAccumulator) {
        m_sdf_accumulator = m_sdf_body->AddAccumulator();
    }

    if (m_patch_body && m_patch_accumulator == kInvalidAccumulator) {
        m_patch_accumulator = m_patch_body->AddAccumulator();
    }
}

void ChSDFContactPair::EmptyAccumulators() {
    if (m_sdf_body && m_sdf_accumulator != kInvalidAccumulator) {
        m_sdf_body->EmptyAccumulator(m_sdf_accumulator);
    }

    if (m_patch_body && m_patch_accumulator != kInvalidAccumulator) {
        m_patch_body->EmptyAccumulator(m_patch_accumulator);
    }
}

ChSDFContactPair::Result ChSDFContactPair::Evaluate() const {
    return EvaluatePatchFrame(m_patch_frame);
}

ChSDFContactPair::Result ChSDFContactPair::EvaluatePatchFrame(const ChFrame<>& patch_frame_local) const {
    Result result;

    if (!IsReady()) {
        return result;
    }

    result.sdf_frame_abs = m_sdf_body->GetFrameRefToAbs().TransformLocalToParent(ChFrameMoving<>(m_sdf_shape_frame));
    result.patch_frame_abs = m_patch_body->GetFrameRefToAbs().TransformLocalToParent(ChFrameMoving<>(patch_frame_local));

    const ChFrameMoving<> patch_frame_shape_moving = result.sdf_frame_abs.TransformParentToLocal(result.patch_frame_abs);
    result.patch_frame_shape = ChFrame<>(patch_frame_shape_moving.GetCoordsys());

    result.relative_kinematics_shape.linear_velocity = patch_frame_shape_moving.GetPosDt();
    result.relative_kinematics_shape.angular_velocity = patch_frame_shape_moving.GetAngVelParent();
    result.relative_kinematics_shape.body_a = ChSDFContactWrenchEvaluator::MakeEffectiveMassProperties(m_sdf_body);
    result.relative_kinematics_shape.body_b = ChSDFContactWrenchEvaluator::MakeEffectiveMassProperties(m_patch_body);
    result.relative_kinematics_shape.fallback_effective_mass = m_pressure_settings.fallback_effective_mass;

    result.patch_result = ChSDFContactWrenchEvaluator::EvaluatePatchLocal(
        *m_sdf_shape, result.patch_frame_shape, m_patch_settings, result.relative_kinematics_shape, m_pressure_settings);

    const ChWrenchd wrench_patch_abs =
        result.patch_frame_abs.TransformWrenchLocalToParent(result.patch_result.wrench_patch);
    const ChWrenchd wrench_shape_abs =
        result.sdf_frame_abs.TransformWrenchLocalToParent(result.patch_result.wrench_shape);

    result.wrench_on_patch_abs = wrench_patch_abs;
    result.wrench_on_sdf_abs = NegateWrench(wrench_shape_abs);
    result.valid = true;

    return result;
}

ChSDFContactPair::Result ChSDFContactPair::EvaluateCandidate(const ChSDFPatchCandidate& candidate) const {
    if (!candidate.valid) {
        return Result();
    }

    return EvaluatePatchFrame(candidate.patch_frame_patch);
}

std::vector<ChSDFPatchCandidate> ChSDFContactPair::GeneratePatchCandidates(
    const ChSDFPatchCandidateSettings& settings) const {
    if (!IsReady()) {
        return {};
    }

    const ChFrameMoving<> sdf_frame_abs =
        m_sdf_body->GetFrameRefToAbs().TransformLocalToParent(ChFrameMoving<>(m_sdf_shape_frame));
    return ChSDFPatchCandidateGenerator::GenerateCandidates(*m_sdf_shape, sdf_frame_abs, m_patch_body->GetFrameRefToAbs(),
                                                            settings);
}

std::vector<ChSDFContactRegion> ChSDFContactPair::BuildContactRegions(
    const ChFrame<>& patch_frame_local,
    const ChSDFContactRegionBuilder::Settings& settings) const {
    if (!IsReady()) {
        return {};
    }

    const ChFrameMoving<> sdf_frame_abs = m_sdf_body->GetFrameRefToAbs().TransformLocalToParent(ChFrameMoving<>(m_sdf_shape_frame));
    const ChFrameMoving<> patch_frame_abs = m_patch_body->GetFrameRefToAbs().TransformLocalToParent(ChFrameMoving<>(patch_frame_local));
    const ChFrameMoving<> patch_frame_shape_moving = sdf_frame_abs.TransformParentToLocal(patch_frame_abs);
    const ChFrame<> patch_frame_shape(patch_frame_shape_moving.GetCoordsys());

    return ChSDFContactRegionBuilder::BuildPatchRegions(*m_sdf_shape, patch_frame_shape, m_patch_settings, settings);
}

std::vector<ChSDFContactRegion> ChSDFContactPair::BuildContactRegions(
    const ChSDFPatchCandidate& candidate,
    const ChSDFContactRegionBuilder::Settings& settings) const {
    if (!candidate.valid) {
        return {};
    }

    return BuildContactRegions(candidate.patch_frame_patch, settings);
}

std::vector<ChSDFContactRegion> ChSDFContactPair::BuildContactRegionsAuto(
    const ChSDFPatchCandidateSettings& candidate_settings,
    const ChSDFContactRegionBuilder::Settings& region_settings,
    std::size_t candidate_index) const {
    const auto candidates = GeneratePatchCandidates(candidate_settings);
    if (candidate_index >= candidates.size()) {
        return {};
    }

    return BuildContactRegions(candidates[candidate_index], region_settings);
}

ChSDFContactPair::Result ChSDFContactPair::EvaluateAuto(const ChSDFPatchCandidateSettings& settings,
                                                        std::size_t candidate_index) const {
    const auto candidates = GeneratePatchCandidates(settings);
    if (candidate_index >= candidates.size()) {
        return Result();
    }

    return EvaluateCandidate(candidates[candidate_index]);
}

bool ChSDFContactPair::Apply(const Result& result) {
    if (!result.valid || !IsReady()) {
        return false;
    }

    if (m_sdf_accumulator == kInvalidAccumulator || m_patch_accumulator == kInvalidAccumulator) {
        return false;
    }

    m_patch_body->AccumulateForce(m_patch_accumulator, result.wrench_on_patch_abs.force, result.patch_frame_abs.GetPos(),
                                  false);
    m_patch_body->AccumulateTorque(m_patch_accumulator, result.wrench_on_patch_abs.torque, false);

    m_sdf_body->AccumulateForce(m_sdf_accumulator, result.wrench_on_sdf_abs.force, result.sdf_frame_abs.GetPos(), false);
    m_sdf_body->AccumulateTorque(m_sdf_accumulator, result.wrench_on_sdf_abs.torque, false);

    return true;
}

ChSDFContactPair::Result ChSDFContactPair::EvaluateAndApply() {
    Result result = Evaluate();
    Apply(result);
    return result;
}

ChSDFContactPair::Result ChSDFContactPair::EvaluateAutoAndApply(const ChSDFPatchCandidateSettings& settings,
                                                                std::size_t candidate_index) {
    Result result = EvaluateAuto(settings, candidate_index);
    Apply(result);
    return result;
}

}  // end namespace chrono
