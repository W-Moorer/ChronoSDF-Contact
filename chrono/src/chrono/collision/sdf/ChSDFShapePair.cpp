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

#include <algorithm>
#include <limits>

#include "chrono/core/ChFrameMoving.h"

namespace chrono {
namespace {

double Clamp01(double value) {
    return std::max(0.0, std::min(value, 1.0));
}

ChVector3d SafeNormalized(const ChVector3d& v, const ChVector3d& fallback = VECT_Z) {
    const double length = v.Length();
    return length > 1.0e-12 ? (v / length) : fallback;
}

double NormalCosine(const ChVector3d& a, const ChVector3d& b) {
    const double la = a.Length();
    const double lb = b.Length();
    if (la <= 1.0e-12 || lb <= 1.0e-12) {
        return -1.0;
    }
    return Vdot(a / la, b / lb);
}

bool PassesThresholds(double area, double penetration, double area_threshold, double penetration_threshold) {
    if (area_threshold > 0 && area < area_threshold) {
        return false;
    }
    if (penetration_threshold > 0 && penetration < penetration_threshold) {
        return false;
    }
    return true;
}

void SuppressRegionResult(ChSDFBrickPairWrenchResult& region_result) {
    region_result.wrench_shape_a = {VNULL, VNULL};
    region_result.wrench_shape_b = {VNULL, VNULL};
    region_result.wrench_world_a = {VNULL, VNULL};
    region_result.wrench_world_b = {VNULL, VNULL};
    region_result.active_samples = 0;
    region_result.active_area = 0;
    region_result.integrated_pressure = 0;
    region_result.mean_pressure = 0;
    region_result.mean_local_stiffness = 0;
    region_result.mean_effective_mass = 0;
    region_result.max_local_stiffness = 0;
    region_result.max_effective_mass = 0;
    region_result.max_penetration = 0;
    region_result.max_pressure = 0;
    region_result.pressure_center_world = VNULL;

    for (auto& sample : region_result.samples) {
        sample.active = false;
        sample.pressure = 0;
        sample.force_world = VNULL;
        sample.force_shape_a = VNULL;
        sample.torque_shape_a = VNULL;
        sample.force_shape_b = VNULL;
        sample.torque_shape_b = VNULL;
        sample.traction_world = VNULL;
    }
}

}  // namespace

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
    const ChSDFNormalPressureSettings& pressure_settings) {
    if (!IsReady()) {
        return {};
    }

    const ChFrame<> shape_a_frame_abs = GetShapeAFrameAbs();
    const ChFrame<> shape_b_frame_abs = GetShapeBFrameAbs();
    const auto regions = StabilizeRegions(BuildContactRegions(pair_settings, region_settings), shape_a_frame_abs, shape_b_frame_abs);

    auto result = ChSDFContactWrenchEvaluator::EvaluateBrickPairRegions(
        regions, GetShapeAFrameAbsMoving(), GetShapeBFrameAbsMoving(),
        ChSDFContactWrenchEvaluator::MakeEffectiveMassProperties(m_body_a),
        ChSDFContactWrenchEvaluator::MakeEffectiveMassProperties(m_body_b), pressure_settings);
    ApplyActivationHysteresis(result);
    UpdateHistory(result);
    return result;
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

void ChSDFShapePair::ClearHistory() {
    m_region_history.clear();
    m_next_persistent_region_id = 1;
}

std::vector<ChSDFBrickPairRegion> ChSDFShapePair::StabilizeRegions(const std::vector<ChSDFBrickPairRegion>& regions,
                                                                   const ChFrame<>& shape_a_frame_abs,
                                                                   const ChFrame<>& shape_b_frame_abs) {
    std::vector<ChSDFBrickPairRegion> stabilized = regions;

    if (!m_stabilization_settings.enable_region_history) {
        for (std::size_t i = 0; i < stabilized.size(); ++i) {
            stabilized[i].persistent_id = i + 1;
            stabilized[i].history_age = 1;
            stabilized[i].matched_history = false;
            stabilized[i].active_history = stabilized[i].HasSamples();
        }
        return stabilized;
    }

    std::vector<char> history_used(m_region_history.size(), 0);
    const double filter_weight = Clamp01(m_stabilization_settings.normal_filter_weight);

    for (auto& region : stabilized) {
        const ChVector3d current_normal = SafeNormalized(region.mean_normal_world);

        int best_match = -1;
        double best_score = std::numeric_limits<double>::infinity();

        for (std::size_t i = 0; i < m_region_history.size(); ++i) {
            if (history_used[i] || m_region_history[i].missed_steps > m_stabilization_settings.max_missed_steps) {
                continue;
            }

            const double distance = (region.centroid_world - m_region_history[i].centroid_world).Length();
            if (m_stabilization_settings.max_match_distance > 0 &&
                distance > m_stabilization_settings.max_match_distance) {
                continue;
            }

            const double cosine = NormalCosine(current_normal, m_region_history[i].filtered_normal_world);
            if (cosine < m_stabilization_settings.min_match_normal_cosine) {
                continue;
            }

            const double distance_score =
                m_stabilization_settings.max_match_distance > 1.0e-12
                    ? distance / m_stabilization_settings.max_match_distance
                    : distance;
            const double score = distance_score + (1.0 - cosine);
            if (score < best_score) {
                best_score = score;
                best_match = static_cast<int>(i);
            }
        }

        if (best_match >= 0) {
            auto& state = m_region_history[static_cast<std::size_t>(best_match)];
            history_used[static_cast<std::size_t>(best_match)] = 1;

            region.persistent_id = state.persistent_id;
            region.history_age = state.age + 1;
            region.matched_history = true;
            region.active_history = state.active;

            if (m_stabilization_settings.enable_normal_filter) {
                for (auto& sample : region.samples) {
                    sample.contact_normal_world =
                        SafeNormalized((1.0 - filter_weight) * state.filtered_normal_world +
                                       filter_weight * sample.contact_normal_world,
                                       state.filtered_normal_world);
                }
                ChSDFContactRegionBuilder::ReparameterizeBrickPairRegion(region, shape_a_frame_abs, shape_b_frame_abs);
            }
        } else {
            region.persistent_id = m_next_persistent_region_id++;
            region.history_age = 1;
            region.matched_history = false;
            region.active_history = false;
        }
    }

    return stabilized;
}

void ChSDFShapePair::ApplyActivationHysteresis(ChSDFShapePairContactResult& result) {
    if (!m_stabilization_settings.enable_region_history) {
        return;
    }

    for (auto& region_result : result.regions) {
        const bool prev_active = region_result.region.active_history;
        const bool activate = PassesThresholds(region_result.active_area, region_result.max_penetration,
                                               m_stabilization_settings.activation_area,
                                               m_stabilization_settings.activation_penetration);
        const bool keep_active =
            prev_active
                ? PassesThresholds(region_result.active_area, region_result.max_penetration,
                                   m_stabilization_settings.deactivation_area,
                                   m_stabilization_settings.deactivation_penetration)
                : activate;

        region_result.region.active_history = keep_active;
        if (!keep_active) {
            SuppressRegionResult(region_result);
        }
    }

    RebuildAggregateResult(result);
}

void ChSDFShapePair::RebuildAggregateResult(ChSDFShapePairContactResult& result) const {
    result.total_regions = result.regions.size();
    result.active_regions = 0;
    result.active_area = 0;
    result.integrated_pressure = 0;
    result.mean_pressure = 0;
    result.mean_local_stiffness = 0;
    result.mean_effective_mass = 0;
    result.max_local_stiffness = 0;
    result.max_effective_mass = 0;
    result.max_penetration = 0;
    result.max_pressure = 0;
    result.wrench_shape_a = {VNULL, VNULL};
    result.wrench_shape_b = {VNULL, VNULL};
    result.wrench_world_a = {VNULL, VNULL};
    result.wrench_world_b = {VNULL, VNULL};

    double local_stiffness_area_sum = 0;
    double effective_mass_area_sum = 0;

    for (const auto& region_result : result.regions) {
        result.wrench_shape_a.force += region_result.wrench_shape_a.force;
        result.wrench_shape_a.torque += region_result.wrench_shape_a.torque;
        result.wrench_shape_b.force += region_result.wrench_shape_b.force;
        result.wrench_shape_b.torque += region_result.wrench_shape_b.torque;
        result.wrench_world_a.force += region_result.wrench_world_a.force;
        result.wrench_world_a.torque += region_result.wrench_world_a.torque;
        result.wrench_world_b.force += region_result.wrench_world_b.force;
        result.wrench_world_b.torque += region_result.wrench_world_b.torque;

        if (!region_result.HasActiveContact()) {
            continue;
        }

        result.active_regions++;
        result.active_area += region_result.active_area;
        result.integrated_pressure += region_result.integrated_pressure;
        local_stiffness_area_sum += region_result.mean_local_stiffness * region_result.active_area;
        effective_mass_area_sum += region_result.mean_effective_mass * region_result.active_area;
        result.max_local_stiffness = std::max(result.max_local_stiffness, region_result.max_local_stiffness);
        result.max_effective_mass = std::max(result.max_effective_mass, region_result.max_effective_mass);
        result.max_penetration = std::max(result.max_penetration, region_result.max_penetration);
        result.max_pressure = std::max(result.max_pressure, region_result.max_pressure);
    }

    if (result.active_area > 0) {
        result.mean_pressure = result.integrated_pressure / result.active_area;
        result.mean_local_stiffness = local_stiffness_area_sum / result.active_area;
        result.mean_effective_mass = effective_mass_area_sum / result.active_area;
    }
}

void ChSDFShapePair::UpdateHistory(const ChSDFShapePairContactResult& result) {
    if (!m_stabilization_settings.enable_region_history) {
        ClearHistory();
        return;
    }

    std::vector<char> matched(m_region_history.size(), 0);
    std::vector<RegionHistoryState> next_history;
    next_history.reserve(result.regions.size() + m_region_history.size());

    for (const auto& region_result : result.regions) {
        RegionHistoryState state;
        bool found_existing = false;

        for (std::size_t i = 0; i < m_region_history.size(); ++i) {
            if (m_region_history[i].persistent_id == region_result.region.persistent_id) {
                state = m_region_history[i];
                matched[i] = 1;
                found_existing = true;
                break;
            }
        }

        if (!found_existing) {
            state.persistent_id =
                region_result.region.persistent_id ? region_result.region.persistent_id : m_next_persistent_region_id++;
        }

        state.centroid_world = region_result.region.centroid_world;
        state.filtered_normal_world = SafeNormalized(region_result.region.mean_normal_world, state.filtered_normal_world);
        state.last_active_area = region_result.active_area;
        state.last_max_penetration = region_result.max_penetration;
        state.age = std::max<std::size_t>(1, region_result.region.history_age);
        state.missed_steps = 0;
        state.active = region_result.region.active_history && region_result.HasActiveContact();
        next_history.push_back(state);
    }

    for (std::size_t i = 0; i < m_region_history.size(); ++i) {
        if (matched[i]) {
            continue;
        }

        auto state = m_region_history[i];
        state.missed_steps++;
        if (state.missed_steps <= m_stabilization_settings.max_missed_steps) {
            next_history.push_back(state);
        }
    }

    std::stable_sort(next_history.begin(), next_history.end(),
                     [](const RegionHistoryState& a, const RegionHistoryState& b) {
                         if (a.active != b.active) {
                             return a.active > b.active;
                         }
                         return a.persistent_id < b.persistent_id;
                     });

    m_region_history = std::move(next_history);
}

}  // end namespace chrono
