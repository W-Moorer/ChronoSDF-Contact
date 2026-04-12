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

#ifndef CH_SDF_CONTACT_PAIR_H
#define CH_SDF_CONTACT_PAIR_H

#include <limits>
#include <memory>

#include "chrono/collision/sdf/ChSDFContactRegion.h"
#include "chrono/collision/sdf/ChSDFPatchCandidate.h"
#include "chrono/collision/sdf/ChSDFContactWrench.h"
#include "chrono/physics/ChBody.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// Minimal rigid-body pair wrapper for distributed SDF contact evaluation.
/// The pair is defined by:
/// - one rigid body carrying an SDF collision shape at a specified local frame
/// - one rigid body carrying a rectangular probe patch at another local frame
///
/// This object does not rely on Bullet contact generation. It explicitly samples the patch against the SDF,
/// evaluates a normal pressure field, and can apply equal-and-opposite wrenches to the two bodies.
class ChApi ChSDFContactPair {
  public:
    static constexpr unsigned int kInvalidAccumulator = std::numeric_limits<unsigned int>::max();

    struct Result {
        bool valid = false;

        ChFrameMoving<> sdf_frame_abs;
        ChFrameMoving<> patch_frame_abs;
        ChFrame<> patch_frame_shape;

        ChSDFPatchKinematics relative_kinematics_shape;
        ChSDFContactWrenchResult patch_result;

        /// Net wrench acting on the probe body, reduced to the probe patch origin and expressed in the absolute frame.
        ChWrenchd wrench_on_patch_abs = {VNULL, VNULL};
        /// Equal-and-opposite wrench acting on the SDF body, reduced to the SDF shape origin and expressed in the absolute frame.
        ChWrenchd wrench_on_sdf_abs = {VNULL, VNULL};
    };

    ChSDFContactPair() = default;

    void SetSDFBody(ChBody* body) { m_sdf_body = body; }
    void SetPatchBody(ChBody* body) { m_patch_body = body; }

    void SetSDFShape(std::shared_ptr<ChCollisionShapeSDF> shape) { m_sdf_shape = shape; }

    /// Set the SDF shape frame in the local reference frame used by the owning body collision model.
    void SetSDFShapeFrame(const ChFrame<>& frame) { m_sdf_shape_frame = frame; }

    /// Set the probe patch frame in the local reference frame used by the probe body collision model.
    void SetPatchFrame(const ChFrame<>& frame) { m_patch_frame = frame; }

    ChBody* GetSDFBody() const { return m_sdf_body; }
    ChBody* GetPatchBody() const { return m_patch_body; }
    std::shared_ptr<ChCollisionShapeSDF> GetSDFShape() const { return m_sdf_shape; }

    const ChFrame<>& GetSDFShapeFrame() const { return m_sdf_shape_frame; }
    const ChFrame<>& GetPatchFrame() const { return m_patch_frame; }

    ChSDFContactPatchSampler::Settings& GetPatchSettings() { return m_patch_settings; }
    const ChSDFContactPatchSampler::Settings& GetPatchSettings() const { return m_patch_settings; }

    ChSDFNormalPressureSettings& GetPressureSettings() { return m_pressure_settings; }
    const ChSDFNormalPressureSettings& GetPressureSettings() const { return m_pressure_settings; }

    /// Return true if all required references are available and the SDF grid is loaded.
    bool IsReady() const;

    /// Create one accumulator on each body if not already available.
    void CreateAccumulators();

    /// Clear the bound accumulators, if any.
    void EmptyAccumulators();

    void SetSDFAccumulatorIndex(unsigned int idx) { m_sdf_accumulator = idx; }
    void SetPatchAccumulatorIndex(unsigned int idx) { m_patch_accumulator = idx; }

    unsigned int GetSDFAccumulatorIndex() const { return m_sdf_accumulator; }
    unsigned int GetPatchAccumulatorIndex() const { return m_patch_accumulator; }

    /// Evaluate the distributed patch contact without applying forces.
    Result Evaluate() const;

    /// Evaluate the distributed contact using an explicitly provided probe patch frame local to the probe body.
    Result EvaluatePatchFrame(const ChFrame<>& patch_frame_local) const;

    /// Evaluate the distributed contact using one generated patch candidate.
    Result EvaluateCandidate(const ChSDFPatchCandidate& candidate) const;

    /// Generate patch candidates for the current pair state.
    std::vector<ChSDFPatchCandidate> GeneratePatchCandidates(const ChSDFPatchCandidateSettings& settings) const;

    /// Build connected contact regions using an explicitly provided probe patch frame local to the probe body.
    std::vector<ChSDFContactRegion> BuildContactRegions(const ChFrame<>& patch_frame_local,
                                                        const ChSDFContactRegionBuilder::Settings& settings) const;

    /// Build connected contact regions using one generated patch candidate.
    std::vector<ChSDFContactRegion> BuildContactRegions(const ChSDFPatchCandidate& candidate,
                                                        const ChSDFContactRegionBuilder::Settings& settings) const;

    /// Build connected contact regions using the first generated patch candidate (or the requested candidate index).
    std::vector<ChSDFContactRegion> BuildContactRegionsAuto(
        const ChSDFPatchCandidateSettings& candidate_settings,
        const ChSDFContactRegionBuilder::Settings& region_settings,
        std::size_t candidate_index = 0) const;

    /// Evaluate the first generated patch candidate (or the candidate at the requested index).
    Result EvaluateAuto(const ChSDFPatchCandidateSettings& settings, std::size_t candidate_index = 0) const;

    /// Apply a previously evaluated result to the currently bound body accumulators.
    bool Apply(const Result& result);

    /// Evaluate the pair and immediately apply the resulting wrenches to the bound accumulators.
    Result EvaluateAndApply();

    /// Generate a patch candidate, evaluate it, and apply the resulting wrenches.
    Result EvaluateAutoAndApply(const ChSDFPatchCandidateSettings& settings, std::size_t candidate_index = 0);

  private:
    ChBody* m_sdf_body = nullptr;
    ChBody* m_patch_body = nullptr;

    std::shared_ptr<ChCollisionShapeSDF> m_sdf_shape;

    ChFrame<> m_sdf_shape_frame;
    ChFrame<> m_patch_frame;

    ChSDFContactPatchSampler::Settings m_patch_settings;
    ChSDFNormalPressureSettings m_pressure_settings;

    unsigned int m_sdf_accumulator = kInvalidAccumulator;
    unsigned int m_patch_accumulator = kInvalidAccumulator;
};

/// @} chrono_collision

}  // end namespace chrono

#endif
