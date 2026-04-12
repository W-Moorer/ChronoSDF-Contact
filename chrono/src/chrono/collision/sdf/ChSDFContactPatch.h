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

#ifndef CH_SDF_CONTACT_PATCH_H
#define CH_SDF_CONTACT_PATCH_H

#include <cstddef>
#include <vector>

#include "chrono/collision/sdf/ChNanoVDBLevelSet.h"
#include "chrono/core/ChFrame.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// One accepted patch sample on a probe plane.
struct ChApi ChSDFContactPatchSample {
    std::size_t iu = 0;
    std::size_t iv = 0;

    ChVector3d point_plane = VNULL;
    ChSDFProbeResult probe;

    double abs_distance = 0;
    double abs_normal_alignment = 0;
};

/// Sample set extracted from a probe plane against a sparse level set.
struct ChApi ChSDFContactPatch {
    ChFrame<> patch_frame;
    ChAABB world_bounds;

    std::size_t total_samples = 0;
    std::size_t accepted_samples = 0;

    std::vector<ChSDFContactPatchSample> samples;
};

/// Rectangular-plane patch sampler that can be used as the first stage of a distributed SDF contact model.
class ChApi ChSDFContactPatchSampler {
  public:
    struct Settings {
        double half_length_u = 0.05;
        double half_length_v = 0.05;
        std::size_t samples_u = 11;
        std::size_t samples_v = 11;
        double max_abs_distance = 0.01;
        double min_abs_normal_alignment = 0.0;
        bool require_zero_crossing = false;
    };

    /// Sample a rectangular patch on the z=0 plane of patch_frame.
    /// Accepted samples satisfy the configured distance and normal filters.
    static ChSDFContactPatch SamplePlanePatch(const ChNanoVDBLevelSet& level_set,
                                              const ChFrame<>& patch_frame,
                                              const Settings& settings);
};

/// @} chrono_collision

}  // end namespace chrono

#endif
