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

#include "chrono/collision/sdf/ChSDFContactPatch.h"

#include <algorithm>
#include <cmath>

namespace chrono {
namespace {

double GridCoordinate(std::size_t index, std::size_t count, double half_length) {
    if (count <= 1) {
        return 0;
    }

    const double alpha = static_cast<double>(index) / static_cast<double>(count - 1);
    return -half_length + 2.0 * half_length * alpha;
}

}  // namespace

ChSDFContactPatch ChSDFContactPatchSampler::SamplePlanePatch(const ChNanoVDBLevelSet& level_set,
                                                             const ChFrame<>& patch_frame,
                                                             const Settings& settings) {
    ChSDFContactPatch patch;
    patch.patch_frame = patch_frame;

    if (!level_set.IsLoaded() || settings.samples_u == 0 || settings.samples_v == 0) {
        return patch;
    }

    const ChVector3d plane_normal = patch_frame.TransformDirectionLocalToParent(ChVector3d(0, 0, 1)).GetNormalized();

    for (std::size_t iu = 0; iu < settings.samples_u; ++iu) {
        const double u = GridCoordinate(iu, settings.samples_u, settings.half_length_u);

        for (std::size_t iv = 0; iv < settings.samples_v; ++iv) {
            const double v = GridCoordinate(iv, settings.samples_v, settings.half_length_v);
            const ChVector3d plane_point(u, v, 0);
            const ChVector3d world_point = patch_frame.TransformPointLocalToParent(plane_point);
            const ChSDFProbeResult probe = level_set.ProbeWorld(world_point);

            ++patch.total_samples;
            if (!probe.valid) {
                continue;
            }

            const double abs_distance = std::abs(probe.distance);
            const double abs_alignment = std::abs(Vdot(probe.normal_world, plane_normal));

            if (abs_distance > settings.max_abs_distance) {
                continue;
            }

            if (abs_alignment < settings.min_abs_normal_alignment) {
                continue;
            }

            if (settings.require_zero_crossing && !probe.zero_crossing) {
                continue;
            }

            ChSDFContactPatchSample sample;
            sample.iu = iu;
            sample.iv = iv;
            sample.point_plane = plane_point;
            sample.probe = probe;
            sample.abs_distance = abs_distance;
            sample.abs_normal_alignment = abs_alignment;

            patch.samples.push_back(sample);
            patch.accepted_samples = patch.samples.size();
            patch.world_bounds += world_point;
        }
    }

    return patch;
}

}  // end namespace chrono
