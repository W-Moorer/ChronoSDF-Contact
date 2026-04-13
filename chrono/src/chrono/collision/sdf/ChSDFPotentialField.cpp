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

#include "chrono/collision/sdf/ChSDFPotentialField.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace chrono {
namespace {

double MinPositiveComponent(const ChVector3d& v) {
    double result = std::numeric_limits<double>::infinity();
    if (v.x() > 0) {
        result = std::min(result, v.x());
    }
    if (v.y() > 0) {
        result = std::min(result, v.y());
    }
    if (v.z() > 0) {
        result = std::min(result, v.z());
    }
    return std::isfinite(result) ? result : 0.0;
}

}  // namespace

double ChSDFPotentialFieldEvaluator::ResolveDepthScale(const ChNanoVDBGridInfo& info,
                                                       const ChSDFPotentialFieldSettings& settings) {
    if (settings.depth_scale > 0) {
        return settings.depth_scale;
    }

    if (info.valid) {
        const double voxel_length = MinPositiveComponent(info.voxel_size);
        if (voxel_length > 0) {
            return voxel_length;
        }
    }

    return 1.0;
}

ChSDFPotentialFieldProbe ChSDFPotentialFieldEvaluator::Evaluate(const ChSDFProbeResult& sdf_probe,
                                                                const ChNanoVDBGridInfo& info,
                                                                const ChSDFPotentialFieldSettings& settings) {
    ChSDFPotentialFieldProbe result;
    result.sdf_probe = sdf_probe;
    result.point_local = sdf_probe.point_world;
    result.gradient_local = sdf_probe.gradient_world;
    result.normal_local = sdf_probe.normal_world;
    result.phi = sdf_probe.distance;

    if (!sdf_probe.valid) {
        return result;
    }

    const double depth_scale = ResolveDepthScale(info, settings);
    const double signed_depth = -sdf_probe.distance;
    result.support_value = settings.support_margin - sdf_probe.distance;

    double depth = settings.clamp_outside_to_zero ? std::max(signed_depth, 0.0) : signed_depth;
    if (settings.depth_cap > 0) {
        if (settings.clamp_outside_to_zero) {
            depth = std::min(depth, settings.depth_cap);
        } else {
            depth = std::max(-settings.depth_cap, std::min(depth, settings.depth_cap));
        }
    }

    result.valid = true;
    result.depth = depth;

    if (depth_scale > 0) {
        result.p0 = settings.modulus * depth / depth_scale;
    }

    bool active_gradient = depth_scale > 0;
    if (settings.clamp_outside_to_zero) {
        active_gradient = active_gradient && signed_depth > 0;
        if (settings.depth_cap > 0) {
            active_gradient = active_gradient && signed_depth < settings.depth_cap;
        }
    } else if (settings.depth_cap > 0) {
        active_gradient = active_gradient && std::abs(signed_depth) < settings.depth_cap;
    }

    if (active_gradient) {
        result.grad_p0_local = sdf_probe.gradient_world * (-settings.modulus / depth_scale);
    }

    return result;
}

}  // end namespace chrono
