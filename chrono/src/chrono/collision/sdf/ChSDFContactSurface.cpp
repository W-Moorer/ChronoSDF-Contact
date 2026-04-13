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

#include "chrono/collision/sdf/ChSDFContactSurface.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "chrono/core/ChFrameMoving.h"
#include "chrono/core/ChMatrix33.h"

namespace chrono {
namespace {

struct ShapeLineInterval {
    bool valid = false;

    double s_min = 0;
    double s_max = 0;
    double surface_s = 0;

    ChVector3d surface_world = VNULL;

    ChSDFPotentialFieldProbe anchor_probe;
    ChSDFPotentialFieldProbe surface_probe;
};

struct ScalarVertex {
    ChVector3d uv = VNULL;
    ChVector3d world = VNULL;
    double value = 0;
};

struct PressureDifferenceSample {
    bool valid = false;

    double s = 0;
    double h = 0;

    ChSDFPotentialFieldProbe probe_a;
    ChSDFPotentialFieldProbe probe_b;
};

struct EqualPressureSolveResult {
    bool valid = false;

    double s = 0;
    double h_value = 0;
    double support_overlap = 0;

    ChVector3d point_world = VNULL;
    ChVector3d normal_world = VNULL;

    ChSDFPotentialFieldProbe probe_a;
    ChSDFPotentialFieldProbe probe_b;
};

ShapeLineInterval BuildShapeLineInterval(const ChCollisionShapeSDF& shape,
                                         const ChFrame<>& shape_frame_abs,
                                         const ChVector3d& origin_world,
                                         const ChVector3d& direction_world,
                                         double half_length,
                                         double step_length);

ChSDFLineSupportInterval IntersectIntervals(const ShapeLineInterval& a, const ShapeLineInterval& b);

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

double ResolvePositive(double value, double fallback) {
    return value > 0 ? value : fallback;
}

ChVector3d SafeNormalized(const ChVector3d& v, const ChVector3d& fallback = VECT_Z) {
    const double length = v.Length();
    return length > 1.0e-12 ? (v / length) : fallback;
}

ChVector3d ProjectOntoTangent(const ChVector3d& direction, const ChVector3d& normal) {
    return direction - normal * Vdot(direction, normal);
}

ChVector3d ChooseAxisTangent(const ChVector3d& normal) {
    const ChVector3d preferred_axes[] = {VECT_X, VECT_Y, VECT_Z};
    ChVector3d best = VNULL;
    double best_length = 0;

    for (const auto& axis : preferred_axes) {
        const ChVector3d tangent = ProjectOntoTangent(axis, normal);
        const double length = tangent.Length();
        if (length > best_length) {
            best = tangent / length;
            best_length = length;
        }
    }

    return best_length > 1.0e-12 ? best : VECT_X;
}

double DistanceToInterval(double value, double interval_min, double interval_max) {
    if (value < interval_min) {
        return interval_min - value;
    }
    if (value > interval_max) {
        return value - interval_max;
    }
    return 0.0;
}

ChSDFPotentialFieldProbe ProbePotentialWorld(const ChCollisionShapeSDF& shape,
                                             const ChFrame<>& shape_frame_abs,
                                             const ChVector3d& point_world) {
    return shape.ProbePotentialLocal(shape_frame_abs.TransformPointParentToLocal(point_world));
}

PressureDifferenceSample SamplePressureDifference(const ChCollisionShapeSDF& shape_a,
                                                 const ChFrame<>& shape_a_frame_abs,
                                                 const ChCollisionShapeSDF& shape_b,
                                                 const ChFrame<>& shape_b_frame_abs,
                                                 const ChVector3d& origin_world,
                                                 const ChVector3d& direction_world,
                                                 double s) {
    PressureDifferenceSample sample;
    sample.s = s;

    const ChVector3d point_world = origin_world + direction_world * s;
    sample.probe_a = ProbePotentialWorld(shape_a, shape_a_frame_abs, point_world);
    sample.probe_b = ProbePotentialWorld(shape_b, shape_b_frame_abs, point_world);
    sample.valid = sample.probe_a.valid && sample.probe_b.valid;
    if (sample.valid) {
        sample.h = sample.probe_a.p0 - sample.probe_b.p0;
    }

    return sample;
}

ChVector3d TransformPotentialGradientWorld(const ChFrame<>& shape_frame_abs, const ChSDFPotentialFieldProbe& probe) {
    return shape_frame_abs.TransformDirectionLocalToParent(probe.grad_p0_local);
}

ChVector3d ResolveContactNormalWorld(const ChFrame<>& shape_a_frame_abs,
                                     const ChFrame<>& shape_b_frame_abs,
                                     const ChSDFPotentialFieldProbe& probe_a,
                                     const ChSDFPotentialFieldProbe& probe_b,
                                     const ChVector3d& fallback_normal_world) {
    const ChVector3d grad_a_world = TransformPotentialGradientWorld(shape_a_frame_abs, probe_a);
    const ChVector3d grad_b_world = TransformPotentialGradientWorld(shape_b_frame_abs, probe_b);
    // Keep the region-contact convention used by the legacy Chrono path:
    // normal points from shape A toward shape B.
    ChVector3d normal_world = grad_b_world - grad_a_world;
    if (normal_world.Length2() <= 1.0e-20) {
        const ChVector3d normal_a_world = shape_a_frame_abs.TransformDirectionLocalToParent(probe_a.normal_local);
        const ChVector3d normal_b_world = shape_b_frame_abs.TransformDirectionLocalToParent(probe_b.normal_local);
        normal_world = normal_a_world - normal_b_world;
    }

    return SafeNormalized(normal_world, fallback_normal_world);
}

EqualPressureSolveResult SolveEqualPressureOnSupportLine(const ChCollisionShapeSDF& shape_a,
                                                         const ChFrame<>& shape_a_frame_abs,
                                                         const ChCollisionShapeSDF& shape_b,
                                                         const ChFrame<>& shape_b_frame_abs,
                                                         const ChVector3d& origin_world,
                                                         const ChVector3d& direction_world,
                                                         double half_length,
                                                         double step_length,
                                                         double equilibrium_tolerance) {
    EqualPressureSolveResult result;

    const ShapeLineInterval interval_a =
        BuildShapeLineInterval(shape_a, shape_a_frame_abs, origin_world, direction_world, half_length, step_length);
    const ShapeLineInterval interval_b =
        BuildShapeLineInterval(shape_b, shape_b_frame_abs, origin_world, direction_world, half_length, step_length);
    const ChSDFLineSupportInterval overlap = IntersectIntervals(interval_a, interval_b);
    if (!overlap.valid) {
        return result;
    }

    result.support_overlap = overlap.overlap_thickness;

    auto finalize_sample = [&](const PressureDifferenceSample& sample) {
        if (!sample.valid) {
            return;
        }
        result.valid = true;
        result.s = sample.s;
        result.h_value = sample.h;
        result.point_world = origin_world + direction_world * sample.s;
        result.probe_a = sample.probe_a;
        result.probe_b = sample.probe_b;
        result.normal_world = ResolveContactNormalWorld(shape_a_frame_abs, shape_b_frame_abs, sample.probe_a, sample.probe_b,
                                                       direction_world);
        if (Vdot(result.normal_world, direction_world) > 0.0) {
            result.normal_world = -result.normal_world;
        }
    };

    if (overlap.overlap_thickness <= 1.0e-12) {
        const PressureDifferenceSample sample =
            SamplePressureDifference(shape_a, shape_a_frame_abs, shape_b, shape_b_frame_abs, origin_world, direction_world,
                                     0.5 * (overlap.s_min + overlap.s_max));
        finalize_sample(sample);
        return result;
    }

    const double safe_step = std::max(step_length, overlap.overlap_thickness / 16.0);
    const int samples_count =
        std::max(4, static_cast<int>(std::ceil(overlap.overlap_thickness / std::max(safe_step, 1.0e-8))));

    PressureDifferenceSample best_sample;
    best_sample.valid = false;
    PressureDifferenceSample previous_sample;
    previous_sample.valid = false;

    auto update_best = [&](const PressureDifferenceSample& sample) {
        if (!sample.valid) {
            return;
        }
        if (!best_sample.valid || std::abs(sample.h) < std::abs(best_sample.h)) {
            best_sample = sample;
        }
    };

    for (int i = 0; i <= samples_count; ++i) {
        const double alpha = static_cast<double>(i) / static_cast<double>(samples_count);
        const double s = overlap.s_min + (overlap.s_max - overlap.s_min) * alpha;
        const PressureDifferenceSample current_sample =
            SamplePressureDifference(shape_a, shape_a_frame_abs, shape_b, shape_b_frame_abs, origin_world, direction_world, s);
        update_best(current_sample);

        if (!current_sample.valid) {
            previous_sample.valid = false;
            continue;
        }

        if (std::abs(current_sample.h) <= equilibrium_tolerance) {
            finalize_sample(current_sample);
            return result;
        }

        if (previous_sample.valid && previous_sample.h * current_sample.h <= 0.0) {
            double low = previous_sample.s;
            double high = current_sample.s;
            PressureDifferenceSample low_sample = previous_sample;
            PressureDifferenceSample high_sample = current_sample;

            for (int iter = 0; iter < 40; ++iter) {
                const double mid = 0.5 * (low + high);
                const PressureDifferenceSample mid_sample = SamplePressureDifference(
                    shape_a, shape_a_frame_abs, shape_b, shape_b_frame_abs, origin_world, direction_world, mid);
                if (!mid_sample.valid) {
                    break;
                }

                update_best(mid_sample);
                if (std::abs(mid_sample.h) <= equilibrium_tolerance || std::abs(high - low) <= 1.0e-8) {
                    finalize_sample(mid_sample);
                    return result;
                }

                if (low_sample.h * mid_sample.h <= 0.0) {
                    high = mid;
                    high_sample = mid_sample;
                } else {
                    low = mid;
                    low_sample = mid_sample;
                }
            }

            finalize_sample(std::abs(low_sample.h) <= std::abs(high_sample.h) ? low_sample : high_sample);
            return result;
        }

        previous_sample = current_sample;
    }

    finalize_sample(best_sample);
    return result;
}

bool FindPhiRootOnLine(const ChCollisionShapeSDF& shape,
                       const ChFrame<>& shape_frame_abs,
                       const ChVector3d& origin_world,
                       const ChVector3d& direction_world,
                       double s_begin,
                       double s_end,
                       const ChSDFPotentialFieldProbe& probe_begin,
                       const ChSDFPotentialFieldProbe& probe_end,
                       double& s_root,
                       ChSDFPotentialFieldProbe& probe_root) {
    if (!probe_begin.valid || !probe_end.valid) {
        return false;
    }

    if (std::abs(probe_begin.phi) <= 1.0e-8) {
        s_root = s_begin;
        probe_root = probe_begin;
        return true;
    }
    if (std::abs(probe_end.phi) <= 1.0e-8) {
        s_root = s_end;
        probe_root = probe_end;
        return true;
    }
    if (probe_begin.phi * probe_end.phi > 0.0) {
        return false;
    }

    double low = s_begin;
    double high = s_end;
    ChSDFPotentialFieldProbe probe_low = probe_begin;

    for (int iter = 0; iter < 32; ++iter) {
        const double mid = 0.5 * (low + high);
        const ChVector3d mid_world = origin_world + direction_world * mid;
        const ChSDFPotentialFieldProbe probe_mid = ProbePotentialWorld(shape, shape_frame_abs, mid_world);
        if (!probe_mid.valid) {
            return false;
        }

        if (std::abs(probe_mid.phi) <= 1.0e-8 || std::abs(high - low) <= 1.0e-8) {
            s_root = mid;
            probe_root = probe_mid;
            return true;
        }

        if (probe_low.phi * probe_mid.phi <= 0.0) {
            high = mid;
        } else {
            low = mid;
            probe_low = probe_mid;
        }
    }

    s_root = 0.5 * (low + high);
    probe_root = ProbePotentialWorld(shape, shape_frame_abs, origin_world + direction_world * s_root);
    return probe_root.valid;
}

std::vector<std::pair<double, double>> ExtractInsideIntervals(const ChCollisionShapeSDF& shape,
                                                              const ChFrame<>& shape_frame_abs,
                                                              const ChVector3d& origin_world,
                                                              const ChVector3d& direction_world,
                                                              double half_length,
                                                              double step_length) {
    std::vector<std::pair<double, double>> intervals;

    const double safe_half_length = std::max(half_length, 1.0e-6);
    const double safe_step = std::max(step_length, 1.0e-6);
    const int steps = std::max(2, static_cast<int>(std::ceil(2.0 * safe_half_length / safe_step)));

    struct LineSample {
        double s = 0;
        ChSDFPotentialFieldProbe probe;
    };

    std::vector<LineSample> samples;
    samples.reserve(static_cast<std::size_t>(steps + 1));

    for (int i = 0; i <= steps; ++i) {
        const double alpha = static_cast<double>(i) / static_cast<double>(steps);
        const double s = -safe_half_length + 2.0 * safe_half_length * alpha;
        samples.push_back(LineSample{s, ProbePotentialWorld(shape, shape_frame_abs, origin_world + direction_world * s)});
    }

    auto is_inside = [](const ChSDFPotentialFieldProbe& probe) { return probe.valid && probe.phi <= 0.0; };

    bool current_inside = is_inside(samples.front().probe);
    double interval_start = samples.front().s;

    for (int i = 1; i <= steps; ++i) {
        const bool next_inside = is_inside(samples[static_cast<std::size_t>(i)].probe);
        if (current_inside != next_inside) {
            double s_root = 0.5 * (samples[static_cast<std::size_t>(i - 1)].s + samples[static_cast<std::size_t>(i)].s);
            ChSDFPotentialFieldProbe probe_root;
            if (!FindPhiRootOnLine(shape, shape_frame_abs, origin_world, direction_world,
                                   samples[static_cast<std::size_t>(i - 1)].s, samples[static_cast<std::size_t>(i)].s,
                                   samples[static_cast<std::size_t>(i - 1)].probe, samples[static_cast<std::size_t>(i)].probe,
                                   s_root, probe_root)) {
                continue;
            }

            if (current_inside) {
                intervals.push_back(std::make_pair(interval_start, s_root));
            } else {
                interval_start = s_root;
            }
        }
        current_inside = next_inside;
    }

    if (current_inside) {
        intervals.push_back(std::make_pair(interval_start, samples.back().s));
    }

    return intervals;
}

ShapeLineInterval BuildShapeLineInterval(const ChCollisionShapeSDF& shape,
                                         const ChFrame<>& shape_frame_abs,
                                         const ChVector3d& origin_world,
                                         const ChVector3d& direction_world,
                                         double half_length,
                                         double step_length) {
    ShapeLineInterval result;
    result.anchor_probe = ProbePotentialWorld(shape, shape_frame_abs, origin_world);

    const auto intervals = ExtractInsideIntervals(shape, shape_frame_abs, origin_world, direction_world, half_length, step_length);
    if (intervals.empty()) {
        if (result.anchor_probe.valid && std::abs(result.anchor_probe.phi) <= 1.0e-8) {
            result.valid = true;
            result.s_min = 0;
            result.s_max = 0;
            result.surface_s = 0;
            result.surface_world = origin_world;
            result.surface_probe = result.anchor_probe;
        }
        return result;
    }

    std::size_t best_index = 0;
    double best_distance = std::numeric_limits<double>::infinity();
    double best_thickness = -1;
    for (std::size_t i = 0; i < intervals.size(); ++i) {
        const double distance = DistanceToInterval(0.0, intervals[i].first, intervals[i].second);
        const double thickness = intervals[i].second - intervals[i].first;
        if (distance < best_distance - 1.0e-12 ||
            (std::abs(distance - best_distance) <= 1.0e-12 && thickness > best_thickness)) {
            best_index = i;
            best_distance = distance;
            best_thickness = thickness;
        }
    }

    result.valid = true;
    result.s_min = intervals[best_index].first;
    result.s_max = intervals[best_index].second;
    result.surface_s = std::abs(result.s_min) <= std::abs(result.s_max) ? result.s_min : result.s_max;
    result.surface_world = origin_world + direction_world * result.surface_s;
    result.surface_probe = ProbePotentialWorld(shape, shape_frame_abs, result.surface_world);
    if (!result.surface_probe.valid) {
        result.surface_probe = result.anchor_probe;
    }

    return result;
}

ChSDFLineSupportInterval IntersectIntervals(const ShapeLineInterval& a, const ShapeLineInterval& b) {
    ChSDFLineSupportInterval overlap;
    if (!a.valid || !b.valid) {
        return overlap;
    }

    overlap.s_min = std::max(a.s_min, b.s_min);
    overlap.s_max = std::min(a.s_max, b.s_max);
    if (overlap.s_max < overlap.s_min) {
        return overlap;
    }

    overlap.valid = true;
    overlap.overlap_thickness = std::max(overlap.s_max - overlap.s_min, 0.0);
    return overlap;
}

double PolygonSignedAreaXY(const std::vector<ChVector3d>& polygon) {
    if (polygon.size() < 3) {
        return 0;
    }

    double twice_area = 0;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        twice_area += a.x() * b.y() - b.x() * a.y();
    }

    return 0.5 * twice_area;
}

ChVector3d PolygonCentroidXY(const std::vector<ChVector3d>& polygon) {
    if (polygon.empty()) {
        return VNULL;
    }

    const double signed_area = PolygonSignedAreaXY(polygon);
    if (std::abs(signed_area) <= 1.0e-12) {
        ChVector3d centroid = VNULL;
        for (const auto& point : polygon) {
            centroid += point;
        }
        return centroid / static_cast<double>(polygon.size());
    }

    double cx = 0;
    double cy = 0;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        const double cross = a.x() * b.y() - b.x() * a.y();
        cx += (a.x() + b.x()) * cross;
        cy += (a.y() + b.y()) * cross;
    }

    const double scale = 1.0 / (6.0 * signed_area);
    return ChVector3d(cx * scale, cy * scale, 0.0);
}

ScalarVertex InterpolateIsoVertex(const ScalarVertex& a, const ScalarVertex& b) {
    const double denominator = a.value - b.value;
    const double t = std::abs(denominator) > 1.0e-12 ? (a.value / denominator) : 0.5;
    ScalarVertex result;
    result.uv = a.uv + (b.uv - a.uv) * t;
    result.world = a.world + (b.world - a.world) * t;
    result.value = 0;
    return result;
}

std::vector<ScalarVertex> ClipTriangleByPositive(const std::array<ScalarVertex, 3>& triangle) {
    std::vector<ScalarVertex> input(triangle.begin(), triangle.end());
    std::vector<ScalarVertex> output;
    output.reserve(4);

    for (std::size_t i = 0; i < input.size(); ++i) {
        const ScalarVertex& current = input[i];
        const ScalarVertex& next = input[(i + 1) % input.size()];
        const bool current_inside = current.value >= 0.0;
        const bool next_inside = next.value >= 0.0;

        if (current_inside && next_inside) {
            output.push_back(next);
        } else if (current_inside && !next_inside) {
            output.push_back(InterpolateIsoVertex(current, next));
        } else if (!current_inside && next_inside) {
            output.push_back(InterpolateIsoVertex(current, next));
            output.push_back(next);
        }
    }

    return output;
}

ChSDFContactSurfacePolygon MakePolygon(const std::vector<ScalarVertex>& clipped_vertices,
                                       const ChFrame<>& chart_frame_world) {
    ChSDFContactSurfacePolygon polygon;
    if (clipped_vertices.size() < 3) {
        return polygon;
    }

    polygon.vertices_uv.reserve(clipped_vertices.size());
    polygon.vertices_world.reserve(clipped_vertices.size());
    for (const auto& vertex : clipped_vertices) {
        polygon.vertices_uv.push_back(vertex.uv);
        polygon.vertices_world.push_back(vertex.world);
    }

    polygon.area_uv = std::abs(PolygonSignedAreaXY(polygon.vertices_uv));
    if (polygon.area_uv <= 1.0e-12) {
        polygon.vertices_uv.clear();
        polygon.vertices_world.clear();
        polygon.area_uv = 0;
        return polygon;
    }

    polygon.centroid_uv = PolygonCentroidXY(polygon.vertices_uv);
    polygon.centroid_world = chart_frame_world.TransformPointLocalToParent(polygon.centroid_uv);
    return polygon;
}

std::size_t NodeIndex(std::size_t iu, std::size_t iv, std::size_t nodes_u) {
    return iv * nodes_u + iu;
}

double TriangleAreaXY(const ChVector3d& a, const ChVector3d& b, const ChVector3d& c) {
    const double twice_area = (b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y());
    return 0.5 * std::abs(twice_area);
}

void AppendPolygonQuadraturePoints(std::vector<ChSDFContactQuadraturePoint>& quadrature_points,
                                   const ChSDFContactSurfacePolygon& polygon,
                                   const ChFrame<>& chart_frame_world,
                                   const ChCollisionShapeSDF& shape_a,
                                   const ChFrame<>& shape_a_frame_abs,
                                   const ChCollisionShapeSDF& shape_b,
                                   const ChFrame<>& shape_b_frame_abs,
                                   const ChVector3d& chart_normal_world,
                                   double search_half_length,
                                   double line_step,
                                   double equilibrium_tolerance) {
    if (polygon.vertices_uv.size() < 3) {
        return;
    }

    const ChVector3d centroid_uv = polygon.centroid_uv;
    for (std::size_t i = 0; i < polygon.vertices_uv.size(); ++i) {
        const ChVector3d& uv0 = polygon.vertices_uv[i];
        const ChVector3d& uv1 = polygon.vertices_uv[(i + 1) % polygon.vertices_uv.size()];

        const double triangle_area_uv = TriangleAreaXY(centroid_uv, uv0, uv1);
        if (triangle_area_uv <= 1.0e-12) {
            continue;
        }

        const ChVector3d quadrature_uv = (centroid_uv + uv0 + uv1) / 3.0;
        const ChVector3d anchor_world = chart_frame_world.TransformPointLocalToParent(quadrature_uv);
        const EqualPressureSolveResult solve = SolveEqualPressureOnSupportLine(
            shape_a, shape_a_frame_abs, shape_b, shape_b_frame_abs, anchor_world, chart_normal_world, search_half_length,
            line_step, equilibrium_tolerance);
        if (!solve.valid) {
            continue;
        }

        ChSDFContactQuadraturePoint quadrature_point;
        quadrature_point.point_world = solve.point_world;
        quadrature_point.point_shape_a = shape_a_frame_abs.TransformPointParentToLocal(solve.point_world);
        quadrature_point.point_shape_b = shape_b_frame_abs.TransformPointParentToLocal(solve.point_world);
        quadrature_point.normal_world = solve.normal_world;
        // Use the clipped chart area directly. The equal-pressure normal is still
        // an approximate lift from the chart, and the projection amplification is
        // numerically unstable in the current SDF chart construction.
        quadrature_point.area_weight = triangle_area_uv;
        quadrature_point.p0_a = solve.probe_a.p0;
        quadrature_point.p0_b = solve.probe_b.p0;
        quadrature_point.p_cap = 0.5 * (quadrature_point.p0_a + quadrature_point.p0_b);
        quadrature_point.h_value = solve.h_value;
        quadrature_point.support_overlap = solve.support_overlap;
        quadrature_point.probe_a = solve.probe_a;
        quadrature_point.probe_b = solve.probe_b;
        quadrature_points.push_back(std::move(quadrature_point));
    }
}

double ResolveBaseSpacing(const ChSDFBrickPairRegion& region,
                          const ChCollisionShapeSDF& shape_a,
                          const ChCollisionShapeSDF& shape_b,
                          const ChSDFRegionChartSettings& settings) {
    if (settings.cell_size > 0) {
        return settings.cell_size;
    }
    if (region.sample_spacing > 0) {
        return region.sample_spacing;
    }

    double spacing = std::numeric_limits<double>::infinity();
    const auto info_a = shape_a.GetGridInfo();
    const auto info_b = shape_b.GetGridInfo();
    if (info_a.valid) {
        spacing = std::min(spacing, MinPositiveComponent(info_a.voxel_size));
    }
    if (info_b.valid) {
        spacing = std::min(spacing, MinPositiveComponent(info_b.voxel_size));
    }

    return std::isfinite(spacing) && spacing > 0 ? spacing : 1.0e-3;
}

ChVector3d BuildPrincipalTangent(const ChSDFBrickPairRegion& region,
                                 const ChVector3d& centroid_world_b,
                                 const ChVector3d& normal_world_b) {
    ChMatrix33<> covariance;
    covariance.setZero();

    for (const auto& sample : region.samples) {
        const ChVector3d projected = ProjectOntoTangent(sample.surface_world_b - centroid_world_b, normal_world_b);
        covariance(0, 0) += projected.x() * projected.x();
        covariance(0, 1) += projected.x() * projected.y();
        covariance(0, 2) += projected.x() * projected.z();
        covariance(1, 0) += projected.y() * projected.x();
        covariance(1, 1) += projected.y() * projected.y();
        covariance(1, 2) += projected.y() * projected.z();
        covariance(2, 0) += projected.z() * projected.x();
        covariance(2, 1) += projected.z() * projected.y();
        covariance(2, 2) += projected.z() * projected.z();
    }

    ChMatrix33<> eigenvectors;
    ChVectorN<double, 3> eigenvalues;
    covariance.SelfAdjointEigenSolve(eigenvectors, eigenvalues);

    for (int column = 2; column >= 0; --column) {
        const ChVector3d raw(eigenvectors(0, column), eigenvectors(1, column), eigenvectors(2, column));
        const ChVector3d tangent = ProjectOntoTangent(raw, normal_world_b);
        const double length = tangent.Length();
        if (length > 1.0e-8) {
            return tangent / length;
        }
    }

    for (const auto& sample : region.samples) {
        const ChVector3d tangent = ProjectOntoTangent(sample.surface_world_b - centroid_world_b, normal_world_b);
        const double length = tangent.Length();
        if (length > 1.0e-8) {
            return tangent / length;
        }
    }

    return ChooseAxisTangent(normal_world_b);
}

std::size_t ClampCellCount(std::size_t requested, std::size_t maximum) {
    if (maximum == 0) {
        return requested;
    }
    return std::max<std::size_t>(1, std::min(requested, maximum));
}

}  // namespace

ChFrame<> ChSDFContactSurfaceBuilder::BuildBSideChartFrame(const ChSDFBrickPairRegion& region,
                                                           const ChFrame<>& shape_b_frame_abs) {
    if (region.samples.empty()) {
        return shape_b_frame_abs;
    }

    ChVector3d centroid_world_b = VNULL;
    ChVector3d mean_normal_world_b = VNULL;
    for (const auto& sample : region.samples) {
        centroid_world_b += sample.surface_world_b;
        mean_normal_world_b += sample.normal_world_b;
    }

    centroid_world_b /= static_cast<double>(region.samples.size());
    const ChVector3d normal_world_b = SafeNormalized(mean_normal_world_b, -region.mean_normal_world);
    const ChVector3d tangent_world = BuildPrincipalTangent(region, centroid_world_b, normal_world_b);

    ChMatrix33<> rotation;
    rotation.SetFromAxisZ(normal_world_b, tangent_world);
    return ChFrame<>(centroid_world_b, rotation);
}

ChSDFContactSurfaceRegion ChSDFContactSurfaceBuilder::BuildRegionSurface(const ChCollisionShapeSDF& shape_a,
                                                                         const ChFrame<>& shape_a_frame_abs,
                                                                         const ChCollisionShapeSDF& shape_b,
                                                                         const ChFrame<>& shape_b_frame_abs,
                                                                         const ChSDFBrickPairRegion& seed_region,
                                                                         const ChSDFRegionChartSettings& chart_settings) {
    ChSDFContactSurfaceRegion surface;
    surface.seed_region = seed_region;
    if (seed_region.samples.empty()) {
        return surface;
    }

    surface.chart_frame_world = BuildBSideChartFrame(seed_region, shape_b_frame_abs);
    surface.chart_frame_shape_a =
        ChFrame<>(shape_a_frame_abs.TransformParentToLocal(ChFrameMoving<>(surface.chart_frame_world)).GetCoordsys());
    surface.chart_frame_shape_b =
        ChFrame<>(shape_b_frame_abs.TransformParentToLocal(ChFrameMoving<>(surface.chart_frame_world)).GetCoordsys());

    ChAABB seed_bounds_uv;
    for (const auto& sample : seed_region.samples) {
        seed_bounds_uv += surface.chart_frame_world.TransformPointParentToLocal(sample.surface_world_b);
    }

    const double base_spacing = ResolveBaseSpacing(seed_region, shape_a, shape_b, chart_settings);
    const double chart_padding = chart_settings.chart_padding >= 0 ? chart_settings.chart_padding : base_spacing;
    const double search_half_length =
        ResolvePositive(chart_settings.normal_search_half_length, 2.5 * std::max(base_spacing, 1.0e-6));
    const double line_step = ResolvePositive(chart_settings.line_search_step, 0.5 * std::max(base_spacing, 1.0e-6));

    double u_min = seed_bounds_uv.min.x();
    double u_max = seed_bounds_uv.max.x();
    double v_min = seed_bounds_uv.min.y();
    double v_max = seed_bounds_uv.max.y();

    if (!std::isfinite(u_min) || !std::isfinite(u_max) || !std::isfinite(v_min) || !std::isfinite(v_max)) {
        u_min = -base_spacing;
        u_max = base_spacing;
        v_min = -base_spacing;
        v_max = base_spacing;
    }

    if (u_max - u_min < base_spacing) {
        const double center = 0.5 * (u_min + u_max);
        u_min = center - 0.5 * base_spacing;
        u_max = center + 0.5 * base_spacing;
    }
    if (v_max - v_min < base_spacing) {
        const double center = 0.5 * (v_min + v_max);
        v_min = center - 0.5 * base_spacing;
        v_max = center + 0.5 * base_spacing;
    }

    u_min -= chart_padding;
    u_max += chart_padding;
    v_min -= chart_padding;
    v_max += chart_padding;

    const std::size_t requested_cells_u =
        std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil((u_max - u_min) / std::max(base_spacing, 1.0e-12))));
    const std::size_t requested_cells_v =
        std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil((v_max - v_min) / std::max(base_spacing, 1.0e-12))));
    surface.cells_u = ClampCellCount(requested_cells_u, chart_settings.max_cells_u);
    surface.cells_v = ClampCellCount(requested_cells_v, chart_settings.max_cells_v);
    surface.nodes_u = surface.cells_u + 1;
    surface.nodes_v = surface.cells_v + 1;
    surface.cell_size_u = (u_max - u_min) / static_cast<double>(surface.cells_u);
    surface.cell_size_v = (v_max - v_min) / static_cast<double>(surface.cells_v);
    surface.cell_size = std::max(surface.cell_size_u, surface.cell_size_v);
    surface.chart_bounds_uv = ChAABB(ChVector3d(u_min, v_min, 0.0), ChVector3d(u_max, v_max, 0.0));

    surface.nodes.reserve(surface.nodes_u * surface.nodes_v);
    const ChVector3d chart_normal_world =
        SafeNormalized(surface.chart_frame_world.TransformDirectionLocalToParent(VECT_Z), VECT_Z);

    for (std::size_t iv = 0; iv < surface.nodes_v; ++iv) {
        for (std::size_t iu = 0; iu < surface.nodes_u; ++iu) {
            const double u = u_min + surface.cell_size_u * static_cast<double>(iu);
            const double v = v_min + surface.cell_size_v * static_cast<double>(iv);

            ChSDFRegionChartNode node;
            node.index = surface.nodes.size();
            node.iu = static_cast<int>(iu);
            node.iv = static_cast<int>(iv);
            node.point_chart = ChVector3d(u, v, 0.0);
            node.anchor_world = surface.chart_frame_world.TransformPointLocalToParent(node.point_chart);

            const ShapeLineInterval interval_a =
                BuildShapeLineInterval(shape_a, shape_a_frame_abs, node.anchor_world, chart_normal_world,
                                       search_half_length, line_step);
            const ShapeLineInterval interval_b =
                BuildShapeLineInterval(shape_b, shape_b_frame_abs, node.anchor_world, chart_normal_world,
                                       search_half_length, line_step);

            node.probe_a = interval_a.anchor_probe;
            node.probe_b = interval_b.anchor_probe;
            node.support_interval = IntersectIntervals(interval_a, interval_b);
            node.support_overlap = node.support_interval.overlap_thickness;
            node.supported = node.support_interval.valid && node.support_overlap > chart_settings.support_epsilon;
            node.point_world_b = interval_b.valid ? interval_b.surface_world : node.anchor_world;

            const ChVector3d normal_world_b =
                interval_b.surface_probe.valid
                    ? shape_b_frame_abs.TransformDirectionLocalToParent(interval_b.surface_probe.normal_local)
                    : shape_b_frame_abs.TransformDirectionLocalToParent(node.probe_b.normal_local);
            node.normal_world_b = SafeNormalized(normal_world_b, chart_normal_world);

            surface.nodes.push_back(node);
        }
    }

    surface.cells.reserve(surface.cells_u * surface.cells_v);
    for (std::size_t iv = 0; iv < surface.cells_v; ++iv) {
        for (std::size_t iu = 0; iu < surface.cells_u; ++iu) {
            ChSDFRegionChartCell cell;
            cell.index = surface.cells.size();
            cell.iu = static_cast<int>(iu);
            cell.iv = static_cast<int>(iv);
            cell.corner_indices = {{NodeIndex(iu, iv, surface.nodes_u), NodeIndex(iu + 1, iv, surface.nodes_u),
                                    NodeIndex(iu + 1, iv + 1, surface.nodes_u), NodeIndex(iu, iv + 1, surface.nodes_u)}};

            const auto& node0 = surface.nodes[cell.corner_indices[0]];
            const auto& node1 = surface.nodes[cell.corner_indices[1]];
            const auto& node2 = surface.nodes[cell.corner_indices[2]];
            const auto& node3 = surface.nodes[cell.corner_indices[3]];

            const ScalarVertex corner0{node0.point_chart, node0.anchor_world, node0.support_overlap - chart_settings.support_epsilon};
            const ScalarVertex corner1{node1.point_chart, node1.anchor_world, node1.support_overlap - chart_settings.support_epsilon};
            const ScalarVertex corner2{node2.point_chart, node2.anchor_world, node2.support_overlap - chart_settings.support_epsilon};
            const ScalarVertex corner3{node3.point_chart, node3.anchor_world, node3.support_overlap - chart_settings.support_epsilon};

            const bool all_inside = corner0.value >= 0.0 && corner1.value >= 0.0 && corner2.value >= 0.0 && corner3.value >= 0.0;
            if (all_inside) {
                ChSDFContactSurfacePolygon polygon;
                polygon.vertices_uv = {corner0.uv, corner1.uv, corner2.uv, corner3.uv};
                polygon.vertices_world = {corner0.world, corner1.world, corner2.world, corner3.world};
                polygon.area_uv = std::abs(PolygonSignedAreaXY(polygon.vertices_uv));
                polygon.centroid_uv = PolygonCentroidXY(polygon.vertices_uv);
                polygon.centroid_world = surface.chart_frame_world.TransformPointLocalToParent(polygon.centroid_uv);
                cell.clipped_polygons.push_back(polygon);
                cell.clipped_area_uv = polygon.area_uv;
            } else {
                const ScalarVertex center{
                    0.25 * (corner0.uv + corner1.uv + corner2.uv + corner3.uv),
                    0.25 * (corner0.world + corner1.world + corner2.world + corner3.world),
                    0.25 * (corner0.value + corner1.value + corner2.value + corner3.value)};

                const std::array<std::array<ScalarVertex, 3>, 4> triangles = {{
                    std::array<ScalarVertex, 3>{{center, corner0, corner1}},
                    std::array<ScalarVertex, 3>{{center, corner1, corner2}},
                    std::array<ScalarVertex, 3>{{center, corner2, corner3}},
                    std::array<ScalarVertex, 3>{{center, corner3, corner0}},
                }};

                for (const auto& triangle : triangles) {
                    const auto clipped = ClipTriangleByPositive(triangle);
                    auto polygon = MakePolygon(clipped, surface.chart_frame_world);
                    if (polygon.area_uv > 0) {
                        cell.clipped_area_uv += polygon.area_uv;
                        cell.clipped_polygons.push_back(std::move(polygon));
                    }
                }
            }

            cell.active = cell.clipped_area_uv > 0;
            if (cell.active) {
                surface.supported_cells++;
                surface.supported_area += cell.clipped_area_uv;
            }

            surface.cells.push_back(std::move(cell));
        }
    }

    for (const auto& cell : surface.cells) {
        if (!cell.active) {
            continue;
        }

        for (const auto& polygon : cell.clipped_polygons) {
            AppendPolygonQuadraturePoints(surface.quadrature_points, polygon, surface.chart_frame_world, shape_a,
                                          shape_a_frame_abs, shape_b, shape_b_frame_abs, chart_normal_world,
                                          search_half_length, line_step, chart_settings.equilibrium_tolerance);
        }
    }

    return surface;
}

std::vector<ChSDFContactSurfaceRegion> ChSDFContactSurfaceBuilder::BuildRegionSurfaces(
    const ChCollisionShapeSDF& shape_a,
    const ChFrame<>& shape_a_frame_abs,
    const ChCollisionShapeSDF& shape_b,
    const ChFrame<>& shape_b_frame_abs,
    const std::vector<ChSDFBrickPairRegion>& seed_regions,
    const ChSDFRegionChartSettings& chart_settings) {
    std::vector<ChSDFContactSurfaceRegion> surfaces;
    surfaces.reserve(seed_regions.size());

    for (const auto& region : seed_regions) {
        surfaces.push_back(
            BuildRegionSurface(shape_a, shape_a_frame_abs, shape_b, shape_b_frame_abs, region, chart_settings));
    }

    return surfaces;
}

}  // end namespace chrono
