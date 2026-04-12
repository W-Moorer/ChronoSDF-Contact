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

#ifndef CH_NANOVDB_LEVEL_SET_H
#define CH_NANOVDB_LEVEL_SET_H

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "chrono/core/ChApiCE.h"
#include "chrono/core/ChVector3.h"
#include "chrono/geometry/ChGeometry.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// Result of probing a sparse signed distance field.
struct ChApi ChSDFProbeResult {
    bool valid = false;
    bool level_set = false;
    bool zero_crossing = false;
    bool inside_world_bbox = false;
    bool inside_index_bbox = false;

    double distance = 0;

    ChVector3d point_world = VNULL;
    ChVector3d point_index = VNULL;
    ChVector3d gradient_world = VNULL;
    ChVector3d gradient_index = VNULL;
    ChVector3d normal_world = VNULL;
};

/// Summary metadata for a NanoVDB level set.
struct ChApi ChNanoVDBGridInfo {
    bool valid = false;
    bool level_set = false;

    std::string grid_name;
    uint64_t active_voxels = 0;

    ChVector3d voxel_size = VNULL;
    ChAABB world_bounds;
    ChAABB index_bounds;
};

/// One sparse narrow-band brick extracted from a NanoVDB leaf node.
struct ChApi ChSDFLeafBrick {
    bool valid = false;

    std::size_t index = 0;
    uint64_t active_voxels = 0;

    double min_value = 0;
    double max_value = 0;
    double min_abs_value = 0;

    ChIntAABB index_bounds;
    ChAABB world_bounds;
    ChVector3d voxel_size = VNULL;
};

/// Reusable NanoVDB-backed sparse level-set query utility.
/// The coordinate system used by ProbeWorld is the coordinate system stored in the NanoVDB grid itself.
class ChApi ChNanoVDBLevelSet {
  public:
    ChNanoVDBLevelSet();
    ~ChNanoVDBLevelSet();

    ChNanoVDBLevelSet(const ChNanoVDBLevelSet&) = delete;
    ChNanoVDBLevelSet& operator=(const ChNanoVDBLevelSet&) = delete;

    ChNanoVDBLevelSet(ChNanoVDBLevelSet&& other) noexcept;
    ChNanoVDBLevelSet& operator=(ChNanoVDBLevelSet&& other) noexcept;

    /// Returns true if this Chrono build was compiled with NanoVDB SDF support enabled.
    static bool HasRuntimeSupport();

    /// Load a float NanoVDB grid from file.
    bool Load(const std::string& filename, const std::string& grid_name = "");

    /// Clear the currently loaded grid.
    void Clear();

    /// Return true if a grid is loaded and ready for queries.
    bool IsLoaded() const;

    /// Return information about the currently loaded grid.
    ChNanoVDBGridInfo GetInfo() const;

    /// Return a cached list of sparse narrow-band bricks, one per NanoVDB leaf node.
    const std::vector<ChSDFLeafBrick>& GetLeafBricks() const;

    /// Probe the sparse level set at a point expressed in the grid's world coordinates.
    ChSDFProbeResult ProbeWorld(const ChVector3d& point_world) const;

    /// Probe the sparse level set at a point expressed in grid index coordinates.
    ChSDFProbeResult ProbeIndex(const ChVector3d& point_index) const;

    /// Return the most recent load/probe error, if any.
    const std::string& GetLastError() const;

  private:
    class Impl;
    std::unique_ptr<Impl> m_impl;
    std::string m_last_error;
};

/// @} chrono_collision

}  // end namespace chrono

#endif
