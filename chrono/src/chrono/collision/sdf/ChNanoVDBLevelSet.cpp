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

#include "chrono/collision/sdf/ChNanoVDBLevelSet.h"

#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

#if CHRONO_HAS_NANOVDB_SDF
#include <nanovdb/GridHandle.h>
#include <nanovdb/NanoVDB.h>
#include <nanovdb/io/IO.h>
#include <nanovdb/math/SampleFromVoxels.h>
#endif

namespace chrono {
namespace {

ChVector3d ToChVector(const ChVector3d& v) {
    return v;
}

ChVector3i ToChIntVector(const ChVector3i& v) {
    return v;
}

#if CHRONO_HAS_NANOVDB_SDF
template <class Vec3T>
ChVector3d ToChVector(const Vec3T& v) {
    return ChVector3d(static_cast<double>(v[0]), static_cast<double>(v[1]), static_cast<double>(v[2]));
}

template <class Vec3T>
ChVector3i ToChIntVector(const Vec3T& v) {
    return ChVector3i(static_cast<int>(v[0]), static_cast<int>(v[1]), static_cast<int>(v[2]));
}

nanovdb::Vec3d ToNanoVector(const ChVector3d& v) {
    return nanovdb::Vec3d(v.x(), v.y(), v.z());
}

ChAABB ToWorldBounds(const nanovdb::NanoGrid<float>& grid, const nanovdb::CoordBBox& bbox) {
    const auto bbox_real = bbox.template asReal<double>();
    const auto pmin = grid.indexToWorldF(bbox_real.min());
    const auto pmax = grid.indexToWorldF(bbox_real.max());
    return ChAABB(ToChVector(pmin), ToChVector(pmax));
}
#endif

}  // namespace

class ChNanoVDBLevelSet::Impl {
  public:
#if CHRONO_HAS_NANOVDB_SDF
    nanovdb::GridHandle<nanovdb::HostBuffer> handle;
    const nanovdb::NanoGrid<float>* grid = nullptr;
    mutable bool bricks_cached = false;
    mutable std::vector<ChSDFLeafBrick> bricks;
#endif
};

ChNanoVDBLevelSet::ChNanoVDBLevelSet() : m_impl(new Impl) {}

ChNanoVDBLevelSet::~ChNanoVDBLevelSet() = default;

ChNanoVDBLevelSet::ChNanoVDBLevelSet(ChNanoVDBLevelSet&& other) noexcept = default;

ChNanoVDBLevelSet& ChNanoVDBLevelSet::operator=(ChNanoVDBLevelSet&& other) noexcept = default;

bool ChNanoVDBLevelSet::HasRuntimeSupport() {
#if CHRONO_HAS_NANOVDB_SDF
    return true;
#else
    return false;
#endif
}

bool ChNanoVDBLevelSet::Load(const std::string& filename, const std::string& grid_name) {
    Clear();

#if CHRONO_HAS_NANOVDB_SDF
    try {
        m_impl->handle =
            grid_name.empty() ? nanovdb::io::readGrid(filename) : nanovdb::io::readGrid(filename, grid_name);
        m_impl->grid = m_impl->handle.grid<float>();
        if (!m_impl->grid) {
            m_last_error = "Requested NanoVDB grid is not a float grid.";
            return false;
        }
        m_last_error.clear();
        return true;
    } catch (const std::exception& e) {
        m_last_error = e.what();
        Clear();
        return false;
    }
#else
    (void)filename;
    (void)grid_name;
    m_last_error = "Chrono was built without CH_ENABLE_NANOVDB_SDF.";
    return false;
#endif
}

void ChNanoVDBLevelSet::Clear() {
#if CHRONO_HAS_NANOVDB_SDF
    m_impl->handle = nanovdb::GridHandle<nanovdb::HostBuffer>();
    m_impl->grid = nullptr;
    m_impl->bricks.clear();
    m_impl->bricks_cached = false;
#endif
}

bool ChNanoVDBLevelSet::IsLoaded() const {
#if CHRONO_HAS_NANOVDB_SDF
    return m_impl->grid != nullptr;
#else
    return false;
#endif
}

ChNanoVDBGridInfo ChNanoVDBLevelSet::GetInfo() const {
    ChNanoVDBGridInfo info;

#if CHRONO_HAS_NANOVDB_SDF
    if (!m_impl->grid) {
        return info;
    }

    const nanovdb::GridMetaData meta(*m_impl->grid);
    info.valid = true;
    info.level_set = m_impl->grid->isLevelSet();
    info.grid_name = meta.shortGridName();
    info.active_voxels = meta.activeVoxelCount();
    info.voxel_size = ToChVector(meta.voxelSize());
    info.world_bounds = ChAABB(ToChVector(meta.worldBBox().min()), ToChVector(meta.worldBBox().max()));
    info.index_bounds = ChAABB(ToChVector(meta.indexBBox().min()), ToChVector(meta.indexBBox().max()));
#endif

    return info;
}

const std::vector<ChSDFLeafBrick>& ChNanoVDBLevelSet::GetLeafBricks() const {
#if CHRONO_HAS_NANOVDB_SDF
    if (!m_impl->grid) {
        m_impl->bricks.clear();
        m_impl->bricks_cached = true;
        return m_impl->bricks;
    }

    if (m_impl->bricks_cached) {
        return m_impl->bricks;
    }

    m_impl->bricks.clear();

    const auto& tree = m_impl->grid->tree();
    using LeafNodeType = typename std::remove_reference<decltype(tree)>::type::LeafNodeType;
    const LeafNodeType* leaf = tree.getFirstLeaf();
    const std::size_t leaf_count = static_cast<std::size_t>(tree.nodeCount(0));
    const ChVector3d voxel_size = ToChVector(m_impl->grid->voxelSize());

    m_impl->bricks.reserve(leaf_count);

    for (std::size_t i = 0; i < leaf_count && leaf; ++i) {
        ChSDFLeafBrick brick;
        brick.index = i;
        brick.valid = true;
        brick.voxel_size = voxel_size;

        const auto bbox = leaf->bbox();
        brick.index_bounds = ChIntAABB(ToChIntVector(bbox.min()), ToChIntVector(bbox.max()));
        brick.world_bounds = ToWorldBounds(*m_impl->grid, bbox);

        brick.min_value = std::numeric_limits<double>::infinity();
        brick.max_value = -std::numeric_limits<double>::infinity();
        brick.min_abs_value = std::numeric_limits<double>::infinity();

        for (auto it = leaf->beginValueOn(); it; ++it) {
            const double value = static_cast<double>(*it);
            brick.active_voxels++;
            brick.min_value = std::min(brick.min_value, value);
            brick.max_value = std::max(brick.max_value, value);
            brick.min_abs_value = std::min(brick.min_abs_value, std::abs(value));
        }

        if (brick.active_voxels == 0) {
            brick.valid = !brick.index_bounds.IsInverted();
            brick.min_value = 0;
            brick.max_value = 0;
            brick.min_abs_value = 0;
        } else {
            brick.valid = !brick.index_bounds.IsInverted();
        }

        m_impl->bricks.push_back(brick);

        const auto* next = reinterpret_cast<const char*>(leaf) + leaf->memUsage();
        leaf = reinterpret_cast<const LeafNodeType*>(next);
    }

    m_impl->bricks_cached = true;
    return m_impl->bricks;
#else
    static const std::vector<ChSDFLeafBrick> empty;
    return empty;
#endif
}

ChSDFProbeResult ChNanoVDBLevelSet::ProbeWorld(const ChVector3d& point_world) const {
    ChSDFProbeResult result;
    result.point_world = point_world;

#if CHRONO_HAS_NANOVDB_SDF
    if (!m_impl->grid) {
        return result;
    }

    const nanovdb::GridMetaData meta(*m_impl->grid);
    const auto accessor = m_impl->grid->getAccessor();
    auto sampler = nanovdb::math::createSampler<1, decltype(accessor), false>(accessor);

    const nanovdb::Vec3d world = ToNanoVector(point_world);
    const auto index = m_impl->grid->worldToIndexF(world);
    const auto grad_index = sampler.gradient(index);
    const auto grad_world = m_impl->grid->indexToWorldGradF(grad_index);

    result.valid = true;
    result.level_set = m_impl->grid->isLevelSet();
    result.point_index = ToChVector(index);
    result.distance = static_cast<double>(sampler(index));
    result.gradient_index = ToChVector(grad_index);
    result.gradient_world = ToChVector(grad_world);
    result.normal_world = result.gradient_world.Length2() > 0 ? result.gradient_world.GetNormalized() : VNULL;
    result.zero_crossing = sampler.zeroCrossing(index);
    result.inside_world_bbox = meta.worldBBox().isInside(world);
    result.inside_index_bbox = meta.indexBBox().template asReal<double>().isInside(index);
#endif

    return result;
}

ChSDFProbeResult ChNanoVDBLevelSet::ProbeIndex(const ChVector3d& point_index) const {
    ChSDFProbeResult result;
    result.point_index = point_index;

#if CHRONO_HAS_NANOVDB_SDF
    if (!m_impl->grid) {
        return result;
    }

    const auto world = m_impl->grid->indexToWorldF(ToNanoVector(point_index));
    result = ProbeWorld(ToChVector(world));
#endif

    return result;
}

const std::string& ChNanoVDBLevelSet::GetLastError() const {
    return m_last_error;
}

}  // end namespace chrono
