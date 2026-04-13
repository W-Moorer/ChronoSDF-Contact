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

#include "chrono/collision/ChCollisionShapeSDF.h"

namespace chrono {

// Register into the object factory, to enable run-time dynamic creation and persistence
CH_FACTORY_REGISTER(ChCollisionShapeSDF)
CH_UPCASTING(ChCollisionShapeSDF, ChCollisionShape)

ChCollisionShapeSDF::ChCollisionShapeSDF() : ChCollisionShape(Type::SDF) {}

ChCollisionShapeSDF::ChCollisionShapeSDF(std::shared_ptr<ChContactMaterial> material) : ChCollisionShape(Type::SDF, material) {}

ChCollisionShapeSDF::ChCollisionShapeSDF(std::shared_ptr<ChContactMaterial> material,
                                         const std::string& filename,
                                         const std::string& grid_name)
    : ChCollisionShape(Type::SDF, material) {
    LoadLevelSet(filename, grid_name);
}

bool ChCollisionShapeSDF::LoadLevelSet(const std::string& filename, const std::string& grid_name) {
    auto level_set = std::make_shared<ChNanoVDBLevelSet>();
    if (!level_set->Load(filename, grid_name)) {
        m_level_set.reset();
        m_level_set_filename = filename;
        m_grid_name = grid_name;
        m_last_error = level_set->GetLastError();
        return false;
    }

    m_level_set = level_set;
    m_level_set_filename = filename;
    m_grid_name = grid_name;
    m_last_error.clear();
    return true;
}

void ChCollisionShapeSDF::SetLevelSet(std::shared_ptr<ChNanoVDBLevelSet> level_set,
                                      const std::string& filename,
                                      const std::string& grid_name) {
    m_level_set = level_set;
    m_level_set_filename = filename;
    m_grid_name = grid_name;
    m_last_error = (!m_level_set || m_level_set->IsLoaded()) ? std::string() : m_level_set->GetLastError();
}

void ChCollisionShapeSDF::ClearLevelSet() {
    m_level_set.reset();
    m_level_set_filename.clear();
    m_grid_name.clear();
    m_last_error.clear();
}

ChNanoVDBGridInfo ChCollisionShapeSDF::GetGridInfo() const {
    return m_level_set ? m_level_set->GetInfo() : ChNanoVDBGridInfo();
}

const std::vector<ChSDFLeafBrick>& ChCollisionShapeSDF::GetLeafBricks() const {
    static const std::vector<ChSDFLeafBrick> empty;
    return m_level_set ? m_level_set->GetLeafBricks() : empty;
}

ChSDFProbeResult ChCollisionShapeSDF::ProbeLocal(const ChVector3d& point_local) const {
    if (!m_level_set) {
        ChSDFProbeResult result;
        result.point_world = point_local;
        return result;
    }

    return m_level_set->ProbeWorld(point_local);
}

ChSDFPotentialFieldProbe ChCollisionShapeSDF::ProbePotentialLocal(const ChVector3d& point_local) const {
    const ChSDFProbeResult sdf_probe = ProbeLocal(point_local);
    return ChSDFPotentialFieldEvaluator::Evaluate(sdf_probe, GetGridInfo(), m_potential_field_settings);
}

ChSDFContactPatch ChCollisionShapeSDF::SamplePatchLocal(
    const ChFrame<>& patch_frame_local,
    const ChSDFContactPatchSampler::Settings& settings) const {
    if (!m_level_set) {
        return ChSDFContactPatch();
    }

    return ChSDFContactPatchSampler::SamplePlanePatch(*m_level_set, patch_frame_local, settings);
}

ChAABB ChCollisionShapeSDF::GetBoundingBox() const {
    const auto info = GetGridInfo();
    return info.valid ? info.world_bounds : ChAABB();
}

void ChCollisionShapeSDF::ArchiveOut(ChArchiveOut& archive_out) {
    archive_out.VersionWrite<ChCollisionShapeSDF>();
    ChCollisionShape::ArchiveOut(archive_out);
    archive_out << CHNVP(m_level_set_filename);
    archive_out << CHNVP(m_grid_name);
    const double potential_field_modulus = m_potential_field_settings.modulus;
    const double potential_field_depth_scale = m_potential_field_settings.depth_scale;
    const double potential_field_depth_cap = m_potential_field_settings.depth_cap;
    const double potential_field_support_margin = m_potential_field_settings.support_margin;
    const bool potential_field_clamp_outside_to_zero = m_potential_field_settings.clamp_outside_to_zero;
    archive_out << CHNVP(potential_field_modulus);
    archive_out << CHNVP(potential_field_depth_scale);
    archive_out << CHNVP(potential_field_depth_cap);
    archive_out << CHNVP(potential_field_support_margin);
    archive_out << CHNVP(potential_field_clamp_outside_to_zero);
}

void ChCollisionShapeSDF::ArchiveIn(ChArchiveIn& archive_in) {
    int version = archive_in.VersionRead<ChCollisionShapeSDF>();
    ChCollisionShape::ArchiveIn(archive_in);
    archive_in >> CHNVP(m_level_set_filename);
    archive_in >> CHNVP(m_grid_name);
    if (version >= 1) {
        double potential_field_modulus = m_potential_field_settings.modulus;
        double potential_field_depth_scale = m_potential_field_settings.depth_scale;
        double potential_field_depth_cap = m_potential_field_settings.depth_cap;
        double potential_field_support_margin = m_potential_field_settings.support_margin;
        bool potential_field_clamp_outside_to_zero = m_potential_field_settings.clamp_outside_to_zero;
        archive_in >> CHNVP(potential_field_modulus);
        archive_in >> CHNVP(potential_field_depth_scale);
        archive_in >> CHNVP(potential_field_depth_cap);
        archive_in >> CHNVP(potential_field_support_margin);
        archive_in >> CHNVP(potential_field_clamp_outside_to_zero);
        m_potential_field_settings.modulus = potential_field_modulus;
        m_potential_field_settings.depth_scale = potential_field_depth_scale;
        m_potential_field_settings.depth_cap = potential_field_depth_cap;
        m_potential_field_settings.support_margin = potential_field_support_margin;
        m_potential_field_settings.clamp_outside_to_zero = potential_field_clamp_outside_to_zero;
    } else {
        m_potential_field_settings = ChSDFPotentialFieldSettings();
    }

    m_last_error.clear();
    m_level_set.reset();
    if (!m_level_set_filename.empty()) {
        LoadLevelSet(m_level_set_filename, m_grid_name);
    }
}

}  // end namespace chrono
