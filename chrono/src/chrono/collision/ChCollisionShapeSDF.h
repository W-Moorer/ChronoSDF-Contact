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

#ifndef CH_COLLISION_SHAPE_SDF_H
#define CH_COLLISION_SHAPE_SDF_H

#include <memory>
#include <string>

#include "chrono/collision/ChCollisionShape.h"
#include "chrono/collision/sdf/ChSDFContactPatch.h"
#include "chrono/collision/sdf/ChSDFPotentialField.h"

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// Collision shape backed by a sparse signed distance field.
/// This shape exposes a reusable Chrono-side SDF resource and patch-sampling interface.
/// Native collision-system population is not yet implemented for this type.
class ChApi ChCollisionShapeSDF : public ChCollisionShape {
  public:
    ChCollisionShapeSDF();
    explicit ChCollisionShapeSDF(std::shared_ptr<ChContactMaterial> material);
    ChCollisionShapeSDF(std::shared_ptr<ChContactMaterial> material,
                        const std::string& filename,
                        const std::string& grid_name = "");

    ~ChCollisionShapeSDF() {}

    /// Load a NanoVDB level set from file.
    bool LoadLevelSet(const std::string& filename, const std::string& grid_name = "");

    /// Attach an already created level set resource to this shape.
    void SetLevelSet(std::shared_ptr<ChNanoVDBLevelSet> level_set,
                     const std::string& filename = "",
                     const std::string& grid_name = "");

    /// Release the currently attached level set resource.
    void ClearLevelSet();

    /// Return the attached level set resource, if any.
    std::shared_ptr<ChNanoVDBLevelSet> GetLevelSet() const { return m_level_set; }

    /// Return true if the shape currently holds a loaded level set.
    bool IsLoaded() const { return m_level_set && m_level_set->IsLoaded(); }

    /// Return the input filename used to load this shape, if available.
    const std::string& GetLevelSetFilename() const { return m_level_set_filename; }

    /// Return the grid name used to load this shape, if available.
    const std::string& GetGridName() const { return m_grid_name; }

    /// Return information about the currently attached level set.
    ChNanoVDBGridInfo GetGridInfo() const;

    /// Return the sparse NanoVDB leaf bricks attached to this shape.
    const std::vector<ChSDFLeafBrick>& GetLeafBricks() const;

    /// Probe the SDF using coordinates expressed in the local frame of this shape.
    ChSDFProbeResult ProbeLocal(const ChVector3d& point_local) const;

    /// Configure the shape-fixed potential pressure field derived from the local SDF depth.
    void SetPotentialFieldSettings(const ChSDFPotentialFieldSettings& settings) { m_potential_field_settings = settings; }
    const ChSDFPotentialFieldSettings& GetPotentialFieldSettings() const { return m_potential_field_settings; }

    /// Probe the shape-fixed potential pressure field using coordinates expressed in the local frame of this shape.
    ChSDFPotentialFieldProbe ProbePotentialLocal(const ChVector3d& point_local) const;

    /// Sample a local patch frame against the SDF carried by this shape.
    ChSDFContactPatch SamplePatchLocal(const ChFrame<>& patch_frame_local,
                                       const ChSDFContactPatchSampler::Settings& settings) const;

    /// Return the most recent load/probe error attached to this shape.
    const std::string& GetLastError() const { return m_last_error; }

    /// Get the shape bounding box in local coordinates.
    virtual ChAABB GetBoundingBox() const override;

    /// Method to allow serialization of transient data to archives.
    virtual void ArchiveOut(ChArchiveOut& archive_out) override;

    /// Method to allow de-serialization of transient data from archives.
    virtual void ArchiveIn(ChArchiveIn& archive_in) override;

  private:
    std::shared_ptr<ChNanoVDBLevelSet> m_level_set;
    std::string m_level_set_filename;
    std::string m_grid_name;
    std::string m_last_error;
    ChSDFPotentialFieldSettings m_potential_field_settings;
};

/// @} chrono_collision

CH_CLASS_VERSION(ChCollisionShapeSDF, 1)

}  // end namespace chrono

#endif
