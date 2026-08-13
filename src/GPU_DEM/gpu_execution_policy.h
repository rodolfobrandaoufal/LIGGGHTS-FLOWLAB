/* ----------------------------------------------------------------------
   GPU_DEM execution-policy and feature-gating helpers.

   The production CPU solver remains the reference path. These helpers make
   unsupported GPU features explicit before any integration with the main
   LIGGGHTS run loop.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_EXECUTION_POLICY_H
#define LMP_GPU_DEM_EXECUTION_POLICY_H

#include "gpu_precision.h"

#include <string>
#include <vector>

namespace LAMMPS_NS {
namespace GPU_DEM {

enum class GpuMode {
  Off,
  Auto,
  Strict
};

enum class GpuContactModel {
  HookeNormal,
  HertzMindlin,
  Cohesive,
  Unknown
};

enum class GpuWallModel {
  None,
  Plane,
  Mesh,
  Unknown
};

struct GpuRequestedFeatures {
  bool atom_style_sphere{true};
  bool nve_sphere{true};
  bool single_mpi_rank{true};
  bool requires_tangential_history{false};
  bool requires_rolling_resistance{false};
  bool requires_cohesion{false};
  bool requires_thermal{false};
  bool requires_restart_history{false};
  bool requires_particle_insertion{false};
  bool requires_multisphere{false};
  GpuContactModel contact_model{GpuContactModel::HookeNormal};
  GpuWallModel wall_model{GpuWallModel::None};
};

struct GpuExecutionPolicy {
  GpuMode mode{GpuMode::Off};
  PrecisionMode precision{PrecisionMode::Mixed};
};

struct GpuFeatureValidation {
  bool gpu_allowed{false};
  bool fallback_required{false};
  std::vector<std::string> unsupported_reasons;

  bool ok() const { return unsupported_reasons.empty(); }
};

const char *gpu_mode_name(GpuMode mode);
const char *gpu_contact_model_name(GpuContactModel model);
const char *gpu_wall_model_name(GpuWallModel model);

GpuFeatureValidation validate_gpu_features(const GpuExecutionPolicy &policy,
                                           const GpuRequestedFeatures &features);

void require_gpu_features_or_throw(const GpuExecutionPolicy &policy,
                                   const GpuRequestedFeatures &features);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
