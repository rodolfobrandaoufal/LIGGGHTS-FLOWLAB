#include "gpu_execution_policy.h"

#include <stdexcept>

namespace LAMMPS_NS {
namespace GPU_DEM {

const char *gpu_mode_name(GpuMode mode)
{
  switch (mode) {
    case GpuMode::Off: return "off";
    case GpuMode::Auto: return "auto";
    case GpuMode::Strict: return "strict";
  }
  return "unknown";
}

const char *gpu_contact_model_name(GpuContactModel model)
{
  switch (model) {
    case GpuContactModel::HookeNormal: return "hooke_normal";
    case GpuContactModel::HertzMindlin: return "hertz_mindlin";
    case GpuContactModel::Cohesive: return "cohesive";
    case GpuContactModel::Unknown: return "unknown";
  }
  return "unknown";
}

const char *gpu_wall_model_name(GpuWallModel model)
{
  switch (model) {
    case GpuWallModel::None: return "none";
    case GpuWallModel::Plane: return "plane";
    case GpuWallModel::Mesh: return "mesh";
    case GpuWallModel::Unknown: return "unknown";
  }
  return "unknown";
}

namespace {

void add_if(bool condition, std::vector<std::string> &reasons, const char *reason)
{
  if (condition) reasons.emplace_back(reason);
}

std::string format_unsupported_message(const GpuFeatureValidation &validation)
{
  std::string message = "GPU_DEM strict mode rejected unsupported feature set";
  for (const std::string &reason : validation.unsupported_reasons) {
    message += "; ";
    message += reason;
  }
  return message;
}

} // namespace

GpuFeatureValidation validate_gpu_features(const GpuExecutionPolicy &policy,
                                           const GpuRequestedFeatures &features)
{
  GpuFeatureValidation validation;

  if (policy.mode == GpuMode::Off) {
    validation.gpu_allowed = false;
    validation.fallback_required = true;
    validation.unsupported_reasons.emplace_back("gpu_mode is off");
    return validation;
  }

  add_if(!features.atom_style_sphere, validation.unsupported_reasons,
         "only atom_style sphere is supported");
  add_if(!features.nve_sphere, validation.unsupported_reasons,
         "only nve/sphere-style integration is supported");
  add_if(!features.single_mpi_rank, validation.unsupported_reasons,
         "multi-rank MPI GPU execution is not implemented");
  add_if(features.requires_cohesion, validation.unsupported_reasons,
         "cohesion is not implemented on GPU");
  add_if(features.requires_thermal, validation.unsupported_reasons,
         "thermal contact is not implemented on GPU");
  add_if(features.requires_restart_history, validation.unsupported_reasons,
         "GPU restart/contact-history serialization is not implemented");
  add_if(features.requires_multisphere, validation.unsupported_reasons,
         "GPU multisphere support is not implemented");

  if (features.contact_model != GpuContactModel::HookeNormal &&
      features.contact_model != GpuContactModel::HertzMindlin)
    validation.unsupported_reasons.emplace_back(
        std::string("unsupported GPU contact model: ") +
        gpu_contact_model_name(features.contact_model));

  if (features.wall_model != GpuWallModel::None &&
      features.wall_model != GpuWallModel::Plane &&
      features.wall_model != GpuWallModel::Mesh)
    validation.unsupported_reasons.emplace_back(
        std::string("unsupported GPU wall model: ") +
        gpu_wall_model_name(features.wall_model));

  validation.gpu_allowed = validation.unsupported_reasons.empty();
  validation.fallback_required = !validation.gpu_allowed;
  return validation;
}

void require_gpu_features_or_throw(const GpuExecutionPolicy &policy,
                                   const GpuRequestedFeatures &features)
{
  const GpuFeatureValidation validation = validate_gpu_features(policy, features);

  if (policy.mode == GpuMode::Strict && !validation.ok())
    throw std::runtime_error(format_unsupported_message(validation));
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
