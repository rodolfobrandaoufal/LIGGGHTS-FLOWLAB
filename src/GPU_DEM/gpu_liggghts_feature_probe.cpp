#include "gpu_liggghts_feature_probe.h"

#include <algorithm>
#include <cctype>
#include <sstream>

namespace LAMMPS_NS {
namespace GPU_DEM {

namespace {

std::string lowercase(std::string text)
{
  std::transform(text.begin(), text.end(), text.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return text;
}

bool contains_token(const std::string &text, const char *token)
{
  return lowercase(text).find(token) != std::string::npos;
}

const char *yes_no(bool value)
{
  return value ? "yes" : "no";
}

} // namespace

GpuContactModel infer_gpu_contact_model(const std::string &pair_style,
                                        bool has_cohesion)
{
  if (has_cohesion) return GpuContactModel::Cohesive;

  const std::string style = lowercase(pair_style);
  if (style.empty() || style == "none") return GpuContactModel::Unknown;
  if (contains_token(style, "hooke")) return GpuContactModel::HookeNormal;
  if (contains_token(style, "hertz") || contains_token(style, "mindlin"))
    return GpuContactModel::HertzMindlin;
  return GpuContactModel::Unknown;
}

GpuWallModel infer_gpu_wall_model(const std::string &wall_style)
{
  const std::string style = lowercase(wall_style);
  if (style.empty() || style == "none") return GpuWallModel::None;
  if (contains_token(style, "mesh")) return GpuWallModel::Mesh;
  if (contains_token(style, "plane") || contains_token(style, "wall/gran"))
    return GpuWallModel::Plane;
  return GpuWallModel::Unknown;
}

GpuRequestedFeatures build_gpu_requested_features(
    const LiggghtsFeatureDescriptor &descriptor)
{
  GpuRequestedFeatures features;
  features.atom_style_sphere = lowercase(descriptor.atom_style) == "sphere";
  features.nve_sphere = contains_token(descriptor.integrator_style, "nve/sphere");
  features.single_mpi_rank = descriptor.mpi_size == 1;
  features.requires_tangential_history = descriptor.has_tangential_history;
  features.requires_rolling_resistance = descriptor.has_rolling_resistance;
  features.requires_cohesion = descriptor.has_cohesion;
  features.requires_thermal = descriptor.has_thermal_contact;
  features.requires_restart_history = descriptor.has_restart_contact_history;
  features.requires_particle_insertion = descriptor.has_particle_insertion_or_deletion;
  features.requires_multisphere = descriptor.has_multisphere;
  features.contact_model = infer_gpu_contact_model(descriptor.pair_style,
                                                   descriptor.has_cohesion);
  features.wall_model = infer_gpu_wall_model(descriptor.wall_style);
  return features;
}

std::string describe_gpu_requested_features(const GpuRequestedFeatures &features)
{
  std::ostringstream out;
  out << "atom_style_sphere=" << yes_no(features.atom_style_sphere)
      << " nve_sphere=" << yes_no(features.nve_sphere)
      << " single_mpi_rank=" << yes_no(features.single_mpi_rank)
      << " contact_model=" << gpu_contact_model_name(features.contact_model)
      << " wall_model=" << gpu_wall_model_name(features.wall_model)
      << " tangential_history=" << yes_no(features.requires_tangential_history)
      << " rolling=" << yes_no(features.requires_rolling_resistance)
      << " cohesion=" << yes_no(features.requires_cohesion)
      << " thermal=" << yes_no(features.requires_thermal)
      << " restart_history=" << yes_no(features.requires_restart_history)
      << " insertion_deletion=" << yes_no(features.requires_particle_insertion)
      << " multisphere=" << yes_no(features.requires_multisphere);
  return out.str();
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
