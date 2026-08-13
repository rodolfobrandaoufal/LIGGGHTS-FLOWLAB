/* ----------------------------------------------------------------------
   LIGGGHTS feature descriptor to GPU_DEM requested-feature mapping.

   This file is intentionally independent of concrete LIGGGHTS classes. The
   production hook can fill this descriptor from Atom, Pair, Fix and MPI
   state without coupling policy code to the full solver.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_LIGGGHTS_FEATURE_PROBE_H
#define LMP_GPU_DEM_LIGGGHTS_FEATURE_PROBE_H

#include "gpu_execution_policy.h"

#include <string>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct LiggghtsFeatureDescriptor {
  std::string atom_style{"sphere"};
  std::string integrator_style{"nve/sphere"};
  std::string pair_style{"gran/hooke"};
  std::string wall_style{"none"};

  int mpi_size{1};
  bool has_tangential_history{false};
  bool has_rolling_resistance{false};
  bool has_cohesion{false};
  bool has_thermal_contact{false};
  bool has_restart_contact_history{false};
  bool has_particle_insertion_or_deletion{false};
  bool has_multisphere{false};
};

GpuContactModel infer_gpu_contact_model(const std::string &pair_style,
                                        bool has_cohesion);

GpuWallModel infer_gpu_wall_model(const std::string &wall_style);

GpuRequestedFeatures build_gpu_requested_features(
    const LiggghtsFeatureDescriptor &descriptor);

std::string describe_gpu_requested_features(const GpuRequestedFeatures &features);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
