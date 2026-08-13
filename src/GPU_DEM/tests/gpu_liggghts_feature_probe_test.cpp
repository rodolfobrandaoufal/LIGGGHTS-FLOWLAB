/* ----------------------------------------------------------------------
   GPU_DEM LIGGGHTS feature-probe regression test.
------------------------------------------------------------------------- */

#include "gpu_liggghts_feature_probe.h"

#include <cstdio>
#include <cstdlib>
#include <string>

using LAMMPS_NS::GPU_DEM::GpuContactModel;
using LAMMPS_NS::GPU_DEM::GpuExecutionPolicy;
using LAMMPS_NS::GPU_DEM::GpuMode;
using LAMMPS_NS::GPU_DEM::GpuWallModel;
using LAMMPS_NS::GPU_DEM::LiggghtsFeatureDescriptor;
using LAMMPS_NS::GPU_DEM::build_gpu_requested_features;
using LAMMPS_NS::GPU_DEM::describe_gpu_requested_features;
using LAMMPS_NS::GPU_DEM::validate_gpu_features;

namespace {

void fail(const char *message)
{
  std::fprintf(stderr, "%s\n", message);
  std::exit(1);
}

void expect_contains(const std::string &text, const char *needle)
{
  if (text.find(needle) == std::string::npos) {
    std::fprintf(stderr, "expected '%s' to contain '%s'\n",
                 text.c_str(), needle);
    std::exit(1);
  }
}

} // namespace

int main()
{
  GpuExecutionPolicy strict;
  strict.mode = GpuMode::Strict;

  LiggghtsFeatureDescriptor supported;
  supported.atom_style = "sphere";
  supported.integrator_style = "fix nve/sphere";
  supported.pair_style = "gran/hooke/history";
  supported.wall_style = "wall/gran/plane";
  supported.has_tangential_history = false;

  const auto supported_features = build_gpu_requested_features(supported);
  if (!supported_features.atom_style_sphere ||
      !supported_features.nve_sphere ||
      !supported_features.single_mpi_rank ||
      supported_features.contact_model != GpuContactModel::HookeNormal ||
      supported_features.wall_model != GpuWallModel::Plane)
    fail("supported descriptor mapped to wrong GPU feature set");

  const auto supported_validation = validate_gpu_features(strict, supported_features);
  if (!supported_validation.gpu_allowed || supported_validation.fallback_required)
    fail("supported descriptor should pass strict GPU validation");

  const std::string supported_summary = describe_gpu_requested_features(supported_features);
  expect_contains(supported_summary, "contact_model=hooke_normal");
  expect_contains(supported_summary, "wall_model=plane");

  LiggghtsFeatureDescriptor blocked;
  blocked.atom_style = "atomic";
  blocked.integrator_style = "fix nve";
  blocked.pair_style = "gran/hertz/history";
  blocked.wall_style = "mesh/surface";
  blocked.mpi_size = 2;
  blocked.has_tangential_history = true;
  blocked.has_cohesion = true;
  blocked.has_particle_insertion_or_deletion = true;
  blocked.has_thermal_contact = true;

  const auto blocked_features = build_gpu_requested_features(blocked);
  const auto blocked_validation = validate_gpu_features(strict, blocked_features);
  if (blocked_validation.gpu_allowed || !blocked_validation.fallback_required)
    fail("blocked descriptor should be rejected by strict GPU validation");
  if (blocked_validation.unsupported_reasons.size() < 4)
    fail("blocked descriptor did not report enough unsupported feature reasons");

  bool saw_atom = false;
  bool saw_integrator = false;
  bool saw_mpi = false;
  bool saw_cohesion = false;
  bool saw_thermal = false;
  for (const std::string &reason : blocked_validation.unsupported_reasons) {
    saw_atom = saw_atom || reason.find("atom_style sphere") != std::string::npos;
    saw_integrator = saw_integrator || reason.find("nve/sphere") != std::string::npos;
    saw_mpi = saw_mpi || reason.find("multi-rank MPI") != std::string::npos;
    saw_cohesion = saw_cohesion || reason.find("cohesion") != std::string::npos;
    saw_thermal = saw_thermal || reason.find("thermal") != std::string::npos;
  }

  if (!saw_atom || !saw_integrator || !saw_mpi ||
      !saw_cohesion || !saw_thermal)
    fail("blocked descriptor missed one or more expected rejection reasons");

  const std::string blocked_summary = describe_gpu_requested_features(blocked_features);
  expect_contains(blocked_summary, "single_mpi_rank=no");
  expect_contains(blocked_summary, "contact_model=cohesive");
  expect_contains(blocked_summary, "wall_model=mesh");

  std::printf("GPU_DEM LIGGGHTS feature-probe test passed\n");
  return 0;
}
