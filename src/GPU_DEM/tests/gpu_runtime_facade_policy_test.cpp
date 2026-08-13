/* ----------------------------------------------------------------------
   GPU_DEM runtime facade policy regression test.
------------------------------------------------------------------------- */

#include "gpu_runtime_facade.h"

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>

using LAMMPS_NS::GPU_DEM::GpuDemDriverConfig;
using LAMMPS_NS::GPU_DEM::GpuExecutionPolicy;
using LAMMPS_NS::GPU_DEM::GpuMode;
using LAMMPS_NS::GPU_DEM::GpuRuntimeFacade;
using LAMMPS_NS::GPU_DEM::LiggghtsFeatureDescriptor;

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

LiggghtsFeatureDescriptor supported_descriptor()
{
  LiggghtsFeatureDescriptor descriptor;
  descriptor.atom_style = "sphere";
  descriptor.integrator_style = "nve/sphere";
  descriptor.pair_style = "gran/hooke";
  descriptor.wall_style = "wall/gran/plane";
  descriptor.mpi_size = 1;
  return descriptor;
}

LiggghtsFeatureDescriptor unsupported_descriptor()
{
  LiggghtsFeatureDescriptor descriptor = supported_descriptor();
  descriptor.pair_style = "gran/hertz/history";
  descriptor.wall_style = "mesh/surface";
  descriptor.mpi_size = 4;
  descriptor.has_tangential_history = true;
  return descriptor;
}

} // namespace

int main()
{
  GpuRuntimeFacade facade;
  GpuDemDriverConfig config;

  GpuExecutionPolicy off;
  off.mode = GpuMode::Off;
  const auto &off_status = facade.configure(off, supported_descriptor(), config, false);
  if (off_status.gpu_active || !off_status.cpu_fallback || facade.driver())
    fail("gpu_mode off should select CPU fallback and create no driver");
  expect_contains(off_status.message, "gpu_mode=off");
  expect_contains(off_status.message, "CPU fallback");

  GpuExecutionPolicy automatic;
  automatic.mode = GpuMode::Auto;
  const auto &auto_status = facade.configure(automatic, unsupported_descriptor(), config, false);
  if (auto_status.gpu_active || !auto_status.cpu_fallback || facade.driver())
    fail("gpu_mode auto should select explicit CPU fallback for unsupported features");
  expect_contains(auto_status.message, "multi-rank MPI");

  GpuExecutionPolicy strict;
  strict.mode = GpuMode::Strict;
  bool strict_threw = false;
  try {
    facade.configure(strict, unsupported_descriptor(), config, false);
  } catch (const std::runtime_error &error) {
    strict_threw = true;
    expect_contains(error.what(), "strict mode rejected");
    expect_contains(error.what(), "multi-rank MPI");
  }
  if (!strict_threw)
    fail("gpu_mode strict should throw for unsupported features");

  const auto &strict_status = facade.configure(strict, supported_descriptor(), config, false);
  if (!strict_status.gpu_active || strict_status.cpu_fallback || facade.driver())
    fail("gpu_mode strict should activate GPU for supported features without creating driver");
  expect_contains(strict_status.message, "GPU_DEM active");

  std::printf("GPU_DEM runtime facade policy test passed\n");
  return 0;
}
