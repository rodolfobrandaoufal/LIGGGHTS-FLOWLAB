/* ----------------------------------------------------------------------
   GPU_DEM runtime facade CUDA smoke test.
------------------------------------------------------------------------- */

#include "gpu_runtime_facade.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using LAMMPS_NS::GPU_DEM::GpuDemDriverConfig;
using LAMMPS_NS::GPU_DEM::GpuExecutionPolicy;
using LAMMPS_NS::GPU_DEM::GpuMode;
using LAMMPS_NS::GPU_DEM::GpuRuntimeFacade;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::LiggghtsFeatureDescriptor;

namespace {

void fail(const char *message)
{
  std::fprintf(stderr, "%s\n", message);
  std::exit(1);
}

HostParticleData make_particles()
{
  HostParticleData host;
  host.resize(1);
  host.position_x[0] = 0.25;
  host.position_y[0] = 0.5;
  host.position_z[0] = 0.75;
  host.velocity_x[0] = 1.0;
  host.velocity_y[0] = 2.0;
  host.velocity_z[0] = 3.0;
  host.omega_x[0] = host.omega_y[0] = host.omega_z[0] = 0.0;
  host.force_x[0] = host.force_y[0] = host.force_z[0] = 0.0;
  host.torque_x[0] = host.torque_y[0] = host.torque_z[0] = 0.0;
  host.radius[0] = 0.1;
  host.mass[0] = 1.0;
  host.density[0] = 2500.0;
  host.type[0] = 1;
  host.mask[0] = 1;
  host.tag[0] = 1;
  host.image_flags[0] = 0;
  return host;
}

void expect_close(double actual, double expected, const char *label)
{
  if (std::fabs(actual - expected) > 1.0e-14) {
    std::fprintf(stderr, "%s mismatch: actual %.17g expected %.17g\n",
                 label, actual, expected);
    std::exit(1);
  }
}

} // namespace

int main()
{
  GpuExecutionPolicy policy;
  policy.mode = GpuMode::Strict;

  LiggghtsFeatureDescriptor descriptor;
  descriptor.atom_style = "sphere";
  descriptor.integrator_style = "nve/sphere";
  descriptor.pair_style = "gran/hooke";
  descriptor.wall_style = "none";
  descriptor.mpi_size = 1;

  GpuDemDriverConfig config;
  config.device_id = 0;
  config.pair_capacity = 1;
  config.num_cells = 1;

  GpuRuntimeFacade facade;
  const auto &status = facade.configure(policy, descriptor, config, true);
  if (!status.gpu_active || status.cpu_fallback || !facade.driver())
    fail("supported strict runtime did not create an active GPU driver");

  HostParticleData input = make_particles();
  facade.driver()->upload_particles(input);
  facade.driver()->synchronize();

  HostParticleData output;
  facade.driver()->download_particles(output);
  facade.driver()->synchronize();

  expect_close(output.position_x[0], input.position_x[0], "position_x");
  expect_close(output.position_y[0], input.position_y[0], "position_y");
  expect_close(output.velocity_z[0], input.velocity_z[0], "velocity_z");

  std::printf("GPU_DEM runtime facade CUDA test passed on device %d (%s)\n",
              facade.driver()->context().device_id(),
              facade.driver()->context().device_properties().name);
  return 0;
}
