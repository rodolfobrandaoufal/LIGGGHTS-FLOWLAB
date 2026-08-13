/* ----------------------------------------------------------------------
   Standalone GPU_DEM gravity/free-fall regression test.
------------------------------------------------------------------------- */

#include "gpu_dem_context.h"
#include "gpu_particle_data.h"
#include "gpu_timestep.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuNeighborList;
using LAMMPS_NS::GPU_DEM::GpuParticleData;
using LAMMPS_NS::GPU_DEM::GpuTimestepParams;
using LAMMPS_NS::GPU_DEM::HookeNormalParams;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::NeighborBuildParams;
using LAMMPS_NS::GPU_DEM::velocity_verlet_hooke_step;

namespace {

void expect_near(double actual, double expected, double tol, const char *label)
{
  if (std::fabs(actual - expected) <= tol) return;
  std::fprintf(stderr, "%s mismatch: actual %.17g expected %.17g tolerance %.3g\n",
               label, actual, expected, tol);
  std::exit(1);
}

int num_cells(const NeighborBuildParams &params)
{
  int cells = 1;
  for (int d = 0; d < 3; ++d) {
    const double length = params.box_max[d] - params.box_min[d];
    const int count = static_cast<int>(std::floor(length / params.cell_size));
    cells *= count > 1 ? count : 1;
  }
  return cells;
}

HostParticleData make_single_particle()
{
  HostParticleData host;
  host.resize(1);
  host.position_x[0] = 0.0;
  host.position_y[0] = 0.0;
  host.position_z[0] = 10.0;
  host.velocity_x[0] = 0.0;
  host.velocity_y[0] = 0.0;
  host.velocity_z[0] = 2.0;
  host.omega_x[0] = host.omega_y[0] = host.omega_z[0] = 0.0;
  host.force_x[0] = host.force_y[0] = host.force_z[0] = 0.0;
  host.torque_x[0] = host.torque_y[0] = host.torque_z[0] = 0.0;
  host.radius[0] = 0.25;
  host.mass[0] = 2.5;
  host.density[0] = 1.0;
  host.type[0] = 1;
  host.mask[0] = 1;
  host.tag[0] = 1;
  host.image_flags[0] = 0;
  return host;
}

} // namespace

int main()
{
  GpuDemContext context(0);

  HostParticleData host = make_single_particle();
  GpuParticleData particles;
  particles.copy_from_host(host, context.stream());

  NeighborBuildParams neighbor_params;
  neighbor_params.box_min[0] = neighbor_params.box_min[1] = neighbor_params.box_min[2] = -100.0;
  neighbor_params.box_max[0] = neighbor_params.box_max[1] = neighbor_params.box_max[2] = 100.0;
  neighbor_params.cell_size = 2.0;
  neighbor_params.skin = 0.0;

  GpuNeighborList neighbors;
  neighbors.allocate(host.size(), 1, num_cells(neighbor_params));

  HookeNormalParams contact_params;
  contact_params.normal_stiffness = 0.0;
  contact_params.normal_damping = 0.0;

  GpuTimestepParams timestep_params;
  timestep_params.dt = 1.0e-4;
  timestep_params.groupbit = 1;
  timestep_params.dimension = 3;
  timestep_params.one_plus_added_mass = 1.0;
  timestep_params.gravity.enabled = true;
  timestep_params.gravity.acceleration[0] = 0.0;
  timestep_params.gravity.acceleration[1] = 0.0;
  timestep_params.gravity.acceleration[2] = -9.81;

  const int steps = 2000;
  for (int step = 0; step < steps; ++step) {
    velocity_verlet_hooke_step(particles, neighbors, neighbor_params, contact_params,
                               nullptr, timestep_params, context.stream());
  }

  HostParticleData result;
  particles.copy_to_host(result, context.stream());
  context.synchronize();

  const double elapsed = steps * timestep_params.dt;
  const double expected_z = host.position_z[0] + host.velocity_z[0] * elapsed +
                            0.5 * timestep_params.gravity.acceleration[2] *
                                elapsed * elapsed;
  const double expected_vz = host.velocity_z[0] +
                             timestep_params.gravity.acceleration[2] * elapsed;
  const double expected_ke = 0.5 * host.mass[0] * expected_vz * expected_vz;
  const double actual_ke = 0.5 * result.mass[0] *
                           result.velocity_z[0] * result.velocity_z[0];

  expect_near(result.position_x[0], 0.0, 1.0e-14, "freefall position_x");
  expect_near(result.position_y[0], 0.0, 1.0e-14, "freefall position_y");
  expect_near(result.position_z[0], expected_z, 2.0e-12, "freefall position_z");
  expect_near(result.velocity_z[0], expected_vz, 2.0e-12, "freefall velocity_z");
  expect_near(actual_ke, expected_ke, 2.0e-12, "freefall kinetic energy");

  std::printf("GPU_DEM gravity free-fall test passed on device %d (%s)\n",
              context.device_id(), context.device_properties().name);
  return 0;
}
