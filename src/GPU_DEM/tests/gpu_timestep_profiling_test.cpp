/* ----------------------------------------------------------------------
   Standalone GPU_DEM timestep profiling regression test.
------------------------------------------------------------------------- */

#include "gpu_dem_context.h"
#include "gpu_particle_data.h"
#include "gpu_timestep.h"
#include "gpu_wall_contact.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuNeighborList;
using LAMMPS_NS::GPU_DEM::GpuParticleData;
using LAMMPS_NS::GPU_DEM::GpuPlaneWalls;
using LAMMPS_NS::GPU_DEM::GpuTimestepParams;
using LAMMPS_NS::GPU_DEM::GpuTimestepTiming;
using LAMMPS_NS::GPU_DEM::HookeNormalParams;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::HostPlaneWalls;
using LAMMPS_NS::GPU_DEM::NeighborBuildParams;
using LAMMPS_NS::GPU_DEM::PlaneWall;
using LAMMPS_NS::GPU_DEM::velocity_verlet_hooke_step;

namespace {

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

void require_nonnegative(float value, const char *label)
{
  if (value >= 0.0f) return;
  std::fprintf(stderr, "%s timing was negative: %.6f ms\n", label, value);
  std::exit(1);
}

HostParticleData make_particles()
{
  HostParticleData host;
  host.resize(3);
  host.position_x = {0.9, 1.95, 3.0};
  host.position_y = {1.0, 1.0, 1.0};
  host.position_z = {1.0, 1.0, 1.0};
  host.velocity_x = {-1.0, -0.8, 0.0};
  host.velocity_y = {0.0, 0.0, 0.0};
  host.velocity_z = {0.0, 0.0, 0.0};
  host.radius = {0.5, 0.5, 0.5};

  for (std::size_t i = 0; i < host.size(); ++i) {
    host.omega_x[i] = host.omega_y[i] = host.omega_z[i] = 0.0;
    host.force_x[i] = host.force_y[i] = host.force_z[i] = 0.0;
    host.torque_x[i] = host.torque_y[i] = host.torque_z[i] = 0.0;
    host.mass[i] = 1.0;
    host.density[i] = 1.0;
    host.type[i] = 1;
    host.mask[i] = 1;
    host.tag[i] = static_cast<long long>(i + 1);
    host.image_flags[i] = 0;
  }
  return host;
}

} // namespace

int main()
{
  GpuDemContext context(0);

  HostParticleData host = make_particles();
  GpuParticleData particles;
  particles.copy_from_host(host, context.stream());

  NeighborBuildParams neighbor_params;
  neighbor_params.box_min[0] = neighbor_params.box_min[1] = neighbor_params.box_min[2] = 0.0;
  neighbor_params.box_max[0] = neighbor_params.box_max[1] = neighbor_params.box_max[2] = 4.0;
  neighbor_params.cell_size = 1.1;
  neighbor_params.skin = 0.0;

  GpuNeighborList neighbors;
  neighbors.allocate(host.size(), 4, num_cells(neighbor_params));

  HostPlaneWalls host_walls;
  PlaneWall lower_x;
  lower_x.normal[0] = 1.0;
  lower_x.normal[1] = 0.0;
  lower_x.normal[2] = 0.0;
  lower_x.offset = 0.0;
  host_walls.walls.push_back(lower_x);

  GpuPlaneWalls walls;
  walls.copy_from_host(host_walls, context.stream());

  HookeNormalParams contact_params;
  contact_params.normal_stiffness = 1000.0;
  contact_params.normal_damping = 2.0;

  GpuTimestepParams timestep_params;
  timestep_params.dt = 1.0e-4;
  timestep_params.groupbit = 1;
  timestep_params.dimension = 3;
  timestep_params.one_plus_added_mass = 1.0;

  GpuTimestepTiming timing;
  velocity_verlet_hooke_step(particles, neighbors, neighbor_params, contact_params,
                             &walls, timestep_params, timing, context.stream());
  context.synchronize();

  require_nonnegative(timing.reset_ms, "reset");
  require_nonnegative(timing.external_force_ms, "external force");
  require_nonnegative(timing.neighbor_build_ms, "neighbor build");
  require_nonnegative(timing.particle_contact_ms, "particle contact");
  require_nonnegative(timing.wall_contact_ms, "wall contact");
  require_nonnegative(timing.initial_integrate_ms, "initial integrate");
  require_nonnegative(timing.final_integrate_ms, "final integrate");
  require_nonnegative(timing.total_ms, "total");

  if (!(timing.total_ms > 0.0f)) {
    std::fprintf(stderr, "total timestep timing was not positive\n");
    return 1;
  }

  std::printf("GPU_DEM timestep profiling test passed on device %d (%s): total=%.6f ms\n",
              context.device_id(), context.device_properties().name, timing.total_ms);
  return 0;
}
