/* ----------------------------------------------------------------------
   Standalone GPU_DEM timestep-driver regression test.
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
using LAMMPS_NS::GPU_DEM::HookeNormalParams;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::HostPlaneWalls;
using LAMMPS_NS::GPU_DEM::NeighborBuildParams;
using LAMMPS_NS::GPU_DEM::PlaneWall;
using LAMMPS_NS::GPU_DEM::velocity_verlet_hooke_step;

namespace {

struct CpuState {
  double x[2]{};
  double v[2]{};
  double f[2]{};
  double radius[2]{0.5, 0.5};
  double mass[2]{1.0, 1.0};
};

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

void add_force_pair(CpuState &state, const HookeNormalParams &params)
{
  const double dx = state.x[1] - state.x[0];
  const double distance = std::fabs(dx);
  if (distance <= 0.0) return;

  const double overlap = state.radius[0] + state.radius[1] - distance;
  if (overlap <= 0.0) return;

  const double nx = dx / distance;
  const double relative_normal_velocity = (state.v[1] - state.v[0]) * nx;
  double normal_force =
      params.normal_stiffness * overlap - params.normal_damping * relative_normal_velocity;
  if (normal_force < 0.0) normal_force = 0.0;

  state.f[0] += -normal_force * nx;
  state.f[1] += normal_force * nx;
}

void add_force_wall(CpuState &state, const HookeNormalParams &params)
{
  for (int i = 0; i < 2; ++i) {
    const double overlap = state.radius[i] - state.x[i];
    if (overlap <= 0.0) continue;

    double normal_force = params.normal_stiffness * overlap - params.normal_damping * state.v[i];
    if (normal_force < 0.0) normal_force = 0.0;
    state.f[i] += normal_force;
  }
}

void compute_cpu_force(CpuState &state, const HookeNormalParams &params)
{
  state.f[0] = 0.0;
  state.f[1] = 0.0;
  add_force_pair(state, params);
  add_force_wall(state, params);
}

void cpu_step(CpuState &state, const HookeNormalParams &params, double dt)
{
  const double dtf = 0.5 * dt;
  compute_cpu_force(state, params);
  state.v[0] += dtf * state.f[0] / state.mass[0];
  state.v[1] += dtf * state.f[1] / state.mass[1];
  state.x[0] += dt * state.v[0];
  state.x[1] += dt * state.v[1];
  compute_cpu_force(state, params);
  state.v[0] += dtf * state.f[0] / state.mass[0];
  state.v[1] += dtf * state.f[1] / state.mass[1];
}

HostParticleData make_particles(const CpuState &state)
{
  HostParticleData host;
  host.resize(2);
  host.position_x = {state.x[0], state.x[1]};
  host.position_y = {1.0, 1.0};
  host.position_z = {1.0, 1.0};
  host.velocity_x = {state.v[0], state.v[1]};
  host.velocity_y = {0.0, 0.0};
  host.velocity_z = {0.0, 0.0};
  host.radius = {state.radius[0], state.radius[1]};

  for (std::size_t i = 0; i < host.size(); ++i) {
    host.omega_x[i] = host.omega_y[i] = host.omega_z[i] = 0.0;
    host.force_x[i] = host.force_y[i] = host.force_z[i] = 0.0;
    host.torque_x[i] = host.torque_y[i] = host.torque_z[i] = 0.0;
    host.mass[i] = state.mass[i];
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

  HookeNormalParams contact_params;
  contact_params.normal_stiffness = 1000.0;
  contact_params.normal_damping = 2.0;

  CpuState cpu;
  cpu.x[0] = 0.9;
  cpu.x[1] = 2.05;
  cpu.v[0] = -1.0;
  cpu.v[1] = -1.2;

  NeighborBuildParams neighbor_params;
  neighbor_params.box_min[0] = neighbor_params.box_min[1] = neighbor_params.box_min[2] = 0.0;
  neighbor_params.box_max[0] = neighbor_params.box_max[1] = neighbor_params.box_max[2] = 4.0;
  neighbor_params.cell_size = 1.1;
  neighbor_params.skin = 0.0;

  HostPlaneWalls host_walls;
  PlaneWall lower_x;
  lower_x.normal[0] = 1.0;
  lower_x.normal[1] = 0.0;
  lower_x.normal[2] = 0.0;
  lower_x.offset = 0.0;
  host_walls.walls.push_back(lower_x);

  GpuPlaneWalls walls;
  walls.copy_from_host(host_walls, context.stream());

  HostParticleData host = make_particles(cpu);
  GpuParticleData particles;
  particles.copy_from_host(host, context.stream());

  GpuNeighborList neighbors;
  neighbors.allocate(host.size(), 2, num_cells(neighbor_params));

  GpuTimestepParams timestep_params;
  timestep_params.dt = 1.0e-4;
  timestep_params.groupbit = 1;
  timestep_params.dimension = 3;
  timestep_params.one_plus_added_mass = 1.0;

  const int steps = 10000;
  for (int step = 0; step < steps; ++step) {
    cpu_step(cpu, contact_params, timestep_params.dt);
    velocity_verlet_hooke_step(particles, neighbors, neighbor_params, contact_params,
                               &walls, timestep_params, context.stream());
  }

  HostParticleData result;
  particles.copy_to_host(result, context.stream());
  context.synchronize();

  expect_near(result.position_x[0], cpu.x[0], 2.0e-12, "driver position_x[0]");
  expect_near(result.position_x[1], cpu.x[1], 2.0e-12, "driver position_x[1]");
  expect_near(result.velocity_x[0], cpu.v[0], 2.0e-12, "driver velocity_x[0]");
  expect_near(result.velocity_x[1], cpu.v[1], 2.0e-12, "driver velocity_x[1]");

  if (!(result.position_x[0] > 0.5)) {
    std::fprintf(stderr, "driver wall contact allowed particle-wall penetration\n");
    return 1;
  }

  std::printf("GPU_DEM timestep driver test passed on device %d (%s)\n",
              context.device_id(), context.device_properties().name);
  return 0;
}
