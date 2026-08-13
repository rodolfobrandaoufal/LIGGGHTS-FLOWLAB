/* ----------------------------------------------------------------------
   Standalone GPU_DEM damped Hookean collision regression test.
------------------------------------------------------------------------- */

#include "gpu_contact_pipeline.h"
#include "gpu_dem_context.h"
#include "gpu_integrator.h"
#include "gpu_neighbor_builder.h"
#include "gpu_particle_data.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuNeighborList;
using LAMMPS_NS::GPU_DEM::GpuParticleData;
using LAMMPS_NS::GPU_DEM::HookeNormalParams;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::NeighborBuildParams;
using LAMMPS_NS::GPU_DEM::build_contact_pairs_uniform_grid;
using LAMMPS_NS::GPU_DEM::compute_hooke_normal_forces;
using LAMMPS_NS::GPU_DEM::nve_sphere_final_integrate;
using LAMMPS_NS::GPU_DEM::nve_sphere_initial_integrate;
using LAMMPS_NS::GPU_DEM::reset_forces;

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

double kinetic_energy(const CpuState &state)
{
  return 0.5 * state.mass[0] * state.v[0] * state.v[0] +
         0.5 * state.mass[1] * state.v[1] * state.v[1];
}

void compute_cpu_force(CpuState &state, const HookeNormalParams &params)
{
  state.f[0] = 0.0;
  state.f[1] = 0.0;

  const double dx = state.x[1] - state.x[0];
  const double distance = std::fabs(dx);
  if (distance <= 0.0) return;

  const double overlap = state.radius[0] + state.radius[1] - distance;
  if (overlap <= 0.0) return;

  const double nx = dx / distance;
  const double normal_relative_velocity = (state.v[1] - state.v[0]) * nx;
  double normal_force =
      params.normal_stiffness * overlap - params.normal_damping * normal_relative_velocity;
  if (normal_force < 0.0) normal_force = 0.0;

  state.f[0] += -normal_force * nx;
  state.f[1] += normal_force * nx;
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

void gpu_step(GpuParticleData &particles, GpuNeighborList &neighbors,
              const NeighborBuildParams &neighbor_params,
              const HookeNormalParams &contact_params, GpuDemContext &context, double dt)
{
  const double dtf = 0.5 * dt;
  reset_forces(particles, context.stream());
  build_contact_pairs_uniform_grid(particles, neighbor_params, neighbors, context.stream());
  compute_hooke_normal_forces(particles, neighbors, contact_params, context.stream());
  nve_sphere_initial_integrate(particles, dt, dtf, 1, 3, 1.0, context.stream());

  reset_forces(particles, context.stream());
  build_contact_pairs_uniform_grid(particles, neighbor_params, neighbors, context.stream());
  compute_hooke_normal_forces(particles, neighbors, contact_params, context.stream());
  nve_sphere_final_integrate(particles, dtf, 1, 3, 1.0, context.stream());
}

} // namespace

int main()
{
  GpuDemContext context(0);

  HookeNormalParams contact_params;
  contact_params.normal_stiffness = 1000.0;
  contact_params.normal_damping = 3.0;

  const double dt = 1.0e-4;
  const int steps = 4000;

  CpuState cpu;
  cpu.x[0] = 1.0;
  cpu.x[1] = 2.2;
  cpu.v[0] = 1.0;
  cpu.v[1] = -1.0;
  const double initial_energy = kinetic_energy(cpu);

  HostParticleData host;
  host.resize(2);
  host.position_x = {cpu.x[0], cpu.x[1]};
  host.position_y = {1.0, 1.0};
  host.position_z = {1.0, 1.0};
  host.velocity_x = {cpu.v[0], cpu.v[1]};
  host.velocity_y = {0.0, 0.0};
  host.velocity_z = {0.0, 0.0};
  host.radius = {cpu.radius[0], cpu.radius[1]};

  for (std::size_t i = 0; i < host.size(); ++i) {
    host.omega_x[i] = host.omega_y[i] = host.omega_z[i] = 0.0;
    host.force_x[i] = host.force_y[i] = host.force_z[i] = 0.0;
    host.torque_x[i] = host.torque_y[i] = host.torque_z[i] = 0.0;
    host.mass[i] = cpu.mass[i];
    host.density[i] = 1.0;
    host.type[i] = 1;
    host.mask[i] = 1;
    host.tag[i] = static_cast<long long>(i + 1);
    host.image_flags[i] = 0;
  }

  NeighborBuildParams neighbor_params;
  neighbor_params.box_min[0] = neighbor_params.box_min[1] = neighbor_params.box_min[2] = 0.0;
  neighbor_params.box_max[0] = neighbor_params.box_max[1] = neighbor_params.box_max[2] = 4.0;
  neighbor_params.cell_size = 1.1;
  neighbor_params.skin = 0.0;

  GpuParticleData particles;
  particles.copy_from_host(host, context.stream());

  GpuNeighborList neighbors;
  neighbors.allocate(host.size(), 2, num_cells(neighbor_params));

  for (int step = 0; step < steps; ++step) {
    cpu_step(cpu, contact_params, dt);
    gpu_step(particles, neighbors, neighbor_params, contact_params, context, dt);
  }

  HostParticleData result;
  particles.copy_to_host(result, context.stream());
  context.synchronize();

  expect_near(result.position_x[0], cpu.x[0], 2.0e-12, "position_x[0]");
  expect_near(result.position_x[1], cpu.x[1], 2.0e-12, "position_x[1]");
  expect_near(result.velocity_x[0], cpu.v[0], 2.0e-12, "velocity_x[0]");
  expect_near(result.velocity_x[1], cpu.v[1], 2.0e-12, "velocity_x[1]");

  const double final_energy = kinetic_energy(cpu);
  if (!(final_energy < initial_energy)) {
    std::fprintf(stderr, "damped collision did not dissipate kinetic energy\n");
    return 1;
  }
  if (!(std::fabs(result.velocity_x[0]) < 1.0 && std::fabs(result.velocity_x[1]) < 1.0)) {
    std::fprintf(stderr, "damped collision did not reduce rebound speed\n");
    return 1;
  }

  std::printf("GPU_DEM damped Hooke collision test passed on device %d (%s)\n",
              context.device_id(), context.device_properties().name);
  return 0;
}
