/* ----------------------------------------------------------------------
   Standalone GPU_DEM two-particle Hertz collision regression test.
------------------------------------------------------------------------- */

#include "gpu_dem_context.h"
#include "gpu_neighbor_builder.h"
#include "gpu_particle_data.h"
#include "gpu_timestep.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuNeighborList;
using LAMMPS_NS::GPU_DEM::GpuParticleData;
using LAMMPS_NS::GPU_DEM::GpuTimestepParams;
using LAMMPS_NS::GPU_DEM::HertzNormalParams;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::NeighborBuildParams;
using LAMMPS_NS::GPU_DEM::velocity_verlet_hertz_step;

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

double effective_radius(const CpuState &state)
{
  return state.radius[0] * state.radius[1] / (state.radius[0] + state.radius[1]);
}

double hertz_force(double overlap, const CpuState &state, const HertzNormalParams &params)
{
  const double sqrt_delta_reff = std::sqrt(effective_radius(state) * overlap);
  const double kn = (4.0 / 3.0) * params.effective_youngs_modulus * sqrt_delta_reff;
  return kn * overlap;
}

double kinetic_energy(const CpuState &state)
{
  return 0.5 * state.mass[0] * state.v[0] * state.v[0] +
         0.5 * state.mass[1] * state.v[1] * state.v[1];
}

double hertz_spring_energy(const CpuState &state, const HertzNormalParams &params)
{
  const double distance = std::fabs(state.x[1] - state.x[0]);
  const double overlap = state.radius[0] + state.radius[1] - distance;
  if (overlap <= 0.0) return 0.0;
  return (8.0 / 15.0) * params.effective_youngs_modulus *
         std::sqrt(effective_radius(state)) * std::pow(overlap, 2.5);
}

void compute_cpu_force(CpuState &state, const HertzNormalParams &params)
{
  state.f[0] = 0.0;
  state.f[1] = 0.0;

  const double dx = state.x[1] - state.x[0];
  const double distance = std::fabs(dx);
  if (distance <= 0.0) return;

  const double overlap = state.radius[0] + state.radius[1] - distance;
  if (overlap <= 0.0) return;

  const double nx = dx / distance;
  const double force = hertz_force(overlap, state, params);
  state.f[0] += -force * nx;
  state.f[1] += force * nx;
}

void cpu_velocity_verlet_step(CpuState &state, const HertzNormalParams &params, double dt)
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

  HertzNormalParams contact_params;
  contact_params.effective_youngs_modulus = 1000.0;
  contact_params.beta_effective = 0.0;
  contact_params.limit_force = true;

  const double dt = 1.0e-5;
  const int steps = 45000;

  CpuState cpu;
  cpu.x[0] = 1.0;
  cpu.x[1] = 2.2;
  cpu.v[0] = 1.0;
  cpu.v[1] = -1.0;
  const double initial_energy = kinetic_energy(cpu) +
                                hertz_spring_energy(cpu, contact_params);

  NeighborBuildParams neighbor_params;
  neighbor_params.box_min[0] = neighbor_params.box_min[1] = neighbor_params.box_min[2] = 0.0;
  neighbor_params.box_max[0] = neighbor_params.box_max[1] = neighbor_params.box_max[2] = 4.0;
  neighbor_params.cell_size = 1.1;
  neighbor_params.skin = 0.0;

  HostParticleData host = make_particles(cpu);
  GpuParticleData particles;
  particles.copy_from_host(host, context.stream());

  GpuNeighborList neighbors;
  neighbors.allocate(host.size(), 2, num_cells(neighbor_params));

  GpuTimestepParams timestep_params;
  timestep_params.dt = dt;
  timestep_params.groupbit = 1;
  timestep_params.dimension = 3;
  timestep_params.one_plus_added_mass = 1.0;

  for (int step = 0; step < steps; ++step) {
    cpu_velocity_verlet_step(cpu, contact_params, dt);
    velocity_verlet_hertz_step(particles, neighbors, neighbor_params, contact_params,
                               nullptr, timestep_params, context.stream());
  }

  HostParticleData result;
  particles.copy_to_host(result, context.stream());
  context.synchronize();

  expect_near(result.position_x[0], cpu.x[0], 2.0e-11, "position_x[0]");
  expect_near(result.position_x[1], cpu.x[1], 2.0e-11, "position_x[1]");
  expect_near(result.velocity_x[0], cpu.v[0], 2.0e-11, "velocity_x[0]");
  expect_near(result.velocity_x[1], cpu.v[1], 2.0e-11, "velocity_x[1]");

  const double final_energy = kinetic_energy(cpu) +
                              hertz_spring_energy(cpu, contact_params);
  expect_near(final_energy, initial_energy, 1.0e-5, "CPU reference total energy");
  if (!(result.velocity_x[0] < -0.99 && result.velocity_x[1] > 0.99)) {
    std::fprintf(stderr, "Hertz collision did not reverse velocities as expected\n");
    return 1;
  }

  std::printf("GPU_DEM Hertz two-particle collision test passed on device %d (%s)\n",
              context.device_id(), context.device_properties().name);
  return 0;
}
