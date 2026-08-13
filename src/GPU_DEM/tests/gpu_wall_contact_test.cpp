/* ----------------------------------------------------------------------
   Standalone GPU_DEM planar wall contact regression test.
------------------------------------------------------------------------- */

#include "gpu_dem_context.h"
#include "gpu_integrator.h"
#include "gpu_particle_data.h"
#include "gpu_wall_contact.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuParticleData;
using LAMMPS_NS::GPU_DEM::GpuPlaneWalls;
using LAMMPS_NS::GPU_DEM::HookeNormalParams;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::HostPlaneWalls;
using LAMMPS_NS::GPU_DEM::PlaneWall;
using LAMMPS_NS::GPU_DEM::compute_hooke_plane_wall_forces;
using LAMMPS_NS::GPU_DEM::nve_sphere_final_integrate;
using LAMMPS_NS::GPU_DEM::nve_sphere_initial_integrate;
using LAMMPS_NS::GPU_DEM::reset_forces;

namespace {

struct CpuState {
  double x{0.9};
  double v{-1.0};
  double f{0.0};
  double radius{0.5};
  double mass{1.0};
};

void expect_near(double actual, double expected, double tol, const char *label)
{
  if (std::fabs(actual - expected) <= tol) return;
  std::fprintf(stderr, "%s mismatch: actual %.17g expected %.17g tolerance %.3g\n",
               label, actual, expected, tol);
  std::exit(1);
}

void compute_cpu_wall_force(CpuState &state, const HookeNormalParams &params)
{
  state.f = 0.0;
  const double distance = state.x;
  const double overlap = state.radius - distance;
  if (overlap <= 0.0) return;

  double normal_force = params.normal_stiffness * overlap - params.normal_damping * state.v;
  if (normal_force < 0.0) normal_force = 0.0;
  state.f = normal_force;
}

void cpu_step(CpuState &state, const HookeNormalParams &params, double dt)
{
  const double dtf = 0.5 * dt;
  compute_cpu_wall_force(state, params);
  state.v += dtf * state.f / state.mass;
  state.x += dt * state.v;
  compute_cpu_wall_force(state, params);
  state.v += dtf * state.f / state.mass;
}

void gpu_step(GpuParticleData &particles, const GpuPlaneWalls &walls,
              const HookeNormalParams &params, GpuDemContext &context, double dt)
{
  const double dtf = 0.5 * dt;
  reset_forces(particles, context.stream());
  compute_hooke_plane_wall_forces(particles, walls, params, context.stream());
  nve_sphere_initial_integrate(particles, dt, dtf, 1, 3, 1.0, context.stream());

  reset_forces(particles, context.stream());
  compute_hooke_plane_wall_forces(particles, walls, params, context.stream());
  nve_sphere_final_integrate(particles, dtf, 1, 3, 1.0, context.stream());
}

HostParticleData make_host_particle(double x, double vx)
{
  HostParticleData host;
  host.resize(1);
  host.position_x[0] = x;
  host.position_y[0] = 1.0;
  host.position_z[0] = 1.0;
  host.velocity_x[0] = vx;
  host.velocity_y[0] = 0.0;
  host.velocity_z[0] = 0.0;
  host.omega_x[0] = host.omega_y[0] = host.omega_z[0] = 0.0;
  host.force_x[0] = host.force_y[0] = host.force_z[0] = 0.0;
  host.torque_x[0] = host.torque_y[0] = host.torque_z[0] = 0.0;
  host.radius[0] = 0.5;
  host.mass[0] = 1.0;
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

  HookeNormalParams params;
  params.normal_stiffness = 1000.0;
  params.normal_damping = 3.0;

  HostPlaneWalls host_walls;
  PlaneWall lower_x;
  lower_x.normal[0] = 1.0;
  lower_x.normal[1] = 0.0;
  lower_x.normal[2] = 0.0;
  lower_x.offset = 0.0;
  host_walls.walls.push_back(lower_x);

  GpuPlaneWalls walls;
  walls.copy_from_host(host_walls, context.stream());

  GpuParticleData static_particle;
  static_particle.copy_from_host(make_host_particle(0.4, 0.0), context.stream());
  reset_forces(static_particle, context.stream());
  compute_hooke_plane_wall_forces(static_particle, walls, params, context.stream());

  HostParticleData static_result;
  static_particle.copy_to_host(static_result, context.stream());
  context.synchronize();
  expect_near(static_result.force_x[0], 100.0, 1.0e-12, "static wall force x");
  expect_near(static_result.force_y[0], 0.0, 1.0e-12, "static wall force y");
  expect_near(static_result.force_z[0], 0.0, 1.0e-12, "static wall force z");

  CpuState cpu;
  GpuParticleData moving_particle;
  moving_particle.copy_from_host(make_host_particle(cpu.x, cpu.v), context.stream());

  const double dt = 1.0e-4;
  const int steps = 8000;
  for (int step = 0; step < steps; ++step) {
    cpu_step(cpu, params, dt);
    gpu_step(moving_particle, walls, params, context, dt);
  }

  HostParticleData result;
  moving_particle.copy_to_host(result, context.stream());
  context.synchronize();

  expect_near(result.position_x[0], cpu.x, 2.0e-12, "wall impact position");
  expect_near(result.velocity_x[0], cpu.v, 2.0e-12, "wall impact velocity");
  if (!(result.velocity_x[0] > 0.0 && result.velocity_x[0] < 1.0)) {
    std::fprintf(stderr, "damped wall impact did not rebound with reduced speed\n");
    return 1;
  }

  std::printf("GPU_DEM planar wall contact test passed on device %d (%s)\n",
              context.device_id(), context.device_properties().name);
  return 0;
}
