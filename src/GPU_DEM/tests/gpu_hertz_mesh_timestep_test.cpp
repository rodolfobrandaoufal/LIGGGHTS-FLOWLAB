/* ----------------------------------------------------------------------
   Standalone GPU_DEM Hertz mesh-wall timestep regression test.
------------------------------------------------------------------------- */

#include "gpu_dem_context.h"
#include "gpu_mesh_wall.h"
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
using LAMMPS_NS::GPU_DEM::GpuTriangleMesh;
using LAMMPS_NS::GPU_DEM::HertzNormalParams;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::HostTriangleMesh;
using LAMMPS_NS::GPU_DEM::NeighborBuildParams;
using LAMMPS_NS::GPU_DEM::TriangleWall;
using LAMMPS_NS::GPU_DEM::velocity_verlet_hertz_mesh_step;

namespace {

struct CpuState {
  double z{0.9};
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

HostTriangleMesh make_floor_mesh()
{
  HostTriangleMesh host;
  TriangleWall triangle;
  triangle.v0[0] = -2.0; triangle.v0[1] = -2.0; triangle.v0[2] = 0.0;
  triangle.v1[0] =  2.0; triangle.v1[1] = -2.0; triangle.v1[2] = 0.0;
  triangle.v2[0] =  0.0; triangle.v2[1] =  2.0; triangle.v2[2] = 0.0;
  triangle.normal[0] = 0.0;
  triangle.normal[1] = 0.0;
  triangle.normal[2] = 1.0;
  host.triangles.push_back(triangle);
  return host;
}

HostParticleData make_particle(const CpuState &state)
{
  HostParticleData host;
  host.resize(1);
  host.position_x[0] = 0.0;
  host.position_y[0] = 0.0;
  host.position_z[0] = state.z;
  host.velocity_x[0] = 0.0;
  host.velocity_y[0] = 0.0;
  host.velocity_z[0] = state.v;
  host.omega_x[0] = host.omega_y[0] = host.omega_z[0] = 0.0;
  host.force_x[0] = host.force_y[0] = host.force_z[0] = 0.0;
  host.torque_x[0] = host.torque_y[0] = host.torque_z[0] = 0.0;
  host.radius[0] = state.radius;
  host.mass[0] = state.mass;
  host.density[0] = 1.0;
  host.type[0] = 1;
  host.mask[0] = 1;
  host.tag[0] = 1;
  host.image_flags[0] = 0;
  return host;
}

void compute_cpu_force(CpuState &state, const HertzNormalParams &params)
{
  state.f = 0.0;
  const double overlap = state.radius - state.z;
  if (overlap <= 0.0) return;

  const double sqrt_delta_reff = std::sqrt(state.radius * overlap);
  const double kn = (4.0 / 3.0) * params.effective_youngs_modulus * sqrt_delta_reff;
  state.f = kn * overlap;
}

void cpu_step(CpuState &state, const HertzNormalParams &params, double dt)
{
  const double dtf = 0.5 * dt;
  compute_cpu_force(state, params);
  state.v += dtf * state.f / state.mass;
  state.z += dt * state.v;
  compute_cpu_force(state, params);
  state.v += dtf * state.f / state.mass;
}

} // namespace

int main()
{
  GpuDemContext context(0);

  HertzNormalParams contact_params;
  contact_params.effective_youngs_modulus = 1000.0;
  contact_params.beta_effective = 0.0;
  contact_params.limit_force = true;

  CpuState cpu;
  GpuParticleData particles;
  particles.copy_from_host(make_particle(cpu), context.stream());

  GpuTriangleMesh mesh;
  mesh.copy_from_host(make_floor_mesh(), context.stream());

  NeighborBuildParams neighbor_params;
  neighbor_params.box_min[0] = neighbor_params.box_min[1] = neighbor_params.box_min[2] = -3.0;
  neighbor_params.box_max[0] = neighbor_params.box_max[1] = neighbor_params.box_max[2] = 3.0;
  neighbor_params.cell_size = 1.1;
  neighbor_params.skin = 0.0;

  GpuNeighborList neighbors;
  neighbors.allocate(1, 1, num_cells(neighbor_params));

  GpuTimestepParams timestep_params;
  timestep_params.dt = 1.0e-5;
  timestep_params.groupbit = 1;
  timestep_params.dimension = 3;
  timestep_params.one_plus_added_mass = 1.0;

  const int steps = 80000;
  for (int step = 0; step < steps; ++step) {
    cpu_step(cpu, contact_params, timestep_params.dt);
    velocity_verlet_hertz_mesh_step(particles, neighbors, neighbor_params,
                                    contact_params, &mesh, timestep_params,
                                    context.stream());
  }

  HostParticleData result;
  particles.copy_to_host(result, context.stream());
  context.synchronize();

  expect_near(result.position_z[0], cpu.z, 2.0e-11, "mesh timestep position_z");
  expect_near(result.velocity_z[0], cpu.v, 2.0e-11, "mesh timestep velocity_z");
  if (!(result.velocity_z[0] > 0.9)) {
    std::fprintf(stderr, "Hertz mesh wall impact did not rebound as expected\n");
    return 1;
  }

  std::printf("GPU_DEM Hertz mesh-wall timestep test passed on device %d (%s)\n",
              context.device_id(), context.device_properties().name);
  return 0;
}
