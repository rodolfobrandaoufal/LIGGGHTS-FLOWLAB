/* ----------------------------------------------------------------------
   AoS bridge plus GPU timestep regression test.

   This test exercises the future LIGGGHTS integration boundary without
   attaching the GPU path to the production run loop.
------------------------------------------------------------------------- */

#include "gpu_dem_context.h"
#include "gpu_liggghts_bridge.h"
#include "gpu_particle_data.h"
#include "gpu_timestep.h"
#include "gpu_wall_contact.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using LAMMPS_NS::GPU_DEM::CpuParticleAoSView;
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
using LAMMPS_NS::GPU_DEM::gather_aos_to_host;
using LAMMPS_NS::GPU_DEM::scatter_dynamic_state_to_aos;
using LAMMPS_NS::GPU_DEM::velocity_verlet_hooke_step;

namespace {

struct OwnedAoS {
  explicit OwnedAoS(int n) : nlocal(n)
  {
    allocate_vector3(x_storage, x);
    allocate_vector3(v_storage, v);
    allocate_vector3(f_storage, f);
    allocate_vector3(omega_storage, omega);
    allocate_vector3(torque_storage, torque);
    radius.resize(nlocal, 0.5);
    rmass.resize(nlocal, 1.0);
    density.resize(nlocal, 1.0);
    type.resize(nlocal, 1);
    mask.resize(nlocal, 1);
    tag.resize(nlocal);
    image.resize(nlocal, 0);
    for (int i = 0; i < nlocal; ++i) tag[i] = i + 1;
  }

  CpuParticleAoSView view()
  {
    CpuParticleAoSView result;
    result.nlocal = static_cast<std::size_t>(nlocal);
    result.x = x.data();
    result.v = v.data();
    result.f = f.data();
    result.omega = omega.data();
    result.torque = torque.data();
    result.radius = radius.data();
    result.rmass = rmass.data();
    result.density = density.data();
    result.type = type.data();
    result.mask = mask.data();
    result.tag = tag.data();
    result.image = image.data();
    return result;
  }

  int nlocal;
  std::vector<double> x_storage, v_storage, f_storage;
  std::vector<double> omega_storage, torque_storage;
  std::vector<double *> x, v, f, omega, torque;
  std::vector<double> radius, rmass, density;
  std::vector<int> type, mask;
  std::vector<LAMMPS_NS::tagint> tag;
  std::vector<LAMMPS_NS::imageint> image;

 private:
  void allocate_vector3(std::vector<double> &storage,
                        std::vector<double *> &ptrs)
  {
    storage.resize(static_cast<std::size_t>(nlocal) * 3, 0.0);
    ptrs.resize(nlocal);
    for (int i = 0; i < nlocal; ++i)
      ptrs[i] = &storage[static_cast<std::size_t>(i) * 3];
  }
};

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

void add_pair_force(CpuState &state, const HookeNormalParams &params)
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

void add_wall_force(CpuState &state, const HookeNormalParams &params)
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
  add_pair_force(state, params);
  add_wall_force(state, params);
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

void fill_aos_from_cpu(const CpuState &state, OwnedAoS &aos)
{
  for (int i = 0; i < 2; ++i) {
    aos.x[i][0] = state.x[i];
    aos.x[i][1] = 1.0;
    aos.x[i][2] = 1.0;
    aos.v[i][0] = state.v[i];
    aos.v[i][1] = 0.0;
    aos.v[i][2] = 0.0;
    aos.radius[i] = state.radius[i];
    aos.rmass[i] = state.mass[i];
  }
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

  OwnedAoS aos(2);
  fill_aos_from_cpu(cpu, aos);

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

  HostParticleData host;
  gather_aos_to_host(aos.view(), host);

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
  scatter_dynamic_state_to_aos(result, aos.view());

  expect_near(aos.x[0][0], cpu.x[0], 2.0e-12, "bridge timestep x[0][0]");
  expect_near(aos.x[1][0], cpu.x[1], 2.0e-12, "bridge timestep x[1][0]");
  expect_near(aos.v[0][0], cpu.v[0], 2.0e-12, "bridge timestep v[0][0]");
  expect_near(aos.v[1][0], cpu.v[1], 2.0e-12, "bridge timestep v[1][0]");

  if (!(aos.x[0][0] > 0.5)) {
    std::fprintf(stderr, "bridge timestep wall contact allowed penetration\n");
    return 1;
  }

  std::printf("GPU_DEM AoS bridge timestep test passed on device %d (%s)\n",
              context.device_id(), context.device_properties().name);
  return 0;
}
