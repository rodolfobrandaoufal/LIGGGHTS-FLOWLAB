/* ----------------------------------------------------------------------
   Standalone GPU_DEM synthetic throughput benchmark.
------------------------------------------------------------------------- */

#include "gpu_dem_context.h"
#include "gpu_particle_data.h"
#include "gpu_timestep.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuNeighborList;
using LAMMPS_NS::GPU_DEM::GpuParticleData;
using LAMMPS_NS::GPU_DEM::GpuTimestepParams;
using LAMMPS_NS::GPU_DEM::GpuTimestepTiming;
using LAMMPS_NS::GPU_DEM::HookeNormalParams;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::NeighborBuildParams;
using LAMMPS_NS::GPU_DEM::velocity_verlet_hooke_step;

namespace {

int parse_int_arg(int argc, char **argv, const char *name, int default_value)
{
  for (int i = 1; i + 1 < argc; ++i) {
    if (std::string(argv[i]) == name) return std::atoi(argv[i + 1]);
  }
  return default_value;
}

double parse_double_arg(int argc, char **argv, const char *name, double default_value)
{
  for (int i = 1; i + 1 < argc; ++i) {
    if (std::string(argv[i]) == name) return std::atof(argv[i + 1]);
  }
  return default_value;
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

std::size_t conservative_pair_capacity(std::size_t n)
{
  return std::max<std::size_t>(1024, n * 32);
}

HostParticleData make_lattice(int nx, int ny, int nz, double spacing, double radius)
{
  HostParticleData host;
  const std::size_t n = static_cast<std::size_t>(nx) * ny * nz;
  host.resize(n);

  std::size_t index = 0;
  for (int iz = 0; iz < nz; ++iz) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int ix = 0; ix < nx; ++ix) {
        host.position_x[index] = 1.0 + ix * spacing;
        host.position_y[index] = 1.0 + iy * spacing;
        host.position_z[index] = 1.0 + iz * spacing;
        host.velocity_x[index] = 0.01 * ((ix % 3) - 1);
        host.velocity_y[index] = 0.01 * ((iy % 3) - 1);
        host.velocity_z[index] = 0.01 * ((iz % 3) - 1);
        host.omega_x[index] = host.omega_y[index] = host.omega_z[index] = 0.0;
        host.force_x[index] = host.force_y[index] = host.force_z[index] = 0.0;
        host.torque_x[index] = host.torque_y[index] = host.torque_z[index] = 0.0;
        host.radius[index] = radius;
        host.mass[index] = 1.0;
        host.density[index] = 1.0;
        host.type[index] = 1;
        host.mask[index] = 1;
        host.tag[index] = static_cast<long long>(index + 1);
        host.image_flags[index] = 0;
        ++index;
      }
    }
  }
  return host;
}

} // namespace

int main(int argc, char **argv)
{
  const int nx = parse_int_arg(argc, argv, "--nx", 16);
  const int ny = parse_int_arg(argc, argv, "--ny", 16);
  const int nz = parse_int_arg(argc, argv, "--nz", 8);
  const int steps = parse_int_arg(argc, argv, "--steps", 20);
  const double spacing = parse_double_arg(argc, argv, "--spacing", 0.95);
  const double radius = parse_double_arg(argc, argv, "--radius", 0.5);
  const double gravity_z = parse_double_arg(argc, argv, "--gravity-z", 0.0);

  if (nx <= 0 || ny <= 0 || nz <= 0 || steps <= 0 || spacing <= 0.0 || radius <= 0.0) {
    std::fprintf(stderr, "invalid benchmark arguments\n");
    return 2;
  }

  GpuDemContext context(0);
  HostParticleData host = make_lattice(nx, ny, nz, spacing, radius);

  GpuParticleData particles;
  particles.copy_from_host(host, context.stream());

  NeighborBuildParams neighbor_params;
  neighbor_params.box_min[0] = neighbor_params.box_min[1] = neighbor_params.box_min[2] = 0.0;
  neighbor_params.box_max[0] = 2.0 + (nx - 1) * spacing;
  neighbor_params.box_max[1] = 2.0 + (ny - 1) * spacing;
  neighbor_params.box_max[2] = 2.0 + (nz - 1) * spacing;
  neighbor_params.cell_size = 2.0 * radius;
  neighbor_params.skin = 0.02;

  GpuNeighborList neighbors;
  neighbors.allocate(host.size(), conservative_pair_capacity(host.size()),
                     num_cells(neighbor_params));

  HookeNormalParams contact_params;
  contact_params.normal_stiffness = 1000.0;
  contact_params.normal_damping = 1.0;

  GpuTimestepParams timestep_params;
  timestep_params.dt = 1.0e-5;
  timestep_params.groupbit = 1;
  timestep_params.dimension = 3;
  timestep_params.one_plus_added_mass = 1.0;
  timestep_params.gravity.enabled = gravity_z != 0.0;
  timestep_params.gravity.acceleration[2] = gravity_z;

  GpuTimestepTiming timing;
  GpuTimestepTiming accumulated;
  for (int step = 0; step < steps; ++step) {
    velocity_verlet_hooke_step(particles, neighbors, neighbor_params, contact_params,
                               nullptr, timestep_params, timing, context.stream());
    accumulated.reset_ms += timing.reset_ms;
    accumulated.external_force_ms += timing.external_force_ms;
    accumulated.neighbor_build_ms += timing.neighbor_build_ms;
    accumulated.particle_contact_ms += timing.particle_contact_ms;
    accumulated.initial_integrate_ms += timing.initial_integrate_ms;
    accumulated.final_integrate_ms += timing.final_integrate_ms;
    accumulated.total_ms += timing.total_ms;
  }
  context.synchronize();

  const int contacts = neighbors.pair_count(context.stream());
  const double total_seconds = accumulated.total_ms / 1000.0;
  const double particle_rate = total_seconds > 0.0 ? host.size() * steps / total_seconds : 0.0;
  const double contact_rate = total_seconds > 0.0 ? contacts * steps / total_seconds : 0.0;

  std::printf(
      "GPU_DEM synthetic benchmark passed on device %d (%s): particles=%zu contacts=%d "
      "steps=%d total=%.6f ms reset=%.6f ms external=%.6f ms neighbor=%.6f ms "
      "contact=%.6f ms initial_integrate=%.6f ms final_integrate=%.6f ms "
      "particles_per_second=%.6e contacts_per_second=%.6e\n",
      context.device_id(), context.device_properties().name, host.size(), contacts, steps,
      accumulated.total_ms, accumulated.reset_ms, accumulated.external_force_ms,
      accumulated.neighbor_build_ms, accumulated.particle_contact_ms,
      accumulated.initial_integrate_ms, accumulated.final_integrate_ms, particle_rate,
      contact_rate);
  return 0;
}
