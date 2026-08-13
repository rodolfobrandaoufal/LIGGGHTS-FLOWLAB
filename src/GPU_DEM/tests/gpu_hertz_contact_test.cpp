/* ----------------------------------------------------------------------
   Standalone GPU_DEM Hertz normal contact-force regression test.
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
using LAMMPS_NS::GPU_DEM::HertzNormalParams;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::NeighborBuildParams;
using LAMMPS_NS::GPU_DEM::build_contact_pairs_uniform_grid;
using LAMMPS_NS::GPU_DEM::compute_hertz_normal_forces;
using LAMMPS_NS::GPU_DEM::reset_forces;

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

double hertz_elastic_force(double overlap, double radi, double radj, double effective_young)
{
  const double reff = radi * radj / (radi + radj);
  const double sqrt_delta_reff = std::sqrt(reff * overlap);
  const double kn = (4.0 / 3.0) * effective_young * sqrt_delta_reff;
  return kn * overlap;
}

} // namespace

int main()
{
  GpuDemContext context(0);

  HostParticleData host;
  host.resize(3);

  host.position_x = {1.0, 1.75, 4.0};
  host.position_y = {1.0, 1.00, 1.0};
  host.position_z = {1.0, 1.00, 1.0};
  host.radius = {0.5, 0.5, 0.5};

  for (std::size_t i = 0; i < host.size(); ++i) {
    host.velocity_x[i] = host.velocity_y[i] = host.velocity_z[i] = 0.0;
    host.omega_x[i] = host.omega_y[i] = host.omega_z[i] = 0.0;
    host.force_x[i] = 100.0 + static_cast<double>(i);
    host.force_y[i] = -50.0;
    host.force_z[i] = 25.0;
    host.torque_x[i] = host.torque_y[i] = host.torque_z[i] = 0.0;
    host.mass[i] = 1.0;
    host.density[i] = 1.0;
    host.type[i] = 1;
    host.mask[i] = 1;
    host.tag[i] = static_cast<long long>(i + 1);
    host.image_flags[i] = 0;
  }

  NeighborBuildParams neighbor_params;
  neighbor_params.box_min[0] = neighbor_params.box_min[1] = neighbor_params.box_min[2] = 0.0;
  neighbor_params.box_max[0] = neighbor_params.box_max[1] = neighbor_params.box_max[2] = 8.0;
  neighbor_params.cell_size = 1.25;
  neighbor_params.skin = 0.0;

  GpuParticleData particles;
  particles.copy_from_host(host, context.stream());
  reset_forces(particles, context.stream());

  GpuNeighborList neighbors;
  neighbors.allocate(host.size(), 8, num_cells(neighbor_params));
  build_contact_pairs_uniform_grid(particles, neighbor_params, neighbors, context.stream());

  HertzNormalParams contact_params;
  contact_params.effective_youngs_modulus = 120.0;
  contact_params.beta_effective = 0.0;
  contact_params.limit_force = true;
  compute_hertz_normal_forces(particles, neighbors, contact_params, context.stream());

  HostParticleData result;
  particles.copy_to_host(result, context.stream());
  context.synchronize();

  const double expected_force = hertz_elastic_force(0.25, 0.5, 0.5,
                                                    contact_params.effective_youngs_modulus);
  expect_near(expected_force, 10.0, 1.0e-12, "expected Hertz force setup");
  expect_near(result.force_x[0], -expected_force, 1.0e-12, "force_x[0]");
  expect_near(result.force_x[1], expected_force, 1.0e-12, "force_x[1]");
  expect_near(result.force_x[2], 0.0, 1.0e-12, "force_x[2]");
  expect_near(result.force_y[0], 0.0, 1.0e-12, "force_y[0]");
  expect_near(result.force_y[1], 0.0, 1.0e-12, "force_y[1]");
  expect_near(result.force_y[2], 0.0, 1.0e-12, "force_y[2]");
  expect_near(result.force_z[0], 0.0, 1.0e-12, "force_z[0]");
  expect_near(result.force_z[1], 0.0, 1.0e-12, "force_z[1]");
  expect_near(result.force_z[2], 0.0, 1.0e-12, "force_z[2]");
  expect_near(result.force_x[0] + result.force_x[1] + result.force_x[2], 0.0,
              1.0e-12, "net force x");

  std::printf("GPU_DEM Hertz normal contact-force test passed on device %d (%s)\n",
              context.device_id(), context.device_properties().name);
  return 0;
}
