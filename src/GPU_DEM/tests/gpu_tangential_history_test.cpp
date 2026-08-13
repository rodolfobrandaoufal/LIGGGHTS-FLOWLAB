/* ----------------------------------------------------------------------
   Standalone GPU_DEM Hertz tangential-history regression test.
------------------------------------------------------------------------- */

#include "gpu_contact_pipeline.h"
#include "gpu_contact_history.h"
#include "gpu_dem_context.h"
#include "gpu_integrator.h"
#include "gpu_neighbor_builder.h"
#include "gpu_particle_data.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using LAMMPS_NS::GPU_DEM::GpuContactHistory;
using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuNeighborList;
using LAMMPS_NS::GPU_DEM::GpuParticleData;
using LAMMPS_NS::GPU_DEM::HertzNormalParams;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::NeighborBuildParams;
using LAMMPS_NS::GPU_DEM::RollingFrictionParams;
using LAMMPS_NS::GPU_DEM::TangentialHistoryParams;
using LAMMPS_NS::GPU_DEM::build_contact_pairs_uniform_grid;
using LAMMPS_NS::GPU_DEM::compute_hertz_history_forces;
using LAMMPS_NS::GPU_DEM::compute_hertz_history_rolling_forces;
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

HostParticleData make_particles(double separation, double tangential_velocity)
{
  HostParticleData host;
  host.resize(2);
  host.position_x = {1.0, 1.0 + separation};
  host.position_y = {1.0, 1.0};
  host.position_z = {1.0, 1.0};
  host.velocity_x = {0.0, 0.0};
  host.velocity_y = {0.0, tangential_velocity};
  host.velocity_z = {0.0, 0.0};
  host.radius = {0.5, 0.5};

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

  NeighborBuildParams neighbor_params;
  neighbor_params.box_min[0] = neighbor_params.box_min[1] = neighbor_params.box_min[2] = 0.0;
  neighbor_params.box_max[0] = neighbor_params.box_max[1] = neighbor_params.box_max[2] = 4.0;
  neighbor_params.cell_size = 1.1;
  neighbor_params.skin = 0.0;

  HertzNormalParams normal_params;
  normal_params.effective_youngs_modulus = 1000.0;
  normal_params.beta_effective = 0.0;
  normal_params.limit_force = true;

  TangentialHistoryParams tangential_params;
  tangential_params.effective_shear_modulus = 500.0;
  tangential_params.friction_coefficient = 10.0;
  tangential_params.tangential_damping = false;
  tangential_params.dt = 1.0e-4;

  GpuParticleData particles;
  particles.copy_from_host(make_particles(0.9, 1.0), context.stream());

  GpuNeighborList neighbors;
  neighbors.allocate(2, 4, num_cells(neighbor_params));

  GpuContactHistory history;
  history.allocate(8);

  history.prepare_for_step(context.stream());
  reset_forces(particles, context.stream());
  build_contact_pairs_uniform_grid(particles, neighbor_params, neighbors, context.stream());
  compute_hertz_history_forces(particles, neighbors, history, normal_params,
                               tangential_params, context.stream());

  HostParticleData first_result;
  particles.copy_to_host(first_result, context.stream());
  context.synchronize();

  const double reff = 0.25;
  const double overlap = 0.1;
  const double kt = 8.0 * tangential_params.effective_shear_modulus *
                    std::sqrt(reff * overlap);
  const double expected_tangential_force = -kt * tangential_params.dt;
  expect_near(first_result.force_y[0], expected_tangential_force, 1.0e-12,
              "first tangential force i");
  expect_near(first_result.force_y[1], -expected_tangential_force, 1.0e-12,
              "first tangential force j");
  expect_near(first_result.torque_z[0], -0.5 * expected_tangential_force, 1.0e-12,
              "first torque i");
  expect_near(first_result.torque_z[1], -0.5 * expected_tangential_force, 1.0e-12,
              "first torque j");

  history.prepare_for_step(context.stream());
  reset_forces(particles, context.stream());
  build_contact_pairs_uniform_grid(particles, neighbor_params, neighbors, context.stream());
  compute_hertz_history_forces(particles, neighbors, history, normal_params,
                               tangential_params, context.stream());

  HostParticleData second_result;
  particles.copy_to_host(second_result, context.stream());
  context.synchronize();
  expect_near(second_result.force_y[0], 2.0 * expected_tangential_force, 1.0e-12,
              "second tangential force i");

  particles.copy_from_host(make_particles(1.2, 0.0), context.stream());
  history.prepare_for_step(context.stream());
  reset_forces(particles, context.stream());
  build_contact_pairs_uniform_grid(particles, neighbor_params, neighbors, context.stream());

  particles.copy_from_host(make_particles(0.9, 0.0), context.stream());
  HostParticleData rolling_particles = make_particles(0.9, 0.0);
  rolling_particles.omega_y[0] = 1.0;
  particles.copy_from_host(rolling_particles, context.stream());
  history.prepare_for_step(context.stream());
  reset_forces(particles, context.stream());
  build_contact_pairs_uniform_grid(particles, neighbor_params, neighbors, context.stream());
  RollingFrictionParams rolling_params;
  rolling_params.coefficient = 0.2;
  rolling_params.torsion_torque = false;
  TangentialHistoryParams rolling_tangential_params = tangential_params;
  rolling_tangential_params.friction_coefficient = 0.0;
  compute_hertz_history_rolling_forces(particles, neighbors, history, normal_params,
                                       rolling_tangential_params, rolling_params,
                                       context.stream());

  HostParticleData reset_result;
  particles.copy_to_host(reset_result, context.stream());
  context.synchronize();
  expect_near(reset_result.force_y[0], 0.0, 1.0e-12,
              "separated contact reset tangential force");
  const double normal_force = (4.0 / 3.0) * normal_params.effective_youngs_modulus *
                              std::sqrt(0.25 * 0.1) * 0.1;
  const double expected_rolling_torque = -rolling_params.coefficient *
                                         normal_force * 0.25;
  expect_near(reset_result.torque_y[0], expected_rolling_torque, 1.0e-12,
              "cdt rolling torque i");
  expect_near(reset_result.torque_y[1], -expected_rolling_torque, 1.0e-12,
              "cdt rolling torque j");

  std::printf("GPU_DEM tangential history test passed on device %d (%s)\n",
              context.device_id(), context.device_properties().name);
  return 0;
}
