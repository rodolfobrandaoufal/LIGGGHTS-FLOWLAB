/* ----------------------------------------------------------------------
   Standalone GPU_DEM neighbor-builder regression test.
------------------------------------------------------------------------- */

#include "gpu_dem_context.h"
#include "gpu_neighbor_builder.h"
#include "gpu_particle_data.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <utility>

using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuNeighborList;
using LAMMPS_NS::GPU_DEM::GpuParticleData;
using LAMMPS_NS::GPU_DEM::HostNeighborPairs;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::NeighborBuildParams;
using LAMMPS_NS::GPU_DEM::build_contact_pairs_uniform_grid;

namespace {

double minimum_image(double delta, double length, int periodic)
{
  if (!periodic) return delta;
  if (delta > 0.5 * length) return delta - length;
  if (delta < -0.5 * length) return delta + length;
  return delta;
}

std::set<std::pair<int, int>> build_reference_pairs(const HostParticleData &host,
                                                    const NeighborBuildParams &params)
{
  std::set<std::pair<int, int>> pairs;
  const int n = static_cast<int>(host.size());
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dx = host.position_x[j] - host.position_x[i];
      double dy = host.position_y[j] - host.position_y[i];
      double dz = host.position_z[j] - host.position_z[i];
      dx = minimum_image(dx, params.box_max[0] - params.box_min[0], params.periodic[0]);
      dy = minimum_image(dy, params.box_max[1] - params.box_min[1], params.periodic[1]);
      dz = minimum_image(dz, params.box_max[2] - params.box_min[2], params.periodic[2]);

      const double cutoff = host.radius[i] + host.radius[j] + params.skin;
      const double rsq = dx * dx + dy * dy + dz * dz;
      if (rsq <= cutoff * cutoff) pairs.insert(std::make_pair(i, j));
    }
  }
  return pairs;
}

std::set<std::pair<int, int>> normalize_gpu_pairs(const HostNeighborPairs &host)
{
  std::set<std::pair<int, int>> pairs;
  for (std::size_t k = 0; k < host.size(); ++k) {
    const int i = std::min(host.first[k], host.second[k]);
    const int j = std::max(host.first[k], host.second[k]);
    pairs.insert(std::make_pair(i, j));
  }
  return pairs;
}

void require_equal(const std::set<std::pair<int, int>> &actual,
                   const std::set<std::pair<int, int>> &expected)
{
  if (actual == expected) return;

  std::fprintf(stderr, "neighbor pair mismatch\nexpected:");
  for (const auto &pair : expected) std::fprintf(stderr, " (%d,%d)", pair.first, pair.second);
  std::fprintf(stderr, "\nactual:");
  for (const auto &pair : actual) std::fprintf(stderr, " (%d,%d)", pair.first, pair.second);
  std::fprintf(stderr, "\n");
  std::exit(1);
}

int num_cells(const NeighborBuildParams &params)
{
  int cells = 1;
  for (int d = 0; d < 3; ++d) {
    const double length = params.box_max[d] - params.box_min[d];
    cells *= std::max(1, static_cast<int>(std::floor(length / params.cell_size)));
  }
  return cells;
}

} // namespace

int main()
{
  GpuDemContext context(0);

  HostParticleData host;
  host.resize(6);

  host.position_x = {1.0, 1.75, 4.0, 9.75, 0.20, 5.0};
  host.position_y = {1.0, 1.00, 4.0, 5.00, 5.00, 8.0};
  host.position_z = {1.0, 1.00, 4.0, 5.00, 5.00, 8.0};
  host.radius = {0.5, 0.5, 0.4, 0.35, 0.35, 0.25};

  for (std::size_t i = 0; i < host.size(); ++i) {
    host.velocity_x[i] = host.velocity_y[i] = host.velocity_z[i] = 0.0;
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

  NeighborBuildParams params;
  params.box_min[0] = params.box_min[1] = params.box_min[2] = 0.0;
  params.box_max[0] = params.box_max[1] = params.box_max[2] = 10.0;
  params.cell_size = 1.25;
  params.skin = 0.05;
  params.periodic[0] = 1;
  params.periodic[1] = 0;
  params.periodic[2] = 0;

  GpuParticleData device_particles;
  device_particles.copy_from_host(host, context.stream());

  GpuNeighborList neighbors;
  neighbors.allocate(host.size(), 16, num_cells(params));
  build_contact_pairs_uniform_grid(device_particles, params, neighbors, context.stream());

  HostNeighborPairs gpu_pairs;
  neighbors.copy_pairs_to_host(gpu_pairs, context.stream());
  if (gpu_pairs.overflow != 0) {
    std::fprintf(stderr, "unexpected neighbor pair overflow\n");
    return 1;
  }

  const auto expected = build_reference_pairs(host, params);
  const auto actual = normalize_gpu_pairs(gpu_pairs);
  require_equal(actual, expected);

  std::printf("GPU_DEM uniform-grid neighbor builder test passed on device %d (%s), pairs=%zu\n",
              context.device_id(), context.device_properties().name, actual.size());
  return 0;
}
