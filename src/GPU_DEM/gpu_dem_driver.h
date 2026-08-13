/* ----------------------------------------------------------------------
   Persistent GPU_DEM execution driver.

   This class owns long-lived CUDA state for the experimental GPU path. It
   intentionally remains outside the production LIGGGHTS timestep loop until
   feature gating and validation are in place.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_DRIVER_H
#define LMP_GPU_DEM_DRIVER_H

#include "gpu_contact_pipeline.h"
#include "gpu_dem_context.h"
#include "gpu_mesh_wall.h"
#include "gpu_neighbor_builder.h"
#include "gpu_particle_data.h"
#include "gpu_timestep.h"
#include "gpu_wall_contact.h"

#include <cstddef>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct GpuDemDriverConfig {
  int device_id{0};
  PrecisionMode precision{PrecisionMode::Mixed};
  std::size_t pair_capacity{0};
  int num_cells{1};
};

class GpuDemDriver {
 public:
  explicit GpuDemDriver(const GpuDemDriverConfig &config = GpuDemDriverConfig());

  GpuDemDriver(const GpuDemDriver &) = delete;
  GpuDemDriver &operator=(const GpuDemDriver &) = delete;

  GpuDemContext &context() { return context_; }
  const GpuDemContext &context() const { return context_; }

  GpuParticleData &particles() { return particles_; }
  const GpuParticleData &particles() const { return particles_; }

  GpuNeighborList &neighbors() { return neighbors_; }
  const GpuNeighborList &neighbors() const { return neighbors_; }

  std::size_t particle_count() const { return particles_.size(); }
  std::size_t pair_capacity() const { return pair_capacity_; }
  int num_cells() const { return num_cells_; }

  void upload_particles(const HostParticleData &host);
  void download_particles(HostParticleData &host) const;
  void set_plane_walls(const HostPlaneWalls &host_walls);
  void set_triangle_mesh(const HostTriangleMesh &host_mesh);

  void step_hooke(const NeighborBuildParams &neighbor_params,
                  const HookeNormalParams &contact_params,
                  const GpuTimestepParams &timestep_params);

  void step_hooke(const NeighborBuildParams &neighbor_params,
                  const HookeNormalParams &contact_params,
                  const GpuTimestepParams &timestep_params,
                  GpuTimestepTiming &timing);

  void step_hertz(const NeighborBuildParams &neighbor_params,
                  const HertzNormalParams &contact_params,
                  const GpuTimestepParams &timestep_params);

  void step_hertz(const NeighborBuildParams &neighbor_params,
                  const HertzNormalParams &contact_params,
                  const GpuTimestepParams &timestep_params,
                  GpuTimestepTiming &timing);

  void step_hertz_history_mesh(const NeighborBuildParams &neighbor_params,
                               const HertzNormalParams &normal_params,
                               const TangentialHistoryParams &tangential_params,
                               const RollingFrictionParams &rolling_params,
                               const GpuTimestepParams &timestep_params);

  void step_hertz_history_mesh(const NeighborBuildParams &neighbor_params,
                               const HertzNormalParams &normal_params,
                               const TangentialHistoryParams &tangential_params,
                               const RollingFrictionParams &rolling_params,
                               const GpuTimestepParams &timestep_params,
                               GpuTimestepTiming &timing);

  void synchronize() const { context_.synchronize(); }

 private:
  void ensure_neighbor_capacity(const NeighborBuildParams &neighbor_params);

  GpuDemContext context_;
  GpuParticleData particles_;
  GpuNeighborList neighbors_;
  GpuPlaneWalls walls_;
  GpuTriangleMesh mesh_;
  GpuContactHistory history_;
  bool has_walls_{false};
  bool has_mesh_{false};
  std::size_t pair_capacity_{0};
  int num_cells_{1};
};

int gpu_dem_num_cells_for_params(const NeighborBuildParams &params);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
