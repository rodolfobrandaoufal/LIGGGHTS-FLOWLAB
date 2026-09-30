/* ----------------------------------------------------------------------
   GPU analytic wall contact scaffold for GPU_DEM.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_WALL_CONTACT_H
#define LMP_GPU_DEM_WALL_CONTACT_H

#include "gpu_contact_pipeline.h"
#include "gpu_particle_data.h"

#include <cuda_runtime.h>

#include <cstddef>
#include <vector>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct PlaneWall {
  double normal[3]{1.0, 0.0, 0.0};
  double offset{0.0};
};

struct HostPlaneWalls {
  std::vector<PlaneWall> walls;
};

class GpuPlaneWalls {
 public:
  GpuPlaneWalls() = default;
  ~GpuPlaneWalls();

  GpuPlaneWalls(const GpuPlaneWalls &) = delete;
  GpuPlaneWalls &operator=(const GpuPlaneWalls &) = delete;

  GpuPlaneWalls(GpuPlaneWalls &&other) noexcept;
  GpuPlaneWalls &operator=(GpuPlaneWalls &&other) noexcept;

  void copy_from_host(const HostPlaneWalls &host, cudaStream_t stream = 0);
  void release() noexcept;

  std::size_t size() const { return size_; }
  const PlaneWall *device_walls() const { return walls_; }

 private:
  std::size_t size_{0};
  std::size_t capacity_{0};
  PlaneWall *walls_{nullptr};
};

void compute_hooke_plane_wall_forces(GpuParticleData &particles,
                                     const GpuPlaneWalls &walls,
                                     const HookeNormalParams &params,
                                     cudaStream_t stream = 0);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
