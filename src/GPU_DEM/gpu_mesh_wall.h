/* ----------------------------------------------------------------------
   GPU triangle-mesh wall contact scaffold for GPU_DEM.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_MESH_WALL_H
#define LMP_GPU_DEM_MESH_WALL_H

#include "gpu_contact_pipeline.h"
#include "gpu_particle_data.h"

#include <cuda_runtime.h>

#include <cstddef>
#include <string>
#include <vector>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct TriangleWall {
  double v0[3]{0.0, 0.0, 0.0};
  double v1[3]{0.0, 0.0, 0.0};
  double v2[3]{0.0, 0.0, 0.0};
  double normal[3]{0.0, 0.0, 1.0};
};

struct HostTriangleMesh {
  std::vector<TriangleWall> triangles;
};

HostTriangleMesh load_ascii_stl_triangle_mesh(const std::string &path, double scale);

class GpuTriangleMesh {
 public:
  GpuTriangleMesh() = default;
  ~GpuTriangleMesh();

  GpuTriangleMesh(const GpuTriangleMesh &) = delete;
  GpuTriangleMesh &operator=(const GpuTriangleMesh &) = delete;

  GpuTriangleMesh(GpuTriangleMesh &&other) noexcept;
  GpuTriangleMesh &operator=(GpuTriangleMesh &&other) noexcept;

  void copy_from_host(const HostTriangleMesh &host, cudaStream_t stream = 0);
  void release() noexcept;

  std::size_t size() const { return size_; }
  const TriangleWall *device_triangles() const { return triangles_; }

 private:
  std::size_t size_{0};
  std::size_t capacity_{0};
  TriangleWall *triangles_{nullptr};
};

void compute_hertz_triangle_wall_forces(GpuParticleData &particles,
                                        const GpuTriangleMesh &mesh,
                                        const HertzNormalParams &params,
                                        cudaStream_t stream = 0);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
