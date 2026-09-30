/* ----------------------------------------------------------------------
   GPU analytic wall contact scaffold for GPU_DEM.
------------------------------------------------------------------------- */

#include "gpu_wall_contact.h"

#include "gpu_error_check.h"
#include "models/normal_hooke.cuh"

#include <stdexcept>
#include <utility>

namespace LAMMPS_NS {
namespace GPU_DEM {

namespace {

constexpr int THREADS_PER_BLOCK = 256;

__global__ void hooke_plane_wall_force_kernel(
    std::size_t n, const double *x, const double *y, const double *z, const double *vx,
    const double *vy, const double *vz, const double *radius, const PlaneWall *walls,
    std::size_t wall_count, HookeNormalParams params, double *fx, double *fy, double *fz)
{
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  double fxi = 0.0;
  double fyi = 0.0;
  double fzi = 0.0;

  for (std::size_t wall_index = 0; wall_index < wall_count; ++wall_index) {
    const PlaneWall wall = walls[wall_index];
    const double nx = wall.normal[0];
    const double ny = wall.normal[1];
    const double nz = wall.normal[2];
    const double distance = x[i] * nx + y[i] * ny + z[i] * nz - wall.offset;
    const double overlap = radius[i] - distance;
    if (overlap <= 0.0) continue;

    const double normal_relative_velocity = vx[i] * nx + vy[i] * ny + vz[i] * nz;
    const double normal_force =
        Models::hooke_normal_force(overlap, normal_relative_velocity, params);

    fxi += normal_force * nx;
    fyi += normal_force * ny;
    fzi += normal_force * nz;
  }

  fx[i] += fxi;
  fy[i] += fyi;
  fz[i] += fzi;
}

} // namespace

GpuPlaneWalls::~GpuPlaneWalls()
{
  release();
}

GpuPlaneWalls::GpuPlaneWalls(GpuPlaneWalls &&other) noexcept
{
  *this = std::move(other);
}

GpuPlaneWalls &GpuPlaneWalls::operator=(GpuPlaneWalls &&other) noexcept
{
  if (this == &other) return *this;
  release();

  size_ = other.size_;
  capacity_ = other.capacity_;
  walls_ = other.walls_;

  other.size_ = 0;
  other.capacity_ = 0;
  other.walls_ = nullptr;
  return *this;
}

void GpuPlaneWalls::copy_from_host(const HostPlaneWalls &host, cudaStream_t stream)
{
  if (host.walls.size() > capacity_) {
    release();
    GPU_DEM_CUDA_CHECK(
        cudaMalloc(reinterpret_cast<void **>(&walls_), host.walls.size() * sizeof(PlaneWall)));
    capacity_ = host.walls.size();
  }

  size_ = host.walls.size();
  if (size_ > 0) {
    GPU_DEM_CUDA_CHECK(cudaMemcpyAsync(walls_, host.walls.data(), size_ * sizeof(PlaneWall),
                                       cudaMemcpyHostToDevice, stream));
  }
}

void GpuPlaneWalls::release() noexcept
{
  if (walls_) {
    cudaFree(walls_);
    walls_ = nullptr;
  }
  size_ = 0;
  capacity_ = 0;
}

void compute_hooke_plane_wall_forces(GpuParticleData &particles, const GpuPlaneWalls &walls,
                                     const HookeNormalParams &params, cudaStream_t stream)
{
  if (!(params.normal_stiffness >= 0.0) || !(params.normal_damping >= 0.0))
    throw std::runtime_error("GPU_DEM Hooke wall parameters must be nonnegative");
  if (walls.size() == 0 || particles.size() == 0) return;

  const int blocks =
      static_cast<int>((particles.size() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  hooke_plane_wall_force_kernel<<<blocks, THREADS_PER_BLOCK, 0, stream>>>(
      particles.size(), particles.position_x(), particles.position_y(), particles.position_z(),
      particles.velocity_x(), particles.velocity_y(), particles.velocity_z(),
      particles.radius(), walls.device_walls(), walls.size(), params, particles.force_x(),
      particles.force_y(), particles.force_z());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
