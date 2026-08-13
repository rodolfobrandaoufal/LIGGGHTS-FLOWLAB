/* ----------------------------------------------------------------------
   Device-side particle insertion helpers for GPU_DEM.
------------------------------------------------------------------------- */

#include "gpu_particle_insertion.h"

#include "gpu_error_check.h"

#include <stdexcept>

namespace LAMMPS_NS {
namespace GPU_DEM {

namespace {

constexpr int THREADS_PER_BLOCK = 256;

__global__ void insert_grid_kernel(std::size_t offset, std::size_t count,
                                   GpuInsertionGridParams params, double *x,
                                   double *y, double *z, double *vx,
                                   double *vy, double *vz, double *omega_x,
                                   double *omega_y, double *omega_z,
                                   double *fx, double *fy, double *fz,
                                   double *torque_x, double *torque_y,
                                   double *torque_z, double *radius,
                                   double *mass, double *density, int *type,
                                   int *mask, long long *tag,
                                   long long *image_flags)
{
  const std::size_t local = static_cast<std::size_t>(blockIdx.x) * blockDim.x +
                            threadIdx.x;
  if (local >= count) return;

  const int nx = params.counts[0];
  const int ny = params.counts[1];
  const int ix = static_cast<int>(local % nx);
  const int iy = static_cast<int>((local / nx) % ny);
  const int iz = static_cast<int>(local / (static_cast<std::size_t>(nx) * ny));
  const std::size_t i = offset + local;

  x[i] = params.origin[0] + ix * params.spacing[0];
  y[i] = params.origin[1] + iy * params.spacing[1];
  z[i] = params.origin[2] + iz * params.spacing[2];
  vx[i] = params.velocity[0];
  vy[i] = params.velocity[1];
  vz[i] = params.velocity[2];
  omega_x[i] = params.omega[0];
  omega_y[i] = params.omega[1];
  omega_z[i] = params.omega[2];
  fx[i] = fy[i] = fz[i] = 0.0;
  torque_x[i] = torque_y[i] = torque_z[i] = 0.0;
  radius[i] = params.radius;
  mass[i] = params.mass;
  density[i] = params.density;
  type[i] = params.type;
  mask[i] = params.mask;
  tag[i] = params.first_tag + static_cast<long long>(local);
  image_flags[i] = 0;
}

} // namespace

std::size_t gpu_insertion_grid_count(const GpuInsertionGridParams &params)
{
  if (params.counts[0] < 0 || params.counts[1] < 0 || params.counts[2] < 0)
    throw std::runtime_error("GPU_DEM insertion grid counts must be nonnegative");
  return static_cast<std::size_t>(params.counts[0]) *
         static_cast<std::size_t>(params.counts[1]) *
         static_cast<std::size_t>(params.counts[2]);
}

void insert_particles_grid(GpuParticleData &particles,
                           const GpuInsertionGridParams &params,
                           cudaStream_t stream)
{
  if (!(params.radius > 0.0) || !(params.mass > 0.0) || !(params.density > 0.0))
    throw std::runtime_error("GPU_DEM insertion requires positive radius, mass, and density");
  if (!(params.spacing[0] > 0.0) || !(params.spacing[1] > 0.0) ||
      !(params.spacing[2] > 0.0))
    throw std::runtime_error("GPU_DEM insertion grid spacing must be positive");

  const std::size_t old_size = particles.size();
  const std::size_t inserted = gpu_insertion_grid_count(params);
  if (inserted == 0) return;

  particles.resize_preserve(old_size + inserted, stream);

  const int blocks =
      static_cast<int>((inserted + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  insert_grid_kernel<<<blocks, THREADS_PER_BLOCK, 0, stream>>>(
      old_size, inserted, params, particles.position_x(), particles.position_y(),
      particles.position_z(), particles.velocity_x(), particles.velocity_y(),
      particles.velocity_z(), particles.omega_x(), particles.omega_y(),
      particles.omega_z(), particles.force_x(), particles.force_y(),
      particles.force_z(), particles.torque_x(), particles.torque_y(),
      particles.torque_z(), particles.radius(), particles.mass(),
      particles.density(), particles.type(), particles.mask(), particles.tag(),
      particles.image_flags());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
