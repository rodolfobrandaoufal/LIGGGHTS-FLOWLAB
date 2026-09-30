/* ----------------------------------------------------------------------
   GPU contact-force pipeline scaffold for GPU_DEM.
------------------------------------------------------------------------- */

#include "gpu_contact_pipeline.h"

#include "gpu_error_check.h"
#include "models/normal_hooke.cuh"

#include <stdexcept>

namespace LAMMPS_NS {
namespace GPU_DEM {

namespace {

constexpr int THREADS_PER_BLOCK = 256;

__global__ void hooke_normal_force_kernel(int pair_count, const int *pair_first,
                                          const int *pair_second, const double *x,
                                          const double *y, const double *z,
                                          const double *vx, const double *vy,
                                          const double *vz, const double *radius,
                                          HookeNormalParams params, double *fx,
                                          double *fy, double *fz)
{
  const int pair_index = blockIdx.x * blockDim.x + threadIdx.x;
  if (pair_index >= pair_count) return;

  const int i = pair_first[pair_index];
  const int j = pair_second[pair_index];

  const double dx = x[j] - x[i];
  const double dy = y[j] - y[i];
  const double dz = z[j] - z[i];
  const double rsq = dx * dx + dy * dy + dz * dz;
  if (rsq <= 0.0) return;

  const double r = sqrt(rsq);
  const double overlap = radius[i] + radius[j] - r;
  if (overlap <= 0.0) return;

  const double nx = dx / r;
  const double ny = dy / r;
  const double nz = dz / r;

  const double dvx = vx[j] - vx[i];
  const double dvy = vy[j] - vy[i];
  const double dvz = vz[j] - vz[i];
  const double normal_relative_velocity = dvx * nx + dvy * ny + dvz * nz;

  const double normal_force =
      Models::hooke_normal_force(overlap, normal_relative_velocity, params);

  const double fix = -normal_force * nx;
  const double fiy = -normal_force * ny;
  const double fiz = -normal_force * nz;

  atomicAdd(&fx[i], fix);
  atomicAdd(&fy[i], fiy);
  atomicAdd(&fz[i], fiz);
  atomicAdd(&fx[j], -fix);
  atomicAdd(&fy[j], -fiy);
  atomicAdd(&fz[j], -fiz);
}

} // namespace

void compute_hooke_normal_forces(GpuParticleData &particles,
                                 const GpuNeighborList &neighbors,
                                 const HookeNormalParams &params,
                                 cudaStream_t stream)
{
  if (!(params.normal_stiffness >= 0.0) || !(params.normal_damping >= 0.0))
    throw std::runtime_error("GPU_DEM Hooke normal parameters must be nonnegative");
  if (neighbors.overflow(stream))
    throw std::runtime_error("GPU_DEM contact-force evaluation refused overflowed pair list");

  const int pairs = neighbors.pair_count(stream);
  if (pairs <= 0) return;
  if (static_cast<std::size_t>(pairs) > neighbors.pair_capacity())
    throw std::runtime_error("GPU_DEM contact-force pair count exceeds pair capacity");

  const int blocks = (pairs + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  hooke_normal_force_kernel<<<blocks, THREADS_PER_BLOCK, 0, stream>>>(
      pairs, neighbors.pair_first(), neighbors.pair_second(), particles.position_x(),
      particles.position_y(), particles.position_z(), particles.velocity_x(),
      particles.velocity_y(), particles.velocity_z(), particles.radius(), params,
      particles.force_x(), particles.force_y(), particles.force_z());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
