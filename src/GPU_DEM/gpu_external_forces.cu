/* ----------------------------------------------------------------------
   GPU external-force kernels for GPU_DEM.
------------------------------------------------------------------------- */

#include "gpu_external_forces.h"

#include "gpu_error_check.h"

namespace LAMMPS_NS {
namespace GPU_DEM {
namespace {

constexpr int THREADS_PER_BLOCK = 256;

__global__ void add_gravity_force_kernel(int n, int groupbit, GravityParams params,
                                         const int *mask, const double *mass,
                                         double *force_x, double *force_y,
                                         double *force_z)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n || !(mask[i] & groupbit)) return;

  force_x[i] += mass[i] * params.acceleration[0];
  force_y[i] += mass[i] * params.acceleration[1];
  force_z[i] += mass[i] * params.acceleration[2];
}

int block_count(std::size_t n)
{
  return static_cast<int>((n + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
}

} // namespace

void add_gravity_force(GpuParticleData &particles, const GravityParams &params,
                       int groupbit, cudaStream_t stream)
{
  if (!params.enabled || particles.size() == 0) return;

  add_gravity_force_kernel<<<block_count(particles.size()), THREADS_PER_BLOCK, 0, stream>>>(
      static_cast<int>(particles.size()), groupbit, params, particles.mask(),
      particles.mass(), particles.force_x(), particles.force_y(), particles.force_z());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
