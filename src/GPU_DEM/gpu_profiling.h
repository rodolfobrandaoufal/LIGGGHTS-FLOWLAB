/* ----------------------------------------------------------------------
   Minimal CUDA event timer for GPU_DEM scaffold tests and early kernels.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_PROFILING_H
#define LMP_GPU_DEM_PROFILING_H

#include "gpu_error_check.h"

namespace LAMMPS_NS {
namespace GPU_DEM {

class GpuEventTimer {
 public:
  GpuEventTimer()
  {
    GPU_DEM_CUDA_CHECK(cudaEventCreate(&start_));
    GPU_DEM_CUDA_CHECK(cudaEventCreate(&stop_));
  }

  ~GpuEventTimer()
  {
    cudaEventDestroy(start_);
    cudaEventDestroy(stop_);
  }

  GpuEventTimer(const GpuEventTimer &) = delete;
  GpuEventTimer &operator=(const GpuEventTimer &) = delete;

  void start(cudaStream_t stream = 0)
  {
    GPU_DEM_CUDA_CHECK(cudaEventRecord(start_, stream));
  }

  float stop(cudaStream_t stream = 0)
  {
    GPU_DEM_CUDA_CHECK(cudaEventRecord(stop_, stream));
    GPU_DEM_CUDA_CHECK(cudaEventSynchronize(stop_));
    float milliseconds = 0.0f;
    GPU_DEM_CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start_, stop_));
    return milliseconds;
  }

 private:
  cudaEvent_t start_{};
  cudaEvent_t stop_{};
};

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
