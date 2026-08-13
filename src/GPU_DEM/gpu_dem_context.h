/* ----------------------------------------------------------------------
   CUDA context wrapper for the GPU-native DEM path.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_CONTEXT_H
#define LMP_GPU_DEM_CONTEXT_H

#include "gpu_precision.h"
#include "gpu_error_check.h"

#include <cuda_runtime.h>

namespace LAMMPS_NS {
namespace GPU_DEM {

class GpuDemContext {
 public:
  explicit GpuDemContext(int device_id = 0,
                         PrecisionMode precision = PrecisionMode::Mixed);
  ~GpuDemContext();

  GpuDemContext(const GpuDemContext &) = delete;
  GpuDemContext &operator=(const GpuDemContext &) = delete;

  GpuDemContext(GpuDemContext &&other) noexcept;
  GpuDemContext &operator=(GpuDemContext &&other) noexcept;

  int device_id() const { return device_id_; }
  PrecisionMode precision_mode() const { return precision_; }
  cudaStream_t stream() const { return stream_; }
  const cudaDeviceProp &device_properties() const { return properties_; }

  void synchronize() const;

 private:
  void release() noexcept;

  int device_id_{-1};
  PrecisionMode precision_{PrecisionMode::Mixed};
  cudaStream_t stream_{nullptr};
  cudaDeviceProp properties_{};
};

int gpu_device_count();

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
