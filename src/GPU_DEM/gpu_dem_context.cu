#include "gpu_dem_context.h"

#include <utility>

namespace LAMMPS_NS {
namespace GPU_DEM {

int gpu_device_count()
{
  int count = 0;
  GPU_DEM_CUDA_CHECK(cudaGetDeviceCount(&count));
  return count;
}

GpuDemContext::GpuDemContext(int device_id, PrecisionMode precision)
    : device_id_(device_id), precision_(precision)
{
  int count = gpu_device_count();
  if (device_id_ < 0 || device_id_ >= count)
    throw std::runtime_error("GPU_DEM requested CUDA device is out of range");

  GPU_DEM_CUDA_CHECK(cudaSetDevice(device_id_));
  GPU_DEM_CUDA_CHECK(cudaGetDeviceProperties(&properties_, device_id_));
  GPU_DEM_CUDA_CHECK(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking));
}

GpuDemContext::~GpuDemContext()
{
  release();
}

GpuDemContext::GpuDemContext(GpuDemContext &&other) noexcept
{
  *this = std::move(other);
}

GpuDemContext &GpuDemContext::operator=(GpuDemContext &&other) noexcept
{
  if (this == &other) return *this;

  release();
  device_id_ = other.device_id_;
  precision_ = other.precision_;
  stream_ = other.stream_;
  properties_ = other.properties_;

  other.device_id_ = -1;
  other.stream_ = nullptr;
  return *this;
}

void GpuDemContext::synchronize() const
{
  if (stream_) GPU_DEM_CUDA_CHECK(cudaStreamSynchronize(stream_));
}

void GpuDemContext::release() noexcept
{
  if (stream_) {
    cudaStreamDestroy(stream_);
    stream_ = nullptr;
  }
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
