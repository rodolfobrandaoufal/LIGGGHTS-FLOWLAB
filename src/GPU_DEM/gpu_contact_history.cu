/* ----------------------------------------------------------------------
   GPU persistent contact-history storage for GPU_DEM.
------------------------------------------------------------------------- */

#include "gpu_contact_history.h"

#include "gpu_error_check.h"

#include <utility>

namespace LAMMPS_NS {
namespace GPU_DEM {
namespace {

constexpr int THREADS_PER_BLOCK = 256;

__global__ void initialize_history_kernel(std::size_t n, unsigned long long *keys,
                                          double *shear_x, double *shear_y,
                                          double *shear_z, int *active,
                                          int *count, int *overflow)
{
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    keys[i] = 0ULL;
    shear_x[i] = 0.0;
    shear_y[i] = 0.0;
    shear_z[i] = 0.0;
    active[i] = 0;
  }
  if (i == 0) {
    *count = 0;
    *overflow = 0;
  }
}

__global__ void prepare_history_kernel(std::size_t n, double *shear_x,
                                       double *shear_y, double *shear_z,
                                       int *active, int *overflow)
{
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  if (active[i] == 0) {
    shear_x[i] = 0.0;
    shear_y[i] = 0.0;
    shear_z[i] = 0.0;
  }
  active[i] = 0;
  if (i == 0) *overflow = 0;
}

int block_count(std::size_t n)
{
  return static_cast<int>((n + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
}

} // namespace

template <class T>
void GpuContactHistory::allocate_array(T *&ptr, std::size_t n)
{
  GPU_DEM_CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&ptr), n * sizeof(T)));
}

template <class T>
void GpuContactHistory::free_array(T *&ptr) noexcept
{
  if (ptr) {
    cudaFree(ptr);
    ptr = nullptr;
  }
}

GpuContactHistory::~GpuContactHistory()
{
  release();
}

GpuContactHistory::GpuContactHistory(GpuContactHistory &&other) noexcept
{
  *this = std::move(other);
}

GpuContactHistory &GpuContactHistory::operator=(GpuContactHistory &&other) noexcept
{
  if (this == &other) return *this;
  release();

  capacity_ = other.capacity_;
  keys_ = other.keys_;
  shear_x_ = other.shear_x_;
  shear_y_ = other.shear_y_;
  shear_z_ = other.shear_z_;
  active_ = other.active_;
  count_ = other.count_;
  overflow_ = other.overflow_;

  other.capacity_ = 0;
  other.keys_ = nullptr;
  other.shear_x_ = nullptr;
  other.shear_y_ = nullptr;
  other.shear_z_ = nullptr;
  other.active_ = nullptr;
  other.count_ = nullptr;
  other.overflow_ = nullptr;
  return *this;
}

void GpuContactHistory::allocate(std::size_t capacity)
{
  if (capacity <= capacity_) return;

  release();
  capacity_ = capacity;
  allocate_array(keys_, capacity_);
  allocate_array(shear_x_, capacity_);
  allocate_array(shear_y_, capacity_);
  allocate_array(shear_z_, capacity_);
  allocate_array(active_, capacity_);
  allocate_array(count_, 1);
  allocate_array(overflow_, 1);

  initialize_history_kernel<<<block_count(capacity_), THREADS_PER_BLOCK>>>(
      capacity_, keys_, shear_x_, shear_y_, shear_z_, active_, count_, overflow_);
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
  GPU_DEM_CUDA_CHECK(cudaDeviceSynchronize());
}

void GpuContactHistory::release() noexcept
{
  free_array(keys_);
  free_array(shear_x_);
  free_array(shear_y_);
  free_array(shear_z_);
  free_array(active_);
  free_array(count_);
  free_array(overflow_);
  capacity_ = 0;
}

void GpuContactHistory::prepare_for_step(cudaStream_t stream)
{
  if (capacity_ == 0) return;
  prepare_history_kernel<<<block_count(capacity_), THREADS_PER_BLOCK, 0, stream>>>(
      capacity_, shear_x_, shear_y_, shear_z_, active_, overflow_);
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
}

int GpuContactHistory::active_count(cudaStream_t stream) const
{
  int count = 0;
  if (count_) {
    GPU_DEM_CUDA_CHECK(cudaMemcpyAsync(&count, count_, sizeof(int),
                                       cudaMemcpyDeviceToHost, stream));
    GPU_DEM_CUDA_CHECK(cudaStreamSynchronize(stream));
  }
  return count;
}

int GpuContactHistory::overflow(cudaStream_t stream) const
{
  int flag = 0;
  if (overflow_) {
    GPU_DEM_CUDA_CHECK(cudaMemcpyAsync(&flag, overflow_, sizeof(int),
                                       cudaMemcpyDeviceToHost, stream));
    GPU_DEM_CUDA_CHECK(cudaStreamSynchronize(stream));
  }
  return flag;
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
