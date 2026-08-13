/* ----------------------------------------------------------------------
   GPU persistent contact-history storage for GPU_DEM.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_CONTACT_HISTORY_H
#define LMP_GPU_DEM_CONTACT_HISTORY_H

#include <cuda_runtime.h>

#include <cstddef>

namespace LAMMPS_NS {
namespace GPU_DEM {

class GpuContactHistory {
 public:
  GpuContactHistory() = default;
  ~GpuContactHistory();

  GpuContactHistory(const GpuContactHistory &) = delete;
  GpuContactHistory &operator=(const GpuContactHistory &) = delete;

  GpuContactHistory(GpuContactHistory &&other) noexcept;
  GpuContactHistory &operator=(GpuContactHistory &&other) noexcept;

  void allocate(std::size_t capacity);
  void release() noexcept;
  void prepare_for_step(cudaStream_t stream = 0);

  std::size_t capacity() const { return capacity_; }
  int active_count(cudaStream_t stream = 0) const;
  int overflow(cudaStream_t stream = 0) const;

  unsigned long long *keys() const { return keys_; }
  double *shear_x() const { return shear_x_; }
  double *shear_y() const { return shear_y_; }
  double *shear_z() const { return shear_z_; }
  int *active() const { return active_; }
  int *count() const { return count_; }
  int *overflow_flag() const { return overflow_; }

 private:
  template <class T>
  static void allocate_array(T *&ptr, std::size_t n);

  template <class T>
  static void free_array(T *&ptr) noexcept;

  std::size_t capacity_{0};
  unsigned long long *keys_{nullptr};
  double *shear_x_{nullptr};
  double *shear_y_{nullptr};
  double *shear_z_{nullptr};
  int *active_{nullptr};
  int *count_{nullptr};
  int *overflow_{nullptr};
};

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
