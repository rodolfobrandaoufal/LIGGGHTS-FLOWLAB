/* ----------------------------------------------------------------------
   Device-resident structure-of-arrays particle storage for GPU_DEM.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_PARTICLE_DATA_H
#define LMP_GPU_DEM_PARTICLE_DATA_H

#include "gpu_error_check.h"

#include <cuda_runtime.h>

#include <cstddef>
#include <vector>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct HostParticleData {
  std::vector<double> position_x, position_y, position_z;
  std::vector<double> velocity_x, velocity_y, velocity_z;
  std::vector<double> omega_x, omega_y, omega_z;
  std::vector<double> force_x, force_y, force_z;
  std::vector<double> torque_x, torque_y, torque_z;
  std::vector<double> radius;
  std::vector<double> mass;
  std::vector<double> density;
  std::vector<int> type;
  std::vector<int> mask;
  std::vector<long long> tag;
  std::vector<long long> image_flags;

  void resize(std::size_t n);
  std::size_t size() const { return position_x.size(); }
  void validate_sizes() const;
};

class GpuParticleData {
 public:
  GpuParticleData() = default;
  ~GpuParticleData();

  GpuParticleData(const GpuParticleData &) = delete;
  GpuParticleData &operator=(const GpuParticleData &) = delete;

  GpuParticleData(GpuParticleData &&other) noexcept;
  GpuParticleData &operator=(GpuParticleData &&other) noexcept;

  void allocate(std::size_t n);
  void resize_preserve(std::size_t n, cudaStream_t stream = 0);
  void release() noexcept;

  void copy_from_host(const HostParticleData &host, cudaStream_t stream = 0);
  void copy_to_host(HostParticleData &host, cudaStream_t stream = 0) const;

  std::size_t size() const { return size_; }
  std::size_t capacity() const { return capacity_; }

  double *position_x() const { return position_x_; }
  double *position_y() const { return position_y_; }
  double *position_z() const { return position_z_; }
  double *velocity_x() const { return velocity_x_; }
  double *velocity_y() const { return velocity_y_; }
  double *velocity_z() const { return velocity_z_; }
  double *omega_x() const { return omega_x_; }
  double *omega_y() const { return omega_y_; }
  double *omega_z() const { return omega_z_; }
  double *force_x() const { return force_x_; }
  double *force_y() const { return force_y_; }
  double *force_z() const { return force_z_; }
  double *torque_x() const { return torque_x_; }
  double *torque_y() const { return torque_y_; }
  double *torque_z() const { return torque_z_; }
  double *radius() const { return radius_; }
  double *mass() const { return mass_; }
  double *density() const { return density_; }
  int *type() const { return type_; }
  int *mask() const { return mask_; }
  long long *tag() const { return tag_; }
  long long *image_flags() const { return image_flags_; }

 private:
  template <class T>
  static void allocate_array(T *&ptr, std::size_t n)
  {
    GPU_DEM_CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&ptr), n * sizeof(T)));
  }

  template <class T>
  static void free_array(T *&ptr) noexcept
  {
    if (ptr) {
      cudaFree(ptr);
      ptr = nullptr;
    }
  }

  template <class T>
  static void copy_to_device(T *dst, const std::vector<T> &src, cudaStream_t stream)
  {
    if (!src.empty())
      GPU_DEM_CUDA_CHECK(cudaMemcpyAsync(dst, src.data(), src.size() * sizeof(T),
                                         cudaMemcpyHostToDevice, stream));
  }

  template <class T>
  static void copy_to_host(std::vector<T> &dst, const T *src, cudaStream_t stream)
  {
    if (!dst.empty())
      GPU_DEM_CUDA_CHECK(cudaMemcpyAsync(dst.data(), src, dst.size() * sizeof(T),
                                         cudaMemcpyDeviceToHost, stream));
  }

  template <class T>
  static void copy_device_to_device(T *dst, const T *src, std::size_t n,
                                    cudaStream_t stream)
  {
    if (n > 0)
      GPU_DEM_CUDA_CHECK(cudaMemcpyAsync(dst, src, n * sizeof(T),
                                         cudaMemcpyDeviceToDevice, stream));
  }

  std::size_t size_{0};
  std::size_t capacity_{0};

  double *position_x_{nullptr}, *position_y_{nullptr}, *position_z_{nullptr};
  double *velocity_x_{nullptr}, *velocity_y_{nullptr}, *velocity_z_{nullptr};
  double *omega_x_{nullptr}, *omega_y_{nullptr}, *omega_z_{nullptr};
  double *force_x_{nullptr}, *force_y_{nullptr}, *force_z_{nullptr};
  double *torque_x_{nullptr}, *torque_y_{nullptr}, *torque_z_{nullptr};
  double *radius_{nullptr}, *mass_{nullptr}, *density_{nullptr};
  int *type_{nullptr};
  int *mask_{nullptr};
  long long *tag_{nullptr}, *image_flags_{nullptr};
};

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
