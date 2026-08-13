/* ----------------------------------------------------------------------
   GPU uniform-grid neighbor candidate construction for GPU_DEM.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_NEIGHBOR_BUILDER_H
#define LMP_GPU_DEM_NEIGHBOR_BUILDER_H

#include "gpu_particle_data.h"

#include <cuda_runtime.h>

#include <cstddef>
#include <vector>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct NeighborBuildParams {
  double box_min[3]{0.0, 0.0, 0.0};
  double box_max[3]{1.0, 1.0, 1.0};
  double cell_size{1.0};
  double skin{0.0};
  int periodic[3]{0, 0, 0};
};

struct HostNeighborPairs {
  std::vector<int> first;
  std::vector<int> second;
  int overflow{0};

  std::size_t size() const { return first.size(); }
};

class GpuNeighborList {
 public:
  GpuNeighborList() = default;
  ~GpuNeighborList();

  GpuNeighborList(const GpuNeighborList &) = delete;
  GpuNeighborList &operator=(const GpuNeighborList &) = delete;

  GpuNeighborList(GpuNeighborList &&other) noexcept;
  GpuNeighborList &operator=(GpuNeighborList &&other) noexcept;

  void allocate(std::size_t particle_capacity, std::size_t pair_capacity, int num_cells);
  void release() noexcept;

  void copy_pairs_to_host(HostNeighborPairs &host, cudaStream_t stream = 0) const;

  std::size_t particle_capacity() const { return particle_capacity_; }
  std::size_t pair_capacity() const { return pair_capacity_; }
  int num_cells() const { return num_cells_; }
  int pair_count(cudaStream_t stream = 0) const;
  int overflow(cudaStream_t stream = 0) const;

  int *cell_keys() const { return cell_keys_; }
  int *particle_ids() const { return particle_ids_; }
  int *cell_head() const { return cell_head_; }
  int *particle_next() const { return particle_next_; }
  int *pair_first() const { return pair_first_; }
  int *pair_second() const { return pair_second_; }
  int *device_pair_count() const { return pair_count_; }
  int *device_overflow() const { return overflow_; }

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

  std::size_t particle_capacity_{0};
  std::size_t pair_capacity_{0};
  int num_cells_{0};

  int *cell_keys_{nullptr};
  int *particle_ids_{nullptr};
  int *cell_head_{nullptr};
  int *particle_next_{nullptr};
  int *pair_first_{nullptr};
  int *pair_second_{nullptr};
  int *pair_count_{nullptr};
  int *overflow_{nullptr};
};

void build_contact_pairs_uniform_grid(const GpuParticleData &particles,
                                      const NeighborBuildParams &params,
                                      GpuNeighborList &neighbors,
                                      cudaStream_t stream = 0);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
