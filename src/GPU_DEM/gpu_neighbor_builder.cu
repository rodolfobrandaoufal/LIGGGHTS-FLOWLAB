/* ----------------------------------------------------------------------
   GPU uniform-grid neighbor candidate construction for GPU_DEM.
------------------------------------------------------------------------- */

#include "gpu_neighbor_builder.h"

#include "gpu_error_check.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace LAMMPS_NS {
namespace GPU_DEM {

namespace {

constexpr int THREADS_PER_BLOCK = 256;

__host__ __device__ int flatten_cell(int ix, int iy, int iz, int nx, int ny)
{
  return ix + nx * (iy + ny * iz);
}

__host__ __device__ int wrap_index(int index, int count)
{
  int wrapped = index % count;
  return wrapped < 0 ? wrapped + count : wrapped;
}

__host__ __device__ int clamp_index(int index, int count)
{
  if (index < 0) return 0;
  if (index >= count) return count - 1;
  return index;
}

int host_max_int(int lhs, int rhs)
{
  return lhs > rhs ? lhs : rhs;
}

std::size_t host_min_size(std::size_t lhs, std::size_t rhs)
{
  return lhs < rhs ? lhs : rhs;
}

void compute_grid_shape(const NeighborBuildParams &params, int shape[3])
{
  if (!(params.cell_size > 0.0))
    throw std::runtime_error("GPU_DEM neighbor build requires positive cell_size");

  for (int d = 0; d < 3; ++d) {
    const double length = params.box_max[d] - params.box_min[d];
    if (!(length > 0.0))
      throw std::runtime_error("GPU_DEM neighbor build requires positive box extent");
    shape[d] = host_max_int(1, static_cast<int>(std::floor(length / params.cell_size)));
  }
}

__global__ void reset_particle_next_kernel(std::size_t n, int *particle_next)
{
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) particle_next[i] = -1;
}

__global__ void assign_cells_kernel(std::size_t n, const double *x, const double *y,
                                    const double *z, NeighborBuildParams params, int nx,
                                    int ny, int nz, int *cell_keys)
{
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  int ix = static_cast<int>(floor((x[i] - params.box_min[0]) / params.cell_size));
  int iy = static_cast<int>(floor((y[i] - params.box_min[1]) / params.cell_size));
  int iz = static_cast<int>(floor((z[i] - params.box_min[2]) / params.cell_size));

  ix = params.periodic[0] ? wrap_index(ix, nx) : clamp_index(ix, nx);
  iy = params.periodic[1] ? wrap_index(iy, ny) : clamp_index(iy, ny);
  iz = params.periodic[2] ? wrap_index(iz, nz) : clamp_index(iz, nz);

  cell_keys[i] = flatten_cell(ix, iy, iz, nx, ny);
}

__global__ void reset_cell_heads_kernel(int num_cells, int *cell_head)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < num_cells) cell_head[i] = -1;
}

__global__ void reset_pair_state_kernel(int *pair_count, int *overflow)
{
  *pair_count = 0;
  *overflow = 0;
}

__global__ void build_cell_linked_list_kernel(std::size_t n, const int *cell_keys,
                                              int *cell_head, int *particle_next)
{
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const int key = cell_keys[i];
  particle_next[i] = atomicExch(&cell_head[key], static_cast<int>(i));
}

__device__ double minimum_image_delta(double delta, double length, int periodic)
{
  if (!periodic) return delta;
  if (delta > 0.5 * length) return delta - length;
  if (delta < -0.5 * length) return delta + length;
  return delta;
}

__global__ void emit_candidate_pairs_kernel(
    std::size_t n, const double *x, const double *y, const double *z, const double *radius,
    NeighborBuildParams params, int nx, int ny, int nz, const int *cell_keys,
    const int *cell_head, const int *particle_next, std::size_t pair_capacity,
    int *pair_first, int *pair_second, int *pair_count, int *overflow)
{
  const std::size_t particle_index = blockIdx.x * blockDim.x + threadIdx.x;
  if (particle_index >= n) return;

  const int i = static_cast<int>(particle_index);
  const int key = cell_keys[i];
  const int iz = key / (nx * ny);
  const int rem = key - iz * nx * ny;
  const int iy = rem / nx;
  const int ix = rem - iy * nx;

  const double box_length_x = params.box_max[0] - params.box_min[0];
  const double box_length_y = params.box_max[1] - params.box_min[1];
  const double box_length_z = params.box_max[2] - params.box_min[2];

  for (int dz = -1; dz <= 1; ++dz) {
    int nz_index = iz + dz;
    if (params.periodic[2])
      nz_index = wrap_index(nz_index, nz);
    else if (nz_index < 0 || nz_index >= nz)
      continue;

    for (int dy = -1; dy <= 1; ++dy) {
      int ny_index = iy + dy;
      if (params.periodic[1])
        ny_index = wrap_index(ny_index, ny);
      else if (ny_index < 0 || ny_index >= ny)
        continue;

      for (int dx = -1; dx <= 1; ++dx) {
        int nx_index = ix + dx;
        if (params.periodic[0])
          nx_index = wrap_index(nx_index, nx);
        else if (nx_index < 0 || nx_index >= nx)
          continue;

        const int neighbor_key = flatten_cell(nx_index, ny_index, nz_index, nx, ny);
        for (int j = cell_head[neighbor_key]; j >= 0; j = particle_next[j]) {
          if (i >= j) continue;

          double dxij = x[j] - x[i];
          double dyij = y[j] - y[i];
          double dzij = z[j] - z[i];
          dxij = minimum_image_delta(dxij, box_length_x, params.periodic[0]);
          dyij = minimum_image_delta(dyij, box_length_y, params.periodic[1]);
          dzij = minimum_image_delta(dzij, box_length_z, params.periodic[2]);

          const double cutoff = radius[i] + radius[j] + params.skin;
          const double rsq = dxij * dxij + dyij * dyij + dzij * dzij;
          if (rsq > cutoff * cutoff) continue;

          const int slot = atomicAdd(pair_count, 1);
          if (static_cast<std::size_t>(slot) < pair_capacity) {
            pair_first[slot] = i;
            pair_second[slot] = j;
          } else {
            *overflow = 1;
          }
        }
      }
    }
  }
}

} // namespace

GpuNeighborList::~GpuNeighborList()
{
  release();
}

GpuNeighborList::GpuNeighborList(GpuNeighborList &&other) noexcept
{
  *this = std::move(other);
}

GpuNeighborList &GpuNeighborList::operator=(GpuNeighborList &&other) noexcept
{
  if (this == &other) return *this;
  release();

  particle_capacity_ = other.particle_capacity_;
  pair_capacity_ = other.pair_capacity_;
  num_cells_ = other.num_cells_;
  cell_keys_ = other.cell_keys_;
  particle_ids_ = other.particle_ids_;
  cell_head_ = other.cell_head_;
  particle_next_ = other.particle_next_;
  pair_first_ = other.pair_first_;
  pair_second_ = other.pair_second_;
  pair_count_ = other.pair_count_;
  overflow_ = other.overflow_;

  other.particle_capacity_ = 0;
  other.pair_capacity_ = 0;
  other.num_cells_ = 0;
  other.cell_keys_ = nullptr;
  other.particle_ids_ = nullptr;
  other.cell_head_ = nullptr;
  other.particle_next_ = nullptr;
  other.pair_first_ = nullptr;
  other.pair_second_ = nullptr;
  other.pair_count_ = nullptr;
  other.overflow_ = nullptr;

  return *this;
}

void GpuNeighborList::allocate(std::size_t particle_capacity, std::size_t pair_capacity,
                               int num_cells)
{
  if (particle_capacity <= particle_capacity_ && pair_capacity <= pair_capacity_ &&
      num_cells <= num_cells_)
    return;

  release();
  particle_capacity_ = particle_capacity;
  pair_capacity_ = pair_capacity;
  num_cells_ = num_cells;

  allocate_array(cell_keys_, particle_capacity_);
  allocate_array(particle_ids_, particle_capacity_);
  allocate_array(cell_head_, num_cells_);
  allocate_array(particle_next_, particle_capacity_);
  allocate_array(pair_first_, pair_capacity_);
  allocate_array(pair_second_, pair_capacity_);
  allocate_array(pair_count_, 1);
  allocate_array(overflow_, 1);
}

void GpuNeighborList::release() noexcept
{
  free_array(cell_keys_);
  free_array(particle_ids_);
  free_array(cell_head_);
  free_array(particle_next_);
  free_array(pair_first_);
  free_array(pair_second_);
  free_array(pair_count_);
  free_array(overflow_);

  particle_capacity_ = 0;
  pair_capacity_ = 0;
  num_cells_ = 0;
}

int GpuNeighborList::pair_count(cudaStream_t stream) const
{
  int count = 0;
  if (pair_count_) {
    GPU_DEM_CUDA_CHECK(
        cudaMemcpyAsync(&count, pair_count_, sizeof(int), cudaMemcpyDeviceToHost, stream));
    GPU_DEM_CUDA_CHECK(cudaStreamSynchronize(stream));
  }
  return count;
}

int GpuNeighborList::overflow(cudaStream_t stream) const
{
  int flag = 0;
  if (overflow_) {
    GPU_DEM_CUDA_CHECK(
        cudaMemcpyAsync(&flag, overflow_, sizeof(int), cudaMemcpyDeviceToHost, stream));
    GPU_DEM_CUDA_CHECK(cudaStreamSynchronize(stream));
  }
  return flag;
}

void GpuNeighborList::copy_pairs_to_host(HostNeighborPairs &host, cudaStream_t stream) const
{
  const int raw_count = pair_count(stream);
  const std::size_t copied_count =
      host_min_size(static_cast<std::size_t>(host_max_int(raw_count, 0)), pair_capacity_);

  host.first.resize(copied_count);
  host.second.resize(copied_count);
  host.overflow = overflow(stream);

  if (copied_count > 0) {
    GPU_DEM_CUDA_CHECK(cudaMemcpyAsync(host.first.data(), pair_first_,
                                       copied_count * sizeof(int), cudaMemcpyDeviceToHost,
                                       stream));
    GPU_DEM_CUDA_CHECK(cudaMemcpyAsync(host.second.data(), pair_second_,
                                       copied_count * sizeof(int), cudaMemcpyDeviceToHost,
                                       stream));
  }
  GPU_DEM_CUDA_CHECK(cudaStreamSynchronize(stream));
}

void build_contact_pairs_uniform_grid(const GpuParticleData &particles,
                                      const NeighborBuildParams &params,
                                      GpuNeighborList &neighbors, cudaStream_t stream)
{
  const std::size_t n = particles.size();
  int shape[3];
  compute_grid_shape(params, shape);
  const int num_cells = shape[0] * shape[1] * shape[2];

  if (neighbors.particle_capacity() < n || neighbors.num_cells() < num_cells)
    throw std::runtime_error("GPU_DEM neighbor list capacity is too small");

  const int particle_blocks = static_cast<int>((n + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  const int cell_blocks = (num_cells + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  reset_particle_next_kernel<<<particle_blocks, THREADS_PER_BLOCK, 0, stream>>>(
      n, neighbors.particle_next());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());

  assign_cells_kernel<<<particle_blocks, THREADS_PER_BLOCK, 0, stream>>>(
      n, particles.position_x(), particles.position_y(), particles.position_z(), params,
      shape[0], shape[1], shape[2], neighbors.cell_keys());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());

  reset_cell_heads_kernel<<<cell_blocks, THREADS_PER_BLOCK, 0, stream>>>(num_cells,
                                                                        neighbors.cell_head());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());

  reset_pair_state_kernel<<<1, 1, 0, stream>>>(neighbors.device_pair_count(),
                                               neighbors.device_overflow());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());

  build_cell_linked_list_kernel<<<particle_blocks, THREADS_PER_BLOCK, 0, stream>>>(
      n, neighbors.cell_keys(), neighbors.cell_head(), neighbors.particle_next());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());

  emit_candidate_pairs_kernel<<<particle_blocks, THREADS_PER_BLOCK, 0, stream>>>(
      n, particles.position_x(), particles.position_y(), particles.position_z(),
      particles.radius(), params, shape[0], shape[1], shape[2], neighbors.cell_keys(),
      neighbors.cell_head(), neighbors.particle_next(), neighbors.pair_capacity(),
      neighbors.pair_first(), neighbors.pair_second(), neighbors.device_pair_count(),
      neighbors.device_overflow());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
