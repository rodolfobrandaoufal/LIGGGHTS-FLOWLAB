/* ----------------------------------------------------------------------
   GPU contact-force pipeline scaffold for GPU_DEM.
------------------------------------------------------------------------- */

#include "gpu_contact_pipeline.h"

#include "gpu_error_check.h"
#include "models/normal_hertz.cuh"
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

__global__ void hertz_normal_force_kernel(int pair_count, const int *pair_first,
                                          const int *pair_second, const double *x,
                                          const double *y, const double *z,
                                          const double *vx, const double *vy,
                                          const double *vz, const double *radius,
                                          const double *mass,
                                          HertzNormalParams params, double *fx,
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
      Models::hertz_normal_force(overlap, normal_relative_velocity,
                                 radius[i], radius[j], mass[i], mass[j], params);

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

__device__ unsigned long long ordered_pair_key(long long tag_i, long long tag_j)
{
  unsigned long long a = static_cast<unsigned long long>(tag_i);
  unsigned long long b = static_cast<unsigned long long>(tag_j);
  if (a > b) {
    const unsigned long long tmp = a;
    a = b;
    b = tmp;
  }
  return (a << 32) ^ b;
}

__device__ int find_or_create_history_slot(
    unsigned long long key, std::size_t capacity, unsigned long long *keys,
    double *shear_x, double *shear_y, double *shear_z, int *active,
    int *count, int *overflow)
{
  const int current_count = *count;
  for (int slot = 0; slot < current_count; ++slot) {
    if (keys[slot] == key) {
      active[slot] = 1;
      return slot;
    }
  }

  const int slot = atomicAdd(count, 1);
  if (static_cast<std::size_t>(slot) >= capacity) {
    *overflow = 1;
    return -1;
  }

  keys[slot] = key;
  shear_x[slot] = 0.0;
  shear_y[slot] = 0.0;
  shear_z[slot] = 0.0;
  active[slot] = 1;
  return slot;
}

__device__ void cross3(double ax, double ay, double az,
                       double bx, double by, double bz,
                       double &cx, double &cy, double &cz)
{
  cx = ay * bz - az * by;
  cy = az * bx - ax * bz;
  cz = ax * by - ay * bx;
}

__device__ void add_cdt_rolling_torque(int i, int j, double nx, double ny,
                                       double nz, double fn, double reff,
                                       const double *omega_x,
                                       const double *omega_y,
                                       const double *omega_z,
                                       RollingFrictionParams params,
                                       double *torque_x, double *torque_y,
                                       double *torque_z)
{
  if (!(params.coefficient > 0.0)) return;

  const double wrx = omega_x[i] - omega_x[j];
  const double wry = omega_y[i] - omega_y[j];
  const double wrz = omega_z[i] - omega_z[j];
  const double wrsq = wrx * wrx + wry * wry + wrz * wrz;
  if (wrsq <= 0.0) return;

  const double scale = params.coefficient * fabs(fn) * reff / sqrt(wrsq);
  double rtx = scale * wrx;
  double rty = scale * wry;
  double rtz = scale * wrz;

  if (!params.torsion_torque) {
    const double normal_part = rtx * nx + rty * ny + rtz * nz;
    rtx -= normal_part * nx;
    rty -= normal_part * ny;
    rtz -= normal_part * nz;
  }

  atomicAdd(&torque_x[i], -rtx);
  atomicAdd(&torque_y[i], -rty);
  atomicAdd(&torque_z[i], -rtz);
  atomicAdd(&torque_x[j], rtx);
  atomicAdd(&torque_y[j], rty);
  atomicAdd(&torque_z[j], rtz);
}

__global__ void hertz_history_force_kernel(
    int pair_count, const int *pair_first, const int *pair_second, const double *x,
    const double *y, const double *z, const double *vx, const double *vy,
    const double *vz, const double *omega_x, const double *omega_y,
    const double *omega_z, const double *radius, const double *mass,
	    const long long *tag, HertzNormalParams normal_params,
	    TangentialHistoryParams tangential_params,
	    RollingFrictionParams rolling_params, std::size_t history_capacity,
	    unsigned long long *history_keys, double *history_shear_x,
    double *history_shear_y, double *history_shear_z, int *history_active,
    int *history_count, int *history_overflow, double *fx, double *fy,
    double *fz, double *torque_x, double *torque_y, double *torque_z)
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

  const unsigned long long key = ordered_pair_key(tag[i], tag[j]);
  const int history_slot = find_or_create_history_slot(
      key, history_capacity, history_keys, history_shear_x, history_shear_y,
      history_shear_z, history_active, history_count, history_overflow);
  if (history_slot < 0) return;

  const double nx = dx / r;
  const double ny = dy / r;
  const double nz = dz / r;

  const double dvx = vx[j] - vx[i];
  const double dvy = vy[j] - vy[i];
  const double dvz = vz[j] - vz[i];
  const double normal_relative_velocity = dvx * nx + dvy * ny + dvz * nz;

  const double fn = Models::hertz_normal_force(overlap, normal_relative_velocity,
                                               radius[i], radius[j], mass[i],
                                               mass[j], normal_params);

  const double reff = Models::effective_radius(radius[i], radius[j]);
  const double meff = Models::effective_mass(mass[i], mass[j]);
  const double sqrt_delta_reff = sqrt(reff * overlap);
  const double kt = 8.0 * tangential_params.effective_shear_modulus * sqrt_delta_reff;
  const double st = kt;
  const double sqrt_five_over_six = 0.91287092917527685576161630466800355659;
  const double gammat = tangential_params.tangential_damping
      ? -2.0 * sqrt_five_over_six * normal_params.beta_effective * sqrt(st * meff)
      : 0.0;

  double wi_cross_n_x, wi_cross_n_y, wi_cross_n_z;
  double wj_cross_n_x, wj_cross_n_y, wj_cross_n_z;
  cross3(omega_x[i], omega_y[i], omega_z[i], nx, ny, nz,
         wi_cross_n_x, wi_cross_n_y, wi_cross_n_z);
  cross3(omega_x[j], omega_y[j], omega_z[j], nx, ny, nz,
         wj_cross_n_x, wj_cross_n_y, wj_cross_n_z);

  const double vrx = dvx - normal_relative_velocity * nx -
                     radius[i] * wi_cross_n_x - radius[j] * wj_cross_n_x;
  const double vry = dvy - normal_relative_velocity * ny -
                     radius[i] * wi_cross_n_y - radius[j] * wj_cross_n_y;
  const double vrz = dvz - normal_relative_velocity * nz -
                     radius[i] * wi_cross_n_z - radius[j] * wj_cross_n_z;

  double sx = history_shear_x[history_slot] + vrx * tangential_params.dt;
  double sy = history_shear_y[history_slot] + vry * tangential_params.dt;
  double sz = history_shear_z[history_slot] + vrz * tangential_params.dt;
  const double shear_normal = sx * nx + sy * ny + sz * nz;
  sx -= shear_normal * nx;
  sy -= shear_normal * ny;
  sz -= shear_normal * nz;

  const double shrsq = sx * sx + sy * sy + sz * sz;
  double ftx = -kt * sx;
  double fty = -kt * sy;
  double ftz = -kt * sz;
  const double ft_friction = tangential_params.friction_coefficient * fabs(fn);
  const double ft_shear_sq = kt * kt * shrsq;
  const double ft_friction_sq = ft_friction * ft_friction;

  if (ft_shear_sq > ft_friction_sq) {
    if (shrsq > 0.0) {
      const double ft_shear = kt * sqrt(shrsq);
      const double ratio = ft_friction / ft_shear;
      ftx *= ratio;
      fty *= ratio;
      ftz *= ratio;
      sx = -ftx / kt;
      sy = -fty / kt;
      sz = -ftz / kt;
    } else {
      ftx = fty = ftz = 0.0;
    }
  } else {
    ftx -= gammat * vrx;
    fty -= gammat * vry;
    ftz -= gammat * vrz;
  }

  history_shear_x[history_slot] = sx;
  history_shear_y[history_slot] = sy;
  history_shear_z[history_slot] = sz;

  const double fnx = -fn * nx;
  const double fny = -fn * ny;
  const double fnz = -fn * nz;
  const double fix = fnx + ftx;
  const double fiy = fny + fty;
  const double fiz = fnz + ftz;

  atomicAdd(&fx[i], fix);
  atomicAdd(&fy[i], fiy);
  atomicAdd(&fz[i], fiz);
  atomicAdd(&fx[j], -fix);
  atomicAdd(&fy[j], -fiy);
  atomicAdd(&fz[j], -fiz);

  double torx, tory, torz;
  cross3(nx, ny, nz, ftx, fty, ftz, torx, tory, torz);
  atomicAdd(&torque_x[i], -radius[i] * torx);
  atomicAdd(&torque_y[i], -radius[i] * tory);
  atomicAdd(&torque_z[i], -radius[i] * torz);
	  atomicAdd(&torque_x[j], -radius[j] * torx);
	  atomicAdd(&torque_y[j], -radius[j] * tory);
	  atomicAdd(&torque_z[j], -radius[j] * torz);

  add_cdt_rolling_torque(i, j, nx, ny, nz, fn, reff, omega_x, omega_y,
                         omega_z, rolling_params, torque_x, torque_y,
                         torque_z);
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

void compute_hertz_normal_forces(GpuParticleData &particles,
                                 const GpuNeighborList &neighbors,
                                 const HertzNormalParams &params,
                                 cudaStream_t stream)
{
  if (!(params.effective_youngs_modulus >= 0.0))
    throw std::runtime_error("GPU_DEM Hertz normal effective Young's modulus must be nonnegative");
  if (neighbors.overflow(stream))
    throw std::runtime_error("GPU_DEM Hertz contact-force evaluation refused overflowed pair list");

  const int pairs = neighbors.pair_count(stream);
  if (pairs <= 0) return;
  if (static_cast<std::size_t>(pairs) > neighbors.pair_capacity())
    throw std::runtime_error("GPU_DEM Hertz contact-force pair count exceeds pair capacity");

  const int blocks = (pairs + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  hertz_normal_force_kernel<<<blocks, THREADS_PER_BLOCK, 0, stream>>>(
      pairs, neighbors.pair_first(), neighbors.pair_second(), particles.position_x(),
      particles.position_y(), particles.position_z(), particles.velocity_x(),
      particles.velocity_y(), particles.velocity_z(), particles.radius(), particles.mass(),
      params, particles.force_x(), particles.force_y(), particles.force_z());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
}

void compute_hertz_history_forces(GpuParticleData &particles,
                                  const GpuNeighborList &neighbors,
                                  GpuContactHistory &history,
                                  const HertzNormalParams &normal_params,
                                  const TangentialHistoryParams &tangential_params,
                                  cudaStream_t stream)
{
  RollingFrictionParams rolling_params;
  compute_hertz_history_rolling_forces(particles, neighbors, history,
                                       normal_params, tangential_params,
                                       rolling_params, stream);
}

void compute_hertz_history_rolling_forces(
    GpuParticleData &particles, const GpuNeighborList &neighbors,
    GpuContactHistory &history, const HertzNormalParams &normal_params,
    const TangentialHistoryParams &tangential_params,
    const RollingFrictionParams &rolling_params, cudaStream_t stream)
{
  if (!(normal_params.effective_youngs_modulus >= 0.0))
    throw std::runtime_error("GPU_DEM Hertz-history effective Young's modulus must be nonnegative");
  if (!(tangential_params.effective_shear_modulus >= 0.0) ||
      !(tangential_params.friction_coefficient >= 0.0) ||
      !(tangential_params.dt > 0.0))
    throw std::runtime_error("GPU_DEM tangential-history parameters are invalid");
  if (!(rolling_params.coefficient >= 0.0))
    throw std::runtime_error("GPU_DEM rolling-friction coefficient must be nonnegative");
  if (neighbors.overflow(stream))
    throw std::runtime_error("GPU_DEM Hertz-history evaluation refused overflowed pair list");

  const int pairs = neighbors.pair_count(stream);
  if (pairs <= 0) return;
  if (history.capacity() == 0)
    throw std::runtime_error("GPU_DEM Hertz-history requires allocated contact history");

  const int blocks = (pairs + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  hertz_history_force_kernel<<<blocks, THREADS_PER_BLOCK, 0, stream>>>(
      pairs, neighbors.pair_first(), neighbors.pair_second(), particles.position_x(),
      particles.position_y(), particles.position_z(), particles.velocity_x(),
      particles.velocity_y(), particles.velocity_z(), particles.omega_x(),
	      particles.omega_y(), particles.omega_z(), particles.radius(), particles.mass(),
	      particles.tag(), normal_params, tangential_params, rolling_params,
	      history.capacity(),
	      history.keys(), history.shear_x(), history.shear_y(), history.shear_z(),
      history.active(), history.count(), history.overflow_flag(), particles.force_x(),
      particles.force_y(), particles.force_z(), particles.torque_x(),
      particles.torque_y(), particles.torque_z());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
  if (history.overflow(stream))
    throw std::runtime_error("GPU_DEM Hertz-history contact history capacity exceeded");
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
