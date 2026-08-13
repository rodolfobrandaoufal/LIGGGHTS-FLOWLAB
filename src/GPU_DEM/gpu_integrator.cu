#include "gpu_integrator.h"

#include "gpu_error_check.h"

namespace LAMMPS_NS {
namespace GPU_DEM {
namespace {

constexpr double SPHERE_INERTIA = 0.4;
constexpr int THREADS_PER_BLOCK = 256;

__device__ double rotational_dtf(double dtf, int dimension)
{
  return dimension == 2 ? dtf / 0.5 : dtf / SPHERE_INERTIA;
}

__global__ void reset_forces_kernel(int n,
                                    double *force_x, double *force_y, double *force_z,
                                    double *torque_x, double *torque_y, double *torque_z)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  force_x[i] = 0.0;
  force_y[i] = 0.0;
  force_z[i] = 0.0;
  torque_x[i] = 0.0;
  torque_y[i] = 0.0;
  torque_z[i] = 0.0;
}

__global__ void nve_sphere_initial_kernel(int n,
                                          double dtv,
                                          double dtf,
                                          int groupbit,
                                          int dimension,
                                          double one_plus_added_mass,
                                          const int *mask,
                                          double *position_x, double *position_y, double *position_z,
                                          double *velocity_x, double *velocity_y, double *velocity_z,
                                          const double *force_x, const double *force_y, const double *force_z,
                                          double *omega_x, double *omega_y, double *omega_z,
                                          const double *torque_x, const double *torque_y, const double *torque_z,
                                          const double *radius,
                                          const double *mass)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n || !(mask[i] & groupbit)) return;

  const double dtfm = dtf / (mass[i] * one_plus_added_mass);
  velocity_x[i] += dtfm * force_x[i];
  velocity_y[i] += dtfm * force_y[i];
  velocity_z[i] += dtfm * force_z[i];

  position_x[i] += dtv * velocity_x[i];
  position_y[i] += dtv * velocity_y[i];
  position_z[i] += dtv * velocity_z[i];

  const double dtirotate = rotational_dtf(dtf, dimension) /
                           (radius[i] * radius[i] * mass[i]);
  omega_x[i] += dtirotate * torque_x[i];
  omega_y[i] += dtirotate * torque_y[i];
  omega_z[i] += dtirotate * torque_z[i];
}

__global__ void nve_sphere_final_kernel(int n,
                                        double dtf,
                                        int groupbit,
                                        int dimension,
                                        double one_plus_added_mass,
                                        const int *mask,
                                        double *velocity_x, double *velocity_y, double *velocity_z,
                                        const double *force_x, const double *force_y, const double *force_z,
                                        double *omega_x, double *omega_y, double *omega_z,
                                        const double *torque_x, const double *torque_y, const double *torque_z,
                                        const double *radius,
                                        const double *mass)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n || !(mask[i] & groupbit)) return;

  const double dtfm = dtf / (mass[i] * one_plus_added_mass);
  velocity_x[i] += dtfm * force_x[i];
  velocity_y[i] += dtfm * force_y[i];
  velocity_z[i] += dtfm * force_z[i];

  const double dtirotate = rotational_dtf(dtf, dimension) /
                           (radius[i] * radius[i] * mass[i]);
  omega_x[i] += dtirotate * torque_x[i];
  omega_y[i] += dtirotate * torque_y[i];
  omega_z[i] += dtirotate * torque_z[i];
}

int block_count(std::size_t n)
{
  return static_cast<int>((n + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
}

} // namespace

void reset_forces(GpuParticleData &particles, cudaStream_t stream)
{
  if (particles.size() == 0) return;
  reset_forces_kernel<<<block_count(particles.size()), THREADS_PER_BLOCK, 0, stream>>>(
      static_cast<int>(particles.size()),
      particles.force_x(), particles.force_y(), particles.force_z(),
      particles.torque_x(), particles.torque_y(), particles.torque_z());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
}

void nve_sphere_initial_integrate(GpuParticleData &particles,
                                  double dtv,
                                  double dtf,
                                  int groupbit,
                                  int dimension,
                                  double one_plus_added_mass,
                                  cudaStream_t stream)
{
  if (particles.size() == 0) return;
  nve_sphere_initial_kernel<<<block_count(particles.size()), THREADS_PER_BLOCK, 0, stream>>>(
      static_cast<int>(particles.size()), dtv, dtf, groupbit, dimension,
      one_plus_added_mass, particles.mask(),
      particles.position_x(), particles.position_y(), particles.position_z(),
      particles.velocity_x(), particles.velocity_y(), particles.velocity_z(),
      particles.force_x(), particles.force_y(), particles.force_z(),
      particles.omega_x(), particles.omega_y(), particles.omega_z(),
      particles.torque_x(), particles.torque_y(), particles.torque_z(),
      particles.radius(), particles.mass());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
}

void nve_sphere_final_integrate(GpuParticleData &particles,
                                double dtf,
                                int groupbit,
                                int dimension,
                                double one_plus_added_mass,
                                cudaStream_t stream)
{
  if (particles.size() == 0) return;
  nve_sphere_final_kernel<<<block_count(particles.size()), THREADS_PER_BLOCK, 0, stream>>>(
      static_cast<int>(particles.size()), dtf, groupbit, dimension,
      one_plus_added_mass, particles.mask(),
      particles.velocity_x(), particles.velocity_y(), particles.velocity_z(),
      particles.force_x(), particles.force_y(), particles.force_z(),
      particles.omega_x(), particles.omega_y(), particles.omega_z(),
      particles.torque_x(), particles.torque_y(), particles.torque_z(),
      particles.radius(), particles.mass());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
