/* ----------------------------------------------------------------------
   GPU kernels for the first no-contact DEM integration scaffold.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_INTEGRATOR_H
#define LMP_GPU_DEM_INTEGRATOR_H

#include "gpu_particle_data.h"

#include <cuda_runtime.h>

namespace LAMMPS_NS {
namespace GPU_DEM {

void reset_forces(GpuParticleData &particles, cudaStream_t stream = 0);

void nve_sphere_initial_integrate(GpuParticleData &particles,
                                  double dtv,
                                  double dtf,
                                  int groupbit,
                                  int dimension,
                                  double one_plus_added_mass,
                                  cudaStream_t stream = 0);

void nve_sphere_final_integrate(GpuParticleData &particles,
                                double dtf,
                                int groupbit,
                                int dimension,
                                double one_plus_added_mass,
                                cudaStream_t stream = 0);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
