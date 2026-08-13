/* ----------------------------------------------------------------------
   Hookean normal contact model for GPU_DEM device kernels.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_MODELS_NORMAL_HOOKE_CUH
#define LMP_GPU_DEM_MODELS_NORMAL_HOOKE_CUH

#include "../gpu_contact_pipeline.h"

namespace LAMMPS_NS {
namespace GPU_DEM {
namespace Models {

__device__ inline double clamp_nonnegative(double value)
{
  return value > 0.0 ? value : 0.0;
}

__device__ inline double hooke_normal_force(double overlap,
                                            double normal_relative_velocity,
                                            const HookeNormalParams &params)
{
  return clamp_nonnegative(params.normal_stiffness * overlap -
                           params.normal_damping * normal_relative_velocity);
}

} // namespace Models
} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
