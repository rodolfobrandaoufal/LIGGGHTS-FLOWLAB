/* ----------------------------------------------------------------------
   Hertz normal contact model for GPU_DEM device kernels.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_MODELS_NORMAL_HERTZ_CUH
#define LMP_GPU_DEM_MODELS_NORMAL_HERTZ_CUH

#include "../gpu_contact_pipeline.h"

#include <math_constants.h>

namespace LAMMPS_NS {
namespace GPU_DEM {
namespace Models {

__device__ inline double effective_radius(double radi, double radj)
{
  return (radi * radj) / (radi + radj);
}

__device__ inline double effective_mass(double mi, double mj)
{
  return (mi * mj) / (mi + mj);
}

__device__ inline double hertz_normal_force(double overlap,
                                            double normal_relative_velocity,
                                            double radi,
                                            double radj,
                                            double mi,
                                            double mj,
                                            const HertzNormalParams &params)
{
  const double reff = effective_radius(radi, radj);
  const double meff = effective_mass(mi, mj);
  const double sqrt_delta_reff = sqrt(reff * overlap);
  const double sn = 2.0 * params.effective_youngs_modulus * sqrt_delta_reff;
  const double kn = (4.0 / 3.0) * params.effective_youngs_modulus * sqrt_delta_reff;
  const double sqrt_five_over_six = 0.91287092917527685576161630466800355659;
  const double gamman = -2.0 * sqrt_five_over_six * params.beta_effective *
                        sqrt(sn * meff);

  double normal_force = kn * overlap - gamman * normal_relative_velocity;
  if (params.limit_force && normal_force < 0.0) normal_force = 0.0;
  return normal_force;
}

} // namespace Models
} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
