/* ----------------------------------------------------------------------
   GPU_DEM precision policy for the GPU-native LIGGGHTS DEM path.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_PRECISION_H
#define LMP_GPU_DEM_PRECISION_H

namespace LAMMPS_NS {
namespace GPU_DEM {

enum class PrecisionMode {
  Double,
  Mixed,
  Single
};

inline const char *precision_mode_name(PrecisionMode mode)
{
  switch (mode) {
    case PrecisionMode::Double: return "double";
    case PrecisionMode::Mixed: return "mixed";
    case PrecisionMode::Single: return "single";
  }
  return "unknown";
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
