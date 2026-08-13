/* ----------------------------------------------------------------------
   Production LIGGGHTS runtime probe for the GPU_DEM CUDA backend.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_RUNTIME_LIGGGHTS_H
#define LMP_GPU_RUNTIME_LIGGGHTS_H

namespace LAMMPS_NS {
namespace GPU_DEM {

struct LiggghtsGpuRuntimeProbe {
  int device_id;
  int device_count;
  char device_name[256];
};

void probe_liggghts_gpu_runtime(int requested_device,
                                int precision_policy,
                                LiggghtsGpuRuntimeProbe &result);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
