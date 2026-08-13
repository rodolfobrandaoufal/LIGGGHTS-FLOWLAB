#include "gpu_runtime_liggghts.h"

#include "gpu_dem_context.h"

#include <cstring>
#include <stdexcept>

namespace LAMMPS_NS {
namespace GPU_DEM {

namespace {

PrecisionMode precision_from_policy(int precision_policy)
{
  if (precision_policy == 0) return PrecisionMode::Double;
  if (precision_policy == 2) return PrecisionMode::Single;
  return PrecisionMode::Mixed;
}

} // namespace

void probe_liggghts_gpu_runtime(int requested_device,
                                int precision_policy,
                                LiggghtsGpuRuntimeProbe &result)
{
  const int count = gpu_device_count();
  if (count <= 0)
    throw std::runtime_error("GPU_DEM found no CUDA-capable devices");

  const int device_id = requested_device < 0 ? 0 : requested_device;
  if (device_id >= count)
    throw std::runtime_error("GPU_DEM requested CUDA device is out of range");

  GpuDemContext context(device_id, precision_from_policy(precision_policy));
  context.synchronize();

  result.device_id = context.device_id();
  result.device_count = count;
  std::strncpy(result.device_name, context.device_properties().name,
               sizeof(result.device_name) - 1);
  result.device_name[sizeof(result.device_name) - 1] = '\0';
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
