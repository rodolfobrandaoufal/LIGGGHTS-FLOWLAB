#include "gpu_runtime_facade.h"

#include <sstream>

namespace LAMMPS_NS {
namespace GPU_DEM {

std::string format_gpu_runtime_message(const GpuExecutionPolicy &policy,
                                       const GpuFeatureValidation &validation)
{
  std::ostringstream out;
  out << "gpu_mode=" << gpu_mode_name(policy.mode);

  if (validation.gpu_allowed) {
    out << " GPU_DEM active";
    return out.str();
  }

  out << " CPU fallback";
  for (const std::string &reason : validation.unsupported_reasons)
    out << "; " << reason;
  return out.str();
}

const GpuRuntimeStatus &GpuRuntimeFacade::configure(
    const GpuExecutionPolicy &policy,
    const LiggghtsFeatureDescriptor &descriptor,
    const GpuDemDriverConfig &driver_config,
    bool create_driver)
{
  reset();

  status_.requested_features = build_gpu_requested_features(descriptor);
  require_gpu_features_or_throw(policy, status_.requested_features);

  status_.validation = validate_gpu_features(policy, status_.requested_features);
  status_.message = format_gpu_runtime_message(policy, status_.validation);

  if (!status_.validation.gpu_allowed) {
    status_.gpu_active = false;
    status_.cpu_fallback = true;
    return status_;
  }

  status_.gpu_active = true;
  status_.cpu_fallback = false;
  if (create_driver)
    driver_.reset(new GpuDemDriver(driver_config));

  return status_;
}

void GpuRuntimeFacade::reset()
{
  driver_.reset();
  status_ = GpuRuntimeStatus();
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
