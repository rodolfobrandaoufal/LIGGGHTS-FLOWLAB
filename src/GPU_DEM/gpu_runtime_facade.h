/* ----------------------------------------------------------------------
   Experimental GPU_DEM runtime facade.

   This class is the narrow decision boundary that will later be called by
   LIGGGHTS. It validates requested features, reports explicit fallback
   reasons, and creates a persistent GPU driver only when the selected policy
   allows GPU execution.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_RUNTIME_FACADE_H
#define LMP_GPU_DEM_RUNTIME_FACADE_H

#include "gpu_dem_driver.h"
#include "gpu_execution_policy.h"
#include "gpu_liggghts_feature_probe.h"

#include <memory>
#include <string>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct GpuRuntimeStatus {
  bool gpu_active{false};
  bool cpu_fallback{false};
  std::string message;
  GpuFeatureValidation validation;
  GpuRequestedFeatures requested_features;
};

class GpuRuntimeFacade {
 public:
  GpuRuntimeFacade() = default;

  GpuRuntimeFacade(const GpuRuntimeFacade &) = delete;
  GpuRuntimeFacade &operator=(const GpuRuntimeFacade &) = delete;

  const GpuRuntimeStatus &status() const { return status_; }
  bool gpu_active() const { return status_.gpu_active; }
  bool cpu_fallback() const { return status_.cpu_fallback; }

  GpuDemDriver *driver() { return driver_.get(); }
  const GpuDemDriver *driver() const { return driver_.get(); }

  const GpuRuntimeStatus &configure(const GpuExecutionPolicy &policy,
                                    const LiggghtsFeatureDescriptor &descriptor,
                                    const GpuDemDriverConfig &driver_config,
                                    bool create_driver = true);

  void reset();

 private:
  GpuRuntimeStatus status_;
  std::unique_ptr<GpuDemDriver> driver_;
};

std::string format_gpu_runtime_message(const GpuExecutionPolicy &policy,
                                       const GpuFeatureValidation &validation);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
