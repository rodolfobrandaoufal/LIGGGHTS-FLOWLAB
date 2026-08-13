/* ----------------------------------------------------------------------
   GPU_DEM execution-policy regression test.
------------------------------------------------------------------------- */

#include "gpu_execution_policy.h"

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>

using LAMMPS_NS::GPU_DEM::GpuContactModel;
using LAMMPS_NS::GPU_DEM::GpuExecutionPolicy;
using LAMMPS_NS::GPU_DEM::GpuMode;
using LAMMPS_NS::GPU_DEM::GpuRequestedFeatures;
using LAMMPS_NS::GPU_DEM::GpuWallModel;
using LAMMPS_NS::GPU_DEM::require_gpu_features_or_throw;
using LAMMPS_NS::GPU_DEM::validate_gpu_features;

namespace {

void fail(const char *message)
{
  std::fprintf(stderr, "%s\n", message);
  std::exit(1);
}

void expect_contains(const std::string &text, const char *needle)
{
  if (text.find(needle) == std::string::npos) {
    std::fprintf(stderr, "expected message to contain '%s', got '%s'\n",
                 needle, text.c_str());
    std::exit(1);
  }
}

} // namespace

int main()
{
  GpuRequestedFeatures supported;
  supported.contact_model = GpuContactModel::HookeNormal;
  supported.wall_model = GpuWallModel::Plane;

  GpuExecutionPolicy off;
  off.mode = GpuMode::Off;
  const auto off_validation = validate_gpu_features(off, supported);
  if (off_validation.gpu_allowed || !off_validation.fallback_required)
    fail("gpu_mode off should force CPU fallback");

  GpuExecutionPolicy strict;
  strict.mode = GpuMode::Strict;
  const auto supported_validation = validate_gpu_features(strict, supported);
  if (!supported_validation.gpu_allowed || supported_validation.fallback_required)
    fail("strict mode should accept the initial supported GPU feature set");
  require_gpu_features_or_throw(strict, supported);

  GpuRequestedFeatures unsupported = supported;
  unsupported.requires_cohesion = true;
  unsupported.requires_thermal = true;
  unsupported.wall_model = GpuWallModel::Mesh;
  unsupported.contact_model = GpuContactModel::HertzMindlin;

  const auto unsupported_validation = validate_gpu_features(strict, unsupported);
  if (unsupported_validation.gpu_allowed || !unsupported_validation.fallback_required)
    fail("strict mode should reject unsupported GPU feature set");
  if (unsupported_validation.unsupported_reasons.size() != 2)
    fail("unsupported feature set should report all requested blockers");

  try {
    require_gpu_features_or_throw(strict, unsupported);
    fail("strict mode did not throw for unsupported feature set");
  } catch (const std::runtime_error &error) {
    const std::string message(error.what());
    expect_contains(message, "cohesion");
    expect_contains(message, "thermal contact");
  }

  GpuExecutionPolicy automatic;
  automatic.mode = GpuMode::Auto;
  const auto auto_validation = validate_gpu_features(automatic, unsupported);
  if (auto_validation.gpu_allowed || !auto_validation.fallback_required)
    fail("auto mode should mark unsupported features for explicit fallback");
  require_gpu_features_or_throw(automatic, unsupported);

  std::printf("GPU_DEM execution-policy test passed\n");
  return 0;
}
