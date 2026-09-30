/* ----------------------------------------------------------------------
   CUDA error checking helpers for GPU_DEM.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_ERROR_CHECK_H
#define LMP_GPU_DEM_ERROR_CHECK_H

#include <cuda_runtime.h>

#include <stdexcept>
#include <sstream>
#include <string>

namespace LAMMPS_NS {
namespace GPU_DEM {

inline void throw_cuda_error(cudaError_t error, const char *call,
                             const char *file, int line)
{
  if (error == cudaSuccess) return;

  std::ostringstream message;
  message << "CUDA call failed: " << call
          << " at " << file << ":" << line
          << " error=" << cudaGetErrorName(error)
          << " (" << cudaGetErrorString(error) << ")";
  throw std::runtime_error(message.str());
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#define GPU_DEM_CUDA_CHECK(call) \
  ::LAMMPS_NS::GPU_DEM::throw_cuda_error((call), #call, __FILE__, __LINE__)

#endif
