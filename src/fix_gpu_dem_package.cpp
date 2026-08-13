#include "fix_gpu_dem_package.h"

#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixGpuDemPackage::FixGpuDemPackage(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp,narg,arg)
{
  if (narg < 3) error->all(FLERR,"Illegal package gpu command");

#ifndef LIGGGHTS_GPU_DEM_RUNTIME
  error->all(FLERR,
             "Package gpu requires the GPU_DEM-linked executable; build and run lmp_hdf5mpi_gpu");
#endif
}

int FixGpuDemPackage::setmask()
{
  return 0;
}
