#ifdef FIX_CLASS

FixStyle(GPU,FixGpuDemPackage)

#else

#ifndef LMP_FIX_GPU_DEM_PACKAGE_H
#define LMP_FIX_GPU_DEM_PACKAGE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGpuDemPackage : public Fix {
 public:
  FixGpuDemPackage(class LAMMPS *, int, char **);
  int setmask();
};

}

#endif
#endif
