/* ----------------------------------------------------------------------
   GPU multisphere template scaffold for GPU_DEM.
------------------------------------------------------------------------- */

#ifndef LMP_GPU_DEM_MULTISPHERE_H
#define LMP_GPU_DEM_MULTISPHERE_H

#include <cstddef>
#include <string>
#include <vector>

namespace LAMMPS_NS {
namespace GPU_DEM {

struct MultisphereChildSphere {
  double offset[3]{0.0, 0.0, 0.0};
  double radius{0.0};
};

struct HostMultisphereTemplate {
  std::vector<MultisphereChildSphere> spheres;
  double density{0.0};
  double mass{0.0};
  double center_of_mass[3]{0.0, 0.0, 0.0};
  double inertia[3]{0.0, 0.0, 0.0};
  double bounding_radius{0.0};
};

HostMultisphereTemplate load_multisphere_template(const std::string &path,
                                                  double density,
                                                  double scale);

} // namespace GPU_DEM
} // namespace LAMMPS_NS

#endif
