/* ----------------------------------------------------------------------
   Standalone GPU_DEM multisphere template loader regression test.
------------------------------------------------------------------------- */

#include "gpu_multisphere.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

using LAMMPS_NS::GPU_DEM::HostMultisphereTemplate;
using LAMMPS_NS::GPU_DEM::load_multisphere_template;

namespace {

void expect_near(double actual, double expected, double tol, const char *label)
{
  if (std::fabs(actual - expected) <= tol) return;
  std::fprintf(stderr, "%s mismatch: actual %.17g expected %.17g tolerance %.3g\n",
               label, actual, expected, tol);
  std::exit(1);
}

} // namespace

int main(int argc, char **argv)
{
  if (argc != 5) {
    std::fprintf(stderr, "usage: %s <template> <density> <scale> <expected-spheres>\n",
                 argv[0]);
    return 2;
  }

  const std::string path = argv[1];
  const double density = std::atof(argv[2]);
  const double scale = std::atof(argv[3]);
  const int expected_spheres = std::atoi(argv[4]);

  HostMultisphereTemplate templ = load_multisphere_template(path, density, scale);
  if (static_cast<int>(templ.spheres.size()) != expected_spheres) {
    std::fprintf(stderr, "sphere count mismatch: actual %zu expected %d\n",
                 templ.spheres.size(), expected_spheres);
    return 1;
  }

  if (!(templ.mass > 0.0 && templ.bounding_radius > 0.0 &&
        templ.inertia[0] > 0.0 && templ.inertia[1] > 0.0 && templ.inertia[2] > 0.0)) {
    std::fprintf(stderr, "invalid multisphere mass properties\n");
    return 1;
  }

  if (expected_spheres == 3) {
    expect_near(templ.spheres[0].radius, 0.00194 * scale, 1.0e-15,
                "sphere 0 radius");
    expect_near(templ.spheres[1].offset[0], -0.0034 * scale, 1.0e-15,
                "sphere 1 offset x");
    expect_near(templ.spheres[2].offset[0], 0.0034 * scale, 1.0e-15,
                "sphere 2 offset x");
  }

  std::printf("GPU_DEM multisphere template test passed: spheres=%zu mass=%.17g "
              "bounding_radius=%.17g path=%s\n",
              templ.spheres.size(), templ.mass, templ.bounding_radius, path.c_str());
  return 0;
}
