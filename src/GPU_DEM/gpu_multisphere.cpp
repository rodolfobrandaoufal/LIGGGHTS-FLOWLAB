/* ----------------------------------------------------------------------
   GPU multisphere template scaffold for GPU_DEM.
------------------------------------------------------------------------- */

#include "gpu_multisphere.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace LAMMPS_NS {
namespace GPU_DEM {
namespace {

constexpr double PI = 3.141592653589793238462643383279502884;

double sphere_volume(double radius)
{
  return (4.0 / 3.0) * PI * radius * radius * radius;
}

double sphere_mass(double radius, double density)
{
  return density * sphere_volume(radius);
}

} // namespace

HostMultisphereTemplate load_multisphere_template(const std::string &path,
                                                  double density,
                                                  double scale)
{
  if (!(density > 0.0)) throw std::runtime_error("GPU_DEM multisphere density must be positive");
  if (!(scale > 0.0)) throw std::runtime_error("GPU_DEM multisphere scale must be positive");

  std::ifstream input(path.c_str());
  if (!input) throw std::runtime_error("GPU_DEM could not open multisphere template: " + path);

  HostMultisphereTemplate result;
  result.density = density;

  std::string line;
  while (std::getline(input, line)) {
    std::istringstream line_stream(line);
    double x, y, z, radius;
    if (!(line_stream >> x >> y >> z >> radius)) continue;
    if (!(radius > 0.0))
      throw std::runtime_error("GPU_DEM multisphere template has nonpositive radius: " + path);

    MultisphereChildSphere sphere;
    sphere.offset[0] = x * scale;
    sphere.offset[1] = y * scale;
    sphere.offset[2] = z * scale;
    sphere.radius = radius * scale;
    result.spheres.push_back(sphere);
  }

  if (result.spheres.empty())
    throw std::runtime_error("GPU_DEM multisphere template contained no spheres: " + path);

  double weighted_center[3]{0.0, 0.0, 0.0};
  for (const MultisphereChildSphere &sphere : result.spheres) {
    const double mass = sphere_mass(sphere.radius, density);
    result.mass += mass;
    weighted_center[0] += mass * sphere.offset[0];
    weighted_center[1] += mass * sphere.offset[1];
    weighted_center[2] += mass * sphere.offset[2];
  }

  for (int d = 0; d < 3; ++d)
    result.center_of_mass[d] = weighted_center[d] / result.mass;

  for (const MultisphereChildSphere &sphere : result.spheres) {
    const double mass = sphere_mass(sphere.radius, density);
    const double dx = sphere.offset[0] - result.center_of_mass[0];
    const double dy = sphere.offset[1] - result.center_of_mass[1];
    const double dz = sphere.offset[2] - result.center_of_mass[2];
    const double own_inertia = 0.4 * mass * sphere.radius * sphere.radius;

    result.inertia[0] += own_inertia + mass * (dy * dy + dz * dz);
    result.inertia[1] += own_inertia + mass * (dx * dx + dz * dz);
    result.inertia[2] += own_inertia + mass * (dx * dx + dy * dy);

    const double distance = std::sqrt(dx * dx + dy * dy + dz * dz) + sphere.radius;
    if (distance > result.bounding_radius) result.bounding_radius = distance;
  }

  return result;
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
