#include "gpu_liggghts_bridge.h"

#include <stdexcept>
#include <string>

namespace LAMMPS_NS {
namespace GPU_DEM {

namespace {

void require_vector3(double **ptr, const char *name)
{
  if (!ptr)
    throw std::runtime_error(std::string("GPU_DEM AoS bridge missing ") + name);
}

double optional_scalar(const double *ptr, std::size_t i, double fallback)
{
  return ptr ? ptr[i] : fallback;
}

int optional_int(const int *ptr, std::size_t i, int fallback)
{
  return ptr ? ptr[i] : fallback;
}

long long optional_tag(const tagint *ptr, std::size_t i)
{
  return ptr ? static_cast<long long>(ptr[i]) : static_cast<long long>(i + 1);
}

long long optional_image(const imageint *ptr, std::size_t i)
{
  return ptr ? static_cast<long long>(ptr[i]) : 0LL;
}

} // namespace

void gather_aos_to_host(const CpuParticleAoSView &aos,
                        HostParticleData &host)
{
  require_vector3(aos.x, "x");
  require_vector3(aos.v, "v");
  require_vector3(aos.f, "f");

  host.resize(aos.nlocal);

  for (std::size_t i = 0; i < aos.nlocal; ++i) {
    host.position_x[i] = aos.x[i][0];
    host.position_y[i] = aos.x[i][1];
    host.position_z[i] = aos.x[i][2];

    host.velocity_x[i] = aos.v[i][0];
    host.velocity_y[i] = aos.v[i][1];
    host.velocity_z[i] = aos.v[i][2];

    host.force_x[i] = aos.f[i][0];
    host.force_y[i] = aos.f[i][1];
    host.force_z[i] = aos.f[i][2];

    host.omega_x[i] = aos.omega ? aos.omega[i][0] : 0.0;
    host.omega_y[i] = aos.omega ? aos.omega[i][1] : 0.0;
    host.omega_z[i] = aos.omega ? aos.omega[i][2] : 0.0;

    host.torque_x[i] = aos.torque ? aos.torque[i][0] : 0.0;
    host.torque_y[i] = aos.torque ? aos.torque[i][1] : 0.0;
    host.torque_z[i] = aos.torque ? aos.torque[i][2] : 0.0;

    host.radius[i] = optional_scalar(aos.radius, i, 0.0);
    host.mass[i] = optional_scalar(aos.rmass, i, 1.0);
    host.density[i] = optional_scalar(aos.density, i, 0.0);

    host.type[i] = optional_int(aos.type, i, 1);
    host.mask[i] = optional_int(aos.mask, i, 1);
    host.tag[i] = optional_tag(aos.tag, i);
    host.image_flags[i] = optional_image(aos.image, i);
  }
}

void scatter_dynamic_state_to_aos(const HostParticleData &host,
                                  const CpuParticleAoSView &aos)
{
  host.validate_sizes();
  if (host.size() != aos.nlocal)
    throw std::runtime_error("GPU_DEM AoS bridge scatter size mismatch");

  require_vector3(aos.x, "x");
  require_vector3(aos.v, "v");
  require_vector3(aos.f, "f");

  for (std::size_t i = 0; i < aos.nlocal; ++i) {
    aos.x[i][0] = host.position_x[i];
    aos.x[i][1] = host.position_y[i];
    aos.x[i][2] = host.position_z[i];

    aos.v[i][0] = host.velocity_x[i];
    aos.v[i][1] = host.velocity_y[i];
    aos.v[i][2] = host.velocity_z[i];

    aos.f[i][0] = host.force_x[i];
    aos.f[i][1] = host.force_y[i];
    aos.f[i][2] = host.force_z[i];

    if (aos.omega) {
      aos.omega[i][0] = host.omega_x[i];
      aos.omega[i][1] = host.omega_y[i];
      aos.omega[i][2] = host.omega_z[i];
    }

    if (aos.torque) {
      aos.torque[i][0] = host.torque_x[i];
      aos.torque[i][1] = host.torque_y[i];
      aos.torque[i][2] = host.torque_z[i];
    }
  }
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
