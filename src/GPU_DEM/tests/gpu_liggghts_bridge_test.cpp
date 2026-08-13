#include "gpu_dem_context.h"
#include "gpu_liggghts_bridge.h"
#include "gpu_particle_data.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using LAMMPS_NS::GPU_DEM::CpuParticleAoSView;
using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuParticleData;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::PrecisionMode;
using LAMMPS_NS::GPU_DEM::gather_aos_to_host;
using LAMMPS_NS::GPU_DEM::scatter_dynamic_state_to_aos;

namespace {

struct OwnedAoS {
  explicit OwnedAoS(int n) : nlocal(n)
  {
    allocate_vector3(x_storage, x);
    allocate_vector3(v_storage, v);
    allocate_vector3(f_storage, f);
    allocate_vector3(omega_storage, omega);
    allocate_vector3(torque_storage, torque);

    radius.resize(nlocal);
    rmass.resize(nlocal);
    density.resize(nlocal);
    type.resize(nlocal);
    mask.resize(nlocal);
    tag.resize(nlocal);
    image.resize(nlocal);
  }

  void fill()
  {
    for (int i = 0; i < nlocal; ++i) {
      const double base = static_cast<double>(i + 1);
      x[i][0] = 0.10 * base;
      x[i][1] = 0.20 * base;
      x[i][2] = 0.30 * base;
      v[i][0] = -0.10 * base;
      v[i][1] = -0.20 * base;
      v[i][2] = -0.30 * base;
      f[i][0] = 1.00 * base;
      f[i][1] = 2.00 * base;
      f[i][2] = 3.00 * base;
      omega[i][0] = 4.00 * base;
      omega[i][1] = 5.00 * base;
      omega[i][2] = 6.00 * base;
      torque[i][0] = 7.00 * base;
      torque[i][1] = 8.00 * base;
      torque[i][2] = 9.00 * base;
      radius[i] = 0.01 * base;
      rmass[i] = 0.05 * base;
      density[i] = 2500.0 + base;
      type[i] = 1 + (i % 2);
      mask[i] = 1 << (i % 3);
      tag[i] = 100 + i;
      image[i] = 200 + i;
    }
  }

  CpuParticleAoSView view()
  {
    CpuParticleAoSView result;
    result.nlocal = static_cast<std::size_t>(nlocal);
    result.x = x.data();
    result.v = v.data();
    result.f = f.data();
    result.omega = omega.data();
    result.torque = torque.data();
    result.radius = radius.data();
    result.rmass = rmass.data();
    result.density = density.data();
    result.type = type.data();
    result.mask = mask.data();
    result.tag = tag.data();
    result.image = image.data();
    return result;
  }

  int nlocal;
  std::vector<double> x_storage, v_storage, f_storage;
  std::vector<double> omega_storage, torque_storage;
  std::vector<double *> x, v, f, omega, torque;
  std::vector<double> radius, rmass, density;
  std::vector<int> type, mask;
  std::vector<LAMMPS_NS::tagint> tag;
  std::vector<LAMMPS_NS::imageint> image;

 private:
  void allocate_vector3(std::vector<double> &storage,
                        std::vector<double *> &ptrs)
  {
    storage.resize(static_cast<std::size_t>(nlocal) * 3);
    ptrs.resize(nlocal);
    for (int i = 0; i < nlocal; ++i)
      ptrs[i] = &storage[static_cast<std::size_t>(i) * 3];
  }
};

void expect_close(double actual, double expected, const char *name)
{
  if (std::fabs(actual - expected) > 1.0e-14)
    throw std::runtime_error(std::string(name) + " mismatch");
}

void check_host_from_aos(const OwnedAoS &aos, const HostParticleData &host)
{
  for (int i = 0; i < aos.nlocal; ++i) {
    expect_close(host.position_x[i], aos.x[i][0], "x0");
    expect_close(host.position_y[i], aos.x[i][1], "x1");
    expect_close(host.position_z[i], aos.x[i][2], "x2");
    expect_close(host.velocity_x[i], aos.v[i][0], "v0");
    expect_close(host.velocity_y[i], aos.v[i][1], "v1");
    expect_close(host.velocity_z[i], aos.v[i][2], "v2");
    expect_close(host.force_x[i], aos.f[i][0], "f0");
    expect_close(host.force_y[i], aos.f[i][1], "f1");
    expect_close(host.force_z[i], aos.f[i][2], "f2");
    expect_close(host.omega_x[i], aos.omega[i][0], "omega0");
    expect_close(host.omega_y[i], aos.omega[i][1], "omega1");
    expect_close(host.omega_z[i], aos.omega[i][2], "omega2");
    expect_close(host.torque_x[i], aos.torque[i][0], "torque0");
    expect_close(host.torque_y[i], aos.torque[i][1], "torque1");
    expect_close(host.torque_z[i], aos.torque[i][2], "torque2");
  }
}

} // namespace

int main()
{
  try {
    GpuDemContext context(0, PrecisionMode::Mixed);

    OwnedAoS aos(8);
    aos.fill();

    HostParticleData host_in;
    gather_aos_to_host(aos.view(), host_in);
    check_host_from_aos(aos, host_in);

    GpuParticleData device_particles;
    device_particles.copy_from_host(host_in, context.stream());
    context.synchronize();

    HostParticleData host_out;
    device_particles.copy_to_host(host_out, context.stream());
    context.synchronize();

    for (std::size_t i = 0; i < host_out.size(); ++i) {
      host_out.position_x[i] += 1.0;
      host_out.velocity_y[i] -= 2.0;
      host_out.force_z[i] += 3.0;
      host_out.omega_x[i] -= 4.0;
      host_out.torque_y[i] += 5.0;
    }

    const double original_radius = aos.radius[3];
    const int original_type = aos.type[3];
    scatter_dynamic_state_to_aos(host_out, aos.view());

    expect_close(aos.x[3][0], host_out.position_x[3], "scatter x");
    expect_close(aos.v[3][1], host_out.velocity_y[3], "scatter v");
    expect_close(aos.f[3][2], host_out.force_z[3], "scatter f");
    expect_close(aos.omega[3][0], host_out.omega_x[3], "scatter omega");
    expect_close(aos.torque[3][1], host_out.torque_y[3], "scatter torque");
    expect_close(aos.radius[3], original_radius, "radius unchanged");
    if (aos.type[3] != original_type)
      throw std::runtime_error("type should not be scattered back");

    std::cout << "GPU_DEM LIGGGHTS AoS bridge test passed on device "
              << context.device_id() << " ("
              << context.device_properties().name << ")\n";
  } catch (const std::exception &error) {
    std::cerr << "GPU_DEM LIGGGHTS AoS bridge test failed: "
              << error.what() << "\n";
    return 1;
  }

  return 0;
}
