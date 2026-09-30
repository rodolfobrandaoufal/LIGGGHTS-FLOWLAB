#include "gpu_dem_context.h"
#include "gpu_particle_data.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuParticleData;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::PrecisionMode;

namespace {

void fill_particles(HostParticleData &particles, int n)
{
  particles.resize(n);
  for (int i = 0; i < n; ++i) {
    const double base = static_cast<double>(i + 1);
    particles.position_x[i] = 0.1 * base;
    particles.position_y[i] = 0.2 * base;
    particles.position_z[i] = 0.3 * base;
    particles.velocity_x[i] = -0.1 * base;
    particles.velocity_y[i] = -0.2 * base;
    particles.velocity_z[i] = -0.3 * base;
    particles.omega_x[i] = 1.0 * base;
    particles.omega_y[i] = 2.0 * base;
    particles.omega_z[i] = 3.0 * base;
    particles.force_x[i] = 4.0 * base;
    particles.force_y[i] = 5.0 * base;
    particles.force_z[i] = 6.0 * base;
    particles.torque_x[i] = 7.0 * base;
    particles.torque_y[i] = 8.0 * base;
    particles.torque_z[i] = 9.0 * base;
    particles.radius[i] = 0.001 * base;
    particles.mass[i] = 0.01 * base;
    particles.density[i] = 2500.0 + base;
    particles.type[i] = 1 + (i % 3);
    particles.mask[i] = (i % 2) ? 1 : 3;
    particles.tag[i] = 1000 + i;
    particles.image_flags[i] = i;
  }
}

template <class T>
void expect_equal(const std::vector<T> &a, const std::vector<T> &b,
                  const char *name)
{
  if (a.size() != b.size())
    throw std::runtime_error(std::string(name) + " size mismatch");
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i] != b[i])
      throw std::runtime_error(std::string(name) + " value mismatch");
  }
}

void compare_particles(const HostParticleData &expected,
                       const HostParticleData &actual)
{
  expect_equal(expected.position_x, actual.position_x, "position_x");
  expect_equal(expected.position_y, actual.position_y, "position_y");
  expect_equal(expected.position_z, actual.position_z, "position_z");
  expect_equal(expected.velocity_x, actual.velocity_x, "velocity_x");
  expect_equal(expected.velocity_y, actual.velocity_y, "velocity_y");
  expect_equal(expected.velocity_z, actual.velocity_z, "velocity_z");
  expect_equal(expected.omega_x, actual.omega_x, "omega_x");
  expect_equal(expected.omega_y, actual.omega_y, "omega_y");
  expect_equal(expected.omega_z, actual.omega_z, "omega_z");
  expect_equal(expected.force_x, actual.force_x, "force_x");
  expect_equal(expected.force_y, actual.force_y, "force_y");
  expect_equal(expected.force_z, actual.force_z, "force_z");
  expect_equal(expected.torque_x, actual.torque_x, "torque_x");
  expect_equal(expected.torque_y, actual.torque_y, "torque_y");
  expect_equal(expected.torque_z, actual.torque_z, "torque_z");
  expect_equal(expected.radius, actual.radius, "radius");
  expect_equal(expected.mass, actual.mass, "mass");
  expect_equal(expected.density, actual.density, "density");
  expect_equal(expected.type, actual.type, "type");
  expect_equal(expected.mask, actual.mask, "mask");
  expect_equal(expected.tag, actual.tag, "tag");
  expect_equal(expected.image_flags, actual.image_flags, "image_flags");
}

} // namespace

int main()
{
  try {
    GpuDemContext context(0, PrecisionMode::Mixed);

    HostParticleData host_in;
    fill_particles(host_in, 16);

    GpuParticleData device_particles;
    device_particles.copy_from_host(host_in, context.stream());
    context.synchronize();

    HostParticleData host_out;
    device_particles.copy_to_host(host_out, context.stream());
    context.synchronize();

    compare_particles(host_in, host_out);

    std::cout << "GPU_DEM particle data transfer test passed on device "
              << context.device_id() << " ("
              << context.device_properties().name << ")\n";
  } catch (const std::exception &error) {
    std::cerr << "GPU_DEM particle data transfer test failed: "
              << error.what() << "\n";
    return 1;
  }

  return 0;
}
