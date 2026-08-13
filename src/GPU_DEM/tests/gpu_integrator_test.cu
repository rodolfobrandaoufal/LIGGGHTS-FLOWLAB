#include "gpu_dem_context.h"
#include "gpu_integrator.h"
#include "gpu_particle_data.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuParticleData;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::PrecisionMode;
using LAMMPS_NS::GPU_DEM::nve_sphere_final_integrate;
using LAMMPS_NS::GPU_DEM::nve_sphere_initial_integrate;
using LAMMPS_NS::GPU_DEM::reset_forces;

namespace {

constexpr double TOL = 1.0e-12;

void expect_near(double actual, double expected, const char *name)
{
  if (std::fabs(actual - expected) > TOL) {
    std::cerr << name << " actual=" << actual
              << " expected=" << expected << "\n";
    throw std::runtime_error("GPU integrator mismatch");
  }
}

HostParticleData make_particles()
{
  HostParticleData p;
  p.resize(2);

  for (int i = 0; i < 2; ++i) {
    p.position_x[i] = 1.0 + i;
    p.position_y[i] = 2.0 + i;
    p.position_z[i] = 3.0 + i;
    p.velocity_x[i] = 0.25 + i;
    p.velocity_y[i] = -0.5 - i;
    p.velocity_z[i] = 0.75 + i;
    p.omega_x[i] = 0.1 + i;
    p.omega_y[i] = 0.2 + i;
    p.omega_z[i] = 0.3 + i;
    p.force_x[i] = 2.0 + i;
    p.force_y[i] = -4.0 - i;
    p.force_z[i] = 6.0 + i;
    p.torque_x[i] = 0.5 + i;
    p.torque_y[i] = -0.25 - i;
    p.torque_z[i] = 0.125 + i;
    p.radius[i] = 0.5 + 0.1 * i;
    p.mass[i] = 2.0 + i;
    p.density[i] = 2500.0;
    p.type[i] = 1;
    p.mask[i] = i == 0 ? 1 : 0;
    p.tag[i] = 100 + i;
    p.image_flags[i] = 0;
  }

  return p;
}

void check_integrated_particle(const HostParticleData &before,
                               const HostParticleData &after)
{
  const double dt = 0.01;
  const double dtf = 0.5 * dt;
  const int i = 0;
  const double mass = before.mass[i];
  const double radius = before.radius[i];
  const double inertia_dtf = (dtf / 0.4) / (radius * radius * mass);

  const double vx_half = before.velocity_x[i] + dtf * before.force_x[i] / mass;
  const double vy_half = before.velocity_y[i] + dtf * before.force_y[i] / mass;
  const double vz_half = before.velocity_z[i] + dtf * before.force_z[i] / mass;

  expect_near(after.position_x[i], before.position_x[i] + dt * vx_half, "position_x");
  expect_near(after.position_y[i], before.position_y[i] + dt * vy_half, "position_y");
  expect_near(after.position_z[i], before.position_z[i] + dt * vz_half, "position_z");

  expect_near(after.velocity_x[i], before.velocity_x[i] + dt * before.force_x[i] / mass, "velocity_x");
  expect_near(after.velocity_y[i], before.velocity_y[i] + dt * before.force_y[i] / mass, "velocity_y");
  expect_near(after.velocity_z[i], before.velocity_z[i] + dt * before.force_z[i] / mass, "velocity_z");

  expect_near(after.omega_x[i], before.omega_x[i] + 2.0 * inertia_dtf * before.torque_x[i], "omega_x");
  expect_near(after.omega_y[i], before.omega_y[i] + 2.0 * inertia_dtf * before.torque_y[i], "omega_y");
  expect_near(after.omega_z[i], before.omega_z[i] + 2.0 * inertia_dtf * before.torque_z[i], "omega_z");
}

void check_ungrouped_particle_unchanged(const HostParticleData &before,
                                        const HostParticleData &after)
{
  const int i = 1;
  expect_near(after.position_x[i], before.position_x[i], "ungrouped position_x");
  expect_near(after.position_y[i], before.position_y[i], "ungrouped position_y");
  expect_near(after.position_z[i], before.position_z[i], "ungrouped position_z");
  expect_near(after.velocity_x[i], before.velocity_x[i], "ungrouped velocity_x");
  expect_near(after.velocity_y[i], before.velocity_y[i], "ungrouped velocity_y");
  expect_near(after.velocity_z[i], before.velocity_z[i], "ungrouped velocity_z");
  expect_near(after.omega_x[i], before.omega_x[i], "ungrouped omega_x");
  expect_near(after.omega_y[i], before.omega_y[i], "ungrouped omega_y");
  expect_near(after.omega_z[i], before.omega_z[i], "ungrouped omega_z");
}

void check_forces_reset(const HostParticleData &after)
{
  for (std::size_t i = 0; i < after.size(); ++i) {
    expect_near(after.force_x[i], 0.0, "force_x reset");
    expect_near(after.force_y[i], 0.0, "force_y reset");
    expect_near(after.force_z[i], 0.0, "force_z reset");
    expect_near(after.torque_x[i], 0.0, "torque_x reset");
    expect_near(after.torque_y[i], 0.0, "torque_y reset");
    expect_near(after.torque_z[i], 0.0, "torque_z reset");
  }
}

} // namespace

int main()
{
  try {
    GpuDemContext context(0, PrecisionMode::Mixed);
    const HostParticleData initial = make_particles();

    GpuParticleData particles;
    particles.copy_from_host(initial, context.stream());
    context.synchronize();

    const double dt = 0.01;
    const double dtf = 0.5 * dt;
    nve_sphere_initial_integrate(particles, dt, dtf, 1, 3, 1.0, context.stream());
    nve_sphere_final_integrate(particles, dtf, 1, 3, 1.0, context.stream());
    context.synchronize();

    HostParticleData integrated;
    particles.copy_to_host(integrated, context.stream());
    context.synchronize();

    check_integrated_particle(initial, integrated);
    check_ungrouped_particle_unchanged(initial, integrated);

    reset_forces(particles, context.stream());
    context.synchronize();

    HostParticleData reset;
    particles.copy_to_host(reset, context.stream());
    context.synchronize();
    check_forces_reset(reset);

    std::cout << "GPU_DEM nve/sphere no-contact integration test passed on device "
              << context.device_id() << " ("
              << context.device_properties().name << ")\n";
  } catch (const std::exception &error) {
    std::cerr << "GPU_DEM nve/sphere no-contact integration test failed: "
              << error.what() << "\n";
    return 1;
  }

  return 0;
}
