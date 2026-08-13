/* ----------------------------------------------------------------------
   Standalone GPU_DEM triangle-mesh wall contact regression test.
------------------------------------------------------------------------- */

#include "gpu_dem_context.h"
#include "gpu_integrator.h"
#include "gpu_mesh_wall.h"
#include "gpu_particle_data.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuParticleData;
using LAMMPS_NS::GPU_DEM::GpuTriangleMesh;
using LAMMPS_NS::GPU_DEM::HertzNormalParams;
using LAMMPS_NS::GPU_DEM::HostParticleData;
using LAMMPS_NS::GPU_DEM::HostTriangleMesh;
using LAMMPS_NS::GPU_DEM::TriangleWall;
using LAMMPS_NS::GPU_DEM::compute_hertz_triangle_wall_forces;
using LAMMPS_NS::GPU_DEM::reset_forces;

namespace {

void expect_near(double actual, double expected, double tol, const char *label)
{
  if (std::fabs(actual - expected) <= tol) return;
  std::fprintf(stderr, "%s mismatch: actual %.17g expected %.17g tolerance %.3g\n",
               label, actual, expected, tol);
  std::exit(1);
}

HostParticleData make_particle(double x, double y, double z)
{
  HostParticleData host;
  host.resize(1);
  host.position_x[0] = x;
  host.position_y[0] = y;
  host.position_z[0] = z;
  host.velocity_x[0] = 0.0;
  host.velocity_y[0] = 0.0;
  host.velocity_z[0] = 0.0;
  host.omega_x[0] = host.omega_y[0] = host.omega_z[0] = 0.0;
  host.force_x[0] = host.force_y[0] = host.force_z[0] = 0.0;
  host.torque_x[0] = host.torque_y[0] = host.torque_z[0] = 0.0;
  host.radius[0] = 0.5;
  host.mass[0] = 2.0;
  host.density[0] = 1.0;
  host.type[0] = 1;
  host.mask[0] = 1;
  host.tag[0] = 1;
  host.image_flags[0] = 0;
  return host;
}

HostTriangleMesh make_single_triangle()
{
  HostTriangleMesh host;
  TriangleWall triangle;
  triangle.v0[0] = -1.0; triangle.v0[1] = -1.0; triangle.v0[2] = 0.0;
  triangle.v1[0] =  1.0; triangle.v1[1] = -1.0; triangle.v1[2] = 0.0;
  triangle.v2[0] =  0.0; triangle.v2[1] =  1.0; triangle.v2[2] = 0.0;
  triangle.normal[0] = 0.0;
  triangle.normal[1] = 0.0;
  triangle.normal[2] = 1.0;
  host.triangles.push_back(triangle);
  return host;
}

double expected_hertz_wall_force(double overlap, double radius, double effective_young)
{
  const double sqrt_delta_reff = std::sqrt(radius * overlap);
  const double kn = (4.0 / 3.0) * effective_young * sqrt_delta_reff;
  return kn * overlap;
}

} // namespace

int main()
{
  GpuDemContext context(0);

  HostTriangleMesh host_mesh = make_single_triangle();
  GpuTriangleMesh mesh;
  mesh.copy_from_host(host_mesh, context.stream());

  GpuParticleData particle;
  particle.copy_from_host(make_particle(0.0, 0.0, 0.4), context.stream());

  HertzNormalParams params;
  params.effective_youngs_modulus = 120.0;
  params.beta_effective = 0.0;
  params.limit_force = true;

  reset_forces(particle, context.stream());
  compute_hertz_triangle_wall_forces(particle, mesh, params, context.stream());

  HostParticleData result;
  particle.copy_to_host(result, context.stream());
  context.synchronize();

  const double expected_force = expected_hertz_wall_force(0.1, 0.5,
                                                         params.effective_youngs_modulus);
  expect_near(result.force_x[0], 0.0, 1.0e-12, "mesh force x");
  expect_near(result.force_y[0], 0.0, 1.0e-12, "mesh force y");
  expect_near(result.force_z[0], expected_force, 1.0e-12, "mesh force z");

  std::printf("GPU_DEM Hertz triangle mesh-wall contact test passed on device %d (%s)\n",
              context.device_id(), context.device_properties().name);
  return 0;
}
