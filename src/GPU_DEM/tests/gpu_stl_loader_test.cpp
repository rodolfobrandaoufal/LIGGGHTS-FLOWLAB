/* ----------------------------------------------------------------------
   Standalone GPU_DEM ASCII STL loader/upload regression test.
------------------------------------------------------------------------- */

#include "gpu_dem_context.h"
#include "gpu_mesh_wall.h"

#include <cstdio>
#include <cstdlib>
#include <string>

using LAMMPS_NS::GPU_DEM::GpuDemContext;
using LAMMPS_NS::GPU_DEM::GpuTriangleMesh;
using LAMMPS_NS::GPU_DEM::HostTriangleMesh;
using LAMMPS_NS::GPU_DEM::load_ascii_stl_triangle_mesh;

namespace {

int parse_int(const char *text)
{
  return std::atoi(text);
}

double parse_double(const char *text)
{
  return std::atof(text);
}

} // namespace

int main(int argc, char **argv)
{
  if (argc != 4) {
    std::fprintf(stderr, "usage: %s <ascii-stl> <scale> <expected-triangles>\n", argv[0]);
    return 2;
  }

  const std::string path = argv[1];
  const double scale = parse_double(argv[2]);
  const int expected_triangles = parse_int(argv[3]);
  if (!(scale > 0.0) || expected_triangles <= 0) {
    std::fprintf(stderr, "invalid scale or expected triangle count\n");
    return 2;
  }

  HostTriangleMesh host = load_ascii_stl_triangle_mesh(path, scale);
  if (static_cast<int>(host.triangles.size()) != expected_triangles) {
    std::fprintf(stderr, "triangle count mismatch for %s: actual %zu expected %d\n",
                 path.c_str(), host.triangles.size(), expected_triangles);
    return 1;
  }

  GpuDemContext context(0);
  GpuTriangleMesh mesh;
  mesh.copy_from_host(host, context.stream());
  context.synchronize();

  if (mesh.size() != host.triangles.size()) {
    std::fprintf(stderr, "GPU mesh size mismatch: actual %zu expected %zu\n",
                 mesh.size(), host.triangles.size());
    return 1;
  }

  std::printf("GPU_DEM STL loader test passed on device %d (%s): triangles=%zu path=%s\n",
              context.device_id(), context.device_properties().name,
              host.triangles.size(), path.c_str());
  return 0;
}
