/* ----------------------------------------------------------------------
   GPU triangle-mesh wall contact scaffold for GPU_DEM.
------------------------------------------------------------------------- */

#include "gpu_mesh_wall.h"

#include "gpu_error_check.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMMPS_NS {
namespace GPU_DEM {
namespace {

constexpr int THREADS_PER_BLOCK = 128;

void compute_host_normal(TriangleWall &triangle)
{
  const double e0x = triangle.v1[0] - triangle.v0[0];
  const double e0y = triangle.v1[1] - triangle.v0[1];
  const double e0z = triangle.v1[2] - triangle.v0[2];
  const double e1x = triangle.v2[0] - triangle.v0[0];
  const double e1y = triangle.v2[1] - triangle.v0[1];
  const double e1z = triangle.v2[2] - triangle.v0[2];
  const double nx = e0y * e1z - e0z * e1y;
  const double ny = e0z * e1x - e0x * e1z;
  const double nz = e0x * e1y - e0y * e1x;
  const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
  if (norm <= 0.0) throw std::runtime_error("GPU_DEM STL contains a degenerate triangle");

  triangle.normal[0] = nx / norm;
  triangle.normal[1] = ny / norm;
  triangle.normal[2] = nz / norm;
}

__device__ double dot3(double ax, double ay, double az, double bx, double by, double bz)
{
  return ax * bx + ay * by + az * bz;
}

__device__ void closest_point_on_triangle(double px, double py, double pz,
                                          const TriangleWall &tri,
                                          double &cx, double &cy, double &cz)
{
  const double ax = tri.v0[0], ay = tri.v0[1], az = tri.v0[2];
  const double bx = tri.v1[0], by = tri.v1[1], bz = tri.v1[2];
  const double cx0 = tri.v2[0], cy0 = tri.v2[1], cz0 = tri.v2[2];

  const double abx = bx - ax, aby = by - ay, abz = bz - az;
  const double acx = cx0 - ax, acy = cy0 - ay, acz = cz0 - az;
  const double apx = px - ax, apy = py - ay, apz = pz - az;
  const double d1 = dot3(abx, aby, abz, apx, apy, apz);
  const double d2 = dot3(acx, acy, acz, apx, apy, apz);
  if (d1 <= 0.0 && d2 <= 0.0) {
    cx = ax; cy = ay; cz = az;
    return;
  }

  const double bpx = px - bx, bpy = py - by, bpz = pz - bz;
  const double d3 = dot3(abx, aby, abz, bpx, bpy, bpz);
  const double d4 = dot3(acx, acy, acz, bpx, bpy, bpz);
  if (d3 >= 0.0 && d4 <= d3) {
    cx = bx; cy = by; cz = bz;
    return;
  }

  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    const double v = d1 / (d1 - d3);
    cx = ax + v * abx;
    cy = ay + v * aby;
    cz = az + v * abz;
    return;
  }

  const double cpx = px - cx0, cpy = py - cy0, cpz = pz - cz0;
  const double d5 = dot3(abx, aby, abz, cpx, cpy, cpz);
  const double d6 = dot3(acx, acy, acz, cpx, cpy, cpz);
  if (d6 >= 0.0 && d5 <= d6) {
    cx = cx0; cy = cy0; cz = cz0;
    return;
  }

  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    const double w = d2 / (d2 - d6);
    cx = ax + w * acx;
    cy = ay + w * acy;
    cz = az + w * acz;
    return;
  }

  const double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    cx = bx + w * (cx0 - bx);
    cy = by + w * (cy0 - by);
    cz = bz + w * (cz0 - bz);
    return;
  }

  const double denom = 1.0 / (va + vb + vc);
  const double v = vb * denom;
  const double w = vc * denom;
  cx = ax + abx * v + acx * w;
  cy = ay + aby * v + acy * w;
  cz = az + abz * v + acz * w;
}

__device__ double hertz_wall_force(double overlap, double normal_relative_velocity,
                                   double radius, double mass,
                                   const HertzNormalParams &params)
{
  const double sqrt_delta_reff = sqrt(radius * overlap);
  const double sn = 2.0 * params.effective_youngs_modulus * sqrt_delta_reff;
  const double kn = (4.0 / 3.0) * params.effective_youngs_modulus * sqrt_delta_reff;
  const double sqrt_five_over_six = 0.91287092917527685576161630466800355659;
  const double gamman = -2.0 * sqrt_five_over_six * params.beta_effective *
                        sqrt(sn * mass);

  double normal_force = kn * overlap - gamman * normal_relative_velocity;
  if (params.limit_force && normal_force < 0.0) normal_force = 0.0;
  return normal_force;
}

__global__ void hertz_triangle_wall_force_kernel(
    std::size_t n, const double *x, const double *y, const double *z,
    const double *vx, const double *vy, const double *vz, const double *radius,
    const double *mass, const TriangleWall *triangles, std::size_t triangle_count,
    HertzNormalParams params, double *fx, double *fy, double *fz)
{
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  double fxi = 0.0;
  double fyi = 0.0;
  double fzi = 0.0;

  for (std::size_t tri_index = 0; tri_index < triangle_count; ++tri_index) {
    const TriangleWall tri = triangles[tri_index];
    double cx, cy, cz;
    closest_point_on_triangle(x[i], y[i], z[i], tri, cx, cy, cz);

    double nx = x[i] - cx;
    double ny = y[i] - cy;
    double nz = z[i] - cz;
    double distance = sqrt(nx * nx + ny * ny + nz * nz);
    if (distance <= 1.0e-30) {
      nx = tri.normal[0];
      ny = tri.normal[1];
      nz = tri.normal[2];
      distance = 0.0;
    } else {
      const double inv_distance = 1.0 / distance;
      nx *= inv_distance;
      ny *= inv_distance;
      nz *= inv_distance;
    }

    const double overlap = radius[i] - distance;
    if (overlap <= 0.0) continue;

    const double normal_relative_velocity = vx[i] * nx + vy[i] * ny + vz[i] * nz;
    const double normal_force =
        hertz_wall_force(overlap, normal_relative_velocity, radius[i], mass[i], params);
    fxi += normal_force * nx;
    fyi += normal_force * ny;
    fzi += normal_force * nz;
  }

  fx[i] += fxi;
  fy[i] += fyi;
  fz[i] += fzi;
}

} // namespace

HostTriangleMesh load_ascii_stl_triangle_mesh(const std::string &path, double scale)
{
  std::ifstream input(path.c_str());
  if (!input) throw std::runtime_error("GPU_DEM could not open ASCII STL: " + path);

  HostTriangleMesh mesh;
  TriangleWall triangle;
  int vertex_index = 0;
  std::string line;

  while (std::getline(input, line)) {
    std::istringstream line_stream(line);
    std::string keyword;
    line_stream >> keyword;
    if (keyword != "vertex") continue;

    double x, y, z;
    if (!(line_stream >> x >> y >> z))
      throw std::runtime_error("GPU_DEM failed to parse STL vertex in: " + path);

    double *vertex = vertex_index == 0 ? triangle.v0 :
                     vertex_index == 1 ? triangle.v1 : triangle.v2;
    vertex[0] = x * scale;
    vertex[1] = y * scale;
    vertex[2] = z * scale;
    ++vertex_index;

    if (vertex_index == 3) {
      compute_host_normal(triangle);
      mesh.triangles.push_back(triangle);
      triangle = TriangleWall();
      vertex_index = 0;
    }
  }

  if (vertex_index != 0)
    throw std::runtime_error("GPU_DEM STL ended with an incomplete triangle: " + path);
  if (mesh.triangles.empty())
    throw std::runtime_error("GPU_DEM STL contained no triangles: " + path);

  return mesh;
}

GpuTriangleMesh::~GpuTriangleMesh()
{
  release();
}

GpuTriangleMesh::GpuTriangleMesh(GpuTriangleMesh &&other) noexcept
{
  *this = std::move(other);
}

GpuTriangleMesh &GpuTriangleMesh::operator=(GpuTriangleMesh &&other) noexcept
{
  if (this == &other) return *this;
  release();

  size_ = other.size_;
  capacity_ = other.capacity_;
  triangles_ = other.triangles_;

  other.size_ = 0;
  other.capacity_ = 0;
  other.triangles_ = nullptr;
  return *this;
}

void GpuTriangleMesh::copy_from_host(const HostTriangleMesh &host, cudaStream_t stream)
{
  if (host.triangles.size() > capacity_) {
    release();
    GPU_DEM_CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&triangles_),
                                  host.triangles.size() * sizeof(TriangleWall)));
    capacity_ = host.triangles.size();
  }

  size_ = host.triangles.size();
  if (size_ > 0) {
    GPU_DEM_CUDA_CHECK(cudaMemcpyAsync(triangles_, host.triangles.data(),
                                       size_ * sizeof(TriangleWall),
                                       cudaMemcpyHostToDevice, stream));
  }
}

void GpuTriangleMesh::release() noexcept
{
  if (triangles_) {
    cudaFree(triangles_);
    triangles_ = nullptr;
  }
  size_ = 0;
  capacity_ = 0;
}

void compute_hertz_triangle_wall_forces(GpuParticleData &particles,
                                        const GpuTriangleMesh &mesh,
                                        const HertzNormalParams &params,
                                        cudaStream_t stream)
{
  if (!(params.effective_youngs_modulus >= 0.0))
    throw std::runtime_error("GPU_DEM Hertz mesh-wall effective Young's modulus must be nonnegative");
  if (mesh.size() == 0 || particles.size() == 0) return;

  const int blocks =
      static_cast<int>((particles.size() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  hertz_triangle_wall_force_kernel<<<blocks, THREADS_PER_BLOCK, 0, stream>>>(
      particles.size(), particles.position_x(), particles.position_y(), particles.position_z(),
      particles.velocity_x(), particles.velocity_y(), particles.velocity_z(),
      particles.radius(), particles.mass(), mesh.device_triangles(), mesh.size(), params,
      particles.force_x(), particles.force_y(), particles.force_z());
  GPU_DEM_CUDA_CHECK(cudaGetLastError());
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
