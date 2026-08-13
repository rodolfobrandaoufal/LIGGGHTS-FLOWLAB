#include "gpu_particle_data.h"

#include <stdexcept>
#include <utility>

namespace LAMMPS_NS {
namespace GPU_DEM {

void HostParticleData::resize(std::size_t n)
{
  position_x.resize(n); position_y.resize(n); position_z.resize(n);
  velocity_x.resize(n); velocity_y.resize(n); velocity_z.resize(n);
  omega_x.resize(n); omega_y.resize(n); omega_z.resize(n);
  force_x.resize(n); force_y.resize(n); force_z.resize(n);
  torque_x.resize(n); torque_y.resize(n); torque_z.resize(n);
  radius.resize(n);
  mass.resize(n);
  density.resize(n);
  type.resize(n);
  mask.resize(n);
  tag.resize(n);
  image_flags.resize(n);
}

void HostParticleData::validate_sizes() const
{
  const std::size_t n = size();
  if (position_y.size() != n || position_z.size() != n ||
      velocity_x.size() != n || velocity_y.size() != n || velocity_z.size() != n ||
      omega_x.size() != n || omega_y.size() != n || omega_z.size() != n ||
      force_x.size() != n || force_y.size() != n || force_z.size() != n ||
      torque_x.size() != n || torque_y.size() != n || torque_z.size() != n ||
      radius.size() != n || mass.size() != n || density.size() != n ||
      type.size() != n || mask.size() != n ||
      tag.size() != n || image_flags.size() != n)
    throw std::runtime_error("GPU_DEM host particle arrays have inconsistent sizes");
}

GpuParticleData::~GpuParticleData()
{
  release();
}

GpuParticleData::GpuParticleData(GpuParticleData &&other) noexcept
{
  *this = std::move(other);
}

GpuParticleData &GpuParticleData::operator=(GpuParticleData &&other) noexcept
{
  if (this == &other) return *this;

  release();
  size_ = other.size_;
  capacity_ = other.capacity_;
  position_x_ = other.position_x_; position_y_ = other.position_y_; position_z_ = other.position_z_;
  velocity_x_ = other.velocity_x_; velocity_y_ = other.velocity_y_; velocity_z_ = other.velocity_z_;
  omega_x_ = other.omega_x_; omega_y_ = other.omega_y_; omega_z_ = other.omega_z_;
  force_x_ = other.force_x_; force_y_ = other.force_y_; force_z_ = other.force_z_;
  torque_x_ = other.torque_x_; torque_y_ = other.torque_y_; torque_z_ = other.torque_z_;
  radius_ = other.radius_; mass_ = other.mass_; density_ = other.density_;
  type_ = other.type_;
  mask_ = other.mask_;
  tag_ = other.tag_;
  image_flags_ = other.image_flags_;

  other.size_ = 0;
  other.capacity_ = 0;
  other.position_x_ = other.position_y_ = other.position_z_ = nullptr;
  other.velocity_x_ = other.velocity_y_ = other.velocity_z_ = nullptr;
  other.omega_x_ = other.omega_y_ = other.omega_z_ = nullptr;
  other.force_x_ = other.force_y_ = other.force_z_ = nullptr;
  other.torque_x_ = other.torque_y_ = other.torque_z_ = nullptr;
  other.radius_ = other.mass_ = other.density_ = nullptr;
  other.type_ = nullptr;
  other.mask_ = nullptr;
  other.tag_ = nullptr;
  other.image_flags_ = nullptr;
  return *this;
}

void GpuParticleData::allocate(std::size_t n)
{
  if (n <= capacity_) {
    size_ = n;
    return;
  }

  release();
  allocate_array(position_x_, n); allocate_array(position_y_, n); allocate_array(position_z_, n);
  allocate_array(velocity_x_, n); allocate_array(velocity_y_, n); allocate_array(velocity_z_, n);
  allocate_array(omega_x_, n); allocate_array(omega_y_, n); allocate_array(omega_z_, n);
  allocate_array(force_x_, n); allocate_array(force_y_, n); allocate_array(force_z_, n);
  allocate_array(torque_x_, n); allocate_array(torque_y_, n); allocate_array(torque_z_, n);
  allocate_array(radius_, n);
  allocate_array(mass_, n);
  allocate_array(density_, n);
  allocate_array(type_, n);
  allocate_array(mask_, n);
  allocate_array(tag_, n);
  allocate_array(image_flags_, n);

  size_ = n;
  capacity_ = n;
}

void GpuParticleData::resize_preserve(std::size_t n, cudaStream_t stream)
{
  if (n <= capacity_) {
    size_ = n;
    return;
  }

  double *new_position_x = nullptr, *new_position_y = nullptr, *new_position_z = nullptr;
  double *new_velocity_x = nullptr, *new_velocity_y = nullptr, *new_velocity_z = nullptr;
  double *new_omega_x = nullptr, *new_omega_y = nullptr, *new_omega_z = nullptr;
  double *new_force_x = nullptr, *new_force_y = nullptr, *new_force_z = nullptr;
  double *new_torque_x = nullptr, *new_torque_y = nullptr, *new_torque_z = nullptr;
  double *new_radius = nullptr, *new_mass = nullptr, *new_density = nullptr;
  int *new_type = nullptr, *new_mask = nullptr;
  long long *new_tag = nullptr, *new_image_flags = nullptr;

  allocate_array(new_position_x, n); allocate_array(new_position_y, n); allocate_array(new_position_z, n);
  allocate_array(new_velocity_x, n); allocate_array(new_velocity_y, n); allocate_array(new_velocity_z, n);
  allocate_array(new_omega_x, n); allocate_array(new_omega_y, n); allocate_array(new_omega_z, n);
  allocate_array(new_force_x, n); allocate_array(new_force_y, n); allocate_array(new_force_z, n);
  allocate_array(new_torque_x, n); allocate_array(new_torque_y, n); allocate_array(new_torque_z, n);
  allocate_array(new_radius, n);
  allocate_array(new_mass, n);
  allocate_array(new_density, n);
  allocate_array(new_type, n);
  allocate_array(new_mask, n);
  allocate_array(new_tag, n);
  allocate_array(new_image_flags, n);

  copy_device_to_device(new_position_x, position_x_, size_, stream);
  copy_device_to_device(new_position_y, position_y_, size_, stream);
  copy_device_to_device(new_position_z, position_z_, size_, stream);
  copy_device_to_device(new_velocity_x, velocity_x_, size_, stream);
  copy_device_to_device(new_velocity_y, velocity_y_, size_, stream);
  copy_device_to_device(new_velocity_z, velocity_z_, size_, stream);
  copy_device_to_device(new_omega_x, omega_x_, size_, stream);
  copy_device_to_device(new_omega_y, omega_y_, size_, stream);
  copy_device_to_device(new_omega_z, omega_z_, size_, stream);
  copy_device_to_device(new_force_x, force_x_, size_, stream);
  copy_device_to_device(new_force_y, force_y_, size_, stream);
  copy_device_to_device(new_force_z, force_z_, size_, stream);
  copy_device_to_device(new_torque_x, torque_x_, size_, stream);
  copy_device_to_device(new_torque_y, torque_y_, size_, stream);
  copy_device_to_device(new_torque_z, torque_z_, size_, stream);
  copy_device_to_device(new_radius, radius_, size_, stream);
  copy_device_to_device(new_mass, mass_, size_, stream);
  copy_device_to_device(new_density, density_, size_, stream);
  copy_device_to_device(new_type, type_, size_, stream);
  copy_device_to_device(new_mask, mask_, size_, stream);
  copy_device_to_device(new_tag, tag_, size_, stream);
  copy_device_to_device(new_image_flags, image_flags_, size_, stream);
  GPU_DEM_CUDA_CHECK(cudaStreamSynchronize(stream));

  release();
  position_x_ = new_position_x; position_y_ = new_position_y; position_z_ = new_position_z;
  velocity_x_ = new_velocity_x; velocity_y_ = new_velocity_y; velocity_z_ = new_velocity_z;
  omega_x_ = new_omega_x; omega_y_ = new_omega_y; omega_z_ = new_omega_z;
  force_x_ = new_force_x; force_y_ = new_force_y; force_z_ = new_force_z;
  torque_x_ = new_torque_x; torque_y_ = new_torque_y; torque_z_ = new_torque_z;
  radius_ = new_radius;
  mass_ = new_mass;
  density_ = new_density;
  type_ = new_type;
  mask_ = new_mask;
  tag_ = new_tag;
  image_flags_ = new_image_flags;
  size_ = n;
  capacity_ = n;
}

void GpuParticleData::release() noexcept
{
  free_array(position_x_); free_array(position_y_); free_array(position_z_);
  free_array(velocity_x_); free_array(velocity_y_); free_array(velocity_z_);
  free_array(omega_x_); free_array(omega_y_); free_array(omega_z_);
  free_array(force_x_); free_array(force_y_); free_array(force_z_);
  free_array(torque_x_); free_array(torque_y_); free_array(torque_z_);
  free_array(radius_);
  free_array(mass_);
  free_array(density_);
  free_array(type_);
  free_array(mask_);
  free_array(tag_);
  free_array(image_flags_);
  size_ = 0;
  capacity_ = 0;
}

void GpuParticleData::copy_from_host(const HostParticleData &host, cudaStream_t stream)
{
  host.validate_sizes();
  allocate(host.size());

  copy_to_device(position_x_, host.position_x, stream);
  copy_to_device(position_y_, host.position_y, stream);
  copy_to_device(position_z_, host.position_z, stream);
  copy_to_device(velocity_x_, host.velocity_x, stream);
  copy_to_device(velocity_y_, host.velocity_y, stream);
  copy_to_device(velocity_z_, host.velocity_z, stream);
  copy_to_device(omega_x_, host.omega_x, stream);
  copy_to_device(omega_y_, host.omega_y, stream);
  copy_to_device(omega_z_, host.omega_z, stream);
  copy_to_device(force_x_, host.force_x, stream);
  copy_to_device(force_y_, host.force_y, stream);
  copy_to_device(force_z_, host.force_z, stream);
  copy_to_device(torque_x_, host.torque_x, stream);
  copy_to_device(torque_y_, host.torque_y, stream);
  copy_to_device(torque_z_, host.torque_z, stream);
  copy_to_device(radius_, host.radius, stream);
  copy_to_device(mass_, host.mass, stream);
  copy_to_device(density_, host.density, stream);
  copy_to_device(type_, host.type, stream);
  copy_to_device(mask_, host.mask, stream);
  copy_to_device(tag_, host.tag, stream);
  copy_to_device(image_flags_, host.image_flags, stream);
}

void GpuParticleData::copy_to_host(HostParticleData &host, cudaStream_t stream) const
{
  host.resize(size_);

  copy_to_host(host.position_x, position_x_, stream);
  copy_to_host(host.position_y, position_y_, stream);
  copy_to_host(host.position_z, position_z_, stream);
  copy_to_host(host.velocity_x, velocity_x_, stream);
  copy_to_host(host.velocity_y, velocity_y_, stream);
  copy_to_host(host.velocity_z, velocity_z_, stream);
  copy_to_host(host.omega_x, omega_x_, stream);
  copy_to_host(host.omega_y, omega_y_, stream);
  copy_to_host(host.omega_z, omega_z_, stream);
  copy_to_host(host.force_x, force_x_, stream);
  copy_to_host(host.force_y, force_y_, stream);
  copy_to_host(host.force_z, force_z_, stream);
  copy_to_host(host.torque_x, torque_x_, stream);
  copy_to_host(host.torque_y, torque_y_, stream);
  copy_to_host(host.torque_z, torque_z_, stream);
  copy_to_host(host.radius, radius_, stream);
  copy_to_host(host.mass, mass_, stream);
  copy_to_host(host.density, density_, stream);
  copy_to_host(host.type, type_, stream);
  copy_to_host(host.mask, mask_, stream);
  copy_to_host(host.tag, tag_, stream);
  copy_to_host(host.image_flags, image_flags_, stream);
}

} // namespace GPU_DEM
} // namespace LAMMPS_NS
