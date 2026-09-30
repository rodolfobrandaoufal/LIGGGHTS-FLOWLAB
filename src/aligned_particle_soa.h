/* ----------------------------------------------------------------------
   Aligned Structure-of-Arrays storage for particle integration kernels.

   This header is intentionally independent from Atom/AtomVec allocation so
   existing MPI exchange, restart, dump, and pair-style code can keep using
   the legacy LIGGGHTS per-atom arrays while selected kernels are ported.
------------------------------------------------------------------------- */

#ifndef LMP_ALIGNED_PARTICLE_SOA_H
#define LMP_ALIGNED_PARTICLE_SOA_H

#include <cstdlib>
#include <new>
#include <vector>

namespace LAMMPS_NS {

template <class T, std::size_t Alignment>
class AlignedAllocator {
 public:
  typedef T value_type;
  typedef T *pointer;
  typedef const T *const_pointer;
  typedef T &reference;
  typedef const T &const_reference;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

  template <class U>
  struct rebind {
    typedef AlignedAllocator<U, Alignment> other;
  };

  AlignedAllocator() throw() {}

  template <class U>
  AlignedAllocator(const AlignedAllocator<U, Alignment> &) throw() {}

  pointer allocate(size_type n, const void * = 0)
  {
    void *ptr = 0;
    if (n == 0) return 0;
    if (posix_memalign(&ptr, Alignment, n * sizeof(T)) != 0)
      throw std::bad_alloc();
    return static_cast<pointer>(ptr);
  }

  void deallocate(pointer p, size_type) throw()
  {
    std::free(p);
  }

  size_type max_size() const throw()
  {
    return static_cast<size_type>(-1) / sizeof(T);
  }

  void construct(pointer p, const T &value)
  {
    new(static_cast<void *>(p)) T(value);
  }

  void destroy(pointer p)
  {
    p->~T();
  }
};

template <class T1, class T2, std::size_t Alignment>
inline bool operator==(const AlignedAllocator<T1, Alignment> &,
                       const AlignedAllocator<T2, Alignment> &)
{
  return true;
}

template <class T1, class T2, std::size_t Alignment>
inline bool operator!=(const AlignedAllocator<T1, Alignment> &,
                       const AlignedAllocator<T2, Alignment> &)
{
  return false;
}

class ParticleSoA {
 public:
  typedef std::vector<double, AlignedAllocator<double, 64> > AlignedDoubleVector;

  ParticleSoA() : n_(0) {}

  void resize(int n)
  {
    if (n <= static_cast<int>(pos_x_.size())) {
      n_ = n;
      return;
    }

    pos_x_.resize(n);
    pos_y_.resize(n);
    pos_z_.resize(n);
    vel_x_.resize(n);
    vel_y_.resize(n);
    vel_z_.resize(n);
    force_x_.resize(n);
    force_y_.resize(n);
    force_z_.resize(n);
    inv_mass_.resize(n);
    n_ = n;
  }

  void load_from_aos(double **x, double **v, double **f, const double *rmass,
                     const double *mass, const int *type, int n)
  {
    resize(n);

    for (int i = 0; i < n; ++i) {
      pos_x_[i] = x[i][0];
      pos_y_[i] = x[i][1];
      pos_z_[i] = x[i][2];
      vel_x_[i] = v[i][0];
      vel_y_[i] = v[i][1];
      vel_z_[i] = v[i][2];
      force_x_[i] = f[i][0];
      force_y_[i] = f[i][1];
      force_z_[i] = f[i][2];
      inv_mass_[i] = rmass ? 1.0 / rmass[i] : 1.0 / mass[type[i]];
    }
  }

  void store_xv_to_aos(double **x, double **v, int n) const
  {
    for (int i = 0; i < n; ++i) {
      x[i][0] = pos_x_[i];
      x[i][1] = pos_y_[i];
      x[i][2] = pos_z_[i];
      v[i][0] = vel_x_[i];
      v[i][1] = vel_y_[i];
      v[i][2] = vel_z_[i];
    }
  }

  void store_v_to_aos(double **v, int n) const
  {
    for (int i = 0; i < n; ++i) {
      v[i][0] = vel_x_[i];
      v[i][1] = vel_y_[i];
      v[i][2] = vel_z_[i];
    }
  }

  void initial_integrate_nve(double dtv, double dtf, const int *mask,
                             int groupbit, int n)
  {
    double * const pos_x = &pos_x_[0];
    double * const pos_y = &pos_y_[0];
    double * const pos_z = &pos_z_[0];
    double * const vel_x = &vel_x_[0];
    double * const vel_y = &vel_y_[0];
    double * const vel_z = &vel_z_[0];
    const double * const force_x = &force_x_[0];
    const double * const force_y = &force_y_[0];
    const double * const force_z = &force_z_[0];
    const double * const inv_mass = &inv_mass_[0];

    for (int i = 0; i < n; ++i) {
      if (mask[i] & groupbit) {
        const double dtfm = dtf * inv_mass[i];
        vel_x[i] += dtfm * force_x[i];
        vel_y[i] += dtfm * force_y[i];
        vel_z[i] += dtfm * force_z[i];
        pos_x[i] += dtv * vel_x[i];
        pos_y[i] += dtv * vel_y[i];
        pos_z[i] += dtv * vel_z[i];
      }
    }
  }

  void final_integrate_nve(double dtf, const int *mask, int groupbit, int n)
  {
    double * const vel_x = &vel_x_[0];
    double * const vel_y = &vel_y_[0];
    double * const vel_z = &vel_z_[0];
    const double * const force_x = &force_x_[0];
    const double * const force_y = &force_y_[0];
    const double * const force_z = &force_z_[0];
    const double * const inv_mass = &inv_mass_[0];

    for (int i = 0; i < n; ++i) {
      if (mask[i] & groupbit) {
        const double dtfm = dtf * inv_mass[i];
        vel_x[i] += dtfm * force_x[i];
        vel_y[i] += dtfm * force_y[i];
        vel_z[i] += dtfm * force_z[i];
      }
    }
  }

  int size() const { return n_; }

 private:
  int n_;
  AlignedDoubleVector pos_x_, pos_y_, pos_z_;
  AlignedDoubleVector vel_x_, vel_y_, vel_z_;
  AlignedDoubleVector force_x_, force_y_, force_z_;
  AlignedDoubleVector inv_mass_;
};

} // namespace LAMMPS_NS

#endif
