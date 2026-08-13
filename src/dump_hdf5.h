/* ----------------------------------------------------------------------
   Parallel HDF5 particle dump.

   Syntax:
     dump ID group-ID hdf5 N file.h5

   Build requires parallel HDF5 and -DLIGGGHTS_HDF5.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(hdf5,DumpHDF5)

#else

#ifndef LMP_DUMP_HDF5_H
#define LMP_DUMP_HDF5_H

#include "dump.h"
#include <vector>

namespace LAMMPS_NS {

class DumpHDF5 : public Dump {
 public:
  DumpHDF5(class LAMMPS *, int, char **);
  virtual ~DumpHDF5();
  virtual void write();

 protected:
  virtual void init_style();
  virtual void write_header(bigint) {}
  virtual void pack(int *) {}
  virtual void write_data(int, double *) {}

 private:
  std::vector<bigint> written_steps_;
  std::vector<long long> written_counts_;

  void pack_local(long long *ids, double *positions, double *velocities,
                  double *forces, double *omegas, double *radii,
                  int nlocal_selected) const;
  void write_xdmf(const char *h5name) const;
  void resolve_filename(char *out, int nout) const;
};

}

#endif
#endif
