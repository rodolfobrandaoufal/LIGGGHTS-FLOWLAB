/* ----------------------------------------------------------------------
   Parallel HDF5 triangular mesh dump.

   Syntax:
     dump ID group-ID mesh/hdf5 N file.h5 [all|mesh-id ...] [properties ...]

   Built-in cell properties:
     id owner area normal aedges acorners index nneighs

   Additional property names are looked up as mesh element containers:
     ScalarContainer<double>      -> /CellData/<name>
     VectorContainer<double,3>    -> /CellData/<name> with Dimensions="M 3"

   Build requires parallel HDF5 and -DLIGGGHTS_HDF5.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(mesh/hdf5,DumpMeshHDF5)
DumpStyle(mesh/gran/HDF5,DumpMeshHDF5)

#else

#ifndef LMP_DUMP_MESH_HDF5_H
#define LMP_DUMP_MESH_HDF5_H

#include "dump.h"
#include <vector>
#include <string>

namespace LAMMPS_NS {

class DumpMeshHDF5 : public Dump {
 public:
  DumpMeshHDF5(class LAMMPS *, int, char **);
  virtual ~DumpMeshHDF5();
  virtual void write();

 protected:
  virtual void init_style();
  virtual void write_header(bigint) {}
  virtual void pack(int *) {}
  virtual void write_data(int, double *) {}

 private:
  enum BuiltinProperty {
    PROP_ID       = 1 << 0,
    PROP_OWNER    = 1 << 1,
    PROP_AREA     = 1 << 2,
    PROP_NORMAL   = 1 << 3,
    PROP_AEDGES   = 1 << 4,
    PROP_ACORNERS = 1 << 5,
    PROP_INDEX    = 1 << 6,
    PROP_NNEIGHS  = 1 << 7
  };

  struct ScalarProperty {
    std::string name;
  };

  struct VectorProperty {
    std::string name;
  };

  std::vector<class FixMeshSurface *> fixes_;
  std::vector<bigint> written_steps_;
  std::vector<long long> written_cells_;
  std::vector<std::string> requested_properties_;
  std::vector<std::string> scalar_properties_;
  std::vector<std::string> vector_properties_;

  int dump_all_meshes_;
  int builtin_mask_;

  void add_mesh_by_id(const char *id);
  void parse_property(const char *name);
  void resolve_filename(char *out, int nout) const;
  void write_xdmf(const char *h5name) const;
};

}

#endif
#endif
