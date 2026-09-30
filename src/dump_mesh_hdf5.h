/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    LIGGGHTS modernization branch

    Parallel HDF5 triangular mesh dump with an XDMF sidecar.

   Syntax:
     dump ID group-ID mesh/hdf5 N file.h5 [all|mesh-id ...] [properties ...]

   Built-in cell properties:
     id owner area normal aedges acorners index nneighs

   Additional property names are looked up as mesh element containers:
     ScalarContainer<double>      -> /CellData/<name>
     VectorContainer<double,3>    -> /CellData/<name> with Dimensions="M 3"

   See doc/dump_hdf5.txt. Build requires parallel HDF5 and -DLIGGGHTS_HDF5.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(mesh/hdf5,DumpMeshHDF5)
DumpStyle(mesh/gran/HDF5,DumpMeshHDF5)

#else

#ifndef LMP_DUMP_MESH_HDF5_H
#define LMP_DUMP_MESH_HDF5_H

#include "dump.h"
#include "dump_hdf5.h"
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
  virtual int modify_param(int, char **);
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
#ifdef LIGGGHTS_HDF5
  DumpHDF5Util::H5Id file_;           // single file; stays open between dumps only with dump_modify flush no
#endif
  std::string open_name_;
  std::vector<DumpHDF5Util::StepEntry> entries_;
  DumpHDF5Util::XdmfSeriesWriter xdmf_;   // rank 0 only
  bool truncate_warned_;
  bool opened_once_;
  double time_offset_;
  bool time_offset_set_;
  std::vector<std::string> requested_properties_;
  std::vector<std::string> scalar_properties_;
  std::vector<std::string> vector_properties_;

  int dump_all_meshes_;
  int builtin_mask_;

  void add_mesh_by_id(const char *id);
  void parse_property(const char *name);
  std::string grid_xml(const DumpHDF5Util::StepEntry &e) const;
  std::string all_grids_xml() const;
  void update_xdmf(const DumpHDF5Util::StepEntry &e, int where);
};

}

#endif
#endif
