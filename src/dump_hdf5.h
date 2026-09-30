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

    Parallel HDF5 particle dump with an XDMF sidecar.

    Syntax:
      dump ID group-ID hdf5 N file.h5

    See doc/dump_hdf5.txt. Build requires parallel HDF5 and -DLIGGGHTS_HDF5.
    This header also declares the HDF5/XDMF helpers shared with
    dump mesh/hdf5 (dump_mesh_hdf5.cpp).
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(hdf5,DumpHDF5)

#else

#ifndef LMP_DUMP_HDF5_H
#define LMP_DUMP_HDF5_H

#include "dump.h"
#include <stdio.h>
#include <string>
#include <vector>

#ifdef LIGGGHTS_HDF5
#include <hdf5.h>
#if !defined(H5_HAVE_PARALLEL)
#error "dump hdf5 / dump mesh/hdf5 (-DLIGGGHTS_HDF5) require an MPI-parallel HDF5 build (H5_HAVE_PARALLEL is not defined). Point the include/library paths to the parallel HDF5 (e.g. /usr/include/hdf5/openmpi, or use h5pcc) or build without -DLIGGGHTS_HDF5."
#endif
#endif

namespace LAMMPS_NS {

class Error;

namespace DumpHDF5Util {

/* one time step as listed in the XDMF temporal collection */
struct StepEntry {
  bigint step;
  double time;
  long long count;      // particles (hdf5) or triangles (mesh/hdf5)
  std::string h5ref;    // HDF5 file name relative to the XDMF file
  bool has_type;        // particle dump: 'type' dataset present
  StepEntry() : step(0), time(0.0), count(0), has_type(false) {}
};

/* incremental XDMF temporal-collection writer (used on one rank only).
   The closing tags are kept at the end of the file; append() seeks to
   the start of the closing tags, writes one <Grid> and re-writes the
   closing tags, so the cost per dump does not depend on the number of
   earlier dumps. rewrite() writes a complete file (rare path). */
class XdmfSeriesWriter {
 public:
  XdmfSeriesWriter();
  ~XdmfSeriesWriter();
  bool rewrite(const std::string &path, const std::string &collection,
               const std::string &grids);
  bool append(const std::string &grid);
  void close();
  bool is_open() const { return fp_ != NULL; }
  const std::string &path() const { return path_; }

 private:
  FILE *fp_;
  long long body_end_;
  std::string path_;
  XdmfSeriesWriter(const XdmfSeriesWriter &);
  XdmfSeriesWriter &operator=(const XdmfSeriesWriter &);
};

/* file name with '*' replaced by the time step (honours dump_modify pad) */
std::string expand_star(const char *pattern, bigint step, int pad);
std::string expand_star_text(const char *pattern, const char *text);
std::string basename_of(const std::string &path);
bool file_exists(const std::string &path);

/* sorted insert; returns 1 if appended at the end, 0 if inserted in the
   middle, 2 if an entry with the same step was replaced */
int insert_entry(std::vector<StepEntry> &entries, const StepEntry &e);

/* append mode: LIGGGHTS does not store the elapsed time in restart files
   (read_restart restarts update->atime at 0). If the file already holds a
   step s_ref <= step, return the offset that makes the time continue from
   it: t(s_ref) + (step - s_ref)*dt - time_now; else 0. */
double continuity_offset(const std::vector<StepEntry> &entries, bigint step,
                         double time_now, double dt);

/* XDMF fragments with explicit NumberType/Precision on every DataItem */
std::string data_item(const char *numtype, int precision, long long rows,
                      int ncol, const std::string &h5ref, bigint step,
                      const char *path);
std::string grid_open(const StepEntry &e);
std::string attribute_xml(const char *name, const char *atype,
                          const char *center, const std::string &item);

#ifdef LIGGGHTS_HDF5

/* RAII owner of one HDF5 identifier */
class H5Id {
 public:
  typedef herr_t (*closer_t)(hid_t);
  H5Id() : id_(-1), closer_(0) {}
  H5Id(hid_t id, closer_t closer) : id_(id), closer_(closer) {}
  ~H5Id() { reset(); }
  void reset(hid_t id = -1, closer_t closer = 0);
  herr_t close();                 // close now, return the HDF5 status
  hid_t get() const { return id_; }
  bool valid() const { return id_ >= 0; }
 private:
  hid_t id_;
  closer_t closer_;
  H5Id(const H5Id &);
  H5Id &operator=(const H5Id &);
};

/* collective error state of one dump write.
   Every HDF5 call is checked with fail(); sync() is called by all ranks
   at stage boundaries (after file open, after each dataset, after
   flush/close). If any rank failed, the message of the lowest failing
   rank is broadcast and all ranks call error->all(). */
class Status {
 public:
  Status(Error *error, MPI_Comm world, const char *style);
  bool ok() const { return !failed_; }
  // record a failure of 'what' on this rank (HDF5 error stack appended)
  void fail(const char *what);
  // check an HDF5 return value (< 0 means failure)
  bool check(long long rv, const char *what) { if (rv < 0) fail(what); return rv >= 0; }
  void sync(const char *file, int line);
 private:
  Error *error_;
  MPI_Comm world_;
  std::string style_;
  bool failed_;
  std::string msg_;
};

/* select this rank's rows [offset, offset+count) of a (rows x ncol) or
   (rows) dataset and write them collectively; handles count == 0 */
void write_dataset(Status &st, hid_t loc, const char *name, hid_t type,
                   int rank, hsize_t rows_global, hsize_t rows_local,
                   hsize_t row_offset, hsize_t ncol, hid_t dxpl,
                   const void *data);

void write_scalar_attribute(Status &st, hid_t loc, const char *name,
                            hid_t type, const void *value);

/* read the list of Step_<n> groups of an open file (all ranks may call
   it). count_dataset is the dataset whose first dimension is the item
   count. Groups without a 'time' attribute (files written before the
   attribute existed) get time = time_base + (step - step_base)*dt and are
   counted in nlegacy. Returns false on an HDF5 error. */
bool read_existing_steps(hid_t file, const std::string &h5ref,
                         const char *count_dataset, double time_base,
                         bigint step_base, double dt,
                         std::vector<StepEntry> &entries, int &nlegacy);

/* serial scan (one rank) of existing files matching a '*' pattern */
void scan_multifile_steps(const char *pattern, const char *count_dataset,
                          double time_base, bigint step_base, double dt,
                          std::vector<StepEntry> &entries, int &nlegacy,
                          int &nbad);

/* open (append) or create/truncate a file for collective access */
hid_t open_parallel_file(Status &st, const std::string &name, bool append,
                         MPI_Comm world, bool &existed);

/* create /Step_<step>; an existing group of that name is deleted first
   (replaced = true) */
void create_step_group(Status &st, hid_t file, bigint step, H5Id &group,
                       bool &replaced);

/* 'timestep' (int64) and 'time' (double) attributes of a step group */
void write_step_attributes(Status &st, hid_t group, bigint step, double time);

/* re-open an existing dump file (read/write) for the next dump */
hid_t reopen_parallel_file(Status &st, const std::string &name, MPI_Comm world);

/* close the step group, then close the file if close_file is set.
   Measured (audit/fixes/hdf5/REPORT.md): close + re-open per dump costs
   ~0.6 ms, H5Fflush on the MPI-IO driver (MPI_File_sync = fsync) ~1.7 ms
   and more, keeping the file open without flushing ~0.35 ms. */
void finish_step(Status &st, H5Id &group, H5Id &file, bool close_file);

#endif

}

class DumpHDF5 : public Dump {
 public:
  DumpHDF5(class LAMMPS *, int, char **);
  virtual ~DumpHDF5();
  virtual void write();

 protected:
  virtual void init_style();
  virtual int modify_param(int, char **);
  virtual void write_header(bigint) {}
  virtual void pack(int *) {}
  virtual void write_data(int, double *) {}

 private:
#ifdef LIGGGHTS_HDF5
  DumpHDF5Util::H5Id file_;           // single file; stays open between dumps only with dump_modify flush no
#endif
  std::string open_name_;             // name of the open file
  std::vector<DumpHDF5Util::StepEntry> entries_;
  DumpHDF5Util::XdmfSeriesWriter xdmf_;   // rank 0 only
  bool truncate_warned_;
  bool opened_once_;                  // single file: created/opened before
  double time_offset_;                // append mode, see continuity_offset()
  bool time_offset_set_;

  std::string grid_xml(const DumpHDF5Util::StepEntry &e) const;
  std::string all_grids_xml() const;
  void update_xdmf(const DumpHDF5Util::StepEntry &e, int where);
  void pack_local(long long *ids, int *types, double *positions,
                  double *velocities, double *forces, double *omegas,
                  double *radii, int nlocal_selected) const;
};

}

#endif
#endif
