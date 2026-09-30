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

    Parallel HDF5 particle dump (dump hdf5) and the HDF5/XDMF helpers
    shared with dump mesh/hdf5. See doc/dump_hdf5.txt.

    MPI note: every HDF5 call that touches file metadata or dataset data
    (H5Fcreate/open/flush/close, H5Gcreate, H5Dcreate, H5Dwrite with a
    collective transfer, H5Acreate/write, H5Ldelete) is made by all ranks
    in the same order. Errors are recorded per rank and reduced at stage
    boundaries (Status::sync) so that all ranks stop together through
    error->all() instead of some ranks entering the next collective call.
------------------------------------------------------------------------- */

#include "dump_hdf5.h"
#include "atom.h"
#include "error.h"
#include "update.h"
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <algorithm>
#include <math.h>

using namespace LAMMPS_NS;
using namespace DumpHDF5Util;

/* ======================================================================
   shared helpers
   ====================================================================== */

static const char XDMF_TAIL[] = "    </Grid>\n  </Domain>\n</Xdmf>\n";

XdmfSeriesWriter::XdmfSeriesWriter() : fp_(NULL), body_end_(0) {}

XdmfSeriesWriter::~XdmfSeriesWriter()
{
  close();
}

void XdmfSeriesWriter::close()
{
  if (fp_) fclose(fp_);
  fp_ = NULL;
  body_end_ = 0;
}

bool XdmfSeriesWriter::rewrite(const std::string &path,
                               const std::string &collection,
                               const std::string &grids)
{
  close();
  path_ = path;
  fp_ = fopen(path.c_str(),"w");
  if (!fp_) return false;
  fprintf(fp_,"<?xml version=\"1.0\" ?>\n");
  fprintf(fp_,"<Xdmf Version=\"3.0\">\n");
  fprintf(fp_,"  <Domain>\n");
  fprintf(fp_,"    <Grid Name=\"%s\" GridType=\"Collection\" CollectionType=\"Temporal\">\n",
          collection.c_str());
  if (!grids.empty()) fwrite(grids.data(),1,grids.size(),fp_);
  body_end_ = static_cast<long long>(ftello(fp_));
  fputs(XDMF_TAIL,fp_);
  if (fflush(fp_) != 0 || ferror(fp_) || body_end_ < 0) return false;
  return true;
}

bool XdmfSeriesWriter::append(const std::string &grid)
{
  if (!fp_) return false;
  // overwrite the closing tags; grid + tail is longer than tail, so no
  // stale bytes can remain after the new end of file
  if (fseeko(fp_,static_cast<off_t>(body_end_),SEEK_SET) != 0) return false;
  fwrite(grid.data(),1,grid.size(),fp_);
  body_end_ = static_cast<long long>(ftello(fp_));
  fputs(XDMF_TAIL,fp_);
  if (fflush(fp_) != 0 || ferror(fp_) || body_end_ < 0) return false;
  return true;
}

std::string DumpHDF5Util::expand_star(const char *pattern, bigint step, int pad)
{
  char num[64];
  if (pad > 0) {
    char fmt[32];
    // BIGINT_FORMAT is "%ld" or "%lld": insert the zero padding after '%'
    snprintf(fmt,sizeof(fmt),"%%0%d%s",pad,&BIGINT_FORMAT[1]);
    snprintf(num,sizeof(num),fmt,step);
  } else {
    snprintf(num,sizeof(num),BIGINT_FORMAT,step);
  }
  return expand_star_text(pattern,num);
}

std::string DumpHDF5Util::expand_star_text(const char *pattern, const char *text)
{
  std::string s(pattern);
  const std::string::size_type pos = s.find('*');
  if (pos == std::string::npos) return s;
  return s.substr(0,pos) + text + s.substr(pos+1);
}

std::string DumpHDF5Util::basename_of(const std::string &path)
{
  const std::string::size_type pos = path.rfind('/');
  return pos == std::string::npos ? path : path.substr(pos+1);
}

bool DumpHDF5Util::file_exists(const std::string &path)
{
  struct stat sb;
  return stat(path.c_str(),&sb) == 0;
}

int DumpHDF5Util::insert_entry(std::vector<StepEntry> &entries, const StepEntry &e)
{
  if (entries.empty() || e.step > entries.back().step) {
    entries.push_back(e);
    return 1;
  }
  size_t i = 0;
  while (i < entries.size() && entries[i].step < e.step) ++i;
  if (i < entries.size() && entries[i].step == e.step) {
    entries[i] = e;
    return 2;
  }
  entries.insert(entries.begin() + i, e);
  return 0;
}

double DumpHDF5Util::continuity_offset(const std::vector<StepEntry> &entries, bigint step,
                                      double time_now, double dt)
{
  const StepEntry *ref = NULL;
  for (size_t i = 0; i < entries.size(); ++i)
    if (entries[i].step <= step) ref = &entries[i];
  if (!ref) return 0.0;
  const double expected = ref->time + static_cast<double>(step - ref->step) * dt;
  const double diff = expected - time_now;
  return (fabs(diff) > 1e-3 * fabs(dt)) ? diff : 0.0;
}

#ifdef LIGGGHTS_HDF5

void H5Id::reset(hid_t id, closer_t closer)
{
  if (id_ >= 0 && closer_) closer_(id_);
  id_ = id;
  closer_ = closer;
}

herr_t H5Id::close()
{
  herr_t rv = 0;
  if (id_ >= 0 && closer_) rv = closer_(id_);
  id_ = -1;
  closer_ = 0;
  return rv;
}

/* ---------------------------------------------------------------------- */

namespace {

struct WalkData {
  std::string text;
  int n;
};

herr_t walk_cb(unsigned, const H5E_error2_t *err, void *data)
{
  WalkData *w = static_cast<WalkData *>(data);
  if (w->n >= 3) return 0;
  if (err && err->desc && err->desc[0]) {
    if (!w->text.empty()) w->text += " <- ";
    if (err->func_name) {
      w->text += err->func_name;
      w->text += "(): ";
    }
    w->text += err->desc;
    w->n++;
  }
  return 0;
}

herr_t collect_step_groups(hid_t, const char *name, const H5L_info_t *, void *data)
{
  if (strncmp(name,"Step_",5) == 0)
    static_cast<std::vector<std::string> *>(data)->push_back(std::string(name));
  return 0;
}

}

Status::Status(Error *error, MPI_Comm world, const char *style) :
  error_(error), world_(world), style_(style), failed_(false)
{
}

void Status::fail(const char *what)
{
  if (failed_) return;   // keep the first failure
  failed_ = true;
  WalkData w;
  w.n = 0;
  H5Ewalk2(H5E_DEFAULT,H5E_WALK_UPWARD,walk_cb,&w);
  H5Eclear2(H5E_DEFAULT);
  msg_ = style_ + ": " + what + " failed";
  if (!w.text.empty()) msg_ += " (HDF5: " + w.text + ")";
}

void Status::sync(const char *file, int line)
{
  int me, nprocs;
  MPI_Comm_rank(world_,&me);
  MPI_Comm_size(world_,&nprocs);
  int mine = failed_ ? me : nprocs;
  int first = nprocs;
  MPI_Allreduce(&mine,&first,1,MPI_INT,MPI_MIN,world_);
  if (first == nprocs) return;
  char buf[1024];
  buf[0] = '\0';
  if (me == first) snprintf(buf,sizeof(buf),"%s",msg_.c_str());
  MPI_Bcast(buf,sizeof(buf),MPI_CHAR,first,world_);
  char full[1100];
  snprintf(full,sizeof(full),"%s [first failing rank %d]",buf,first);
  error_->all(file,line,full);
}

/* ---------------------------------------------------------------------- */

void DumpHDF5Util::write_dataset(Status &st, hid_t loc, const char *name, hid_t type,
                                 int rank, hsize_t rows_global, hsize_t rows_local,
                                 hsize_t row_offset, hsize_t ncol, hid_t dxpl,
                                 const void *data)
{
  const hsize_t dims[2] = {rows_global,ncol};
  const hsize_t count[2] = {rows_local,ncol};
  const hsize_t offset[2] = {row_offset,0};
  const hsize_t one[2] = {1,1};
  const hsize_t max_chunk_rows = 1048576;
  const hsize_t chunk_rows = rows_global == 0 ? 1 :
    (rows_global < max_chunk_rows ? rows_global : max_chunk_rows);
  const hsize_t chunk[2] = {chunk_rows,ncol};
  const bool empty = (rows_local == 0);
  const long long dummy[3] = {0,0,0};

  std::string what;

  // local (non-collective) preparation
  H5Id filespace(H5Screate_simple(rank,dims,NULL),H5Sclose);
  what = std::string("H5Screate_simple(file) for ") + name;
  st.check(filespace.get(),what.c_str());
  H5Id memspace(H5Screate_simple(rank,empty ? one : count,NULL),H5Sclose);
  what = std::string("H5Screate_simple(memory) for ") + name;
  st.check(memspace.get(),what.c_str());
  if (st.ok()) {
    herr_t rv;
    if (empty) {
      rv = H5Sselect_none(filespace.get());
      if (rv >= 0) rv = H5Sselect_none(memspace.get());
    } else {
      rv = H5Sselect_hyperslab(filespace.get(),H5S_SELECT_SET,offset,NULL,count,NULL);
    }
    what = std::string("H5Sselect for ") + name;
    st.check(rv,what.c_str());
  }
  H5Id dcpl(H5Pcreate(H5P_DATASET_CREATE),H5Pclose);
  what = std::string("H5Pcreate(dataset create) for ") + name;
  st.check(dcpl.get(),what.c_str());
  if (st.ok()) {
    what = std::string("H5Pset_chunk for ") + name;
    st.check(H5Pset_chunk(dcpl.get(),rank,chunk),what.c_str());
  }
  st.sync(FLERR);

  // collective creation and write
  H5Id dset(H5Dcreate2(loc,name,type,filespace.get(),H5P_DEFAULT,dcpl.get(),H5P_DEFAULT),H5Dclose);
  what = std::string("H5Dcreate ") + name;
  if (st.check(dset.get(),what.c_str())) {
    what = std::string("H5Dwrite ") + name;
    st.check(H5Dwrite(dset.get(),type,memspace.get(),filespace.get(),dxpl,
                      (empty || !data) ? static_cast<const void *>(dummy) : data),
             what.c_str());
  }
  what = std::string("H5Dclose ") + name;
  st.check(dset.close(),what.c_str());
  st.sync(FLERR);
}

void DumpHDF5Util::write_scalar_attribute(Status &st, hid_t loc, const char *name,
                                          hid_t type, const void *value)
{
  const hsize_t dims[1] = {1};
  std::string what = std::string("H5Acreate ") + name;
  H5Id space(H5Screate_simple(1,dims,NULL),H5Sclose);
  st.check(space.get(),what.c_str());
  st.sync(FLERR);
  H5Id attr(H5Acreate2(loc,name,type,space.get(),H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
  if (st.check(attr.get(),what.c_str())) {
    what = std::string("H5Awrite ") + name;
    st.check(H5Awrite(attr.get(),type,value),what.c_str());
  }
  what = std::string("H5Aclose ") + name;
  st.check(attr.close(),what.c_str());
  st.sync(FLERR);
}

/* ----------------------------------------------------------------------
   list Step_<n> groups of an open file. Returns false on an HDF5 error.
   time_base/step_base/dt: extrapolation for legacy groups without a
   'time' attribute (time = time_base + (step - step_base)*dt).
------------------------------------------------------------------------- */

static bool read_steps_impl(hid_t file, const std::string &h5ref,
                            const char *count_dataset, double time_base,
                            bigint step_base, double dt,
                            std::vector<StepEntry> &entries, int &nlegacy)
{
  std::vector<std::string> names;
  hsize_t idx = 0;
  if (H5Literate(file,H5_INDEX_NAME,H5_ITER_NATIVE,&idx,collect_step_groups,&names) < 0)
    return false;

  for (size_t i = 0; i < names.size(); ++i) {
    H5Id group(H5Gopen2(file,names[i].c_str(),H5P_DEFAULT),H5Gclose);
    if (!group.valid()) return false;

    StepEntry e;
    e.h5ref = h5ref;
    long long step_ll = atoll(names[i].c_str() + 5);
    if (H5Aexists(group.get(),"timestep") > 0) {
      H5Id attr(H5Aopen(group.get(),"timestep",H5P_DEFAULT),H5Aclose);
      if (!attr.valid() || H5Aread(attr.get(),H5T_NATIVE_LLONG,&step_ll) < 0) return false;
    }
    e.step = static_cast<bigint>(step_ll);

    double t = 0.0;
    if (H5Aexists(group.get(),"time") > 0) {
      H5Id attr(H5Aopen(group.get(),"time",H5P_DEFAULT),H5Aclose);
      if (!attr.valid() || H5Aread(attr.get(),H5T_NATIVE_DOUBLE,&t) < 0) return false;
    } else {
      t = time_base + static_cast<double>(e.step - step_base) * dt;
      nlegacy++;
    }
    e.time = t;

    H5Id dset(H5Dopen2(group.get(),count_dataset,H5P_DEFAULT),H5Dclose);
    if (!dset.valid()) return false;
    H5Id space(H5Dget_space(dset.get()),H5Sclose);
    if (!space.valid()) return false;
    hsize_t dims[2] = {0,0};
    if (H5Sget_simple_extent_dims(space.get(),dims,NULL) < 0) return false;
    e.count = static_cast<long long>(dims[0]);
    e.has_type = H5Lexists(group.get(),"type",H5P_DEFAULT) > 0;

    insert_entry(entries,e);
  }
  return true;
}

bool DumpHDF5Util::read_existing_steps(hid_t file, const std::string &h5ref,
                                       const char *count_dataset,
                                       double time_base, bigint step_base, double dt,
                                       std::vector<StepEntry> &entries, int &nlegacy)
{
  return read_steps_impl(file,h5ref,count_dataset,time_base,step_base,dt,entries,nlegacy);
}

/* ----------------------------------------------------------------------
   serial scan (call on one rank) of existing files matching a '*' pattern
   (the '*' must be in the file name, not in a directory component).
   Files that cannot be read are skipped and counted in nbad.
------------------------------------------------------------------------- */

void DumpHDF5Util::scan_multifile_steps(const char *pattern, const char *count_dataset,
                                        double time_base, bigint step_base, double dt,
                                        std::vector<StepEntry> &entries,
                                        int &nlegacy, int &nbad)
{
  const std::string pat(pattern);
  const std::string::size_type slash = pat.rfind('/');
  const std::string dir = (slash == std::string::npos) ? std::string(".") : pat.substr(0,slash);
  const std::string base = (slash == std::string::npos) ? pat : pat.substr(slash+1);
  const std::string::size_type star = base.find('*');
  if (star == std::string::npos) return;
  const std::string prefix = base.substr(0,star);
  const std::string suffix = base.substr(star+1);

  DIR *d = opendir(dir.c_str());
  if (!d) return;
  struct dirent *de;
  std::vector<std::string> matches;
  while ((de = readdir(d)) != NULL) {
    const std::string name(de->d_name);
    if (name.size() <= prefix.size() + suffix.size()) continue;
    if (name.compare(0,prefix.size(),prefix) != 0) continue;
    if (name.compare(name.size()-suffix.size(),suffix.size(),suffix) != 0) continue;
    const std::string mid = name.substr(prefix.size(),name.size()-prefix.size()-suffix.size());
    bool digits = !mid.empty();
    for (size_t i = 0; i < mid.size(); ++i)
      if (mid[i] < '0' || mid[i] > '9') digits = false;
    if (digits) matches.push_back(name);
  }
  closedir(d);
  std::sort(matches.begin(),matches.end());

  for (size_t i = 0; i < matches.size(); ++i) {
    const std::string path = (slash == std::string::npos) ? matches[i] : dir + "/" + matches[i];
    H5Id f(H5Fopen(path.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose);
    if (!f.valid() ||
        !read_steps_impl(f.get(),matches[i],count_dataset,time_base,step_base,dt,entries,nlegacy)) {
      H5Eclear2(H5E_DEFAULT);
      nbad++;
    }
  }
}

/* ---------------------------------------------------------------------- */

hid_t DumpHDF5Util::open_parallel_file(Status &st, const std::string &name, bool append,
                                       MPI_Comm world, bool &existed)
{
  int me;
  MPI_Comm_rank(world,&me);
  int ex = 0;
  if (me == 0) ex = file_exists(name) ? 1 : 0;
  MPI_Bcast(&ex,1,MPI_INT,0,world);
  existed = (ex != 0);

  H5Id fapl(H5Pcreate(H5P_FILE_ACCESS),H5Pclose);
  if (st.check(fapl.get(),"H5Pcreate(file access)"))
    st.check(H5Pset_fapl_mpio(fapl.get(),world,MPI_INFO_NULL),"H5Pset_fapl_mpio");
  st.sync(FLERR);

  hid_t file;
  std::string what;
  if (append && existed) {
    file = H5Fopen(name.c_str(),H5F_ACC_RDWR,fapl.get());
    what = "H5Fopen(append) of '" + name + "'";
  } else {
    file = H5Fcreate(name.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,fapl.get());
    what = "H5Fcreate of '" + name + "'";
  }
  st.check(file,what.c_str());
  return file;
}

hid_t DumpHDF5Util::reopen_parallel_file(Status &st, const std::string &name, MPI_Comm world)
{
  H5Id fapl(H5Pcreate(H5P_FILE_ACCESS),H5Pclose);
  if (st.check(fapl.get(),"H5Pcreate(file access)"))
    st.check(H5Pset_fapl_mpio(fapl.get(),world,MPI_INFO_NULL),"H5Pset_fapl_mpio");
  st.sync(FLERR);
  const hid_t file = H5Fopen(name.c_str(),H5F_ACC_RDWR,fapl.get());
  const std::string what = "H5Fopen (re-open for the next dump) of '" + name + "'";
  st.check(file,what.c_str());
  return file;
}

void DumpHDF5Util::create_step_group(Status &st, hid_t file, bigint step,
                                     H5Id &group, bool &replaced)
{
  char name[64];
  snprintf(name,sizeof(name),"Step_" BIGINT_FORMAT,step);
  replaced = false;
  const htri_t exists = H5Lexists(file,name,H5P_DEFAULT);
  st.check(exists,"H5Lexists(step group)");
  st.sync(FLERR);
  if (exists > 0) {
    replaced = true;
    st.check(H5Ldelete(file,name,H5P_DEFAULT),"H5Ldelete(existing step group)");
    st.sync(FLERR);
  }
  group.reset(H5Gcreate2(file,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
  st.check(group.get(),"H5Gcreate(step group)");
  st.sync(FLERR);
}

void DumpHDF5Util::write_step_attributes(Status &st, hid_t group, bigint step, double time)
{
  const long long step_ll = static_cast<long long>(step);
  write_scalar_attribute(st,group,"timestep",H5T_NATIVE_LLONG,&step_ll);
  write_scalar_attribute(st,group,"time",H5T_NATIVE_DOUBLE,&time);
}

void DumpHDF5Util::finish_step(Status &st, H5Id &group, H5Id &file, bool close_file)
{
  st.check(group.close(),"H5Gclose(step group)");
  st.sync(FLERR);
  if (close_file) {
    st.check(file.close(),"H5Fclose");
    st.sync(FLERR);
  }
}

#endif // LIGGGHTS_HDF5

/* ----------------------------------------------------------------------
   XDMF DataItem with explicit number type (F-22): readers default to
   Float/Precision 4 when NumberType/Precision are missing.
------------------------------------------------------------------------- */

std::string DumpHDF5Util::data_item(const char *numtype, int precision,
                                    long long rows, int ncol,
                                    const std::string &h5ref, bigint step,
                                    const char *path)
{
  char buf[1024];
  char dims[64];
  if (ncol > 1) snprintf(dims,sizeof(dims),"%lld %d",rows,ncol);
  else snprintf(dims,sizeof(dims),"%lld",rows);
  snprintf(buf,sizeof(buf),
           "<DataItem Format=\"HDF\" NumberType=\"%s\" Precision=\"%d\" Dimensions=\"%s\">%s:/Step_"
           BIGINT_FORMAT "/%s</DataItem>",
           numtype,precision,dims,h5ref.c_str(),step,path);
  return std::string(buf);
}

std::string DumpHDF5Util::grid_open(const StepEntry &e)
{
  char buf[256];
  snprintf(buf,sizeof(buf),
           "      <Grid Name=\"Step_" BIGINT_FORMAT "\" GridType=\"Uniform\">\n"
           "        <Time Value=\"%.17g\" />\n",
           e.step,e.time);
  return std::string(buf);
}

std::string DumpHDF5Util::attribute_xml(const char *name, const char *atype,
                                        const char *center, const std::string &item)
{
  return std::string("        <Attribute Name=\"") + name + "\" AttributeType=\"" + atype +
         "\" Center=\"" + center + "\">\n          " + item + "\n        </Attribute>\n";
}

/* ======================================================================
   DumpHDF5
   ====================================================================== */

DumpHDF5::DumpHDF5(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg),
  truncate_warned_(false),
  opened_once_(false),
  time_offset_(0.0),
  time_offset_set_(false)
{
  if (narg != 5)
    error->all(FLERR,"Illegal dump hdf5 command: expected 'dump ID group hdf5 N file.h5'");
  if (strchr(filename,'%'))
    error->all(FLERR,"Dump hdf5 does not support '%' in the file name "
               "(all ranks write one file collectively)");
  if (compressed)
    error->all(FLERR,"Dump hdf5 does not support gzip (.gz) file names");

  binary = 1;
  buffer_allow = 0;
  buffer_flag = 0;
  size_one = 15; // id + type + xyz + vxyz + fxyz + omegaxyz + radius

#ifdef LIGGGHTS_HDF5
  // errors are reported through error->all() with the HDF5 error stack
  // text (Status::fail), so switch off HDF5's own stderr printing
  H5Eset_auto2(H5E_DEFAULT,NULL,NULL);
#endif
}

/* ---------------------------------------------------------------------- */

DumpHDF5::~DumpHDF5()
{
#ifdef LIGGGHTS_HDF5
  // collective: all ranks destroy their dumps in the same order
  if (file_.valid() && file_.close() < 0 && me == 0)
    error->warning(FLERR,"Dump hdf5: H5Fclose failed when closing the dump file");
#endif
  xdmf_.close();
}

/* ---------------------------------------------------------------------- */

void DumpHDF5::init_style()
{
#ifndef LIGGGHTS_HDF5
  error->all(FLERR,"Dump hdf5 requires a parallel HDF5 build with -DLIGGGHTS_HDF5");
#else
  if (!atom->tag_enable)
    error->all(FLERR,"Dump hdf5 requires atom IDs");
  if (!atom->radius_flag)
    error->all(FLERR,"Dump hdf5 requires per-atom radius");
  if (!atom->omega_flag)
    error->all(FLERR,"Dump hdf5 requires per-atom omega");
  if (format_user)
    error->all(FLERR,"dump_modify format is not supported by dump hdf5 "
               "(data are stored in binary double/integer datasets)");
#endif
}

/* ----------------------------------------------------------------------
   dump_modify keywords not handled by Dump::modify_params() end up here.
   Reject them instead of letting SortBuffer accept 'sort' silently.
------------------------------------------------------------------------- */

int DumpHDF5::modify_param(int /*narg*/, char **arg)
{
  char msg[512];
  snprintf(msg,sizeof(msg),
           "dump_modify keyword '%s' is not supported by dump %s "
           "(supported: append, every, first, flush, pad)",arg[0],style);
  error->all(FLERR,msg);
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpHDF5::pack_local(long long *ids, int *types, double *positions,
                          double *velocities, double *forces, double *omegas,
                          double *radii, int nlocal_selected) const
{
  tagint *tag = atom->tag;
  int *type = atom->type;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double *radius = atom->radius;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  size_t m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    ids[m] = static_cast<long long>(tag[i]);
    types[m] = type[i];
    for (int k = 0; k < 3; ++k) {
      positions[3*m+k] = x[i][k];
      velocities[3*m+k] = v[i][k];
      forces[3*m+k] = f[i][k];
      omegas[3*m+k] = omega[i][k];
    }
    radii[m] = radius[i];
    m++;
  }

  if (m != static_cast<size_t>(nlocal_selected))
    error->one(FLERR,"Internal dump hdf5 packing error");
}

/* ---------------------------------------------------------------------- */

std::string DumpHDF5::grid_xml(const StepEntry &e) const
{
  const long long n = e.count;
  std::string s = grid_open(e);
  char buf[128];
  snprintf(buf,sizeof(buf),
           "        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%lld\" />\n",n);
  s += buf;
  s += "        <Geometry GeometryType=\"XYZ\">\n          ";
  s += data_item("Float",8,n,3,e.h5ref,e.step,"position");
  s += "\n        </Geometry>\n";
  s += attribute_xml("id","Scalar","Node",data_item("Int",8,n,1,e.h5ref,e.step,"id"));
  if (e.has_type)
    s += attribute_xml("type","Scalar","Node",data_item("Int",4,n,1,e.h5ref,e.step,"type"));
  s += attribute_xml("velocity","Vector","Node",data_item("Float",8,n,3,e.h5ref,e.step,"velocity"));
  s += attribute_xml("force","Vector","Node",data_item("Float",8,n,3,e.h5ref,e.step,"force"));
  s += attribute_xml("omega","Vector","Node",data_item("Float",8,n,3,e.h5ref,e.step,"omega"));
  s += attribute_xml("radius","Scalar","Node",data_item("Float",8,n,1,e.h5ref,e.step,"radius"));
  s += "      </Grid>\n";
  return s;
}

std::string DumpHDF5::all_grids_xml() const
{
  std::string s;
  for (size_t i = 0; i < entries_.size(); ++i) s += grid_xml(entries_[i]);
  return s;
}

/* ----------------------------------------------------------------------
   rank 0: bring the XDMF file(s) up to date after entry e was recorded
   where = return value of insert_entry()
------------------------------------------------------------------------- */

void DumpHDF5::update_xdmf(const StepEntry &e, int where)
{
  bool ok;
  if (where == 1 && xdmf_.is_open()) ok = xdmf_.append(grid_xml(e));
  else ok = xdmf_.rewrite(xdmf_.path(),"LIGGGHTS",all_grids_xml());
  if (!ok) {
    std::string msg = "Dump hdf5: cannot write XDMF file '" + xdmf_.path() + "'";
    error->one(FLERR,msg.c_str());
  }

  if (multifile) {
    // self-contained XDMF next to each per-step file (lists only its step)
    XdmfSeriesWriter one;
    const std::string path = open_name_ + ".xdmf";
    if (!one.rewrite(path,"LIGGGHTS",grid_xml(e))) {
      std::string msg = "Dump hdf5: cannot write XDMF file '" + path + "'";
      error->one(FLERR,msg.c_str());
    }
    one.close();
  }
}

/* ---------------------------------------------------------------------- */

void DumpHDF5::write()
{
#ifndef LIGGGHTS_HDF5
  error->all(FLERR,"Dump hdf5 requires a parallel HDF5 build with -DLIGGGHTS_HDF5");
#else
  Status st(error,world,"dump hdf5");

  const int nlocal_selected = count();
  long long local_count_ll = nlocal_selected;
  long long global_count_ll = 0;
  long long offset_ll = 0;

  MPI_Allreduce(&local_count_ll,&global_count_ll,1,MPI_LONG_LONG,MPI_SUM,world);
  MPI_Exscan(&local_count_ll,&offset_ll,1,MPI_LONG_LONG,MPI_SUM,world);
  if (me == 0) offset_ll = 0;

  const size_t n = static_cast<size_t>(nlocal_selected);
  std::vector<long long> ids(n);
  std::vector<int> types(n);
  std::vector<double> positions(3*n);
  std::vector<double> velocities(3*n);
  std::vector<double> forces(3*n);
  std::vector<double> omegas(3*n);
  std::vector<double> radii(n);
  if (n > 0)
    pack_local(&ids[0],&types[0],&positions[0],&velocities[0],&forces[0],
               &omegas[0],&radii[0],nlocal_selected);

  const bigint step = update->ntimestep;
  // elapsed simulation time as LIGGGHTS computes it (thermo keyword 'time')
  const double time_now = update->atime + (update->ntimestep - update->atimestep)*update->dt;

  // ---- open the file

  H5Id local_file;                // multifile: one file per dump
  H5Id *file = &file_;
  if (multifile) {
    open_name_ = expand_star(filename,step,padflag);
    bool existed;
    local_file.reset(open_parallel_file(st,open_name_,false,world,existed),H5Fclose);
    st.sync(FLERR);
    file = &local_file;
    if (me == 0 && !xdmf_.is_open()) {
      // master XDMF of the series (created on the first dump)
      const std::string master = expand_star_text(filename,"series") + ".xdmf";
      int nlegacy = 0, nbad = 0;
      if (append_flag)
        scan_multifile_steps(filename,"id",0.0,0,update->dt,
                             entries_,nlegacy,nbad);
      if (nbad) error->warning(FLERR,"Dump hdf5 append: some existing files matching "
                               "the '*' pattern could not be read and are not listed in the XDMF");
      if (nlegacy) error->warning(FLERR,"Dump hdf5 append: existing steps without a 'time' "
                                  "attribute get time = step*dt (current dt)");
      if (!xdmf_.rewrite(master,"LIGGGHTS",all_grids_xml()))
        error->one(FLERR,"Dump hdf5: cannot write the XDMF series file");
    }
  } else if (!file_.valid() && opened_once_) {
    // 'dump_modify flush yes' (default): the file was closed after the
    // previous dump so that it is complete on disk between dumps
    file_.reset(reopen_parallel_file(st,open_name_,world),H5Fclose);
    st.sync(FLERR);
  } else if (!file_.valid()) {
    opened_once_ = true;
    open_name_ = filename;
    bool existed;
    file_.reset(open_parallel_file(st,open_name_,append_flag != 0,world,existed),H5Fclose);
    st.sync(FLERR);
    if (existed && !append_flag && me == 0 && !truncate_warned_) {
      std::string msg = "Dump hdf5: existing file '" + open_name_ +
        "' is truncated (use 'dump_modify <ID> append yes' to keep its steps)";
      error->warning(FLERR,msg.c_str());
      truncate_warned_ = true;
    }
    entries_.clear();
    if (existed && append_flag) {
      int nlegacy = 0;
      st.check(read_existing_steps(file_.get(),basename_of(open_name_),"id",0.0,
                                   0,update->dt,entries_,nlegacy) ? 0 : -1,
               "reading the existing Step_ groups for append");
      st.sync(FLERR);
      if (nlegacy && me == 0)
        error->warning(FLERR,"Dump hdf5 append: existing steps without a 'time' "
                       "attribute get time = step*dt (current dt)");
    }
    if (me == 0 && !xdmf_.rewrite(open_name_ + ".xdmf","LIGGGHTS",all_grids_xml())) {
      std::string msg = "Dump hdf5: cannot write XDMF file '" + open_name_ + ".xdmf'";
      error->one(FLERR,msg.c_str());
    }
  }

  // ---- append: keep the time axis continuous across a restart
  //      (computed once, when the file is first opened)
  if (append_flag && !time_offset_set_) {
    double off = 0.0;
    if (me == 0) off = continuity_offset(entries_,step,time_now,update->dt);
    MPI_Bcast(&off,1,MPI_DOUBLE,0,world);
    time_offset_ = off;
    time_offset_set_ = true;
    if (off != 0.0 && me == 0) {
      char msg[512];
      snprintf(msg,sizeof(msg),"Dump hdf5 append: elapsed time restarted after "
               "read_restart; the stored time is continued from the steps already in the "
               "file (offset %.17g)",off);
      error->warning(FLERR,msg);
    }
  }
  const double time = time_now + time_offset_;

  // ---- step group

  H5Id group;
  bool replaced = false;
  create_step_group(st,file->get(),step,group,replaced);
  if (replaced && me == 0) {
    char msg[256];
    snprintf(msg,sizeof(msg),"Dump hdf5: step group Step_" BIGINT_FORMAT
             " already exists in the file and is replaced",step);
    error->warning(FLERR,msg);
  }

  // ---- datasets

  H5Id dxpl(H5Pcreate(H5P_DATASET_XFER),H5Pclose);
  if (st.check(dxpl.get(),"H5Pcreate(dataset transfer)"))
    st.check(H5Pset_dxpl_mpio(dxpl.get(),H5FD_MPIO_COLLECTIVE),"H5Pset_dxpl_mpio");
  st.sync(FLERR);

  const hsize_t rows = static_cast<hsize_t>(global_count_ll);
  const hsize_t lrows = static_cast<hsize_t>(nlocal_selected);
  const hsize_t off = static_cast<hsize_t>(offset_ll);
  const hid_t g = group.get();
  const hid_t x = dxpl.get();

  write_dataset(st,g,"id",H5T_NATIVE_LLONG,1,rows,lrows,off,1,x,n ? &ids[0] : NULL);
  write_dataset(st,g,"type",H5T_NATIVE_INT,1,rows,lrows,off,1,x,n ? &types[0] : NULL);
  write_dataset(st,g,"position",H5T_NATIVE_DOUBLE,2,rows,lrows,off,3,x,n ? &positions[0] : NULL);
  write_dataset(st,g,"velocity",H5T_NATIVE_DOUBLE,2,rows,lrows,off,3,x,n ? &velocities[0] : NULL);
  write_dataset(st,g,"force",H5T_NATIVE_DOUBLE,2,rows,lrows,off,3,x,n ? &forces[0] : NULL);
  write_dataset(st,g,"omega",H5T_NATIVE_DOUBLE,2,rows,lrows,off,3,x,n ? &omegas[0] : NULL);
  write_dataset(st,g,"radius",H5T_NATIVE_DOUBLE,1,rows,lrows,off,1,x,n ? &radii[0] : NULL);
  write_step_attributes(st,g,step,time);

  st.check(dxpl.close(),"H5Pclose(dataset transfer)");
  // close after the dump unless 'dump_modify flush no' (then kept open until
  // the dump is destroyed); multifile: always close
  finish_step(st,group,*file,multifile != 0 || flush_flag != 0);

  // ---- XDMF (rank 0)

  StepEntry e;
  e.step = step;
  e.time = time;
  e.count = global_count_ll;
  e.h5ref = basename_of(open_name_);
  e.has_type = true;
  const int where = insert_entry(entries_,e);
  if (me == 0) update_xdmf(e,where);
#endif
}
