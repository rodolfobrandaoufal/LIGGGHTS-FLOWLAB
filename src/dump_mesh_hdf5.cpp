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

    Parallel HDF5 triangular mesh dump (dump mesh/hdf5). Uses the HDF5 /
    XDMF helpers of dump_hdf5.cpp; see the MPI note there and
    doc/dump_hdf5.txt.
------------------------------------------------------------------------- */

#include "dump_mesh_hdf5.h"
#include "error.h"
#include "update.h"
#include "modify.h"
#include "fix_mesh_surface.h"
#include "tri_mesh.h"
#include "container.h"
#include "comm.h"
#include <mpi.h>
#include <stdio.h>
#include <string.h>

#include <limits.h>

using namespace LAMMPS_NS;
using namespace DumpHDF5Util;

/* ---------------------------------------------------------------------- */

DumpMeshHDF5::DumpMeshHDF5(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg),
  truncate_warned_(false),
  opened_once_(false),
  time_offset_(0.0),
  time_offset_set_(false),
  dump_all_meshes_(0),
  builtin_mask_(0)
{
  if (narg < 5)
    error->all(FLERR,"Illegal dump mesh/hdf5 command");
  if (strchr(filename,'%'))
    error->all(FLERR,"Dump mesh/hdf5 does not support '%' in the file name "
               "(all ranks write one file collectively)");
  if (compressed)
    error->all(FLERR,"Dump mesh/hdf5 does not support gzip (.gz) file names");
#ifdef LIGGGHTS_HDF5
  H5Eset_auto2(H5E_DEFAULT,NULL,NULL);
#endif

  binary = 1;
  buffer_allow = 0;
  buffer_flag = 0;
  format = NULL;
  format_default = NULL;
  format_user = NULL;
  size_one = 0;

  int iarg = 5;
  bool parsing_meshes = true;

  while (iarg < narg) {
    if (parsing_meshes) {
      if (strcmp(arg[iarg],"all") == 0) {
        dump_all_meshes_ = 1;
        ++iarg;
        continue;
      }

      const int ifix = modify->find_fix(arg[iarg]);
      FixMeshSurface *fms = (ifix < 0) ? 0 : dynamic_cast<FixMeshSurface*>(modify->fix[ifix]);
      if (fms) {
        if (dump_all_meshes_)
          error->all(FLERR,"dump mesh/hdf5 cannot combine 'all' with explicit mesh ids");
        fixes_.push_back(fms);
        ++iarg;
        continue;
      }

      parsing_meshes = false;
    }

    parse_property(arg[iarg++]);
  }

  if (!dump_all_meshes_ && fixes_.empty())
    dump_all_meshes_ = 1;
}

/* ---------------------------------------------------------------------- */

DumpMeshHDF5::~DumpMeshHDF5()
{
#ifdef LIGGGHTS_HDF5
  if (file_.valid() && file_.close() < 0 && me == 0)
    error->warning(FLERR,"Dump mesh/hdf5: H5Fclose failed when closing the dump file");
#endif
  xdmf_.close();
}

/* ---------------------------------------------------------------------- */

int DumpMeshHDF5::modify_param(int /*narg*/, char **arg)
{
  char msg[512];
  snprintf(msg,sizeof(msg),
           "dump_modify keyword '%s' is not supported by dump %s "
           "(supported: append, every, first, flush, pad)",arg[0],style);
  error->all(FLERR,msg);
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpMeshHDF5::parse_property(const char *name)
{
  if (strcmp(name,"id") == 0)
    builtin_mask_ |= PROP_ID;
  else if (strcmp(name,"owner") == 0)
    builtin_mask_ |= PROP_OWNER;
  else if (strcmp(name,"area") == 0)
    builtin_mask_ |= PROP_AREA;
  else if (strcmp(name,"normal") == 0)
    builtin_mask_ |= PROP_NORMAL;
  else if (strcmp(name,"aedges") == 0)
    builtin_mask_ |= PROP_AEDGES;
  else if (strcmp(name,"acorners") == 0)
    builtin_mask_ |= PROP_ACORNERS;
  else if (strcmp(name,"index") == 0)
    builtin_mask_ |= PROP_INDEX;
  else if (strcmp(name,"nneighs") == 0)
    builtin_mask_ |= PROP_NNEIGHS;
  else
    requested_properties_.push_back(std::string(name));
}

/* ---------------------------------------------------------------------- */

void DumpMeshHDF5::init_style()
{
#ifndef LIGGGHTS_HDF5
  error->all(FLERR,"Dump mesh/hdf5 requires a parallel HDF5 build with -DLIGGGHTS_HDF5");
#else
  if (format_user)
    error->all(FLERR,"dump_modify format is not supported by dump mesh/hdf5");
  if (dump_all_meshes_) {
    fixes_.clear();
    const int nmesh = modify->n_fixes_style("mesh/surface");
    for (int i = 0; i < nmesh; ++i) {
      FixMeshSurface *fms =
        static_cast<FixMeshSurface*>(modify->find_fix_style("mesh/surface",i));
      if (fms) fixes_.push_back(fms);
    }
  }

  if (fixes_.empty())
    error->all(FLERR,"Dump mesh/hdf5 found no fix mesh/surface instances");

  for (size_t i = 0; i < fixes_.size(); ++i)
    if (!fixes_[i]->triMesh())
      error->all(FLERR,"Dump mesh/hdf5 requires triangular mesh/surface fixes");

  scalar_properties_.clear();
  vector_properties_.clear();
  for (size_t ip = 0; ip < requested_properties_.size(); ++ip) {
    bool found_scalar = false;
    bool found_vector = false;
    for (size_t im = 0; im < fixes_.size(); ++im) {
      TriMesh *mesh = fixes_[im]->triMesh();
      if (mesh->prop().getElementProperty<ScalarContainer<double> >(requested_properties_[ip].c_str()))
        found_scalar = true;
      if (mesh->prop().getElementProperty<VectorContainer<double,3> >(requested_properties_[ip].c_str()))
        found_vector = true;
    }
    if (found_scalar) scalar_properties_.push_back(requested_properties_[ip]);
    if (found_vector) vector_properties_.push_back(requested_properties_[ip]);
    if (!found_scalar && !found_vector)
      error->all(FLERR,"Dump mesh/hdf5 requested unknown mesh element property");
  }
#endif
}


/* ---------------------------------------------------------------------- */

std::string DumpMeshHDF5::grid_xml(const StepEntry &e) const
{
  const long long cells = e.count;
  const long long vertices = 3 * cells;
  std::string s = grid_open(e);
  char buf[128];
  snprintf(buf,sizeof(buf),
           "        <Topology TopologyType=\"Triangle\" NumberOfElements=\"%lld\">\n          ",cells);
  s += buf;
  s += data_item("Int",4,cells,3,e.h5ref,e.step,"Topology/Connectivity");
  s += "\n        </Topology>\n        <Geometry GeometryType=\"XYZ\">\n          ";
  s += data_item("Float",8,vertices,3,e.h5ref,e.step,"Geometry/Vertices");
  s += "\n        </Geometry>\n";

  if (builtin_mask_ & PROP_ID)
    s += attribute_xml("id","Scalar","Cell",data_item("Int",4,cells,1,e.h5ref,e.step,"CellData/id"));
  if (builtin_mask_ & PROP_OWNER)
    s += attribute_xml("owner","Scalar","Cell",data_item("Int",4,cells,1,e.h5ref,e.step,"CellData/owner"));
  if (builtin_mask_ & PROP_AREA)
    s += attribute_xml("area","Scalar","Cell",data_item("Float",8,cells,1,e.h5ref,e.step,"CellData/area"));
  if (builtin_mask_ & PROP_NORMAL)
    s += attribute_xml("normal","Vector","Cell",data_item("Float",8,cells,3,e.h5ref,e.step,"CellData/normal"));
  if (builtin_mask_ & PROP_AEDGES)
    s += attribute_xml("active_edges","Scalar","Cell",data_item("Int",4,cells,1,e.h5ref,e.step,"CellData/active_edges"));
  if (builtin_mask_ & PROP_ACORNERS)
    s += attribute_xml("active_corners","Scalar","Cell",data_item("Int",4,cells,1,e.h5ref,e.step,"CellData/active_corners"));
  if (builtin_mask_ & PROP_INDEX)
    s += attribute_xml("index","Scalar","Cell",data_item("Int",4,cells,1,e.h5ref,e.step,"CellData/index"));
  if (builtin_mask_ & PROP_NNEIGHS)
    s += attribute_xml("nneighs","Scalar","Cell",data_item("Int",4,cells,1,e.h5ref,e.step,"CellData/nneighs"));

  for (size_t ip = 0; ip < scalar_properties_.size(); ++ip) {
    const std::string path = "CellData/" + scalar_properties_[ip];
    s += attribute_xml(scalar_properties_[ip].c_str(),"Scalar","Cell",
                       data_item("Float",8,cells,1,e.h5ref,e.step,path.c_str()));
  }
  for (size_t ip = 0; ip < vector_properties_.size(); ++ip) {
    const std::string path = "CellData/" + vector_properties_[ip];
    s += attribute_xml(vector_properties_[ip].c_str(),"Vector","Cell",
                       data_item("Float",8,cells,3,e.h5ref,e.step,path.c_str()));
  }

  s += "      </Grid>\n";
  return s;
}

std::string DumpMeshHDF5::all_grids_xml() const
{
  std::string s;
  for (size_t i = 0; i < entries_.size(); ++i) s += grid_xml(entries_[i]);
  return s;
}

void DumpMeshHDF5::update_xdmf(const StepEntry &e, int where)
{
  bool ok;
  if (where == 1 && xdmf_.is_open()) ok = xdmf_.append(grid_xml(e));
  else ok = xdmf_.rewrite(xdmf_.path(),"LIGGGHTS_Mesh",all_grids_xml());
  if (!ok) {
    std::string msg = "Dump mesh/hdf5: cannot write XDMF file '" + xdmf_.path() + "'";
    error->one(FLERR,msg.c_str());
  }
  if (multifile) {
    XdmfSeriesWriter one;
    const std::string path = open_name_ + ".xdmf";
    if (!one.rewrite(path,"LIGGGHTS_Mesh",grid_xml(e))) {
      std::string msg = "Dump mesh/hdf5: cannot write XDMF file '" + path + "'";
      error->one(FLERR,msg.c_str());
    }
    one.close();
  }
}

/* ---------------------------------------------------------------------- */

void DumpMeshHDF5::write()
{
#ifndef LIGGGHTS_HDF5
  error->all(FLERR,"Dump mesh/hdf5 requires a parallel HDF5 build with -DLIGGGHTS_HDF5");
#else
  Status st(error,world,"dump mesh/hdf5");

  long long local_cells = 0;
  for (size_t im = 0; im < fixes_.size(); ++im) {
    TriMesh *mesh = fixes_[im]->triMesh();
    if (!mesh->isParallel() && comm->me != 0) continue;
    local_cells += mesh->sizeLocal();
  }

  long long global_cells = 0;
  long long cell_offset = 0;
  MPI_Allreduce(&local_cells,&global_cells,1,MPI_LONG_LONG,MPI_SUM,world);
  MPI_Exscan(&local_cells,&cell_offset,1,MPI_LONG_LONG,MPI_SUM,world);
  if (comm->me == 0) cell_offset = 0;

  // connectivity is stored as 32-bit int (layout unchanged)
  if (3 * global_cells > static_cast<long long>(INT_MAX))
    error->all(FLERR,"Dump mesh/hdf5: more than INT_MAX mesh vertices");

  const size_t nc = static_cast<size_t>(local_cells);
  const long long vertex_offset = 3 * cell_offset;

  std::vector<double> vertices(9 * nc);
  std::vector<int> connectivity(3 * nc);

  std::vector<int> id_data, owner_data, aedges_data, acorners_data, index_data, nneighs_data;
  std::vector<double> area_data, normal_data;
  if (builtin_mask_ & PROP_ID) id_data.resize(nc);
  if (builtin_mask_ & PROP_OWNER) owner_data.resize(nc);
  if (builtin_mask_ & PROP_AREA) area_data.resize(nc);
  if (builtin_mask_ & PROP_NORMAL) normal_data.resize(3 * nc);
  if (builtin_mask_ & PROP_AEDGES) aedges_data.resize(nc);
  if (builtin_mask_ & PROP_ACORNERS) acorners_data.resize(nc);
  if (builtin_mask_ & PROP_INDEX) index_data.resize(nc);
  if (builtin_mask_ & PROP_NNEIGHS) nneighs_data.resize(nc);

  std::vector<std::vector<double> > scalar_data(scalar_properties_.size());
  std::vector<std::vector<double> > vector_data(vector_properties_.size());
  for (size_t i = 0; i < scalar_data.size(); ++i)
    scalar_data[i].resize(nc,0.0);
  for (size_t i = 0; i < vector_data.size(); ++i)
    vector_data[i].resize(3 * nc,0.0);

  size_t c = 0;
  for (size_t im = 0; im < fixes_.size(); ++im) {
    TriMesh *mesh = fixes_[im]->triMesh();
    if (!mesh->isParallel() && comm->me != 0) continue;
    const int nlocal = mesh->sizeLocal();

    std::vector<ScalarContainer<double>*> scalar_cont(scalar_properties_.size(),0);
    std::vector<VectorContainer<double,3>*> vector_cont(vector_properties_.size(),0);
    for (size_t ip = 0; ip < scalar_properties_.size(); ++ip) {
      scalar_cont[ip] =
        mesh->prop().getElementProperty<ScalarContainer<double> >(scalar_properties_[ip].c_str());
    }
    for (size_t ip = 0; ip < vector_properties_.size(); ++ip) {
      vector_cont[ip] =
        mesh->prop().getElementProperty<VectorContainer<double,3> >(vector_properties_[ip].c_str());
    }

    for (int itri = 0; itri < nlocal; ++itri, ++c) {
      for (int j = 0; j < 3; ++j) {
        double node[3];
        mesh->node(itri,j,node);
        const size_t v = 3*c + j;
        vertices[3*v+0] = node[0];
        vertices[3*v+1] = node[1];
        vertices[3*v+2] = node[2];
        connectivity[3*c+j] = static_cast<int>(vertex_offset + static_cast<long long>(v));
      }

      if (builtin_mask_ & PROP_ID) id_data[c] = mesh->id(itri);
      if (builtin_mask_ & PROP_OWNER) owner_data[c] = comm->me;
      if (builtin_mask_ & PROP_AREA) area_data[c] = mesh->areaElem(itri);
      if (builtin_mask_ & PROP_NORMAL) {
        double normal[3];
        mesh->surfaceNorm(itri,normal);
        normal_data[3*c+0] = normal[0];
        normal_data[3*c+1] = normal[1];
        normal_data[3*c+2] = normal[2];
      }
      if (builtin_mask_ & PROP_AEDGES) aedges_data[c] = mesh->n_active_edges(itri);
      if (builtin_mask_ & PROP_ACORNERS) acorners_data[c] = mesh->n_active_corners(itri);
      if (builtin_mask_ & PROP_INDEX) index_data[c] = itri;
      if (builtin_mask_ & PROP_NNEIGHS) nneighs_data[c] = mesh->nNeighs(itri);

      for (size_t ip = 0; ip < scalar_cont.size(); ++ip)
        if (scalar_cont[ip]) scalar_data[ip][c] = scalar_cont[ip]->get(itri);
      for (size_t ip = 0; ip < vector_cont.size(); ++ip)
        if (vector_cont[ip]) {
          double value[3];
          vector_cont[ip]->get(itri,value);
          vector_data[ip][3*c+0] = value[0];
          vector_data[ip][3*c+1] = value[1];
          vector_data[ip][3*c+2] = value[2];
        }
    }
  }

  const bigint step = update->ntimestep;
  // elapsed simulation time as LIGGGHTS computes it (thermo keyword 'time')
  const double time_now = update->atime + (update->ntimestep - update->atimestep)*update->dt;

  // ---- open the file (same policy as dump hdf5)

  H5Id local_file;
  H5Id *file = &file_;
  if (multifile) {
    open_name_ = expand_star(filename,step,padflag);
    bool existed;
    local_file.reset(open_parallel_file(st,open_name_,false,world,existed),H5Fclose);
    st.sync(FLERR);
    file = &local_file;
    if (comm->me == 0 && !xdmf_.is_open()) {
      const std::string master = expand_star_text(filename,"series") + ".xdmf";
      int nlegacy = 0, nbad = 0;
      if (append_flag)
        scan_multifile_steps(filename,"Topology/Connectivity",0.0,
                             0,update->dt,entries_,nlegacy,nbad);
      if (nbad) error->warning(FLERR,"Dump mesh/hdf5 append: some existing files matching "
                               "the '*' pattern could not be read and are not listed in the XDMF");
      if (nlegacy) error->warning(FLERR,"Dump mesh/hdf5 append: existing steps without a 'time' "
                                  "attribute get time = step*dt (current dt)");
      if (!xdmf_.rewrite(master,"LIGGGHTS_Mesh",all_grids_xml()))
        error->one(FLERR,"Dump mesh/hdf5: cannot write the XDMF series file");
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
    if (existed && !append_flag && comm->me == 0 && !truncate_warned_) {
      std::string msg = "Dump mesh/hdf5: existing file '" + open_name_ +
        "' is truncated (use 'dump_modify <ID> append yes' to keep its steps)";
      error->warning(FLERR,msg.c_str());
      truncate_warned_ = true;
    }
    entries_.clear();
    if (existed && append_flag) {
      int nlegacy = 0;
      st.check(read_existing_steps(file_.get(),basename_of(open_name_),"Topology/Connectivity",
                                   0.0,0,update->dt,
                                   entries_,nlegacy) ? 0 : -1,
               "reading the existing Step_ groups for append");
      st.sync(FLERR);
      if (nlegacy && comm->me == 0)
        error->warning(FLERR,"Dump mesh/hdf5 append: existing steps without a 'time' "
                       "attribute get time = step*dt (current dt)");
    }
    if (comm->me == 0 && !xdmf_.rewrite(open_name_ + ".xdmf","LIGGGHTS_Mesh",all_grids_xml())) {
      std::string msg = "Dump mesh/hdf5: cannot write XDMF file '" + open_name_ + ".xdmf'";
      error->one(FLERR,msg.c_str());
    }
  }

  // ---- append: keep the time axis continuous across a restart
  //      (computed once, when the file is first opened)
  if (append_flag && !time_offset_set_) {
    double off = 0.0;
    if (comm->me == 0) off = continuity_offset(entries_,step,time_now,update->dt);
    MPI_Bcast(&off,1,MPI_DOUBLE,0,world);
    time_offset_ = off;
    time_offset_set_ = true;
    if (off != 0.0 && comm->me == 0) {
      char msg[512];
      snprintf(msg,sizeof(msg),"Dump mesh/hdf5 append: elapsed time restarted after "
               "read_restart; the stored time is continued from the steps already in the "
               "file (offset %.17g)",off);
      error->warning(FLERR,msg);
    }
  }
  const double time = time_now + time_offset_;

  // ---- step group and sub-groups

  H5Id step_group;
  bool replaced = false;
  create_step_group(st,file->get(),step,step_group,replaced);
  if (replaced && comm->me == 0) {
    char msg[256];
    snprintf(msg,sizeof(msg),"Dump mesh/hdf5: step group Step_" BIGINT_FORMAT
             " already exists in the file and is replaced",step);
    error->warning(FLERR,msg);
  }

  H5Id geometry_group(H5Gcreate2(step_group.get(),"Geometry",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
  st.check(geometry_group.get(),"H5Gcreate(Geometry)");
  H5Id topology_group(H5Gcreate2(step_group.get(),"Topology",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
  st.check(topology_group.get(),"H5Gcreate(Topology)");
  H5Id celldata_group(H5Gcreate2(step_group.get(),"CellData",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
  st.check(celldata_group.get(),"H5Gcreate(CellData)");
  H5Id dxpl(H5Pcreate(H5P_DATASET_XFER),H5Pclose);
  if (st.check(dxpl.get(),"H5Pcreate(dataset transfer)"))
    st.check(H5Pset_dxpl_mpio(dxpl.get(),H5FD_MPIO_COLLECTIVE),"H5Pset_dxpl_mpio");
  st.sync(FLERR);

  const hid_t x = dxpl.get();
  const hsize_t rows = static_cast<hsize_t>(global_cells);
  const hsize_t lrows = static_cast<hsize_t>(local_cells);
  const hsize_t off = static_cast<hsize_t>(cell_offset);

  write_dataset(st,geometry_group.get(),"Vertices",H5T_NATIVE_DOUBLE,2,3*rows,3*lrows,3*off,3,x,
                nc ? &vertices[0] : NULL);
  write_dataset(st,topology_group.get(),"Connectivity",H5T_NATIVE_INT,2,rows,lrows,off,3,x,
                nc ? &connectivity[0] : NULL);

  const hid_t cd = celldata_group.get();
  if (builtin_mask_ & PROP_ID)
    write_dataset(st,cd,"id",H5T_NATIVE_INT,1,rows,lrows,off,1,x,nc ? &id_data[0] : NULL);
  if (builtin_mask_ & PROP_OWNER)
    write_dataset(st,cd,"owner",H5T_NATIVE_INT,1,rows,lrows,off,1,x,nc ? &owner_data[0] : NULL);
  if (builtin_mask_ & PROP_AREA)
    write_dataset(st,cd,"area",H5T_NATIVE_DOUBLE,1,rows,lrows,off,1,x,nc ? &area_data[0] : NULL);
  if (builtin_mask_ & PROP_NORMAL)
    write_dataset(st,cd,"normal",H5T_NATIVE_DOUBLE,2,rows,lrows,off,3,x,nc ? &normal_data[0] : NULL);
  if (builtin_mask_ & PROP_AEDGES)
    write_dataset(st,cd,"active_edges",H5T_NATIVE_INT,1,rows,lrows,off,1,x,nc ? &aedges_data[0] : NULL);
  if (builtin_mask_ & PROP_ACORNERS)
    write_dataset(st,cd,"active_corners",H5T_NATIVE_INT,1,rows,lrows,off,1,x,nc ? &acorners_data[0] : NULL);
  if (builtin_mask_ & PROP_INDEX)
    write_dataset(st,cd,"index",H5T_NATIVE_INT,1,rows,lrows,off,1,x,nc ? &index_data[0] : NULL);
  if (builtin_mask_ & PROP_NNEIGHS)
    write_dataset(st,cd,"nneighs",H5T_NATIVE_INT,1,rows,lrows,off,1,x,nc ? &nneighs_data[0] : NULL);
  for (size_t ip = 0; ip < scalar_properties_.size(); ++ip)
    write_dataset(st,cd,scalar_properties_[ip].c_str(),H5T_NATIVE_DOUBLE,1,rows,lrows,off,1,x,
                  nc ? &scalar_data[ip][0] : NULL);
  for (size_t ip = 0; ip < vector_properties_.size(); ++ip)
    write_dataset(st,cd,vector_properties_[ip].c_str(),H5T_NATIVE_DOUBLE,2,rows,lrows,off,3,x,
                  nc ? &vector_data[ip][0] : NULL);

  write_step_attributes(st,step_group.get(),step,time);

  st.check(dxpl.close(),"H5Pclose(dataset transfer)");
  st.check(celldata_group.close(),"H5Gclose(CellData)");
  st.check(topology_group.close(),"H5Gclose(Topology)");
  st.check(geometry_group.close(),"H5Gclose(Geometry)");
  finish_step(st,step_group,*file,multifile != 0 || flush_flag != 0);

  // ---- XDMF (rank 0)

  StepEntry e;
  e.step = step;
  e.time = time;
  e.count = global_cells;
  e.h5ref = basename_of(open_name_);
  const int where = insert_entry(entries_,e);
  if (comm->me == 0) update_xdmf(e,where);
#endif
}
