/* ----------------------------------------------------------------------
   Parallel HDF5 triangular mesh dump implementation.
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

#ifdef LIGGGHTS_HDF5
#include <hdf5.h>
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpMeshHDF5::DumpMeshHDF5(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg),
  dump_all_meshes_(0),
  builtin_mask_(0)
{
  if (narg < 5)
    error->all(FLERR,"Illegal dump mesh/hdf5 command");

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

void DumpMeshHDF5::resolve_filename(char *out, int nout) const
{
  const char *star = strchr(filename,'*');
  if (!star) {
    snprintf(out,nout,"%s",filename);
    return;
  }

  const int prefix = static_cast<int>(star - filename);
  snprintf(out,nout,"%.*s" BIGINT_FORMAT "%s",
           prefix,filename,update->ntimestep,star+1);
}

/* ---------------------------------------------------------------------- */

#ifdef LIGGGHTS_HDF5
static void create_group_if_needed(hid_t parent, const char *name)
{
  hid_t group = H5Gcreate(parent,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if (group >= 0) H5Gclose(group);
}

static bool empty_selection(int rank, const hsize_t *count)
{
  for (int i = 0; i < rank; ++i)
    if (count[i] == 0) return true;
  return false;
}

static hid_t create_selected_memspace(int rank, const hsize_t *count, bool empty)
{
  if (!empty) return H5Screate_simple(rank,count,NULL);

  hsize_t one[3] = {1,1,1};
  hid_t memspace = H5Screate_simple(rank,one,NULL);
  H5Sselect_none(memspace);
  return memspace;
}

static void write_double_dataset(hid_t group, const char *name, int rank,
                                 const hsize_t *dims, const hsize_t *count,
                                 const hsize_t *offset, const hsize_t *chunk,
                                 hid_t dxpl, const double *data)
{
  const double dummy = 0.0;
  const bool empty = empty_selection(rank,count);
  hid_t filespace = H5Screate_simple(rank,dims,NULL);
  hid_t memspace = create_selected_memspace(rank,count,empty);
  if (empty)
    H5Sselect_none(filespace);
  else
    H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl,rank,chunk);
  hid_t dset = H5Dcreate(group,name,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,dcpl,H5P_DEFAULT);
  H5Dwrite(dset,H5T_NATIVE_DOUBLE,memspace,filespace,dxpl,empty ? &dummy : data);
  H5Dclose(dset);
  H5Pclose(dcpl);
  H5Sclose(memspace);
  H5Sclose(filespace);
}

static void write_int_dataset(hid_t group, const char *name, int rank,
                              const hsize_t *dims, const hsize_t *count,
                              const hsize_t *offset, const hsize_t *chunk,
                              hid_t dxpl, const int *data)
{
  const int dummy = 0;
  const bool empty = empty_selection(rank,count);
  hid_t filespace = H5Screate_simple(rank,dims,NULL);
  hid_t memspace = create_selected_memspace(rank,count,empty);
  if (empty)
    H5Sselect_none(filespace);
  else
    H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl,rank,chunk);
  hid_t dset = H5Dcreate(group,name,H5T_NATIVE_INT,filespace,H5P_DEFAULT,dcpl,H5P_DEFAULT);
  H5Dwrite(dset,H5T_NATIVE_INT,memspace,filespace,dxpl,empty ? &dummy : data);
  H5Dclose(dset);
  H5Pclose(dcpl);
  H5Sclose(memspace);
  H5Sclose(filespace);
}
#endif

/* ---------------------------------------------------------------------- */

void DumpMeshHDF5::write()
{
#ifndef LIGGGHTS_HDF5
  error->all(FLERR,"Dump mesh/hdf5 requires a parallel HDF5 build with -DLIGGGHTS_HDF5");
#else
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

  const long long local_vertices = 3 * local_cells;
  const long long global_vertices = 3 * global_cells;
  const long long vertex_offset = 3 * cell_offset;

  std::vector<double> vertices(static_cast<size_t>(3 * local_vertices));
  std::vector<int> connectivity(static_cast<size_t>(3 * local_cells));

  std::vector<int> id_data, owner_data, aedges_data, acorners_data, index_data, nneighs_data;
  std::vector<double> area_data, normal_data;
  if (builtin_mask_ & PROP_ID) id_data.resize(static_cast<size_t>(local_cells));
  if (builtin_mask_ & PROP_OWNER) owner_data.resize(static_cast<size_t>(local_cells));
  if (builtin_mask_ & PROP_AREA) area_data.resize(static_cast<size_t>(local_cells));
  if (builtin_mask_ & PROP_NORMAL) normal_data.resize(static_cast<size_t>(3 * local_cells));
  if (builtin_mask_ & PROP_AEDGES) aedges_data.resize(static_cast<size_t>(local_cells));
  if (builtin_mask_ & PROP_ACORNERS) acorners_data.resize(static_cast<size_t>(local_cells));
  if (builtin_mask_ & PROP_INDEX) index_data.resize(static_cast<size_t>(local_cells));
  if (builtin_mask_ & PROP_NNEIGHS) nneighs_data.resize(static_cast<size_t>(local_cells));

  std::vector<std::vector<double> > scalar_data(scalar_properties_.size());
  std::vector<std::vector<double> > vector_data(vector_properties_.size());
  for (size_t i = 0; i < scalar_data.size(); ++i)
    scalar_data[i].resize(static_cast<size_t>(local_cells),0.0);
  for (size_t i = 0; i < vector_data.size(); ++i)
    vector_data[i].resize(static_cast<size_t>(3 * local_cells),0.0);

  long long c = 0;
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
        const long long v = 3*c + j;
        vertices[3*v+0] = node[0];
        vertices[3*v+1] = node[1];
        vertices[3*v+2] = node[2];
        connectivity[3*c+j] = static_cast<int>(vertex_offset + v);
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

  char h5name[1024];
  resolve_filename(h5name,sizeof(h5name));

  hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist,world,MPI_INFO_NULL);
  hid_t file = -1;
  if (multifile || written_steps_.empty())
    file = H5Fcreate(h5name,H5F_ACC_TRUNC,H5P_DEFAULT,plist);
  else
    file = H5Fopen(h5name,H5F_ACC_RDWR,plist);
  H5Pclose(plist);
  if (file < 0) error->all(FLERR,"Cannot open parallel HDF5 mesh dump file");

  char step_group_name[128];
  snprintf(step_group_name,sizeof(step_group_name),"/Step_" BIGINT_FORMAT,update->ntimestep);
  hid_t step_group = H5Gcreate(file,step_group_name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if (step_group < 0) error->all(FLERR,"Cannot create HDF5 mesh timestep group");
  create_group_if_needed(step_group,"Geometry");
  create_group_if_needed(step_group,"Topology");
  create_group_if_needed(step_group,"CellData");
  hid_t geometry_group = H5Gopen(step_group,"Geometry",H5P_DEFAULT);
  hid_t topology_group = H5Gopen(step_group,"Topology",H5P_DEFAULT);
  hid_t celldata_group = H5Gopen(step_group,"CellData",H5P_DEFAULT);

  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl,H5FD_MPIO_COLLECTIVE);

  hsize_t vertex_dims[2] = {static_cast<hsize_t>(global_vertices),3};
  hsize_t vertex_count[2] = {static_cast<hsize_t>(local_vertices),3};
  hsize_t vertex_off[2] = {static_cast<hsize_t>(vertex_offset),0};
  hsize_t vertex_chunk[2] = {global_vertices > 1048576 ? 1048576 : (global_vertices > 0 ? static_cast<hsize_t>(global_vertices) : 1),3};
  write_double_dataset(geometry_group,"Vertices",2,vertex_dims,vertex_count,vertex_off,vertex_chunk,dxpl,
                       vertices.empty() ? 0 : &vertices[0]);

  hsize_t cell_dims[2] = {static_cast<hsize_t>(global_cells),3};
  hsize_t cell_count[2] = {static_cast<hsize_t>(local_cells),3};
  hsize_t cell_off[2] = {static_cast<hsize_t>(cell_offset),0};
  hsize_t cell_chunk[2] = {global_cells > 1048576 ? 1048576 : (global_cells > 0 ? static_cast<hsize_t>(global_cells) : 1),3};
  write_int_dataset(topology_group,"Connectivity",2,cell_dims,cell_count,cell_off,cell_chunk,dxpl,
                    connectivity.empty() ? 0 : &connectivity[0]);

  hsize_t scalar_dims[1] = {static_cast<hsize_t>(global_cells)};
  hsize_t scalar_count[1] = {static_cast<hsize_t>(local_cells)};
  hsize_t scalar_off[1] = {static_cast<hsize_t>(cell_offset)};
  hsize_t scalar_chunk[1] = {global_cells > 1048576 ? 1048576 : (global_cells > 0 ? static_cast<hsize_t>(global_cells) : 1)};

  if (builtin_mask_ & PROP_ID)
    write_int_dataset(celldata_group,"id",1,scalar_dims,scalar_count,scalar_off,scalar_chunk,dxpl,id_data.empty() ? 0 : &id_data[0]);
  if (builtin_mask_ & PROP_OWNER)
    write_int_dataset(celldata_group,"owner",1,scalar_dims,scalar_count,scalar_off,scalar_chunk,dxpl,owner_data.empty() ? 0 : &owner_data[0]);
  if (builtin_mask_ & PROP_AREA)
    write_double_dataset(celldata_group,"area",1,scalar_dims,scalar_count,scalar_off,scalar_chunk,dxpl,area_data.empty() ? 0 : &area_data[0]);
  if (builtin_mask_ & PROP_NORMAL)
    write_double_dataset(celldata_group,"normal",2,cell_dims,cell_count,cell_off,cell_chunk,dxpl,normal_data.empty() ? 0 : &normal_data[0]);
  if (builtin_mask_ & PROP_AEDGES)
    write_int_dataset(celldata_group,"active_edges",1,scalar_dims,scalar_count,scalar_off,scalar_chunk,dxpl,aedges_data.empty() ? 0 : &aedges_data[0]);
  if (builtin_mask_ & PROP_ACORNERS)
    write_int_dataset(celldata_group,"active_corners",1,scalar_dims,scalar_count,scalar_off,scalar_chunk,dxpl,acorners_data.empty() ? 0 : &acorners_data[0]);
  if (builtin_mask_ & PROP_INDEX)
    write_int_dataset(celldata_group,"index",1,scalar_dims,scalar_count,scalar_off,scalar_chunk,dxpl,index_data.empty() ? 0 : &index_data[0]);
  if (builtin_mask_ & PROP_NNEIGHS)
    write_int_dataset(celldata_group,"nneighs",1,scalar_dims,scalar_count,scalar_off,scalar_chunk,dxpl,nneighs_data.empty() ? 0 : &nneighs_data[0]);

  for (size_t ip = 0; ip < scalar_properties_.size(); ++ip)
    write_double_dataset(celldata_group,scalar_properties_[ip].c_str(),1,scalar_dims,scalar_count,scalar_off,scalar_chunk,dxpl,
                         scalar_data[ip].empty() ? 0 : &scalar_data[ip][0]);
  for (size_t ip = 0; ip < vector_properties_.size(); ++ip)
    write_double_dataset(celldata_group,vector_properties_[ip].c_str(),2,cell_dims,cell_count,cell_off,cell_chunk,dxpl,
                         vector_data[ip].empty() ? 0 : &vector_data[ip][0]);

  H5Pclose(dxpl);
  H5Gclose(celldata_group);
  H5Gclose(topology_group);
  H5Gclose(geometry_group);
  H5Gclose(step_group);
  H5Fclose(file);

  written_steps_.push_back(update->ntimestep);
  written_cells_.push_back(global_cells);
  if (comm->me == 0) write_xdmf(h5name);
#endif
}

/* ---------------------------------------------------------------------- */

void DumpMeshHDF5::write_xdmf(const char *h5name) const
{
  char xdmfname[1200];
  snprintf(xdmfname,sizeof(xdmfname),"%s.xdmf",h5name);

  const char *h5ref = strrchr(h5name,'/');
  h5ref = h5ref ? h5ref + 1 : h5name;

  FILE *xmf = fopen(xdmfname,"w");
  if (!xmf) error->one(FLERR,"Cannot open XDMF sidecar file for dump mesh/hdf5");

  fprintf(xmf,"<?xml version=\"1.0\" ?>\n");
  fprintf(xmf,"<Xdmf Version=\"3.0\">\n");
  fprintf(xmf,"  <Domain>\n");
  fprintf(xmf,"    <Grid Name=\"LIGGGHTS_Mesh\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");

  for (size_t i = 0; i < written_steps_.size(); ++i) {
    const long long step = static_cast<long long>(written_steps_[i]);
    const long long cells = written_cells_[i];
    const long long vertices = 3 * cells;
    fprintf(xmf,"      <Grid Name=\"Step_%lld\" GridType=\"Uniform\">\n",step);
    fprintf(xmf,"        <Time Value=\"%lld\" />\n",step);
    fprintf(xmf,"        <Topology TopologyType=\"Triangle\" NumberOfElements=\"%lld\">\n",cells);
    fprintf(xmf,"          <DataItem Format=\"HDF\" Dimensions=\"%lld 3\">%s:/Step_%lld/Topology/Connectivity</DataItem>\n",cells,h5ref,step);
    fprintf(xmf,"        </Topology>\n");
    fprintf(xmf,"        <Geometry GeometryType=\"XYZ\">\n");
    fprintf(xmf,"          <DataItem Format=\"HDF\" Dimensions=\"%lld 3\">%s:/Step_%lld/Geometry/Vertices</DataItem>\n",vertices,h5ref,step);
    fprintf(xmf,"        </Geometry>\n");

    if (builtin_mask_ & PROP_ID)
      fprintf(xmf,"        <Attribute Name=\"id\" AttributeType=\"Scalar\" Center=\"Cell\"><DataItem Format=\"HDF\" Dimensions=\"%lld\">%s:/Step_%lld/CellData/id</DataItem></Attribute>\n",cells,h5ref,step);
    if (builtin_mask_ & PROP_OWNER)
      fprintf(xmf,"        <Attribute Name=\"owner\" AttributeType=\"Scalar\" Center=\"Cell\"><DataItem Format=\"HDF\" Dimensions=\"%lld\">%s:/Step_%lld/CellData/owner</DataItem></Attribute>\n",cells,h5ref,step);
    if (builtin_mask_ & PROP_AREA)
      fprintf(xmf,"        <Attribute Name=\"area\" AttributeType=\"Scalar\" Center=\"Cell\"><DataItem Format=\"HDF\" Dimensions=\"%lld\">%s:/Step_%lld/CellData/area</DataItem></Attribute>\n",cells,h5ref,step);
    if (builtin_mask_ & PROP_NORMAL)
      fprintf(xmf,"        <Attribute Name=\"normal\" AttributeType=\"Vector\" Center=\"Cell\"><DataItem Format=\"HDF\" Dimensions=\"%lld 3\">%s:/Step_%lld/CellData/normal</DataItem></Attribute>\n",cells,h5ref,step);
    if (builtin_mask_ & PROP_AEDGES)
      fprintf(xmf,"        <Attribute Name=\"active_edges\" AttributeType=\"Scalar\" Center=\"Cell\"><DataItem Format=\"HDF\" Dimensions=\"%lld\">%s:/Step_%lld/CellData/active_edges</DataItem></Attribute>\n",cells,h5ref,step);
    if (builtin_mask_ & PROP_ACORNERS)
      fprintf(xmf,"        <Attribute Name=\"active_corners\" AttributeType=\"Scalar\" Center=\"Cell\"><DataItem Format=\"HDF\" Dimensions=\"%lld\">%s:/Step_%lld/CellData/active_corners</DataItem></Attribute>\n",cells,h5ref,step);
    if (builtin_mask_ & PROP_INDEX)
      fprintf(xmf,"        <Attribute Name=\"index\" AttributeType=\"Scalar\" Center=\"Cell\"><DataItem Format=\"HDF\" Dimensions=\"%lld\">%s:/Step_%lld/CellData/index</DataItem></Attribute>\n",cells,h5ref,step);
    if (builtin_mask_ & PROP_NNEIGHS)
      fprintf(xmf,"        <Attribute Name=\"nneighs\" AttributeType=\"Scalar\" Center=\"Cell\"><DataItem Format=\"HDF\" Dimensions=\"%lld\">%s:/Step_%lld/CellData/nneighs</DataItem></Attribute>\n",cells,h5ref,step);

    for (size_t ip = 0; ip < scalar_properties_.size(); ++ip)
      fprintf(xmf,"        <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\"><DataItem Format=\"HDF\" Dimensions=\"%lld\">%s:/Step_%lld/CellData/%s</DataItem></Attribute>\n",
              scalar_properties_[ip].c_str(),cells,h5ref,step,scalar_properties_[ip].c_str());
    for (size_t ip = 0; ip < vector_properties_.size(); ++ip)
      fprintf(xmf,"        <Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"Cell\"><DataItem Format=\"HDF\" Dimensions=\"%lld 3\">%s:/Step_%lld/CellData/%s</DataItem></Attribute>\n",
              vector_properties_[ip].c_str(),cells,h5ref,step,vector_properties_[ip].c_str());

    fprintf(xmf,"      </Grid>\n");
  }

  fprintf(xmf,"    </Grid>\n");
  fprintf(xmf,"  </Domain>\n");
  fprintf(xmf,"</Xdmf>\n");
  fclose(xmf);
}
