/* ----------------------------------------------------------------------
   Parallel HDF5 particle dump implementation.
------------------------------------------------------------------------- */

#include "dump_hdf5.h"
#include "atom.h"
#include "error.h"
#include "update.h"
#include <mpi.h>
#include <stdio.h>
#include <string.h>

#ifdef LIGGGHTS_HDF5
#include <hdf5.h>
#include <vector>
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpHDF5::DumpHDF5(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg)
{
  if (narg != 5)
    error->all(FLERR,"Illegal dump hdf5 command: expected 'dump ID group hdf5 N file.h5'");

  binary = 1;
  buffer_allow = 0;
  buffer_flag = 0;
  size_one = 14; // id + xyz + vxyz + fxyz + omegaxyz + radius
}

/* ---------------------------------------------------------------------- */

DumpHDF5::~DumpHDF5()
{
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
#endif
}

/* ---------------------------------------------------------------------- */

void DumpHDF5::resolve_filename(char *out, int nout) const
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

void DumpHDF5::pack_local(long long *ids, double *positions, double *velocities,
                          double *forces, double *omegas, double *radii,
                          int nlocal_selected) const
{
  tagint *tag = atom->tag;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    ids[m] = static_cast<long long>(tag[i]);
    positions[3*m+0] = x[i][0];
    positions[3*m+1] = x[i][1];
    positions[3*m+2] = x[i][2];
    velocities[3*m+0] = v[i][0];
    velocities[3*m+1] = v[i][1];
    velocities[3*m+2] = v[i][2];
    forces[3*m+0] = f[i][0];
    forces[3*m+1] = f[i][1];
    forces[3*m+2] = f[i][2];
    omegas[3*m+0] = omega[i][0];
    omegas[3*m+1] = omega[i][1];
    omegas[3*m+2] = omega[i][2];
    radii[m] = radius[i];
    m++;
  }

  if (m != nlocal_selected)
    error->one(FLERR,"Internal dump hdf5 packing error");
}

/* ---------------------------------------------------------------------- */

void DumpHDF5::write()
{
#ifndef LIGGGHTS_HDF5
  error->all(FLERR,"Dump hdf5 requires a parallel HDF5 build with -DLIGGGHTS_HDF5");
#else
  const int nlocal_selected = count();
  long long local_count_ll = nlocal_selected;
  long long global_count_ll = 0;
  long long offset_ll = 0;

  MPI_Allreduce(&local_count_ll,&global_count_ll,1,MPI_LONG_LONG,MPI_SUM,world);
  MPI_Exscan(&local_count_ll,&offset_ll,1,MPI_LONG_LONG,MPI_SUM,world);
  if (me == 0) offset_ll = 0;

  std::vector<long long> ids(nlocal_selected);
  std::vector<double> positions(3*nlocal_selected);
  std::vector<double> velocities(3*nlocal_selected);
  std::vector<double> forces(3*nlocal_selected);
  std::vector<double> omegas(3*nlocal_selected);
  std::vector<double> radii(nlocal_selected);
  pack_local(ids.empty() ? 0 : &ids[0],
             positions.empty() ? 0 : &positions[0],
             velocities.empty() ? 0 : &velocities[0],
             forces.empty() ? 0 : &forces[0],
             omegas.empty() ? 0 : &omegas[0],
             radii.empty() ? 0 : &radii[0],
             nlocal_selected);

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
  if (file < 0) error->all(FLERR,"Cannot open parallel HDF5 dump file");

  char step_group_name[128];
  snprintf(step_group_name,sizeof(step_group_name),"/Step_" BIGINT_FORMAT,update->ntimestep);
  hid_t group = H5Gcreate(file,step_group_name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if (group < 0) error->all(FLERR,"Cannot create HDF5 timestep group");

  hsize_t scalar_dims[1] = {static_cast<hsize_t>(global_count_ll)};
  hsize_t vector_dims[2] = {static_cast<hsize_t>(global_count_ll),3};
  hsize_t scalar_count[1] = {static_cast<hsize_t>(nlocal_selected)};
  hsize_t vector_count[2] = {static_cast<hsize_t>(nlocal_selected),3};
  hsize_t scalar_offset[1] = {static_cast<hsize_t>(offset_ll)};
  hsize_t vector_offset[2] = {static_cast<hsize_t>(offset_ll),0};

  const hsize_t chunk_rows =
    static_cast<hsize_t>(global_count_ll < 1048576 ? global_count_ll : 1048576);
  hsize_t scalar_chunk[1] = {chunk_rows > 0 ? chunk_rows : 1};
  hsize_t vector_chunk[2] = {chunk_rows > 0 ? chunk_rows : 1,3};

  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl,H5FD_MPIO_COLLECTIVE);

  hid_t filespace = H5Screate_simple(1,scalar_dims,NULL);
  hid_t memspace = H5Screate_simple(1,scalar_count,NULL);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,scalar_offset,NULL,scalar_count,NULL);
  hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl,1,scalar_chunk);
  hid_t dset = H5Dcreate(group,"id",H5T_NATIVE_LLONG,filespace,H5P_DEFAULT,dcpl,H5P_DEFAULT);
  H5Dwrite(dset,H5T_NATIVE_LLONG,memspace,filespace,dxpl,ids.empty() ? 0 : &ids[0]);
  H5Dclose(dset);
  H5Pclose(dcpl);
  H5Sclose(memspace);
  H5Sclose(filespace);

  filespace = H5Screate_simple(2,vector_dims,NULL);
  memspace = H5Screate_simple(2,vector_count,NULL);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,vector_offset,NULL,vector_count,NULL);
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl,2,vector_chunk);
  dset = H5Dcreate(group,"position",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,dcpl,H5P_DEFAULT);
  H5Dwrite(dset,H5T_NATIVE_DOUBLE,memspace,filespace,dxpl,positions.empty() ? 0 : &positions[0]);
  H5Dclose(dset);
  H5Pclose(dcpl);
  H5Sclose(memspace);
  H5Sclose(filespace);

  filespace = H5Screate_simple(2,vector_dims,NULL);
  memspace = H5Screate_simple(2,vector_count,NULL);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,vector_offset,NULL,vector_count,NULL);
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl,2,vector_chunk);
  dset = H5Dcreate(group,"velocity",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,dcpl,H5P_DEFAULT);
  H5Dwrite(dset,H5T_NATIVE_DOUBLE,memspace,filespace,dxpl,velocities.empty() ? 0 : &velocities[0]);
  H5Dclose(dset);
  H5Pclose(dcpl);
  H5Sclose(memspace);
  H5Sclose(filespace);

  filespace = H5Screate_simple(2,vector_dims,NULL);
  memspace = H5Screate_simple(2,vector_count,NULL);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,vector_offset,NULL,vector_count,NULL);
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl,2,vector_chunk);
  dset = H5Dcreate(group,"force",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,dcpl,H5P_DEFAULT);
  H5Dwrite(dset,H5T_NATIVE_DOUBLE,memspace,filespace,dxpl,forces.empty() ? 0 : &forces[0]);
  H5Dclose(dset);
  H5Pclose(dcpl);
  H5Sclose(memspace);
  H5Sclose(filespace);

  filespace = H5Screate_simple(2,vector_dims,NULL);
  memspace = H5Screate_simple(2,vector_count,NULL);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,vector_offset,NULL,vector_count,NULL);
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl,2,vector_chunk);
  dset = H5Dcreate(group,"omega",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,dcpl,H5P_DEFAULT);
  H5Dwrite(dset,H5T_NATIVE_DOUBLE,memspace,filespace,dxpl,omegas.empty() ? 0 : &omegas[0]);
  H5Dclose(dset);
  H5Pclose(dcpl);
  H5Sclose(memspace);
  H5Sclose(filespace);

  filespace = H5Screate_simple(1,scalar_dims,NULL);
  memspace = H5Screate_simple(1,scalar_count,NULL);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,scalar_offset,NULL,scalar_count,NULL);
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl,1,scalar_chunk);
  dset = H5Dcreate(group,"radius",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,dcpl,H5P_DEFAULT);
  H5Dwrite(dset,H5T_NATIVE_DOUBLE,memspace,filespace,dxpl,radii.empty() ? 0 : &radii[0]);
  H5Dclose(dset);
  H5Pclose(dcpl);
  H5Sclose(memspace);
  H5Sclose(filespace);

  hsize_t attr_dims[1] = {1};
  hid_t attr_space = H5Screate_simple(1,attr_dims,NULL);
  hid_t attr = H5Acreate(group,"timestep",H5T_NATIVE_LLONG,attr_space,H5P_DEFAULT,H5P_DEFAULT);
  long long step_ll = static_cast<long long>(update->ntimestep);
  H5Awrite(attr,H5T_NATIVE_LLONG,&step_ll);
  H5Aclose(attr);
  H5Sclose(attr_space);

  H5Pclose(dxpl);
  H5Gclose(group);
  H5Fclose(file);

  written_steps_.push_back(update->ntimestep);
  written_counts_.push_back(global_count_ll);
  if (me == 0) write_xdmf(h5name);
#endif
}

/* ---------------------------------------------------------------------- */

void DumpHDF5::write_xdmf(const char *h5name) const
{
  char xdmfname[1200];
  snprintf(xdmfname,sizeof(xdmfname),"%s.xdmf",h5name);

  const char *h5ref = strrchr(h5name,'/');
  h5ref = h5ref ? h5ref + 1 : h5name;

  FILE *xmf = fopen(xdmfname,"w");
  if (!xmf) error->one(FLERR,"Cannot open XDMF sidecar file for dump hdf5");

  fprintf(xmf,"<?xml version=\"1.0\" ?>\n");
  fprintf(xmf,"<Xdmf Version=\"3.0\">\n");
  fprintf(xmf,"  <Domain>\n");
  fprintf(xmf,"    <Grid Name=\"LIGGGHTS\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");

  for (size_t i = 0; i < written_steps_.size(); ++i) {
    const long long step = static_cast<long long>(written_steps_[i]);
    const long long count = written_counts_[i];
    fprintf(xmf,"      <Grid Name=\"Step_%lld\" GridType=\"Uniform\">\n",step);
    fprintf(xmf,"        <Time Value=\"%lld\" />\n",step);
    fprintf(xmf,"        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%lld\" />\n",count);
    fprintf(xmf,"        <Geometry GeometryType=\"XYZ\">\n");
    fprintf(xmf,"          <DataItem Format=\"HDF\" Dimensions=\"%lld 3\">%s:/Step_%lld/position</DataItem>\n",count,h5ref,step);
    fprintf(xmf,"        </Geometry>\n");
    fprintf(xmf,"        <Attribute Name=\"id\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf,"          <DataItem Format=\"HDF\" Dimensions=\"%lld\">%s:/Step_%lld/id</DataItem>\n",count,h5ref,step);
    fprintf(xmf,"        </Attribute>\n");
    fprintf(xmf,"        <Attribute Name=\"velocity\" AttributeType=\"Vector\" Center=\"Node\">\n");
    fprintf(xmf,"          <DataItem Format=\"HDF\" Dimensions=\"%lld 3\">%s:/Step_%lld/velocity</DataItem>\n",count,h5ref,step);
    fprintf(xmf,"        </Attribute>\n");
    fprintf(xmf,"        <Attribute Name=\"force\" AttributeType=\"Vector\" Center=\"Node\">\n");
    fprintf(xmf,"          <DataItem Format=\"HDF\" Dimensions=\"%lld 3\">%s:/Step_%lld/force</DataItem>\n",count,h5ref,step);
    fprintf(xmf,"        </Attribute>\n");
    fprintf(xmf,"        <Attribute Name=\"omega\" AttributeType=\"Vector\" Center=\"Node\">\n");
    fprintf(xmf,"          <DataItem Format=\"HDF\" Dimensions=\"%lld 3\">%s:/Step_%lld/omega</DataItem>\n",count,h5ref,step);
    fprintf(xmf,"        </Attribute>\n");
    fprintf(xmf,"        <Attribute Name=\"radius\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf,"          <DataItem Format=\"HDF\" Dimensions=\"%lld\">%s:/Step_%lld/radius</DataItem>\n",count,h5ref,step);
    fprintf(xmf,"        </Attribute>\n");
    fprintf(xmf,"      </Grid>\n");
  }

  fprintf(xmf,"    </Grid>\n");
  fprintf(xmf,"  </Domain>\n");
  fprintf(xmf,"</Xdmf>\n");
  fclose(xmf);
}
