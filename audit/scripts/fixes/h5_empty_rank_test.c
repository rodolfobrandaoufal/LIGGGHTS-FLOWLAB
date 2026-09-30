/* Replicates the dataset-write pattern of src/dump_hdf5.cpp (hyperslab with count 0 on empty ranks,
   chunk dims = max(global,1)) and of src/dump_mesh_hdf5.cpp (H5Sselect_none on empty ranks).
   argv[1] = number of particles on rank 0 (others get 0 unless argv[2] given). */
#include <hdf5.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
static int write_like_dump_hdf5(hid_t g,const char*name,long long nloc,long long nglob,long long off,int use_none){
  hsize_t dims[1]={(hsize_t)nglob}, cnt[1]={(hsize_t)nloc}, of[1]={(hsize_t)off};
  hsize_t chunk[1]={nglob<1048576?(hsize_t)nglob:1048576}; if(chunk[0]==0)chunk[0]=1;
  double *buf = nloc? malloc(sizeof(double)*nloc):NULL; for(long long i=0;i<nloc;i++)buf[i]=off+i;
  hid_t fs=H5Screate_simple(1,dims,NULL); hid_t ms;
  herr_t e1;
  if(use_none && nloc==0){ hsize_t one[1]={1}; ms=H5Screate_simple(1,one,NULL); H5Sselect_none(ms); e1=H5Sselect_none(fs);}  
  else { ms=H5Screate_simple(1,cnt,NULL); e1=H5Sselect_hyperslab(fs,H5S_SELECT_SET,of,NULL,cnt,NULL);} 
  hid_t dcpl=H5Pcreate(H5P_DATASET_CREATE); herr_t e2=H5Pset_chunk(dcpl,1,chunk);
  hid_t ds=H5Dcreate(g,name,H5T_NATIVE_DOUBLE,fs,H5P_DEFAULT,dcpl,H5P_DEFAULT);
  hid_t dx=H5Pcreate(H5P_DATASET_XFER); H5Pset_dxpl_mpio(dx,H5FD_MPIO_COLLECTIVE);
  double dummy=0; herr_t e3=H5Dwrite(ds,H5T_NATIVE_DOUBLE,ms,fs,dx,buf?buf:&dummy);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  printf("[rank %d] %s nloc=%lld nglob=%lld: select=%d set_chunk=%d dcreate=%lld dwrite=%d\n",rank,name,nloc,nglob,(int)e1,(int)e2,(long long)ds,(int)e3);
  if(ds>=0)H5Dclose(ds); H5Pclose(dcpl);H5Pclose(dx);H5Sclose(ms);H5Sclose(fs); free(buf); return (ds<0||e3<0);
}
int main(int argc,char**argv){
  MPI_Init(&argc,&argv); int rank,size; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
  H5Eset_auto(H5E_DEFAULT,NULL,NULL);
  long long n0 = argc>1? atoll(argv[1]):10;
  long long nloc = rank==0? n0:0, nglob=0, off=0;
  MPI_Allreduce(&nloc,&nglob,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  MPI_Exscan(&nloc,&off,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD); if(rank==0)off=0;
  hid_t fa=H5Pcreate(H5P_FILE_ACCESS); H5Pset_fapl_mpio(fa,MPI_COMM_WORLD,MPI_INFO_NULL);
  const char* fn = argc>2? argv[2] : "h5_empty_rank_test.h5";
  hid_t f=H5Fcreate(fn,H5F_ACC_TRUNC,H5P_DEFAULT,fa); H5Pclose(fa);
  int bad=0;
  bad|=write_like_dump_hdf5(f,"particle_style",nloc,nglob,off,0);
  bad|=write_like_dump_hdf5(f,"mesh_style",nloc,nglob,off,1);
  H5Fclose(f);
  unsigned maj,min,rel; H5get_libversion(&maj,&min,&rel);
  if(rank==0)printf("HDF5 %u.%u.%u, ranks=%d, any_failure_rank0=%d\n",maj,min,rel,size,bad);
  MPI_Finalize(); return 0;
}
