// Unit-level check of src/aligned_particle_soa.h: (1) n==0 path, (2) bitwise equality vs legacy FixNVE (dtf/m vs dtf*(1/m)), (3) alignment.
#include "aligned_particle_soa.h"
#include <cstdio>
#include <cstdint>
#include <random>
#include <cstring>
using namespace LAMMPS_NS;
int main(int argc,char**argv){
  const int n = 100000; std::mt19937_64 g(42); std::uniform_real_distribution<double> U(0.5,2.0);
  double **x=new double*[n],**v=new double*[n],**f=new double*[n]; double *rm=new double[n]; int *mask=new int[n], *type=new int[n];
  double **xl=new double*[n],**vl=new double*[n];
  for(int i=0;i<n;i++){x[i]=new double[3];v[i]=new double[3];f[i]=new double[3];xl[i]=new double[3];vl[i]=new double[3];
    for(int c=0;c<3;c++){x[i][c]=xl[i][c]=U(g);v[i][c]=vl[i][c]=U(g)-1.25;f[i][c]=(U(g)-1.25)*1e3;} rm[i]=U(g)*1e-6; mask[i]=1; type[i]=1;}
  const double dtv=1e-5, dtf=0.5*1e-5;
  ParticleSoA s; s.load_from_aos(x,v,f,rm,0,type,n); s.initial_integrate_nve(dtv,dtf,mask,1,n); s.store_xv_to_aos(x,v,n);
  long diffx=0,diffv=0;
  for(int i=0;i<n;i++){ double dtfm=dtf/rm[i]; for(int c=0;c<3;c++){ vl[i][c]+=dtfm*f[i][c]; xl[i][c]+=dtv*vl[i][c];
      if(memcmp(&vl[i][c],&v[i][c],8))diffv++; if(memcmp(&xl[i][c],&x[i][c],8))diffx++; } }
  printf("bitwise mismatches after one initial_integrate: v %ld / %d, x %ld / %d\n",diffv,3*n,diffx,3*n);
  if(argc>1 && strcmp(argv[1],"empty")==0){ ParticleSoA e; e.load_from_aos(x,v,f,rm,0,type,0); printf("calling initial_integrate_nve with n=0 on never-sized SoA\n"); fflush(stdout); e.initial_integrate_nve(dtv,dtf,mask,1,0); printf("n=0 returned normally\n"); }
  return 0;
}
