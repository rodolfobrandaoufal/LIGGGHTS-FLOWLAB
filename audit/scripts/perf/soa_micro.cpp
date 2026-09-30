// SoA vs AoS fix-nve micro-benchmark: times each phase of the LIGGGHTS_USE_SOA_NVE path
// (src/aligned_particle_soa.h, used verbatim) against the in-place AoS loop of src/fix_nve.cpp.
// Build: g++ -O3 -march=native -fno-fast-math -I../../../src soa_micro.cpp -o soa_micro
#include "aligned_particle_soa.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
using namespace LAMMPS_NS;
static double **alloc2(int n) { double *d = new double[3*n]; double **p = new double*[n]; for (int i=0;i<n;i++) p[i]=d+3*i; return p; }
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
int main(int argc, char **argv) {
  int n = atoi(argv[1]); int iters = atoi(argv[2]);
  double **x = alloc2(n), **v = alloc2(n), **f = alloc2(n);
  std::vector<double> rmass(n); std::vector<int> mask(n, 1), type(n, 1); double mass[2] = {1, 1};
  for (int i=0;i<n;i++){ rmass[i]=1e-5*(1+0.1*(i%7)); for(int k=0;k<3;k++){x[i][k]=i*1e-3+k; v[i][k]=1e-3*k; f[i][k]=1e-4*(k+1);} }
  const double dtv = 1e-5, dtf = 0.5e-5; const int groupbit = 1;
  ParticleSoA soa;
  double t_load=0,t_kern=0,t_store=0,t_aos=0, t_load2=0,t_kern2=0,t_store2=0,t_aos2=0;
  std::vector<double> sl, sa;
  for (int it=0; it<iters; it++) {
    double t0=now(); soa.load_from_aos(x,v,f,rmass.data(),mass,type.data(),n);
    double t1=now(); soa.initial_integrate_nve(dtv,dtf,mask.data(),groupbit,n);
    double t2=now(); soa.store_xv_to_aos(x,v,n);
    double t3=now(); soa.load_from_aos(x,v,f,rmass.data(),mass,type.data(),n);
    double t4=now(); soa.final_integrate_nve(dtf,mask.data(),groupbit,n);
    double t5=now(); soa.store_v_to_aos(v,n);
    double t6=now();
    // AoS reference (fix_nve.cpp rmass branch, initial + final)
    for (int i=0;i<n;i++) if (mask[i]&groupbit){ double dtfm=dtf/rmass[i]; v[i][0]+=dtfm*f[i][0]; v[i][1]+=dtfm*f[i][1]; v[i][2]+=dtfm*f[i][2]; x[i][0]+=dtv*v[i][0]; x[i][1]+=dtv*v[i][1]; x[i][2]+=dtv*v[i][2]; }
    double t7=now();
    for (int i=0;i<n;i++) if (mask[i]&groupbit){ double dtfm=dtf/rmass[i]; v[i][0]+=dtfm*f[i][0]; v[i][1]+=dtfm*f[i][1]; v[i][2]+=dtfm*f[i][2]; }
    double t8=now();
    if (it==0) continue; // warm-up
    t_load+=t1-t0; t_kern+=t2-t1; t_store+=t3-t2; t_load2+=t4-t3; t_kern2+=t5-t4; t_store2+=t6-t5; t_aos+=t7-t6; t_aos2+=t8-t7;
    sl.push_back((t6-t0)); sa.push_back(t8-t6);
  }
  int m = iters-1; double ns = 1e9/(double)m/n;
  std::sort(sl.begin(), sl.end()); std::sort(sa.begin(), sa.end());
  printf("n=%d iters=%d [ns/atom/step]\n", n, m);
  printf(" SoA initial: load %.3f kernel %.3f store %.3f | final: load %.3f kernel %.3f store %.3f | SoA total %.3f (median %.3f)\n",
    t_load*ns,t_kern*ns,t_store*ns,t_load2*ns,t_kern2*ns,t_store2*ns,(t_load+t_kern+t_store+t_load2+t_kern2+t_store2)*ns, sl[m/2]*1e9/n);
  printf(" AoS initial %.3f final %.3f | AoS total %.3f (median %.3f)\n", t_aos*ns, t_aos2*ns, (t_aos+t_aos2)*ns, sa[m/2]*1e9/n);
  printf(" ratio SoA/AoS (median) %.2f ; copy share of SoA %.1f%%\n", sl[m/2]/sa[m/2], 100*(t_load+t_store+t_load2+t_store2)/(t_load+t_kern+t_store+t_load2+t_kern2+t_store2));
  return (int)(x[n/2][0]*0);
}
