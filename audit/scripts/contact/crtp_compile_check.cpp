// Does src/contact_model_crtp_api.h compile and can it be called with the
// double** matrices the PropertyRegistry hands out?  (It is included nowhere.)
#include "../../../src/contact_model_crtp_api.h"
using namespace LIGGGHTS::ContactModels::CRTP;
int main() {
  double row[2] = {0, 1e7}; double *Y[2] = {row, row};   // registry-style double**
  double brow[2] = {0, -0.2}; double *B[2] = {brow, brow};
  double x0[3]={0,0,0}, x1[3]={0.0019,0,0}; double *x[2]={x0,x1};
  double v0[3]={0,0,0}, v1[3]={0,0,0}; double *v[2]={v0,v1};
  double f0[3]={0,0,0}, f1[3]={0,0,0}; double *f[2]={f0,f1};
  double rad[2]={1e-3,1e-3}, m[2]={1e-5,1e-5}; int type[2]={1,1}, mask[2]={1,1};
  AtomSoAView a{x,v,f,rad,m,type,mask,2};
  int il[1]={0}, nn[1]={1}, nb[1]={1}; int *fn[1]={nb};
  HalfNeighborListView l{1,il,nn,fn};
#ifdef TRY_DOUBLE_PTR
  compute_with_plugin<CONTACT_PLUGIN_HERTZ_MINDLIN>(a,l,1,1.0,Y,Y,B,0.9128709291752769);
#else
  compute_with_plugin<CONTACT_PLUGIN_HERTZ_MINDLIN>(a,l,1,1.0,(const double**)Y,(const double**)Y,(const double**)B,0.9128709291752769);
#endif
  __builtin_printf("f_i=%g f_j=%g (j<nlocal so j updated)\n", f[0][0], f[1][0]);
  return 0;
}
