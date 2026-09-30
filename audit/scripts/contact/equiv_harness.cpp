// Equivalence harness for the "non-semantic" contact-model micro-optimizations
// (LIGGGHTS_MODIFICATION_REPORT.txt section 8).
//
// Each OLD_* function is a verbatim transcription of `git show HEAD:<file>`,
// each NEW_* function a verbatim transcription of the working tree.
// Build and run: see run_equiv.sh in this directory.
//
// For every pair we count: (a) branch-selection mismatches, (b) outputs that
// are not bitwise identical (NaN == NaN treated as equal, since NaN payloads
// are not semantically meaningful), (c) max ULP distance when both finite.

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <random>
#include <vector>
#include <limits>

static const double NaN = std::numeric_limits<double>::quiet_NaN();
static const double INF = std::numeric_limits<double>::infinity();
static const double DMIN = std::numeric_limits<double>::denorm_min();

static int64_t ulp_dist(double a, double b) {
  if (std::isnan(a) && std::isnan(b)) return 0;
  if (std::isnan(a) || std::isnan(b)) return INT64_MAX;
  if (a == b) return 0;                 // also +0 == -0
  int64_t ia, ib; std::memcpy(&ia, &a, 8); std::memcpy(&ib, &b, 8);
  if (ia < 0) ia = INT64_MIN - ia;      // map to monotone integer line
  if (ib < 0) ib = INT64_MIN - ib;
  int64_t d = ia - ib; return d < 0 ? -d : d;
}
static bool same(double a, double b) {
  if (std::isnan(a) && std::isnan(b)) return true;
  return std::memcmp(&a, &b, 8) == 0;
}

/* ------------------------------------------------------------------ */
/* 1. tangential_model_no_history.h : gamma selection                   */
/* ------------------------------------------------------------------ */
struct NH { double gamma; int branch; };  // branch 1 = rescaled (Ft_friction/vrel)

__attribute__((noinline)) NH OLD_no_history(double xmu, double Fn, double gammat,
                                             double vtr1, double vtr2, double vtr3) {
  const double vrel = sqrt(vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3);
  const double Ft_friction = xmu * fabs(Fn);
  double gamma = 0.0; int br;
  if (Ft_friction < gammat*vrel) { gamma = Ft_friction/vrel; br = 1; }
  else { gamma = gammat; br = 0; }
  return {gamma, br};
}
__attribute__((noinline)) NH NEW_no_history(double xmu, double Fn, double gammat,
                                             double vtr1, double vtr2, double vtr3) {
  const double vrelsq = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  const double Ft_friction = xmu * fabs(Fn);
  double gamma = gammat; int br = 0;
  if (gammat > 0.0 && vrelsq > 0.0 &&
      Ft_friction*Ft_friction < gammat*gammat*vrelsq) { gamma = Ft_friction/sqrt(vrelsq); br = 1; }
  return {gamma, br};
}

/* ------------------------------------------------------------------ */
/* 2. tangential_model_history.h : Coulomb sliding check                */
/* ------------------------------------------------------------------ */
struct TH { double Ft[3]; double shear[3]; int branch; }; // 0 stick, 1 slide, 2 slide&shrmag==0

__attribute__((noinline)) TH OLD_history(double kt, double xmu, double Fn, double gammat,
                                          const double sh[3], const double vtr[3]) {
  TH r; double shear[3] = {sh[0], sh[1], sh[2]};
  const double shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);
  double Ft1 = -(kt*shear[0]), Ft2 = -(kt*shear[1]), Ft3 = -(kt*shear[2]);
  const double Ft_shear = kt * shrmag;
  const double Ft_friction = xmu * fabs(Fn);
  if (Ft_shear > Ft_friction) {
    if (shrmag != 0.0) {
      const double ratio = Ft_friction / Ft_shear;
      Ft1 *= ratio; Ft2 *= ratio; Ft3 *= ratio;
      shear[0] = -Ft1/kt; shear[1] = -Ft2/kt; shear[2] = -Ft3/kt; r.branch = 1;
    } else { Ft1 = Ft2 = Ft3 = 0.0; r.branch = 2; }
  } else {
    Ft1 -= gammat*vtr[0]; Ft2 -= gammat*vtr[1]; Ft3 -= gammat*vtr[2]; r.branch = 0;
  }
  r.Ft[0]=Ft1; r.Ft[1]=Ft2; r.Ft[2]=Ft3; r.shear[0]=shear[0]; r.shear[1]=shear[1]; r.shear[2]=shear[2];
  return r;
}
__attribute__((noinline)) TH NEW_history(double kt, double xmu, double Fn, double gammat,
                                          const double sh[3], const double vtr[3]) {
  TH r; double shear[3] = {sh[0], sh[1], sh[2]};
  const double shrsq = shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2];
  double Ft1 = -(kt*shear[0]), Ft2 = -(kt*shear[1]), Ft3 = -(kt*shear[2]);
  const double Ft_friction = xmu * fabs(Fn);
  const double Ft_shear_sq = kt * kt * shrsq;
  const double Ft_friction_sq = Ft_friction * Ft_friction;
  if (Ft_shear_sq > Ft_friction_sq) {
    if (shrsq != 0.0) {
      const double shrmag = sqrt(shrsq);
      const double Ft_shear = kt * shrmag;
      const double ratio = Ft_friction / Ft_shear;
      Ft1 *= ratio; Ft2 *= ratio; Ft3 *= ratio;
      shear[0] = -Ft1/kt; shear[1] = -Ft2/kt; shear[2] = -Ft3/kt; r.branch = 1;
    } else { Ft1 = Ft2 = Ft3 = 0.0; r.branch = 2; }
  } else {
    Ft1 -= gammat*vtr[0]; Ft2 -= gammat*vtr[1]; Ft3 -= gammat*vtr[2]; r.branch = 0;
  }
  r.Ft[0]=Ft1; r.Ft[1]=Ft2; r.Ft[2]=Ft3; r.shear[0]=shear[0]; r.shear[1]=shear[1]; r.shear[2]=shear[2];
  return r;
}

/* ------------------------------------------------------------------ */
/* 3. normal_model_luding.h : shared gamma                              */
/* ------------------------------------------------------------------ */
__attribute__((noinline)) void OLD_luding(double meff, double kn, double crl, bool td, double &gn, double &gt) {
  gn = sqrt(4.*meff*kn/(1.+(M_PI/crl)*(M_PI/crl)));
  gt = sqrt(4.*meff*kn/(1.+(M_PI/crl)*(M_PI/crl)));
  if (!td) gt = 0.0;
}
__attribute__((noinline)) void NEW_luding(double meff, double kn, double crl, bool td, double &gn, double &gt) {
  const double coeffRestLog = crl;
  const double coeffRestLogTerm = M_PI/coeffRestLog;
  const double gamma = sqrt(4.*meff*kn/(1.+coeffRestLogTerm*coeffRestLogTerm));
  gn = gamma; gt = td ? gamma : 0.0;
}

/* ------------------------------------------------------------------ */
/* 4. rolling_model_cdt.h (wall branch) and epsd/luding superquadric     */
/* ------------------------------------------------------------------ */
__attribute__((noinline)) int OLD_cdt(double w1, double w2, double w3, double out[3]) {
  const double wrmag = sqrt(w1*w1+w2*w2+w3*w3);
  if (wrmag > 0.) { out[0]=w1/wrmag; out[1]=w2/wrmag; out[2]=w3/wrmag; return 1; }
  out[0]=out[1]=out[2]=0; return 0;
}
__attribute__((noinline)) int NEW_cdt(double w1, double w2, double w3, double out[3]) {
  const double wrsq = w1*w1+w2*w2+w3*w3;
  if (wrsq > 0.) { const double wrmag = sqrt(wrsq); out[0]=w1/wrmag; out[1]=w2/wrmag; out[2]=w3/wrmag; return 1; }
  out[0]=out[1]=out[2]=0; return 0;
}
__attribute__((noinline)) int OLD_epsd(double w1, double w2, double w3, double out[3]) {
  const double omega_mag = sqrt(w1*w1 + w2*w2 + w3*w3);
  if (omega_mag != 0.0) { out[0]=w1/omega_mag; out[1]=w2/omega_mag; out[2]=w3/omega_mag; return 1; }
  out[0]=out[1]=out[2]=0; return 0;
}
__attribute__((noinline)) int NEW_epsd(double w1, double w2, double w3, double out[3]) {
  const double omega_sq = w1*w1 + w2*w2 + w3*w3;
  if (omega_sq != 0.0) { const double omega_mag = sqrt(omega_sq); out[0]=w1/omega_mag; out[1]=w2/omega_mag; out[2]=w3/omega_mag; return 1; }
  out[0]=out[1]=out[2]=0; return 0;
}

/* ------------------------------------------------------------------ */
struct Stat { long n=0, branch_mismatch=0, not_bitwise=0; int64_t max_ulp=0; };
static void acc_val(Stat &s, double a, double b) {
  if (!same(a,b)) { s.not_bitwise++; int64_t d = ulp_dist(a,b); if (d > s.max_ulp) s.max_ulp = d; }
}
static void report(const char *name, const Stat &s) {
  std::printf("%-44s n=%9ld branch_mismatch=%7ld non_bitwise_outputs=%7ld max_ulp=%lld\n",
              name, s.n, s.branch_mismatch, s.not_bitwise, (long long)s.max_ulp);
}

static std::mt19937_64 rng(20260929);
static double U(double a, double b) { return std::uniform_real_distribution<double>(a,b)(rng); }
static double logu(double lo, double hi) { return std::pow(10.0, U(lo,hi)); }
static double sgn() { return (rng() & 1) ? 1.0 : -1.0; }

int main() {
  const double special[] = {0.0, -0.0, DMIN, 1e-310, 1e-160, 1e-155, 1e-100, 1e-20, 1.0, 1e20, 1e150, 1e154, 1e155,
                            1e300, INF, -INF, NaN, -1.0};
  const int nsp = sizeof(special)/sizeof(special[0]);

  /* ---------- no_history: exhaustive special-value grid --------------- */
  {
    Stat s; long printed = 0;
    for (int a=0;a<nsp;a++) for (int b=0;b<nsp;b++) for (int c=0;c<nsp;c++) for (int d=0;d<nsp;d++) {
      const double Fn = special[a], gammat = special[b], v1 = special[c], v2 = special[d]*0.5;
      NH o = OLD_no_history(0.5, Fn, gammat, v1, v2, 0.0);
      NH n = NEW_no_history(0.5, Fn, gammat, v1, v2, 0.0);
      s.n++; if (o.branch != n.branch) s.branch_mismatch++; acc_val(s, o.gamma, n.gamma);
      if ((o.branch != n.branch || !same(o.gamma,n.gamma)) && printed < 40) {
        printed++;
        std::printf("  no_hist diff: Fn=%g gammat=%g v=(%g,%g,0): OLD br=%d g=%.17g  NEW br=%d g=%.17g\n",
                    Fn, gammat, v1, v2, o.branch, o.gamma, n.branch, n.gamma);
      }
    }
    report("no_history special-value grid (18^4)", s);
  }
  /* mandatory checkpoints */
  {
    std::printf("  checkpoint vrel==0  gammat=5 Fn=1  : OLD g=%.17g  NEW g=%.17g\n",
                OLD_no_history(0.5,1,5,0,0,0).gamma, NEW_no_history(0.5,1,5,0,0,0).gamma);
    std::printf("  checkpoint vrel==0  gammat=5 Fn=0  : OLD g=%.17g  NEW g=%.17g\n",
                OLD_no_history(0.5,0,5,0,0,0).gamma, NEW_no_history(0.5,0,5,0,0,0).gamma);
    std::printf("  checkpoint gammat==0 v=1 Fn=1      : OLD g=%.17g  NEW g=%.17g\n",
                OLD_no_history(0.5,1,0,1,0,0).gamma, NEW_no_history(0.5,1,0,1,0,0).gamma);
    std::printf("  checkpoint gammat==0 v=0 Fn=0      : OLD g=%.17g  NEW g=%.17g\n",
                OLD_no_history(0.5,0,0,0,0,0).gamma, NEW_no_history(0.5,0,0,0,0,0).gamma);
    std::printf("  checkpoint Fn==0 (Ft_friction=0) v=1 gammat=5 : OLD g=%.17g NEW g=%.17g\n",
                OLD_no_history(0.5,0,5,1,0,0).gamma, NEW_no_history(0.5,0,5,1,0,0).gamma);
    std::printf("  checkpoint gammat<0 (-5) v=1 Fn=1  : OLD g=%.17g  NEW g=%.17g\n",
                OLD_no_history(0.5,1,-5,1,0,0).gamma, NEW_no_history(0.5,1,-5,1,0,0).gamma);
    std::printf("  checkpoint gammat=inf v=0 Fn=1     : OLD g=%.17g  NEW g=%.17g\n",
                OLD_no_history(0.5,1,INF,0,0,0).gamma, NEW_no_history(0.5,1,INF,0,0,0).gamma);
    std::printf("  checkpoint vtr=1e-160 gammat=1 Fn=1e-170 : OLD br=%d g=%.17g  NEW br=%d g=%.17g\n",
                OLD_no_history(1,1e-170,1,1e-160,0,0).branch, OLD_no_history(1,1e-170,1,1e-160,0,0).gamma,
                NEW_no_history(1,1e-170,1,1e-160,0,0).branch, NEW_no_history(1,1e-170,1,1e-160,0,0).gamma);
    std::printf("  checkpoint gammat*vrel=1e200, Fn=1e190 (overflow of squares): OLD br=%d NEW br=%d\n",
                OLD_no_history(1,1e190,1e100,1e100,0,0).branch, NEW_no_history(1,1e190,1e100,1e100,0,0).branch);
  }
  /* ---------- no_history: random physical range ---------------------- */
  {
    Stat s;
    for (long k=0;k<20000000;k++) {
      const double xmu = U(0,1), Fn = sgn()*logu(-9,4), gammat = logu(-6,3);
      const double v1 = sgn()*logu(-12,2), v2 = sgn()*logu(-12,2), v3 = sgn()*logu(-12,2);
      NH o = OLD_no_history(xmu,Fn,gammat,v1,v2,v3), n = NEW_no_history(xmu,Fn,gammat,v1,v2,v3);
      s.n++; if (o.branch != n.branch) s.branch_mismatch++; acc_val(s,o.gamma,n.gamma);
    }
    report("no_history random physical (2e7)", s);
  }
  /* ---------- no_history: constructed near-ties ---------------------- */
  {
    Stat s;
    for (long k=0;k<5000000;k++) {
      const double gammat = logu(-6,3);
      const double v1 = sgn()*logu(-8,1), v2 = sgn()*logu(-8,1), v3 = sgn()*logu(-8,1);
      const double vrel = sqrt(v1*v1+v2*v2+v3*v3);
      double Ft = gammat*vrel; int steps = (int)(rng()%9) - 4;          // +-4 ulp around the tie
      for (int q=0;q<std::abs(steps);q++) Ft = std::nextafter(Ft, steps>0 ? INF : 0.0);
      NH o = OLD_no_history(1.0,Ft,gammat,v1,v2,v3), n = NEW_no_history(1.0,Ft,gammat,v1,v2,v3);
      s.n++; if (o.branch != n.branch) s.branch_mismatch++; acc_val(s,o.gamma,n.gamma);
    }
    report("no_history near-tie (+-4ulp, 5e6)", s);
  }
  /* ---------- history: random physical, special, near-tie ------------ */
  {
    Stat s, sN;
    for (long k=0;k<20000000;k++) {
      const double kt = logu(1,9), xmu = U(0,1), Fn = sgn()*logu(-6,4), gammat = logu(-4,3);
      double sh[3] = {sgn()*logu(-14,-2), sgn()*logu(-14,-2), sgn()*logu(-14,-2)};
      double vt[3] = {sgn()*logu(-8,1), sgn()*logu(-8,1), sgn()*logu(-8,1)};
      TH o = OLD_history(kt,xmu,Fn,gammat,sh,vt), n = NEW_history(kt,xmu,Fn,gammat,sh,vt);
      s.n++; if (o.branch != n.branch) s.branch_mismatch++;
      for (int c=0;c<3;c++) { acc_val(s,o.Ft[c],n.Ft[c]); acc_val(s,o.shear[c],n.shear[c]); }
    }
    report("history random physical (2e7)", s);
    for (long k=0;k<5000000;k++) {
      const double kt = logu(1,9), gammat = logu(-4,3);
      double sh[3] = {sgn()*logu(-10,-3), sgn()*logu(-10,-3), sgn()*logu(-10,-3)};
      double vt[3] = {1e-3, -2e-3, 0.5e-3};
      double Fs = kt*sqrt(sh[0]*sh[0]+sh[1]*sh[1]+sh[2]*sh[2]);
      int steps = (int)(rng()%9) - 4;
      for (int q=0;q<std::abs(steps);q++) Fs = std::nextafter(Fs, steps>0 ? INF : 0.0);
      TH o = OLD_history(kt,1.0,Fs,gammat,sh,vt), n = NEW_history(kt,1.0,Fs,gammat,sh,vt);
      sN.n++; if (o.branch != n.branch) sN.branch_mismatch++;
      for (int c=0;c<3;c++) { acc_val(sN,o.Ft[c],n.Ft[c]); acc_val(sN,o.shear[c],n.shear[c]); }
    }
    report("history near-tie (+-4ulp, 5e6)", sN);
    Stat sS; long printed=0;
    for (int a=0;a<nsp;a++) for (int b=0;b<nsp;b++) for (int c=0;c<nsp;c++) {
      const double kt = special[a], Fn = special[b];
      double sh[3] = {special[c], 0.5*special[c], 0.0}; double vt[3] = {1.0, 0.0, 0.0};
      TH o = OLD_history(kt,0.5,Fn,1.0,sh,vt), n = NEW_history(kt,0.5,Fn,1.0,sh,vt);
      sS.n++; bool bm = o.branch != n.branch; if (bm) sS.branch_mismatch++;
      bool vd=false; for (int q=0;q<3;q++) { if(!same(o.Ft[q],n.Ft[q])||!same(o.shear[q],n.shear[q])) vd=true;
        acc_val(sS,o.Ft[q],n.Ft[q]); acc_val(sS,o.shear[q],n.shear[q]); }
      if ((bm||vd) && printed < 25) { printed++;
        std::printf("  history diff: kt=%g Fn=%g shear0=%g : OLD br=%d Ft0=%.6g  NEW br=%d Ft0=%.6g\n",
                    kt, Fn, special[c], o.branch, o.Ft[0], n.branch, n.Ft[0]); }
    }
    report("history special-value grid (18^3)", sS);
  }
  /* ---------- luding ---------------------------------------------------- */
  {
    Stat s;
    const double crl_sp[] = {0.0, -0.0, log(0.05), log(0.5), log(0.999999), -INF, NaN, -1e-300};
    for (double crl : crl_sp) for (int td=0; td<2; td++) {
      double a,b,c,d; OLD_luding(1e-3, 1e5, crl, td, a, b); NEW_luding(1e-3, 1e5, crl, td, c, d);
      s.n++; acc_val(s,a,c); acc_val(s,b,d);
    }
    for (long k=0;k<5000000;k++) {
      double a,b,c,d; const double m=logu(-9,1), kn=logu(0,9), crl=log(U(0.05,1.0)); const bool td = rng()&1;
      OLD_luding(m,kn,crl,td,a,b); NEW_luding(m,kn,crl,td,c,d); s.n++; acc_val(s,a,c); acc_val(s,b,d);
    }
    report("luding gamma (5e6 + specials)", s);
  }
  /* ---------- rolling --------------------------------------------------- */
  {
    Stat sc, se;
    auto run = [&](double w1,double w2,double w3){
      double o[3], n[3];
      int bo = OLD_cdt(w1,w2,w3,o), bn = NEW_cdt(w1,w2,w3,n); sc.n++; if (bo!=bn) sc.branch_mismatch++;
      for (int q=0;q<3;q++) acc_val(sc,o[q],n[q]);
      bo = OLD_epsd(w1,w2,w3,o); bn = NEW_epsd(w1,w2,w3,n); se.n++; if (bo!=bn) se.branch_mismatch++;
      for (int q=0;q<3;q++) acc_val(se,o[q],n[q]);
    };
    for (int a=0;a<nsp;a++) for (int b=0;b<nsp;b++) for (int c=0;c<nsp;c++) run(special[a],special[b],special[c]);
    for (long k=0;k<5000000;k++) run(sgn()*logu(-200,200), sgn()*logu(-200,200), sgn()*logu(-200,200));
    report("rolling cdt wrmag>0 vs wrsq>0", sc);
    report("rolling epsd/luding omega!=0 vs omega_sq!=0", se);
  }
  return 0;
}
