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
    Contributing author for this file:

    LIGGGHTS modernization branch (roadmap item B6, findings C-18/V-01,
    C-23/V-02)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Restitution mappings used by the opt-in keyword 'correctRestitution'
   of the normal models hooke, hertz and luding.

   Linear spring-dashpot and Hertz/Tsuji dashpot with the attractive
   normal force clipped to zero ('limitForce on'):

     Write the damping ratio zeta = gamma_n/(2 sqrt(m* k_n)) for hooke,
     and zeta = -beta for hertz (gamma_n = -2 sqrt(5/6) beta sqrt(S_n m*)).
     With the clipped force, the realised restitution coefficient is

       e_clip(zeta) = G(zeta)^2,
       G(zeta) = exp(-phi(zeta)),
       phi(zeta) = zeta arccos(zeta)/sqrt(1-zeta^2)      (zeta < 1)
                 = 1                                    (zeta = 1)
                 = zeta arccosh(zeta)/sqrt(zeta^2-1)    (zeta > 1)

     (Schwager & Poeschel, Phys. Rev. E 78 (2008) 051304: the particles
     separate when the force vanishes, before the overlap is zero; one
     factor G comes from the compression phase, one from the restitution
     phase). For the linear model this is exact. For the Hertz law with
     the delta^(1/4) Tsuji dashpot, the same expression reproduces an
     independent integration of the clipped ODE
       x'' = -max(0, x^(3/2) + sqrt(5) zeta x^(1/4) x')
     to a relative 1e-9 for zeta in [0.03, 5] (audit/fixes/phaseB/normal).
     The corrected damping ratio solves e_clip(zeta) = e_input, i.e.
     phi(zeta) = -ln(e_input)/2. phi is strictly increasing, so a
     bisection converges to machine precision.

   Luding (2008) hysteretic law (k1 loading, k2(deltaMax) unloading,
   plasticity limit deltaMaxLim, kc = 0, f0 = 0, force clipped at 0):

     In units k1 = 1, m* = 1, deltaMaxLim = 1, time 1/sqrt(k1/m*), the
     contact is fully determined by kappa = k2max/k1, the damping ratio
     zeta = gamma_n/(2 sqrt(m* k1)) and the normalised impact velocity
     s = v0/(sqrt(k1/m*) deltaMaxLim). The realised restitution is
     velocity-dependent by design (deltaMax/deltaMaxLim grows with v0),
     so zeta(s) is tabulated per type pair (log-spaced s, RK4 integration
     of the continuous law) and looked up once per contact from the
     approach velocity at first touch. If the hysteresis alone already
     dissipates more than requested (e_hys(s) < e_input), zeta = 0.
------------------------------------------------------------------------- */

#ifndef LIGGGHTS_RESTITUTION_MAPPING_H
#define LIGGGHTS_RESTITUTION_MAPPING_H

#include <cmath>
#include <vector>

namespace LIGGGHTS {
namespace ContactModels {
namespace RestitutionMapping {

  // phi(zeta) as defined above, zeta >= 0
  inline double phiClipped(const double zeta)
  {
    if (zeta <= 0.0) return 0.0;
    if (zeta < 1.0)
    {
      const double th = acos(zeta);
      return th > 0.0 ? zeta*th/sin(th) : 1.0;
    }
    if (zeta > 1.0)
    {
      const double u = acosh(zeta);
      return u > 0.0 ? zeta*u/sinh(u) : 1.0;
    }
    return 1.0;
  }

  // restitution of the clipped linear / Hertz-Tsuji dashpot
  inline double eClipped(const double zeta)
  {
    const double g = exp(-phiClipped(zeta));
    return g*g;
  }

  // damping ratio for which the clipped law gives restitution e (0 < e <= 1)
  inline double zetaClipped(const double e)
  {
    if (!(e < 1.0)) return 0.0;
    if (!(e > 0.0)) return HUGE_VAL;
    const double target = -0.5*log(e);
    double lo = 0.0, hi = 1.0;
    while (phiClipped(hi) < target) { lo = hi; hi *= 2.0; }
    for (int it = 0; it < 200; ++it)
    {
      const double mid = 0.5*(lo+hi);
      if (mid <= lo || mid >= hi) break;
      if (phiClipped(mid) < target) lo = mid; else hi = mid;
    }
    return 0.5*(lo+hi);
  }

  // restitution coefficient belonging to the Tsuji beta (beta <= 0):
  // inverse of beta = ln e / sqrt(ln^2 e + pi^2)
  inline double eFromBeta(const double beta)
  {
    if (!(beta < 0.0)) return 1.0;
    if (!(beta > -1.0)) return 0.0;
    return exp(M_PI*beta/sqrt(1.0-beta*beta));
  }

  /* ------------------------------------------------------------------
     Luding law, normalised as described above. Returns the realised
     restitution coefficient of a head-on collision with approach speed s.
  ------------------------------------------------------------------ */

  inline double ludingHysForce(const double x, const double dm, const double kappa)
  {
    double f;
    if (dm >= 1.0)
    {
      f = kappa*(x-1.0) + 1.0;
      if (f < 0.0) f = 0.0;
    }
    else
    {
      const double k2 = 1.0 + (kappa-1.0)*dm;
      const double fTmp = k2*(x-dm) + dm;
      if (fTmp >= x) f = x;
      else f = fTmp > 0.0 ? fTmp : 0.0;
    }
    return f;
  }

  inline double ludingAcc(const double x, const double v, const double dm,
                          const double kappa, const double zeta)
  {
    const double dmx = x > dm ? x : dm;
    const double F = ludingHysForce(x, dmx, kappa) + 2.0*zeta*v;
    return F > 0.0 ? -F : 0.0;
  }

  inline double ludingRestitution(const double kappa, const double zeta, const double s)
  {
    if (!(s > 0.0)) return eClipped(zeta);
    // time step: the fastest oscillation is sqrt(kappa); the dashpot adds
    // a rate 2 zeta, so shrink the step for strong damping
    const double w2 = kappa > 1.0 ? kappa : 1.0;
    const double h = 2e-3/sqrt(w2) / (1.0 + zeta);
    double x = 0.0, v = s, dm = 0.0;
    const long nmax = 50000000L;
    for (long n = 0; n < nmax; ++n)
    {
      const double k1x = v,                 k1v = ludingAcc(x, v, dm, kappa, zeta);
      const double k2x = v + 0.5*h*k1v,     k2v = ludingAcc(x+0.5*h*k1x, v+0.5*h*k1v, dm, kappa, zeta);
      const double k3x = v + 0.5*h*k2v,     k3v = ludingAcc(x+0.5*h*k2x, v+0.5*h*k2v, dm, kappa, zeta);
      const double k4x = v + h*k3v,         k4v = ludingAcc(x+h*k3x, v+h*k3v, dm, kappa, zeta);
      const double xn = x + h/6.0*(k1x + 2.0*k2x + 2.0*k3x + k4x);
      const double vn = v + h/6.0*(k1v + 2.0*k2v + 2.0*k3v + k4v);
      const double dmn = xn > dm ? xn : dm;
      if (vn < 0.0)
      {
        // separation: the clipped force vanishes (and stays zero) or x < 0
        const double Fn = ludingHysForce(xn, dmn, kappa) + 2.0*zeta*vn;
        if (Fn <= 0.0 || xn <= 0.0)
        {
          const double F0 = ludingHysForce(x, dm, kappa) + 2.0*zeta*v;
          double vsep = vn;
          if (Fn <= 0.0 && F0 > 0.0) vsep = v + (vn-v)*F0/(F0-Fn);
          return -vsep/s;
        }
      }
      x = xn; v = vn; dm = dmn;
    }
    return -v/s;
  }

  // damping ratio giving restitution e at normalised impact speed s;
  // 0 if the hysteresis alone already gives e_hys <= e
  inline double ludingZeta(const double kappa, const double e, const double s)
  {
    if (!(e < 1.0)) return 0.0;
    if (kappa <= 1.0 || !(s > 0.0)) return zetaClipped(e);
    const double e0 = ludingRestitution(kappa, 0.0, s);
    if (e0 <= e) return 0.0;
    const double lt = log(e);
    // bracket: e decreases monotonically with zeta
    double lo = 0.0, flo = log(e0) - lt;
    double hi = zetaClipped(e), fhi = log(ludingRestitution(kappa, hi, s)) - lt;
    while (fhi > 0.0) { lo = hi; flo = fhi; hi *= 2.0; fhi = log(ludingRestitution(kappa, hi, s)) - lt; }
    // Illinois false position
    int side = 0;
    double z = hi;
    for (int it = 0; it < 100; ++it)
    {
      z = (lo*fhi - hi*flo)/(fhi - flo);
      const double fz = log(ludingRestitution(kappa, z, s)) - lt;
      if (fabs(fz) < 1e-9 || hi - lo < 1e-12*hi) break;
      if (fz > 0.0) { lo = z; flo = fz; if (side == -1) fhi *= 0.5; side = -1; }
      else          { hi = z; fhi = fz; if (side == +1) flo *= 0.5; side = +1; }
    }
    return z;
  }

  // zeta(s) tabulated on a log grid for one (kappa, e) pair
  class LudingZetaTable
  {
  public:
    LudingZetaTable() : kappa_(-1.0), e_(-1.0) {}

    bool matches(const double kappa, const double e) const
    { return kappa == kappa_ && e == e_; }

    void build(const double kappa, const double e)
    {
      kappa_ = kappa; e_ = e;
      table_.assign(NPTS, 0.0);
      if (kappa <= 1.0)
      {
        const double z = zetaClipped(e);
        for (int i = 0; i < NPTS; ++i) table_[i] = z;
        return;
      }
      for (int i = 0; i < NPTS; ++i)
        table_[i] = ludingZeta(kappa, e, pow(10.0, LOG10_SMIN + i/double(PER_DECADE)));
    }

    double zeta(const double s) const
    {
      if (!(s > 0.0)) return table_[0];
      const double u = (log10(s) - LOG10_SMIN)*PER_DECADE;
      if (u <= 0.0) return table_[0];
      if (u >= NPTS-1) return table_[NPTS-1];
      const int i = int(u);
      const double f = u - i;
      return (1.0-f)*table_[i] + f*table_[i+1];
    }

  private:
    enum { PER_DECADE = 64, NPTS = 8*PER_DECADE + 1, LOG10_SMIN = -4 };
    double kappa_, e_;
    std::vector<double> table_;
  };

}
}
}

#endif
