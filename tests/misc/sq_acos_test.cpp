/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine, released by DCS Computing GmbH.
    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses .

    Contributing author: LIGGGHTS modernization branch
------------------------------------------------------------------------- */

// Degenerate-input test of Superquadric::pre_initial_estimate (phase E,
// misc agent; signfma report: acos(cos_phi) with cos_phi > 1 by ulps).
// Points on the local x-z plane (y = 0) and on/near the local axes of an
// unrotated superquadric: cos_phi = x/sin_theta/r is +-1 exactly in real
// arithmetic, but can exceed 1 after rounding, giving acos = NaN and an
// undefined (out of range) grid index. Prints the number of bad points and
// returns 1 if there are any.
//
// build: g++ -O2 -march=native -DSUPERQUADRIC_ACTIVE_FLAG
//        -DNONSPHERICAL_ACTIVE_FLAG -I<src> sq_acos_test.cpp
//        <src>/superquadric.cpp <src>/math_extra_liggghts_superquadric.cpp
//        <src>/math_extra_liggghts_nonspherical.cpp   (mpicxx: needs mpi.h)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "superquadric.h"
#include "error.h"

// the math helpers report errors through LAMMPS_NS::Error; this test links
// without the LAMMPS library, so provide the one symbol they need
void LAMMPS_NS::Error::one(const char *file, int line, const char *str)
{
  fprintf(stderr, "ERROR: %s (%s:%d)\n", str, file, line);
  abort();
}

int main()
{
  double center[3] = {0., 0., 0.};
  double quat[4] = {1., 0., 0., 0.};
  double shape[3] = {1.e-3, 1.e-3, 1.e-3};       // identical semi-axes
  double blockiness[2] = {2., 2.};
  Superquadric sq(center, quat, shape, blockiness);

  const int nphi = 40, ntheta = 20;
  long nbad = 0, ntot = 0;
  // deterministic sweep: x-z plane (y = +0 and y = -0) and the x axis
  for (int k = 1; k < 200000; k++) {
    const double t = 1.e-3 * (0.05 + 0.9 * k / 200000.);   // polar angle-ish
    for (int s = 0; s < 4; s++) {
      double p[3];
      p[0] = ((s & 1) ? -1. : 1.) * 1.e-3 * sin(3.0 * t / 1.e-3 + 0.1) * (1. + k * 1.e-6);
      p[1] = (s & 2) ? -0.0 : 0.0;
      p[2] = 1.e-3 * cos(3.0 * t / 1.e-3 + 0.1);
      int iphi = -1, itheta = -1;
      sq.pre_initial_estimate(p, nphi, ntheta, &iphi, &itheta);
      ntot++;
      if (iphi < 0 || iphi >= nphi - 1 || itheta < 0 || itheta >= ntheta - 1)
        nbad++;
    }
  }
  // points exactly on the local x axis
  const double xs[] = {1.e-3, 7.3e-4, 3.1e-9, 1.0, 123.456};
  for (double x : xs) {
    double p[3] = {x, 0., 0.};
    int iphi = -1, itheta = -1;
    sq.pre_initial_estimate(p, nphi, ntheta, &iphi, &itheta);
    ntot++;
    if (iphi < 0 || iphi >= nphi - 1 || itheta < 0 || itheta >= ntheta - 1)
      nbad++;
  }
  printf("sq_acos_test: %ld of %ld degenerate points gave an invalid (NaN) grid index\n", nbad, ntot);
  return nbad ? 1 : 0;
}
