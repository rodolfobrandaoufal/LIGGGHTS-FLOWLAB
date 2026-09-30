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
    LIGGGHTS modernization branch (size-class multi-level grid for
    polydisperse granular neighbor lists, roadmap C4 / finding S-09)

    Idea of per-collection bins and collection-pair stencils as in
    LAMMPS "neighbor multi" (Shire, Hanley & Stratford, Comp. Part. Mech.
    8 (2021) 653); implementation written for LIGGGHTS granular lists.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "neigh_multi_level_grid.h"
#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "mpi_liggghts.h"

using namespace LAMMPS_NS;

#define MAXCLASS 8
#define SMALL 1.0e-6
#define CUT2BIN_RATIO 100
#define STENCIL_GROWTH 1.02    // headroom when a class stencil must grow

/* ---------------------------------------------------------------------- */

MultiLevelGrid::MultiLevelGrid(LAMMPS *lmp) : Pointers(lmp),
  nclass(0), edge(NULL), rmin(0.), rmax(0.),
  binsize(NULL), bininv(NULL), nbin(NULL), mbin(NULL), mbinlo(NULL),
  mbins(NULL), maxhead(NULL), binhead(NULL),
  count(NULL), rcur(NULL), rstencil(NULL),
  bins(NULL), aclass(NULL), maxatom(0),
  used(NULL), nused(0), maxused(0),
  nstencil(NULL), maxstencil(NULL), stencil(NULL),
  nstencil_rebuild(0),
  dimension(3), cdf(1.), skin(0.), cutmax(0.),
  bboxlo(NULL), bboxhi(NULL), grids_ready(0)
{
}

/* ---------------------------------------------------------------------- */

MultiLevelGrid::~MultiLevelGrid()
{
  destroy_classes();
  memory->destroy(bins);
  memory->destroy(aclass);
  memory->destroy(used);
}

/* ---------------------------------------------------------------------- */

void MultiLevelGrid::destroy_classes()
{
  if (nclass == 0) return;
  for (int c = 0; c < nclass; c++) memory->destroy(binhead[c]);
  for (int p = 0; p < nclass*nclass; p++) memory->destroy(stencil[p]);
  delete [] edge;
  delete [] binsize; delete [] bininv;
  delete [] nbin; delete [] mbin; delete [] mbinlo;
  delete [] mbins; delete [] maxhead; delete [] binhead;
  delete [] count; delete [] rcur; delete [] rstencil;
  delete [] nstencil; delete [] maxstencil; delete [] stencil;
  edge = NULL; binhead = NULL; stencil = NULL;
  nclass = 0;
  grids_ready = 0;
}

/* ---------------------------------------------------------------------- */

void MultiLevelGrid::allocate_classes()
{
  edge = new double[nclass > 1 ? nclass-1 : 1];
  binsize = new double[nclass][3];
  bininv = new double[nclass][3];
  nbin = new int[nclass][3];
  mbin = new int[nclass][3];
  mbinlo = new int[nclass][3];
  mbins = new int[nclass];
  maxhead = new int[nclass];
  binhead = new int*[nclass];
  count = new int[nclass];
  rcur = new double[nclass];
  rstencil = new double[nclass];
  nstencil = new int[nclass*nclass];
  maxstencil = new int[nclass*nclass];
  stencil = new int*[nclass*nclass];
  for (int c = 0; c < nclass; c++) {
    mbins[c] = maxhead[c] = 0;
    binhead[c] = NULL;
    count[c] = 0;
    rcur[c] = rstencil[c] = 0.;
  }
  for (int p = 0; p < nclass*nclass; p++) {
    nstencil[p] = maxstencil[p] = 0;
    stencil[p] = NULL;
  }
}

/* ----------------------------------------------------------------------
   define the radius classes
   automatic: one class per factor of ~2 in radius, at most MAXCLASS;
   a (nearly) monodisperse system gets a single class
------------------------------------------------------------------------- */

void MultiLevelGrid::init(int nclass_user, double rmin_in, double rmax_in,
                          double cdf_in, double skin_in, double cutmax_in,
                          int dimension_in)
{
  cdf = cdf_in;
  skin = skin_in;
  cutmax = cutmax_in;
  dimension = dimension_in;

  double lo = rmin_in, hi = rmax_in;
  int valid = (lo > 0. && hi > 0. && lo <= hi);
  if (!valid) lo = hi = (hi > 0. ? hi : 0.);

  int n;
  if (!valid || hi <= lo) n = 1;
  else if (nclass_user > 0) n = nclass_user;
  else {
    n = static_cast<int>(ceil(log(hi/lo)/log(2.) - 1.e-9));
    if (n < 1) n = 1;
  }
  if (n > MAXCLASS) n = MAXCLASS;

  if (n != nclass || lo != rmin || hi != rmax) {
    destroy_classes();
    nclass = n;
    rmin = lo;
    rmax = hi;
    allocate_classes();
    // geometric edges between rmin and rmax
    for (int c = 0; c < nclass-1; c++)
      edge[c] = rmin*pow(rmax/rmin,static_cast<double>(c+1)/nclass);
  }
  grids_ready = 0;
}

/* ----------------------------------------------------------------------
   closest distance between central bin (0,0,0) and bin (i,j,k) of class c
------------------------------------------------------------------------- */

double MultiLevelGrid::bin_distance(int c, int i, int j, int k) const
{
  double delx,dely,delz;

  if (i > 0) delx = (i-1)*binsize[c][0];
  else if (i == 0) delx = 0.0;
  else delx = (i+1)*binsize[c][0];

  if (j > 0) dely = (j-1)*binsize[c][1];
  else if (j == 0) dely = 0.0;
  else dely = (j+1)*binsize[c][1];

  if (k > 0) delz = (k-1)*binsize[c][2];
  else if (k == 0) delz = 0.0;
  else delz = (k+1)*binsize[c][2];

  return (delx*delx + dely*dely + delz*delz);
}

/* ----------------------------------------------------------------------
   stencil for owned atoms of class ci searching atoms of class cj,
   expressed as bin offsets in the grid of class cj
   cut covers every pair ci-cj; it is capped at cutmax (= cutneighmax),
   the search range style bin uses, which keeps the stencil of an owned
   atom inside the local bin range (comm cutghost >= cutneighmax)
------------------------------------------------------------------------- */

void MultiLevelGrid::create_stencil(int ci, int cj)
{
  double cut = ((rstencil[ci]+rstencil[cj])*cdf + skin)*(1.+SMALL);
  if (cut > cutmax) cut = cutmax;
  const double cutsq = cut*cut;

  int sx = static_cast<int> (cut*bininv[cj][0]);
  if (sx*binsize[cj][0] < cut) sx++;
  int sy = static_cast<int> (cut*bininv[cj][1]);
  if (sy*binsize[cj][1] < cut) sy++;
  int sz = static_cast<int> (cut*bininv[cj][2]);
  if (sz*binsize[cj][2] < cut) sz++;
  if (dimension == 2) sz = 0;

  const int p = ci*nclass+cj;
  const int smax = (2*sx+1)*(2*sy+1)*(2*sz+1);
  if (smax > maxstencil[p]) {
    maxstencil[p] = smax;
    memory->destroy(stencil[p]);
    memory->create(stencil[p],smax,"neigh_mlg:stencil");
  }

  const int mx = mbin[cj][0], my = mbin[cj][1];
  int n = 0;
  int *s = stencil[p];
  for (int k = -sz; k <= sz; k++)
    for (int j = -sy; j <= sy; j++)
      for (int i = -sx; i <= sx; i++)
        if (bin_distance(cj,i,j,k) < cutsq)
          s[n++] = k*my*mx + j*mx + i;
  nstencil[p] = n;
}

/* ----------------------------------------------------------------------
   per-class bin grids (same construction as Neighbor::setup_bins)
   plus stencils for all class pairs
   collective: the nominal class radius is the global max radius of the
   atoms currently in that class (class edge if the class is empty)
------------------------------------------------------------------------- */

void MultiLevelGrid::setup(double *bboxlo_in, double *bboxhi_in,
                           const double *bsubboxlo, const double *bsubboxhi)
{
  bboxlo = bboxlo_in;
  bboxhi = bboxhi_in;

  // nominal radius per class from current owned atoms

  double *rloc = new double[nclass];
  double *rnom = new double[nclass];
  for (int c = 0; c < nclass; c++) rloc[c] = 0.;
  double *radius = atom->radius;
  const int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    const int c = classof(radius[i]);
    if (radius[i] > rloc[c]) rloc[c] = radius[i];
  }
  MPI_Allreduce(rloc,rnom,nclass,MPI_DOUBLE,MPI_MAX,world);
  for (int c = 0; c < nclass; c++) {
    if (rnom[c] <= 0.) rnom[c] = (c < nclass-1) ? edge[c] : rmax;
    if (rnom[c] <= 0.) rnom[c] = 0.5*cutmax;
  }

  double bbox[3];
  bbox[0] = bboxhi[0] - bboxlo[0];
  bbox[1] = bboxhi[1] - bboxlo[1];
  bbox[2] = bboxhi[2] - bboxlo[2];

  for (int c = 0; c < nclass; c++) {

    // bin size = full class-class cutoff 2*r_c*cdf + skin
    // (style bin uses 1/2 of cutneighmax; granular cutoffs are ~1 diameter,
    // so half-cutoff bins hold ~0.15 atoms and the 5x5x5 stencil is mostly
    // empty bins - measured 23-36% slower builds, see C4 report)
    // a coarser class is then scanned by any finer class with a 3x3x3
    // stencil, since cut(c,c) >= cut(c',c) for c' < c

    double bs = 2.*rnom[c]*cdf + skin;
    if (bs > cutmax && cutmax > 0.) bs = cutmax;
    if (bs <= 0.) bs = bbox[0];
    const double inv = 1.0/bs;

    if (bbox[0]*inv > MAXSMALLINT || bbox[1]*inv > MAXSMALLINT ||
        bbox[2]*inv > MAXSMALLINT)
      error->all(FLERR,"Domain too large for neighbor multi bins");

    for (int d = 0; d < 3; d++) {
      nbin[c][d] = static_cast<int> (bbox[d]*inv);
      if (dimension == 2 && d == 2) nbin[c][d] = 1;
      if (nbin[c][d] == 0) nbin[c][d] = 1;
      binsize[c][d] = bbox[d]/nbin[c][d];
      bininv[c][d] = 1.0/binsize[c][d];
    }
    // a class grid finer than CUT2BIN_RATIO would mean a huge number of bins
    // (flat non-periodic system) - fall back to a coarser class bin

    for (int d = 0; d < dimension; d++)
      if (bs*bininv[c][d] > CUT2BIN_RATIO)
        error->all(FLERR,"Cannot use neighbor multi bins - box size << cutoff");

    int hi[3];
    for (int d = 0; d < 3; d++) {
      if (dimension == 2 && d == 2) { mbinlo[c][d] = hi[d] = 0; continue; }
      double coord = bsubboxlo[d] - SMALL*bbox[d];
      mbinlo[c][d] = static_cast<int> ((coord-bboxlo[d])*bininv[c][d]);
      if (coord < bboxlo[d]) mbinlo[c][d] = mbinlo[c][d] - 1;
      coord = bsubboxhi[d] + SMALL*bbox[d];
      hi[d] = static_cast<int> ((coord-bboxlo[d])*bininv[c][d]);
      mbinlo[c][d] -= 1;
      hi[d] += 1;
    }
    for (int d = 0; d < 3; d++) mbin[c][d] = hi[d] - mbinlo[c][d] + 1;

    bigint bbin = ((bigint) mbin[c][0]) * ((bigint) mbin[c][1]) *
      ((bigint) mbin[c][2]);
    if (bbin > MAXSMALLINT) error->one(FLERR,"Too many neighbor multi bins");
    mbins[c] = bbin;
    if (mbins[c] > maxhead[c]) {
      maxhead[c] = mbins[c];
      memory->destroy(binhead[c]);
      memory->create(binhead[c],maxhead[c],"neigh_mlg:binhead");
    }
    for (int b = 0; b < mbins[c]; b++) binhead[c][b] = -1;
    count[c] = 0;
    nused = 0;

    rstencil[c] = rnom[c];
  }

  for (int ci = 0; ci < nclass; ci++)
    for (int cj = 0; cj < nclass; cj++)
      create_stencil(ci,cj);

  delete [] rloc;
  delete [] rnom;
  grids_ready = 1;
}

/* ----------------------------------------------------------------------
   sort owned and ghost atoms into their class grids
   binned in reverse order so each bin list is in forward order
   if a class holds atoms larger than its stencils assume, the stencils
   of that class are rebuilt (local operation, no communication)
------------------------------------------------------------------------- */

void MultiLevelGrid::bin_atoms(int bitmask)
{
  if (!grids_ready)
    error->one(FLERR,"Neighbor multi bins used before setup");

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  // clear the bins filled at the previous binning (O(atoms), not O(bins))

  for (int u = 0; u < nused; u++)
    binhead[used[2*u]][used[2*u+1]] = -1;
  nused = 0;

  if (atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(bins);
    memory->destroy(aclass);
    memory->create(bins,maxatom,"neigh_mlg:bins");
    memory->create(aclass,maxatom,"neigh_mlg:aclass");
  }
  if (nall > maxused) {
    maxused = atom->nmax > nall ? atom->nmax : nall;
    memory->destroy(used);
    memory->create(used,2*maxused,"neigh_mlg:used");
  }

  for (int c = 0; c < nclass; c++) {
    count[c] = 0;
    rcur[c] = 0.;
  }

  double **x = atom->x;
  double *radius = atom->radius;
  int *mask = atom->mask;
  const int nfirst = bitmask ? atom->nfirst : nlocal;

  for (int i = nall-1; i >= 0; i--) {
    if (i >= nlocal) {
      if (bitmask && !(mask[i] & bitmask)) continue;
    } else if (i >= nfirst) continue;
    const double r = radius[i];
    const int c = classof(r);
    aclass[i] = c;
    const int ibin = coord2bin(c,x[i]);
    if (binhead[c][ibin] < 0) {
      used[2*nused] = c;
      used[2*nused+1] = ibin;
      nused++;
    }
    bins[i] = binhead[c][ibin];
    binhead[c][ibin] = i;
    count[c]++;
    if (r > rcur[c]) rcur[c] = r;
  }

  // grow stencils where the actual radii exceed the assumed ones

  int grow = 0;
  for (int c = 0; c < nclass; c++)
    if (rcur[c] > rstencil[c]) {
      rstencil[c] = rcur[c]*STENCIL_GROWTH;
      grow = 1;
    }
  if (grow) {
    nstencil_rebuild++;
    for (int ci = 0; ci < nclass; ci++)
      for (int cj = 0; cj < nclass; cj++)
        create_stencil(ci,cj);
  }
}

/* ---------------------------------------------------------------------- */

void MultiLevelGrid::print_summary()
{
  if (comm->me != 0) return;
  char line[512];
  sprintf(line,"neighbor multi (granular): %d radius class(es), "
          "r = %g .. %g\n",nclass,rmin,rmax);
  if (screen) fputs(line,screen);
  if (logfile) fputs(line,logfile);
  for (int c = 0; c < nclass; c++) {
    const double lo = (c == 0) ? 0. : edge[c-1];
    if (c < nclass-1)
      sprintf(line,"  class %d: %g <= r < %g, bin size %g, own-class stencil %d bins\n",
              c,lo,edge[c],grids_ready ? binsize[c][0] : 0.,
              grids_ready ? nstencil[c*nclass+c] : 0);
    else
      sprintf(line,"  class %d: r >= %g, bin size %g, own-class stencil %d bins\n",
              c,lo,grids_ready ? binsize[c][0] : 0.,
              grids_ready ? nstencil[c*nclass+c] : 0);
    if (screen) fputs(line,screen);
    if (logfile) fputs(line,logfile);
  }
}

/* ---------------------------------------------------------------------- */

bigint MultiLevelGrid::memory_usage()
{
  bigint bytes = 0;
  bytes += 2*memory->usage(bins,maxatom);
  bytes += memory->usage(used,2*maxused);
  for (int c = 0; c < nclass; c++) bytes += memory->usage(binhead[c],maxhead[c]);
  for (int p = 0; p < nclass*nclass; p++)
    bytes += memory->usage(stencil[p],maxstencil[p]);
  return bytes;
}
