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

#ifndef LMP_NEIGHBOR_MULTI_LEVEL_GRID_H
#define LMP_NEIGHBOR_MULTI_LEVEL_GRID_H

#include "pointers.h"

namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
   MultiLevelGrid: size-class binning for granular neighbor lists
   ("neighbor <skin> multi" with a granular pair style)

   - atoms are sorted into nclass radius classes (geometric class edges
     between the smallest and largest radius known at init)
   - every class has its own bin grid, bin size = 1/2 of the class-class
     cutoff (r_c*cdf + skin/2), same rule as style bin uses for cutneighmax
   - for every class pair (ci,cj) there is a stencil on the grid of cj,
     sized by (rmax_ci + rmax_cj)*cdf + skin, i.e. a small atom does not
     scan a big-atom-sized neighbourhood of small atoms and vice versa
   - class membership is decided from the current radius at every build;
     stencils grow automatically if a class holds larger atoms than
     assumed, so the class edges only affect efficiency, never the result
   - the set of pairs found is the same as for style bin: the pair test
     in the builders is identical and every stencil covers the pair cutoff
------------------------------------------------------------------------- */

class MultiLevelGrid : protected Pointers
{
  public:

    MultiLevelGrid(class LAMMPS *);
    ~MultiLevelGrid();

    // define radius classes; nclass_user <= 0 means automatic
    void init(int nclass_user, double rmin, double rmax,
              double cdf, double skin, double cutmax, int dimension);

    // per-class grids and stencils; collective (MPI) call
    void setup(double *bboxlo, double *bboxhi,
               const double *bsubboxlo, const double *bsubboxhi);

    // classify and bin owned + ghost atoms; refresh stencils if needed
    void bin_atoms(int includegroup_bitmask);

    // human-readable summary (rank 0 prints it)
    void print_summary();

    bigint memory_usage();

    inline int classof(double r) const
    {
        int c = 0;
        while (c < nclass-1 && r >= edge[c]) c++;
        return c;
    }

    inline int coord2bin(int c, const double *x) const
    {
        int ix,iy,iz;
        const double *inv = bininv[c];
        const int *nb = nbin[c];
        const int *lo = mbinlo[c];

        if (x[0] >= bboxhi[0])
          ix = static_cast<int> ((x[0]-bboxhi[0])*inv[0]) + nb[0];
        else if (x[0] >= bboxlo[0]) {
          ix = static_cast<int> ((x[0]-bboxlo[0])*inv[0]);
          ix = (ix < nb[0]-1) ? ix : nb[0]-1;
        } else
          ix = static_cast<int> ((x[0]-bboxlo[0])*inv[0]) - 1;

        if (x[1] >= bboxhi[1])
          iy = static_cast<int> ((x[1]-bboxhi[1])*inv[1]) + nb[1];
        else if (x[1] >= bboxlo[1]) {
          iy = static_cast<int> ((x[1]-bboxlo[1])*inv[1]);
          iy = (iy < nb[1]-1) ? iy : nb[1]-1;
        } else
          iy = static_cast<int> ((x[1]-bboxlo[1])*inv[1]) - 1;

        if (dimension == 2) iz = 0;
        else if (x[2] >= bboxhi[2])
          iz = static_cast<int> ((x[2]-bboxhi[2])*inv[2]) + nb[2];
        else if (x[2] >= bboxlo[2]) {
          iz = static_cast<int> ((x[2]-bboxlo[2])*inv[2]);
          iz = (iz < nb[2]-1) ? iz : nb[2]-1;
        } else
          iz = static_cast<int> ((x[2]-bboxlo[2])*inv[2]) - 1;

        return (iz-lo[2])*mbin[c][1]*mbin[c][0] + (iy-lo[1])*mbin[c][0] + (ix-lo[0]);
    }

    int nclass;              // # of radius classes
    double *edge;            // class c holds edge[c-1] <= r < edge[c]
    double rmin,rmax;        // radius range used to place the edges

    // per-class grid
    double (*binsize)[3],(*bininv)[3];
    int (*nbin)[3],(*mbin)[3],(*mbinlo)[3];
    int *mbins,*maxhead;
    int **binhead;

    int *count;              // # of binned atoms (owned+ghost) per class
    double *rcur;            // max radius per class at last binning
    double *rstencil;        // radius the class stencils were built for

    int *bins;               // next atom in same bin of same class
    int *aclass;             // class of each owned/ghost atom
    int maxatom;

    int *used;               // (class,bin) of every bin filled at last
    int nused,maxused;       // binning, so only those need clearing

    int *nstencil;           // class-pair stencils, index ci*nclass+cj
    int *maxstencil;         // (bins of the grid of class cj)
    int **stencil;

    bigint nstencil_rebuild; // # of stencil refreshes due to grown radii

  private:

    int dimension;
    double cdf,skin,cutmax;
    double *bboxlo,*bboxhi;
    int grids_ready;

    void allocate_classes();
    void destroy_classes();
    void create_stencil(int ci, int cj);
    double bin_distance(int c, int i, int j, int k) const;
};

}

#endif
