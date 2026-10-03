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
    LIGGGHTS modernization branch

    OpenMP threaded kernel of pair gran (roadmap C2) with deterministic
    (thread-count independent) force summation (roadmap B8). Per-thread
    accumulation mode after the LAMMPS OPENMP package
    (pair_gran_hooke_history_omp.cpp, thr_omp.cpp; Axel Kohlmeyer).

    The kernel is a member of PairStyles::Granular<> but is defined and
    explicitly instantiated in this translation unit only. lammps.cpp, which
    instantiates all contact models, therefore compiles the serial kernel
    exactly as without OpenMP (GCC's unit-wide inlining budget, and with it
    FMA contraction, would otherwise change).
------------------------------------------------------------------------- */

#ifdef LIGGGHTS_OMP

#include <omp.h>
#include <mpi.h>
#include <cmath>
#include <string.h>
#include <string>
#include <vector>
#include "pair.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "fix_insert_stream.h"
#include "granular_pair_style.h"
#include "contact_models.h"
#include "style_surface_model.h"
#include "style_normal_model.h"
#include "style_tangential_model.h"
#include "style_cohesion_model.h"
#include "style_rolling_model.h"
#include "pair_gran_base.h"
#include "thr_granular.h"

namespace LIGGGHTS {

namespace ThrGranular {

struct PairState {
  struct Scratch {
    ContactModels::SurfacesIntersectData *sidata;
    ContactModels::ForceData *i_forces;
    ContactModels::ForceData *j_forces;
  };
  std::vector<Scratch> scratch;   // thread-private per-contact data
  PairTranspose tr;               // per-atom contribution lists
  std::vector<double> slot;       // 12 per pair slot: dF_i, dT_i, dF_j, dT_j
  std::vector<char> flag;         // 1 = pair slot has a force update
  std::vector<double> delta;      // 3 per pair slot (energy/virial steps only)
  std::vector<double> fbuf;       // per-thread mode: nthreads x natom x 6

  // blocked deterministic mode: thread t owns the contiguous list rows
  // [rowstart[t], rowstart[t+1]); only pairs whose j belongs to another
  // thread's rows ("cross pairs") are buffered
  struct Blocked {
    int64_t key_ncalls;
    int key_inum, key_nlocal, key_nall, key_nthr;
    const void *key_list;
    bool ok;                      // every local j has a list row; newton off
    std::vector<int> offset;      // inum+1: first pair slot of row ii
    std::vector<int> rowof;       // nlocal: list row of atom k (-1: none)
    std::vector<int> rowblock;    // inum: thread block of row ii
    std::vector<int> rowstart;    // nthr+1
    std::vector<int> cidx;        // per pair slot: cross index or -1
    std::vector<int> c_ii, c_jj;  // cross pairs in row order
    std::vector<int> cstart;      // nthr+1: cross pairs of the rows of thread t
    std::vector<int> in_start;    // 2*nthr+1: incoming j-side contributions of block t,
    std::vector<int> in_list;     //   group 2t = from earlier rows, 2t+1 = from later rows
    Blocked() : key_ncalls(-1), key_inum(-1), key_nlocal(-1), key_nall(-1),
      key_nthr(-1), key_list(0), ok(false) {}
  } blk;
  std::vector<double> cbuf;       // 12 per cross pair
  std::vector<char> cflag;

  void build_blocked(int64_t ncalls, int inum, const int *ilist, const int *numneigh,
                     int * const *firstneigh, int nlocal, int nall, int newton,
                     const void *list, int nthr)
  {
    Blocked & b = blk;
    b.key_ncalls = ncalls; b.key_inum = inum; b.key_nlocal = nlocal;
    b.key_nall = nall; b.key_nthr = nthr; b.key_list = list;
    b.ok = !newton;
    b.offset.resize(inum+1);
    b.offset[0] = 0;
    for (int ii = 0; ii < inum; ii++) b.offset[ii+1] = b.offset[ii] + numneigh[ilist[ii]];
    const int npair = b.offset[inum];
    b.rowof.assign(nlocal, -1);
    for (int ii = 0; ii < inum; ii++) b.rowof[ilist[ii]] = ii;
    // contiguous row blocks with about the same number of pairs
    b.rowstart.assign(nthr+1, inum);
    b.rowstart[0] = 0;
    b.rowblock.resize(inum);
    for (int ii = 0, t = 0; ii < inum; ii++) {
      while (t+1 < nthr && (int64_t)b.offset[ii] * nthr >= (int64_t)npair * (t+1)) b.rowstart[++t] = ii;
      b.rowblock[ii] = t;
    }
    for (int t = 1; t <= nthr; t++) if (b.rowstart[t] < b.rowstart[t-1]) b.rowstart[t] = b.rowstart[t-1];
    b.cidx.assign(npair, -1);
    b.c_ii.clear(); b.c_jj.clear();
    b.cstart.assign(nthr+1, 0);
    for (int ii = 0; ii < inum && b.ok; ii++) {
      const int t = b.rowblock[ii];
      const int *jlist = firstneigh[ilist[ii]];
      const int jnum = numneigh[ilist[ii]];
      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;
        if (j >= nlocal) continue;               // ghost: no j-side force (newton off)
        if (b.rowof[j] < 0) { b.ok = false; break; }
        if (b.rowblock[b.rowof[j]] == t) continue;
        b.cidx[b.offset[ii]+jj] = (int)b.c_ii.size();
        b.c_ii.push_back(ii);
        b.c_jj.push_back(jj);
        b.cstart[t+1]++;
      }
    }
    if (!b.ok) return;
    for (int t = 0; t < nthr; t++) b.cstart[t+1] += b.cstart[t];
    // incoming j-side contributions per block, in row order
    const int nc = b.c_ii.size();
    std::vector<int> grp(nc);
    b.in_start.assign(2*nthr+1, 0);
    for (int c = 0; c < nc; c++) {
      const int ii = b.c_ii[c];
      const int j = firstneigh[ilist[ii]][b.c_jj[c]] & NEIGHMASK;
      const int u = b.rowblock[b.rowof[j]];
      grp[c] = b.rowblock[ii] < u ? 2*u : 2*u+1;
      b.in_start[grp[c]+1]++;
    }
    for (int g = 0; g < 2*nthr; g++) b.in_start[g+1] += b.in_start[g];
    b.in_list.resize(nc);
    std::vector<int> fillp(b.in_start.begin(), b.in_start.end()-1);
    for (int c = 0; c < nc; c++) b.in_list[fillp[grp[c]]++] = c;
  }

  void free_scratch()
  {
    for (size_t t = 0; t < scratch.size(); t++) {
      if (scratch[t].sidata) aligned_free(scratch[t].sidata);
      if (scratch[t].i_forces) aligned_free(scratch[t].i_forces);
      if (scratch[t].j_forces) aligned_free(scratch[t].j_forces);
    }
    scratch.clear();
  }
  ~PairState() { free_scratch(); }
  void invalidate() { tr.invalidate(); blk.key_ncalls = -1; }
};

void pair_state_free(PairState *state)
{
  delete state;
}

}

namespace PairStyles {

using namespace ContactModels;
using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   Each list row ii (atom i) is handled by one thread, so the contact
   history of pair (i,j), stored in row i, is written by one thread only
   (with newton off every pair is stored once per rank). Per-contact scratch
   (sidata, i_forces, j_forces) is thread-private.

   deterministic mode: the pair loop stores the contribution of each pair in
   a slot buffer; a second loop over atoms adds them to f/torque in the order
   of the serial loop (ThrGranular::PairTranspose). The result does not
   depend on the thread count and equals the serial kernel bitwise.
   per-thread mode ("package omp N deterministic no"): per-thread f/torque
   arrays, static schedule, reduction in thread order (LAMMPS OPENMP).
   A step with an energy/virial tally always uses the deterministic path and
   tallies serially in pair order.
------------------------------------------------------------------------- */

template<typename ContactModel>
bool Granular<ContactModel>::compute_force_thr(PairGran * pg, int eflag, int vflag, int addflag)
{
  const ThrGranular::Config cfg = ThrGranular::config(lmp);

  // Neighbor::init() restarts neighbor->ncalls in every run: the setup
  // pass (lists always rebuilt) invalidates the cached per-atom lists
  if (update->setupflag && thr_state_) thr_state_->invalidate();

  // finding S-17: with the velocity predictor one thread also runs this
  // kernel, so that results are bitwise identical for any thread count
  if (cfg.nthreads <= 1 && !velocity_predictor_dv(addflag)) return false;

  // conditions under which the serial kernel runs
  std::string why;
  if (update->setupflag) why = "setup";
  else {
    why = ThrGranular::unsafe_reason(pg);
    if (!why.empty()) ;
    else if (pg->cpl() && addflag) why = "compute pair/gran/local";
    else if (pg->storeContactForces()) why = "per-contact force storage (fix contactproperty/atom)";
    else if (pg->storeContactForcesStress()) why = "per-contact force storage for stress";
    else if (pg->storeSumDelta()) why = "multicontact delta storage";
    else if (pg->store_sum_normal_force()) why = "sum of normal forces per atom";
    else if (atom->superquadric_flag) why = "superquadric particles";
    else if (atom->shapetype_flag) why = "non-spherical (convex) particles";
    else if (modify->n_fixes_style("insert/stream/predefined") > 0) why = "fix insert/stream/predefined";
    else if (synchronized_verlet_) why = "synchronized_verlet";
    if (!why.empty()) ThrGranular::fallback_warning(lmp, pg, "pair gran", why);
  }
  if (!why.empty()) return false;

  if (!thr_state_) thr_state_ = new ThrGranular::PairState();
  ThrGranular::PairState & st = *thr_state_;

  if (eflag || vflag)
    pg->ev_setup(eflag, vflag);
  else
    pg->evflag = pg->vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  // contiguous per-atom 3-vectors: index the data block directly (pair_gran_base.h)
  const double * const x0 = x ? x[0] : NULL;
  double * const v0 = v ? v[0] : NULL;
  double * const f0 = f ? f[0] : NULL;
  double * const om0 = omega ? omega[0] : NULL;
  double * const t0 = torque ? torque[0] : NULL;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  const int newton_pair = force->newton_pair;
  const int inum = pg->list->inum;
  int * const ilist = pg->list->ilist;
  int * const numneigh = pg->list->numneigh;
  int ** const firstneigh = pg->list->firstneigh;
  int ** const first_contact_flag = pg->listgranhistory ? pg->listgranhistory->firstneigh : NULL;
  double ** const first_contact_hist = pg->listgranhistory ? pg->listgranhistory->firstdouble : NULL;

  const int dnum = pg->dnum();
  const int freeze_group_bit = pg->freeze_group_bit();
  const int computeflag = pg->computeflag();
  const int shearupdate = pg->shearupdate();
  const int evflag = pg->evflag;
  const bool det = cfg.deterministic || evflag;
  const int nthr = cfg.nthreads;
  const int chunk = cfg.chunk;
  const bool sphere_flag = atom->sphere_flag;
  const bool shapetype_flag = atom->shapetype_flag;
  const double contactDistanceMultiplier = neighbor->contactDistanceFactor*neighbor->contactDistanceFactor;
  double ** const vpred_dv = velocity_predictor_dv(addflag);   // finding S-17

  // blocked deterministic mode (default); the slot mode below is used for
  // energy/virial steps, newton on, or lists without a row for every atom
  bool blocked = false;
  if (det && !evflag && !newton_pair && chunk == 0) {
    ThrGranular::PairState::Blocked & b = st.blk;
    if (b.key_ncalls != neighbor->ncalls || b.key_inum != inum || b.key_nlocal != nlocal ||
        b.key_nall != nall || b.key_nthr != nthr || b.key_list != (const void*)pg->list)
      st.build_blocked(neighbor->ncalls, inum, ilist, numneigh, firstneigh, nlocal, nall,
                       newton_pair, pg->list, nthr);
    blocked = b.ok;
  }

  // pair slots and per-atom contribution lists (deterministic slot mode)
  const bool slotmode = det && !blocked;
  ThrGranular::PairTranspose & tr = st.tr;
  if (slotmode && !tr.valid(neighbor->ncalls, inum, nlocal, nall, newton_pair, pg->list))
    tr.build(neighbor->ncalls, inum, ilist, numneigh, firstneigh, nlocal, nall, newton_pair, pg->list, NEIGHMASK);
  const int npair = slotmode ? tr.offset[inum] : 0;
  if (slotmode) {
    if ((int)st.flag.size() < npair) st.flag.resize(npair);
    if (st.slot.size() < 12*(size_t)npair) st.slot.resize(12*(size_t)npair);
    if (evflag && st.delta.size() < 3*(size_t)npair) st.delta.resize(3*(size_t)npair);
  }
  const int * const offset = slotmode ? &tr.offset[0] : NULL;
  char * const pflag = slotmode && npair > 0 ? &st.flag[0] : NULL;
  double * const slot = slotmode && npair > 0 ? &st.slot[0] : NULL;
  double * const pdelta = slotmode && evflag && npair > 0 ? &st.delta[0] : NULL;

  const int ncross = blocked ? (int)st.blk.c_ii.size() : 0;
  if (blocked) {
    if ((int)st.cflag.size() < ncross) st.cflag.resize(ncross);
    if (st.cbuf.size() < 12*(size_t)ncross) st.cbuf.resize(12*(size_t)ncross);
  }

  // per-thread force arrays (per-thread mode)
  const int natom = newton_pair ? nall : nlocal;
  if (!det && st.fbuf.size() < (size_t)nthr*6*natom) st.fbuf.resize((size_t)nthr*6*natom);
  double * const fbuf = !det && natom > 0 ? &st.fbuf[0] : NULL;

  if ((int)st.scratch.size() != nthr) {
    st.free_scratch();
    ThrGranular::PairState::Scratch empty = { NULL, NULL, NULL };
    st.scratch.assign(nthr, empty);
  }

  #pragma omp parallel num_threads(nthr)
  {
    const int tid = omp_get_thread_num();
    ThrGranular::PairState::Scratch & scr = st.scratch[tid];
    if (!scr.sidata) {
      scr.sidata = aligned_malloc<SurfacesIntersectData>(32);
      scr.i_forces = aligned_malloc<ForceData>(32);
      scr.j_forces = aligned_malloc<ForceData>(32);
    }
    // same initial state as the serial kernel
    memset((void*)scr.sidata, 0, sizeof(SurfacesIntersectData));
    memset((void*)scr.i_forces, 0, sizeof(ForceData));
    memset((void*)scr.j_forces, 0, sizeof(ForceData));
    scr.sidata->area_ratio = 1.0;

    SurfacesIntersectData & sidata = *scr.sidata;
    ForceData & i_forces = *scr.i_forces;
    VelocityPredictorScratch vpred;   // finding S-17, per thread
    ForceData & j_forces = *scr.j_forces;
    sidata.is_wall = false;
    sidata.computeflag = computeflag;
    sidata.shearupdate = shearupdate;

    double * const myf = fbuf ? fbuf + (size_t)tid*6*natom : NULL;
    if (myf) memset(myf, 0, sizeof(double)*6*natom);

    // sinks of a pair contribution
    enum { M_SLOT, M_THREAD, M_CROSS, M_BLOCK };
    ThrGranular::PairState::Blocked & bk = st.blk;
    const int * const bk_offset = blocked ? &bk.offset[0] : NULL;
    const int * const bk_cidx = blocked && !bk.cidx.empty() ? &bk.cidx[0] : NULL;
    double * const cbuf = ncross > 0 ? &st.cbuf[0] : NULL;
    char * const cflag = ncross > 0 ? &st.cflag[0] : NULL;
    int cur_cross = -1;

    // pairs jj0..jj1-1 of list row ii (jj1 < 0: whole row), with the
    // arithmetic of compute_force_serial()
    auto row = [&](const int ii, const int jj0, const int jj1, const int mode)
    {
      const int i = ilist[ii];
      const double xtmp = x0[3*i];
      const double ytmp = x0[3*i+1];
      const double ztmp = x0[3*i+2];
      const double radi = radius[i];
      int * const contact_flags = first_contact_flag ? first_contact_flag[i] : NULL;
      double * const all_contact_hist = first_contact_hist ? first_contact_hist[i] : NULL;
      int * const jlist = firstneigh[i];
      const int jnum = jj1 < 0 ? numneigh[i] : jj1;
      const int base = mode == M_SLOT ? offset[ii] : (mode == M_BLOCK ? bk_offset[ii] : 0);

      sidata.i = i;
      sidata.radi = radi;

      for (int jj = jj0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;

        if (mode == M_BLOCK && bk_cidx) {
          // cross pair: evaluated in the first phase, add its i-side now
          const int c = bk_cidx[base + jj];
          if (c >= 0) {
            if (computeflag && cflag[c]) {
              const double relax_i = pg->relax(i);
              const double * const sc = cbuf + 12*(size_t)c;
              for (int coord = 0; coord < 3; coord++) {
                f[i][coord] += relax_i*sc[coord];
                torque[i][coord] += relax_i*sc[3+coord];
              }
            }
            continue;
          }
        }

        const double xj = x0[3*j];
        const double yj = x0[3*j+1];
        const double zj = x0[3*j+2];
        const double delx = xtmp - xj;
        const double dely = ytmp - yj;
        const double delz = ztmp - zj;
        const double rsq = delx * delx + dely * dely + delz * delz;
        const double radj = radius[j];

        sidata.radj = radj;
        const double radsum = radi + radj;

        sidata.j = j;
        sidata.delta[0] = delx;
        sidata.delta[1] = dely;
        sidata.delta[2] = delz;
        sidata.rsq = rsq;
        sidata.radsum = radsum;
        sidata.contact_flags = contact_flags ? &contact_flags[jj] : NULL;
        sidata.contact_history = all_contact_hist ? &all_contact_hist[dnum*jj] : NULL;

        i_forces.reset();
        j_forces.reset();

        #ifdef SUPERQUADRIC_ACTIVE_FLAG
        if (rmass) {
          sidata.mi = rmass[i];
          sidata.mj = rmass[j];
        } else {
          sidata.mi = mass[type[i]];
          sidata.mj = mass[type[j]];
        }
        sidata.omega_i = om0+3*i;
        sidata.omega_j = om0+3*j;
        #endif

        sidata.v_i     = v0+3*i;
        sidata.v_j     = v0+3*j;
        const int itype = type[i];
        const int jtype = type[j];
        sidata.itype = itype;
        sidata.jtype = jtype;

        if (rsq < radsum * radsum && cmodel.checkSurfaceIntersect(sidata)) {
          const double r = sqrt(rsq);
          const double rinv = 1.0 / r;

          const double enx_sphere = delx * rinv;
          const double eny_sphere = dely * rinv;
          const double enz_sphere = delz * rinv;

          double mi, mj;

          if (rmass) {
            mi = rmass[i];
            mj = rmass[j];
          } else {
            mi = mass[itype];
            mj = mass[jtype];
          }
          if (pg->fr_pair()) {
            const double * mass_rigid = pg->mr_pair();
            if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
            if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
          }

          double meff = mi * mj / (mi + mj);
          if (mask[i] & freeze_group_bit)
            meff = mj;
          if (mask[j] & freeze_group_bit)
            meff = mi;

          sidata.r = r;
          sidata.rinv = rinv;
          sidata.meff = meff;
          sidata.mi = mi;
          sidata.mj = mj;

          if(sphere_flag) {
              sidata.en[0]   = enx_sphere;
              sidata.en[1]   = eny_sphere;
              sidata.en[2]   = enz_sphere;
          }
          sidata.omega_i = om0+3*i;
          sidata.omega_j = om0+3*j;

          // finding S-17 (opt-in), as in the serial kernel
          if (vpred_dv) {
            if (vpred_full_)
              velocity_predictor_full(sidata, vpred_dv[i], vpred_dv[j], vpred);
            else {
              velocity_predictor_normal(v0+3*i, vpred_dv[i], sidata.en, vpred.vi);
              velocity_predictor_normal(v0+3*j, vpred_dv[j], sidata.en, vpred.vj);
              sidata.v_i = vpred.vi;
              sidata.v_j = vpred.vj;
            }
          }

          cmodel.surfacesIntersect(sidata, i_forces, j_forces);

          cmodel.endSurfacesIntersect(sidata, 0, i_forces, j_forces);

          sidata.has_force_update = true;

        } else if(rsq < contactDistanceMultiplier * radsum * radsum && !shapetype_flag) {
          sidata.has_force_update = false;
          cmodel.surfacesClose(sidata, i_forces, j_forces);
        } else {
          sidata.has_force_update = false;
          if (sidata.contact_flags && *sidata.contact_flags)
            reset_separated_pair(sidata, dnum);  // finding X-05 (pair_gran_base.h)
        }

        if (mode == M_BLOCK) {
          if (sidata.has_force_update && computeflag) {
            const double relax_i = pg->relax(i);
            force_update(relax_i, f0+3*i, t0+3*i, i_forces);
            if (j < nlocal) {
              const double relax_j = pg->relax(j);
              force_update(relax_j, f0+3*j, t0+3*j, j_forces);
            }
          }
        } else if (mode == M_CROSS) {
          double * const sc = cbuf + 12*(size_t)cur_cross;
          if (sidata.has_force_update) {
            for (int coord = 0; coord < 3; coord++) {
              sc[coord]   = i_forces.delta_F[coord];
              sc[3+coord] = i_forces.delta_torque[coord];
              sc[6+coord] = j_forces.delta_F[coord];
              sc[9+coord] = j_forces.delta_torque[coord];
            }
            cflag[cur_cross] = 1;
          } else
            cflag[cur_cross] = 0;
        } else if (mode == M_SLOT) {
          const int p = base + jj;
          if (sidata.has_force_update) {
            double * const s = slot + 12*(size_t)p;
            s[0] = i_forces.delta_F[0];
            s[1] = i_forces.delta_F[1];
            s[2] = i_forces.delta_F[2];
            s[3] = i_forces.delta_torque[0];
            s[4] = i_forces.delta_torque[1];
            s[5] = i_forces.delta_torque[2];
            s[6] = j_forces.delta_F[0];
            s[7] = j_forces.delta_F[1];
            s[8] = j_forces.delta_F[2];
            s[9] = j_forces.delta_torque[0];
            s[10] = j_forces.delta_torque[1];
            s[11] = j_forces.delta_torque[2];
            if (pdelta) {
              pdelta[3*(size_t)p]   = sidata.delta[0];
              pdelta[3*(size_t)p+1] = sidata.delta[1];
              pdelta[3*(size_t)p+2] = sidata.delta[2];
            }
            pflag[p] = 1;
          } else
            pflag[p] = 0;
        } else if (sidata.has_force_update && computeflag) {
          const double relax_i = pg->relax(i);
          force_update(relax_i, &myf[6*i], &myf[6*i+3], i_forces);
          if (newton_pair || j < nlocal) {
            const double relax_j = pg->relax(j);
            force_update(relax_j, &myf[6*j], &myf[6*j+3], j_forces);
          }
        }
      }
    };

    // deterministic results do not depend on the schedule: dynamic chunks
    // only on request ("package omp N chunk M"); default static blocks
    // (contiguous rows, better cache locality)
    if (blocked) {
      // blocked deterministic mode. Thread tid owns rows [rowstart[tid],
      // rowstart[tid+1]) and the atoms of these rows. The contributions to an
      // atom are then added in the serial order: first the j-side of cross
      // pairs from earlier rows, then the pairs of the own rows (cross pairs
      // there contribute their buffered i-side), then the j-side of cross
      // pairs from later rows.
      if (omp_get_num_threads() != nthr)
        error->one(FLERR, "pair gran (OpenMP): got fewer threads than requested (OMP_DYNAMIC?)");
      // phase 1: evaluate the cross pairs of the own rows
      for (int c = bk.cstart[tid]; c < bk.cstart[tid+1]; c++) {
        cur_cross = c;
        row(bk.c_ii[c], bk.c_jj[c], bk.c_jj[c]+1, M_CROSS);
      }
      #pragma omp barrier
      auto add_jside = [&](const int g)
      {
        if (!computeflag) return;
        for (int e = bk.in_start[g]; e < bk.in_start[g+1]; e++) {
          const int c = bk.in_list[e];
          if (!cflag[c]) continue;
          const int j = firstneigh[ilist[bk.c_ii[c]]][bk.c_jj[c]] & NEIGHMASK;
          const double relax_j = pg->relax(j);
          const double * const sc = cbuf + 12*(size_t)c + 6;
          for (int coord = 0; coord < 3; coord++) {
            f[j][coord] += relax_j*sc[coord];
            torque[j][coord] += relax_j*sc[3+coord];
          }
        }
      };
      // phase 2
      add_jside(2*tid);
      for (int ii = bk.rowstart[tid]; ii < bk.rowstart[tid+1]; ii++)
        row(ii, 0, -1, M_BLOCK);
      add_jside(2*tid+1);
    } else if (det && chunk > 0) {
      #pragma omp for schedule(dynamic,chunk)
      for (int ii = 0; ii < inum; ii++)
        row(ii, 0, -1, M_SLOT);
    } else {
      #pragma omp for schedule(static)
      for (int ii = 0; ii < inum; ii++)
        row(ii, 0, -1, det ? M_SLOT : M_THREAD);
    }
    // implicit barrier: all contributions are stored

    if (computeflag && slotmode) {
      // add the contributions per atom in the order of the serial loop
      const int * const tstart = &tr.start[0];
      const int * const tentry = tr.entry.empty() ? NULL : &tr.entry[0];
      #pragma omp for schedule(static)
      for (int k = 0; k < tr.natom; k++) {
        const int e0 = tstart[k];
        const int e1 = tstart[k+1];
        if (e0 == e1) continue;
        const double relax_k = pg->relax(k);
        double * const fk = f[k];
        double * const tk = torque[k];
        for (int e = e0; e < e1; e++) {
          const int en = tentry[e];
          if (en < 0) {
            const int ii = -en - 1;
            const int p1 = offset[ii+1];
            for (int p = offset[ii]; p < p1; p++) {
              if (!pflag[p]) continue;
              const double * const s = slot + 12*(size_t)p;
              for (int coord = 0; coord < 3; coord++) {
                fk[coord] += relax_k*s[coord];
                tk[coord] += relax_k*s[3+coord];
              }
            }
          } else if (pflag[en]) {
            const double * const s = slot + 12*(size_t)en + 6;
            for (int coord = 0; coord < 3; coord++) {
              fk[coord] += relax_k*s[coord];
              tk[coord] += relax_k*s[3+coord];
            }
          }
        }
      }
    } else if (computeflag && !det) {
      // reduce the per-thread arrays in thread order
      #pragma omp for schedule(static)
      for (int k = 0; k < natom; k++) {
        for (int t = 0; t < nthr; t++) {
          const double * const b = fbuf + ((size_t)t*natom + k)*6;
          f[k][0] += b[0];
          f[k][1] += b[1];
          f[k][2] += b[2];
          torque[k][0] += b[3];
          torque[k][1] += b[4];
          torque[k][2] += b[5];
        }
      }
    }
  }

  // energy/virial tally in the order of the serial loop
  if (evflag) {
    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      const int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];
      for (int jj = 0; jj < jnum; jj++) {
        const int p = offset[ii] + jj;
        if (!pflag[p]) continue;
        const int j = jlist[jj] & NEIGHMASK;
        const double * const s = slot + 12*(size_t)p;
        const double * const d = pdelta + 3*(size_t)p;
        pg->ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, s[0], s[1], s[2], d[0], d[1], d[2]);
      }
    }
  }

  cmodel.endPass(*st.scratch[0].sidata, *st.scratch[0].i_forces, *st.scratch[0].j_forces);
  return true;
}

/* ----------------------------------------------------------------------
   explicit instantiation for every contact model of the pair factory
   (the same list as granular_styles.h)
------------------------------------------------------------------------- */

#define GRAN_MODEL(MODEL,TANGENTIAL,COHESION,ROLLING,SURFACE) \
  template bool Granular<ContactModel<GranStyle<MODEL, TANGENTIAL, COHESION, ROLLING, SURFACE> > >:: \
    compute_force_thr(PairGran *, int, int, int);
#include "style_contact_model.h"
#undef GRAN_MODEL

#ifndef LIGGGHTS_NO_CONTACT_MODEL_FALLBACK
template bool Granular<ContactModel<GranStyle<NORMAL_OFF, TANGENTIAL_OFF, COHESION_OFF, ROLLING_OFF, SURFACE_DEFAULT> > >::
  compute_force_thr(PairGran *, int, int, int);
#endif

}

}

#endif
