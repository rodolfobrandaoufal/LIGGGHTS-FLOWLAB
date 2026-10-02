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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)
    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef PAIR_GRAN_BASE_H_
#define PAIR_GRAN_BASE_H_

#include <vector>
#include "contact_interface.h"
#include "math_extra_liggghts.h"

#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include "math_extra_liggghts_nonspherical.h"
#include "math_const.h"
#endif

#include "pair_gran.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "update.h"
#include "error.h"
#include "fix_contact_property_atom.h"
#include "os_specific.h"
#include "fix_insert_stream_predefined.h"
#include "fix_property_atom.h"
#include "modify.h"
#include "velocity_predictor.h"

#include "granular_pair_style.h"

#ifdef LIGGGHTS_OMP
#include <omp.h>
#include "thr_granular.h"
#endif

namespace LIGGGHTS {
using namespace ContactModels;

namespace PairStyles {

using namespace LAMMPS_NS;

template<typename ContactModel>
class Granular : private Pointers, public IGranularPairStyle {
  SurfacesIntersectData * aligned_sidata;
  ForceData * aligned_i_forces;
  ForceData * aligned_j_forces;
  ContactModel cmodel;
#ifdef LIGGGHTS_OMP
  // OpenMP kernel state (thr_granular.h, pair_gran_omp.cpp)
  LIGGGHTS::ThrGranular::PairState * thr_state_ = NULL;
#endif

  // LIGGGHTS modernization branch, finding X-05: the contact history of a
  // pair that is beyond the contact distance is reset at once (as fix
  // wall/gran does for walls and the neighbor build does at the next
  // rebuild). 'history_clear_legacy on' keeps the record until the next
  // neighbor build (behaviour up to phase D).
  bool history_clear_legacy_;
  // multicontact: atom i uses its per-contact expanded radius, as atom j does.
  // 'multicontact_radius_legacy on' restores the unexpanded radius of i.
  bool multicontact_radius_legacy_;
  // one-time warning of the corrected default (see reset_separated_pair())
  std::vector<std::vector<int*> > reset_slots_;
  bigint reset_ncalls_;
  int history_changed_;
  bool history_clear_warned_;
  bool history_warning_pending_;
  int me_;
  bool multicontact_warned_;
  FixPropertyAtom * fix_vpred_;   // finding S-17, NULL unless velocity_predictor normal|full
  bool vpred_full_;

  inline void force_update(double relax,double *const f, double *const torque,
      const ForceData & forces)
  {
    for (int coord = 0; coord < 3; coord++)
    {
      f[coord] += relax*forces.delta_F[coord];
      torque[coord] += relax*forces.delta_torque[coord];
    }
  }

public:
  Granular(class LAMMPS * lmp, PairGran* parent, const int64_t hash) : Pointers(lmp),
    aligned_sidata(aligned_malloc<SurfacesIntersectData>(32)),
    aligned_i_forces(aligned_malloc<ForceData>(32)),
    aligned_j_forces(aligned_malloc<ForceData>(32)),
    cmodel(lmp, parent,false /*is_wall*/, hash),
    history_clear_legacy_(false),
    multicontact_radius_legacy_(false),
    reset_ncalls_(-1),
    history_changed_(0),
    history_clear_warned_(false),
    history_warning_pending_(false),
    me_(comm->me),
    multicontact_warned_(false),
    fix_vpred_(NULL),
    vpred_full_(false)
  {
  }

  virtual ~Granular() {
    // a pending X-05 warning of the last run (comm is already deleted here)
    if (history_warning_pending_ && me_ == 0)
      print_history_warning();
    aligned_free(aligned_sidata);
    aligned_free(aligned_i_forces);
    aligned_free(aligned_j_forces);
#ifdef LIGGGHTS_OMP
    LIGGGHTS::ThrGranular::pair_state_free(thr_state_);
#endif
  }

  int64_t hashcode()
  { return cmodel.hashcode(); }

  virtual void settings(int nargs, char ** args, IContactHistorySetup *hsetup) {
    Settings settings(lmp);
    cmodel.registerSettings(settings);
    settings.registerOnOff("history_clear_legacy", history_clear_legacy_, false);
    settings.registerOnOff("multicontact_radius_legacy", multicontact_radius_legacy_, false);
    bool success = settings.parseArguments(nargs, args);
    cmodel.postSettings(hsetup);

#ifdef LIGGGHTS_DEBUG
    if(comm->me == 0) {
      fprintf(screen, "==== PAIR SETTINGS ====\n");
      settings.print_all(screen);
      fprintf(screen, "==== PAIR SETTINGS ====\n");

      fprintf(logfile, "==== PAIR SETTINGS ====\n");
      settings.print_all(logfile);
      fprintf(logfile, "==== PAIR SETTINGS ====\n");
    }
#endif

    if(!success) {
      error->all(FLERR,settings.error_message.c_str());
    }
  }

  virtual void init_granular() {
    cmodel.connectToProperties(force->registry);
    // finding S-17 (velocity_predictor.h)
    fix_vpred_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(
        VELOCITY_PREDICTOR_FULL,"property/atom","vector",6,0,"pair gran",false));
    vpred_full_ = (fix_vpred_ != NULL);
    if (!fix_vpred_)
      fix_vpred_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(
          VELOCITY_PREDICTOR_NORMAL,"property/atom","vector",3,0,"pair gran",false));
    if (vpred_full_ && !(cmodel.contact_match("surface","default") &&
                         (cmodel.contact_match("tangential","history") ||
                          cmodel.contact_match("tangential","no_history"))))
      error->all(FLERR,"pair gran: 'velocity_predictor full' supports surface default with tangential "
                 "history or no_history only; use 'velocity_predictor normal'");

#ifdef LIGGGHTS_DEBUG
    if(comm->me == 0) {
      fprintf(screen, "==== PAIR GLOBAL PROPERTIES ====\n");
      force->registry.print_all(screen);
      fprintf(screen, "==== PAIR GLOBAL PROPERTIES ====\n");

      fprintf(logfile, "==== PAIR GLOBAL PROPERTIES ====\n");
      force->registry.print_all(logfile);
      fprintf(logfile, "==== PAIR GLOBAL PROPERTIES ====\n");
    }
#endif
  }

  virtual void write_restart_settings(FILE * fp)
  {
    int64_t hashcode = cmodel.hashcode();
    fwrite(&hashcode, sizeof(int64_t), 1, fp);
  }

  virtual void read_restart_settings(FILE * fp, int64_t hashcode)
  {
    int me = comm->me;
    int64_t selected = -1;
    if(me == 0){
      size_t dummy = fread(&selected, sizeof(int64_t), 1, fp);
      UNUSED(dummy);
      // sanity check
      if(hashcode != -1) { // backward compability
          if(hashcode != cmodel.hashcode())
              error->one(FLERR,"wrong pair style loaded!");
      } else {
          if(selected != cmodel.hashcode())
              error->one(FLERR,"wrong pair style loaded!");
      }
    }
  }

  inline bool contact_match(const std::string mtype, const std::string model) {
    return cmodel.contact_match(mtype, model);
  }

  int get_history_offset(const std::string hname)
  {
    return cmodel.get_history_offset(hname);
  }

  double stressStrainExponent()
  {
    return cmodel.stressStrainExponent();
  }

  /* ----------------------------------------------------------------------
     finding X-05: called for a pair beyond the contact distance (no
     surfacesIntersect, no surfacesClose) whose contact flag is set. Resets
     its contact flag and history as the next neighbor build would, so that a
     pair touching again before that build starts from zero history. Only for
     a force evaluation that updates the history (not in setup, not for
     compute pair/gran/local). Each pair is evaluated by one thread, so the
     record itself needs no synchronisation.
     For the one-time warning, a reset record that held non-zero history is
     remembered (per thread) until the next neighbor build; if the pair
     touches again before that build (its flag is set again), the corrected
     default can have changed the result (history_pass_end()). This is a
     superset: a stale value that the models overwrite before reading it
     (e.g. kt_old without a stored spring) leaves the result unchanged.
  ------------------------------------------------------------------------- */
  void reset_separated_pair(SurfacesIntersectData & sidata, const int dnum)
  {
    int * const cflag = sidata.contact_flags;
    if (!cflag || !*cflag || history_clear_legacy_ || !sidata.computeflag || !sidata.shearupdate)
      return;
    *cflag = 0;
    bool nonzero = false;
    double * const hist = sidata.contact_history;
    if (hist)
      for (int d = 0; d < dnum; d++) {
        if (hist[d] != 0.0) nonzero = true;
        hist[d] = 0.0;
      }
    if (nonzero && !history_clear_warned_) {
#ifdef LIGGGHTS_OMP
      const size_t tid = omp_in_parallel() ? (size_t)omp_get_thread_num() : 0;
#else
      const size_t tid = 0;
#endif
      if (tid < reset_slots_.size())
        reset_slots_[tid].push_back(cflag);
    }
  }

  // before a force evaluation: forget the reset records of an old neighbor
  // list (their addresses are no longer valid), one record list per thread
  void history_pass_begin(PairGran * pg)
  {
    if (history_clear_legacy_ || history_clear_warned_ || pg->dnum() == 0)
      return;
    if (reset_ncalls_ != neighbor->ncalls) {
      for (size_t t = 0; t < reset_slots_.size(); t++)
        reset_slots_[t].clear();
      reset_ncalls_ = neighbor->ncalls;
    }
#ifdef LIGGGHTS_OMP
    size_t nthr = (size_t)omp_get_max_threads();
    const size_t nthr_pkg = (size_t)LIGGGHTS::ThrGranular::config(lmp).nthreads;
    if (nthr_pkg > nthr) nthr = nthr_pkg;
#else
    const size_t nthr = 1;
#endif
    if (reset_slots_.size() < nthr)
      reset_slots_.resize(nthr);
  }

  // after a force evaluation: did a pair touch again after its history was
  // reset? The flag is reduced over all procs only at the last step of a run
  // (the same step on all procs) and only until the first warning, so that
  // the per-step cost is negligible. The warning itself is printed by rank 0
  // at the setup of the next run or at the end of the input, i.e. outside the
  // thermo output of the run.
  void history_pass_end(PairGran * pg)
  {
    if (history_clear_legacy_ || history_clear_warned_ || pg->dnum() == 0)
      return;
    for (size_t t = 0; t < reset_slots_.size() && !history_changed_; t++)
      for (size_t k = 0; k < reset_slots_[t].size(); k++)
        if (*reset_slots_[t][k]) { history_changed_ = 1; break; }
    if (update->ntimestep == update->laststep && !update->setupflag) {
      int any = 0;
      MPI_Allreduce(&history_changed_, &any, 1, MPI_INT, MPI_MAX, world);
      if (any) {
        history_clear_warned_ = true;
        history_warning_pending_ = true;
        for (size_t t = 0; t < reset_slots_.size(); t++)
          std::vector<int*>().swap(reset_slots_[t]);
      }
    }
  }

  void print_history_warning()
  {
    history_warning_pending_ = false;
    error->warning(FLERR, "pair gran: in the previous run a pair touched again within one neighbor "
                   "interval after separating; its contact history was reset at the separation "
                   "(finding X-05), up to now it was kept until the next neighbor build, so results "
                   "can differ from earlier versions. Use 'history_clear_legacy on' for the old behaviour");
  }

  // one-time warning for the multicontact radius (rank 0)
  void multicontact_warning(PairGran * pg)
  {
    if (!multicontact_radius_legacy_ && !multicontact_warned_ && pg->storeSumDelta()) {
      multicontact_warned_ = true;
      if (comm->me == 0)
        error->warning(FLERR, "pair gran: surface multicontact now uses the expanded radius of both particles "
                       "of a contact (the radius of atom i was unexpanded up to now, so results depended on "
                       "the pair orientation). Use 'multicontact_radius_legacy on' for the old behaviour");
    }
  }

  /* ----------------------------------------------------------------------
     finding S-17: the stored increments dv = dt/2 a(n) of 'fix nve/sphere
     ... velocity_predictor yes', or NULL. Only during the force evaluation
     of a time step: in setup and for compute pair/gran/local (addflag) the
     velocities are already synchronised.
  ------------------------------------------------------------------------- */
  inline double ** velocity_predictor_dv(const int addflag) const
  {
    if (!fix_vpred_ || addflag || update->setupflag) return NULL;
    return fix_vpred_->array_atom;
  }

  virtual void compute_force(PairGran * pg, int eflag, int vflag, int addflag)
  {
    if (history_warning_pending_ && update->setupflag) {
      if (comm->me == 0) print_history_warning();
      history_warning_pending_ = false;
    }
    multicontact_warning(pg);
    history_pass_begin(pg);
#ifdef LIGGGHTS_OMP
    if (!compute_force_thr(pg, eflag, vflag, addflag))
#endif
      compute_force_serial(pg, eflag, vflag, addflag);
    history_pass_end(pg);
  }

#ifdef LIGGGHTS_OMP
  // OpenMP kernel (roadmap C2/B8). Defined in pair_gran_omp.cpp and
  // explicitly instantiated there, so that this translation unit (which
  // instantiates all contact models) and hence the serial kernel's code
  // generation stay unchanged. Returns false if the serial kernel must run.
  bool compute_force_thr(PairGran * pg, int eflag, int vflag, int addflag);
#endif

  // the serial kernel
  void compute_force_serial(PairGran * pg, int eflag, int vflag, int addflag)
  {
    // finding S-17: a separate instantiation with the velocity predictor, so
    // that the code generated for the default kernel (and hence its floating
    // point results, e.g. FMA contraction) is unchanged
    if (!velocity_predictor_dv(addflag))
      compute_force_serial_t<0>(pg, eflag, vflag, addflag);
    else if (vpred_full_)
      compute_force_serial_t<2>(pg, eflag, vflag, addflag);
    else
      compute_force_serial_t<1>(pg, eflag, vflag, addflag);
  }

  // VPRED: 0 no predictor, 1 'normal', 2 'full'
  template<int VPRED>
  void compute_force_serial_t(PairGran * pg, int eflag, int vflag, int addflag)
  {
    if (eflag || vflag)
      pg->ev_setup(eflag, vflag);
    else
      pg->evflag = pg->vflag_fdotr = 0;

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double **omega = atom->omega;
    double **torque = atom->torque;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    int *mask = atom->mask;
    int *tag = atom->tag;
    int nlocal = atom->nlocal;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    int superquadric_flag = atom->superquadric_flag;
#endif // SUPERQUADRIC_ACTIVE_FLAG
    const int newton_pair = force->newton_pair;
    int inum = pg->list->inum;
    int * ilist = pg->list->ilist;
    int * numneigh = pg->list->numneigh;

    int ** firstneigh = pg->list->firstneigh;
    int ** first_contact_flag = pg->listgranhistory ? pg->listgranhistory->firstneigh : NULL;
    double ** first_contact_hist = pg->listgranhistory ? pg->listgranhistory->firstdouble : NULL;

    const int dnum = pg->dnum();
    const bool store_contact_forces = pg->storeContactForces();
    const bool store_contact_forces_stress = pg->storeContactForcesStress();
    const int freeze_group_bit = pg->freeze_group_bit();

    const double contactDistanceMultiplier = neighbor->contactDistanceFactor*neighbor->contactDistanceFactor;

    // finding S-17: predicted normal velocities (velocity_predictor.h)
    double ** const vpred_dv = VPRED ? velocity_predictor_dv(addflag) : NULL;
    VelocityPredictorScratch vpred;

    // fix insert/stream/predefined
    // check if inserted
    // use most_recent_ins_step of this fix for this
    std::vector<FixInsertStreamPredefined*> fix_insert;
    int nfix_insert = modify->n_fixes_style("insert/stream/predefined");
    for (int i = 0; i < nfix_insert; i++)
    {
        FixInsertStreamPredefined * fix = static_cast<FixInsertStreamPredefined*>(modify->find_fix_style("insert/stream/predefined", i));
        if (fix->has_inserted())
            fix_insert.push_back(fix);
    }

    // clear data, just to be safe
    memset(aligned_sidata, 0, sizeof(SurfacesIntersectData));
    memset(aligned_i_forces, 0, sizeof(ForceData));
    memset(aligned_j_forces, 0, sizeof(ForceData));
    aligned_sidata->area_ratio = 1.0;

    SurfacesIntersectData & sidata = *aligned_sidata;
    ForceData & i_forces = *aligned_i_forces;
    ForceData & j_forces = *aligned_j_forces;
    sidata.is_wall = false;
    sidata.computeflag = pg->computeflag();
    sidata.shearupdate = pg->shearupdate();

    cmodel.beginPass(sidata, i_forces, j_forces);

    // loop over neighbors of my atoms

    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      double radi = radius[i];
      int * const contact_flags = first_contact_flag ? first_contact_flag[i] : NULL;
      double * const all_contact_hist = first_contact_hist ? first_contact_hist[i] : NULL;
      int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      sidata.i = i;
      #ifdef SUPERQUADRIC_ACTIVE_FLAG
          if(superquadric_flag) {
            sidata.radi = cbrt(0.75 * atom->volume[i] / M_PI);
          } else
            sidata.radi = radi;
      #else
          sidata.radi = radi;
      #endif

      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;

        const double xj = x[j][0];
        const double yj = x[j][1];
        const double zj = x[j][2];
        const double delx = xtmp - xj;
        const double dely = ytmp - yj;
        const double delz = ztmp - zj;
        const double rsq = delx * delx + dely * dely + delz * delz;
        double radj = radius[j];

        // In case of multicontact models use the computed delta_ij and delta_ji to expand the radius (on a per contact basis)
        if (pg->storeSumDelta()) {
            FixContactPropertyAtom* mcFix = pg->fix_store_multicontact_delta();
            const int cj = mcFix->has_partner(i, tag[j]);
            radi = radius[i];
            if (cj != -1)
            {
                const double * const dataI = mcFix->contacthistory(i, cj);
                radi += dataI[3];
            }
            // the expanded radius of i, as for j (sidata.radj) below
            if (!multicontact_radius_legacy_)
                sidata.radi = radi;
            const int ci = mcFix->has_partner(j, tag[i]);
            if (ci != -1)
            {
                const double * const dataJ = mcFix->contacthistory(j, ci);
                radj += dataJ[3];
            }
        }

#ifdef SUPERQUADRIC_ACTIVE_FLAG
        if(superquadric_flag) {
          sidata.radj = cbrt(0.75 * atom->volume[j] / M_PI);
        } else
          sidata.radj = radj;
#else
        sidata.radj = radj;
#endif
        const double radsum = radi + radj;

        sidata.j = j;
        sidata.delta[0] = delx;
        sidata.delta[1] = dely;
        sidata.delta[2] = delz;
        sidata.rsq = rsq;
        sidata.radsum = radsum;
        sidata.contact_flags = contact_flags ? &contact_flags[jj] : NULL;
        sidata.contact_history = all_contact_hist ? &all_contact_hist[dnum*jj] : NULL;

        if (!fix_insert.empty())
        {
            std::vector<FixInsertStreamPredefined*>::iterator it = fix_insert.begin();
            for (; it != fix_insert.end(); it++)
            {
                (*it)->copy_history(i, j, sidata.contact_history);
            }
        }

        i_forces.reset();
        j_forces.reset();

        // rsq < radsum * radsum is broad phase check with bounding spheres
        // cmodel.checkSurfaceIntersect() is narrow phase check
        
        #ifdef SUPERQUADRIC_ACTIVE_FLAG
        if (rmass) {
          sidata.mi = rmass[i];
          sidata.mj = rmass[j];
        } else {
          sidata.mi = mass[type[i]];
          sidata.mj = mass[type[j]];
        }
        sidata.omega_i = omega[i];
        sidata.omega_j = omega[j];
        #endif

        sidata.v_i     = v[i];
        sidata.v_j     = v[j];
        const int itype = type[i];
        const int jtype = type[j];
        sidata.itype = itype;
        sidata.jtype = jtype;

        if (rsq < radsum * radsum && cmodel.checkSurfaceIntersect(sidata)) {
          const double r = sqrt(rsq);
          const double rinv = 1.0 / r;

          // unit normal vector for case of spherical particles
          // for non-spherical, this is done by surface model
          const double enx_sphere = delx * rinv;
          const double eny_sphere = dely * rinv;
          const double enz_sphere = delz * rinv;

          // meff = effective mass of pair of particles
          // if I or J part of rigid body, use body mass
          // if I or J is frozen, meff is other particle
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

          // copy collision data to struct (compiler can figure out a better way to
          // interleave these stores with the double calculations above.
          sidata.r = r;
          sidata.rinv = rinv;
          sidata.meff = meff;
          sidata.mi = mi;
          sidata.mj = mj;
          
          if(atom->sphere_flag) {
              sidata.en[0]   = enx_sphere;
              sidata.en[1]   = eny_sphere;
              sidata.en[2]   = enz_sphere;
          }
          sidata.omega_i = omega[i];
          sidata.omega_j = omega[j];

          if (VPRED == 2)
            velocity_predictor_full(sidata, vpred_dv[i], vpred_dv[j], vpred);
          else if (VPRED == 1) {
            velocity_predictor_normal(v[i], vpred_dv[i], sidata.en, vpred.vi);
            velocity_predictor_normal(v[j], vpred_dv[j], sidata.en, vpred.vj);
            sidata.v_i = vpred.vi;
            sidata.v_j = vpred.vj;
          }

          cmodel.surfacesIntersect(sidata, i_forces, j_forces);

          cmodel.endSurfacesIntersect(sidata, 0, i_forces, j_forces);

          // if there is a surface touch, there will always be a force
          sidata.has_force_update = true;

        // surfacesClose is not supported for convex particles
        } else if(rsq < contactDistanceMultiplier * radsum * radsum && !atom->shapetype_flag) {
          // apply force update only if selected contact models have requested it
          sidata.has_force_update = false;
          cmodel.surfacesClose(sidata, i_forces, j_forces);
        } else {
          sidata.has_force_update = false;
          if (sidata.contact_flags && *sidata.contact_flags)
            reset_separated_pair(sidata, dnum);   // finding X-05
        }

        if(sidata.has_force_update) {
          if (sidata.computeflag) {

            const double relax_i = pg->relax(i);
            force_update(relax_i,f[i],torque[i],i_forces);

            if(newton_pair || j < nlocal) {
              const double relax_j = pg->relax(j);
              force_update(relax_j,f[j],torque[j],j_forces);
            }

            // summation of f.n to compute a simplistic pressure
            if (pg->store_sum_normal_force())
            {
                double * const iforce = pg->get_sum_normal_force_ptr(i);
                const double fDotN = vectorDot3D(i_forces.delta_F, sidata.en);
                *iforce += fDotN;
                if (j < nlocal || newton_pair)
                {
                    double * const jforce = pg->get_sum_normal_force_ptr(j);
                    *jforce += fDotN;
                }
            }
          }

          if (pg->cpl() && addflag)
            pg->cpl_add_pair(sidata, i_forces);

          if (pg->evflag)
            pg->ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0,i_forces.delta_F[0],i_forces.delta_F[1],i_forces.delta_F[2],sidata.delta[0],sidata.delta[1],sidata.delta[2]);

          if (store_contact_forces && 0 == update->ntimestep % pg->storeContactForcesEvery())
          {
            double forces_torques_i[6],forces_torques_j[6];

            if(pg->fix_contact_forces()->has_partner(i,atom->tag[j]) == -1)
            {
                vectorCopy3D(i_forces.delta_F,&(forces_torques_i[0]));
                vectorCopy3D(i_forces.delta_torque,&(forces_torques_i[3]));
                pg->fix_contact_forces()->add_partner(i,atom->tag[j],forces_torques_i);
            }
            if(pg->fix_contact_forces()->has_partner(j,atom->tag[i]) == -1)
            {
                vectorCopy3D(j_forces.delta_F,&(forces_torques_j[0]));
                vectorCopy3D(j_forces.delta_torque,&(forces_torques_j[3]));
                pg->fix_contact_forces()->add_partner(j,atom->tag[i],forces_torques_j);
            }
          }

          if (store_contact_forces_stress)
          {
            double forces_pos[4];

            if(pg->fix_contact_forces_stress()->has_partner(i,atom->tag[j]) == -1)
            {
                vectorCopy3D(i_forces.delta_F,&(forces_pos[0]));
                forces_pos[3] = (double) j;
                pg->fix_contact_forces_stress()->add_partner(i,atom->tag[j],forces_pos);
            }
            if(pg->fix_contact_forces_stress()->has_partner(j,atom->tag[i]) == -1)
            {
                vectorCopy3D(j_forces.delta_F,&(forces_pos[0]));
                forces_pos[3] = (double) i;
                pg->fix_contact_forces_stress()->add_partner(j,atom->tag[i],forces_pos);
            }
          }
        }
      }
    }

    cmodel.endPass(sidata, i_forces, j_forces);

    if (pg->cpl() && addflag)
        pg->cpl_pair_finalize();

    if(store_contact_forces)
        pg->fix_contact_forces()->do_forward_comm();
    if(store_contact_forces_stress)
        pg->fix_contact_forces_stress()->do_forward_comm();
  }
};

}

}
#endif /* PAIR_GRAN_BASE_H_ */
