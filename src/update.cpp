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
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <cstddef>
#include <exception>
#include "update.h"
#include "atom.h"
#include "comm.h"
#include "integrate.h"
#include "min.h"
#include "style_integrate.h"
#include "style_minimize.h"
#include "neighbor.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "region.h"
#include "compute.h"
#include "output.h"
#include "memory.h"
#include "error.h"

#ifdef LIGGGHTS_GPU_DEM_RUNTIME
#include "GPU_DEM/gpu_runtime_liggghts.h"
#endif

using namespace LAMMPS_NS;

namespace {

int contains_text(const char *text, const char *needle)
{
  return text && needle && strstr(text,needle) != NULL;
}

const char *safe_text(const char *text)
{
  return text ? text : "none";
}

const char *safe_string_text(const std::string &text)
{
  return text.empty() ? "none" : text.c_str();
}

int model_is_enabled(const std::string &text)
{
  return !text.empty() && text != "none" && text != "off";
}

int fix_style_contains(Modify *modify, const char *needle)
{
  if (!modify || !needle) return 0;
  for (int i = 0; i < modify->nfix; ++i) {
    if (modify->fix[i] && contains_text(modify->fix[i]->style,needle))
      return 1;
  }
  return 0;
}

int atom_style_is_gpu_dem_sphere(const char *style)
{
  return style && (strcmp(style,"sphere") == 0 || strcmp(style,"granular") == 0);
}

void write_gpu_dem_line(FILE *screen, FILE *logfile, const char *line)
{
  if (screen) fprintf(screen,"%s\n",line);
  if (logfile) fprintf(logfile,"%s\n",line);
}

void write_gpu_dem_blocker(FILE *screen, FILE *logfile, const char *reason,
                           int &nblockers)
{
  if (screen) fprintf(screen,"GPU_DEM unsupported: %s\n",reason);
  if (logfile) fprintf(logfile,"GPU_DEM unsupported: %s\n",reason);
  ++nblockers;
}

void write_gpu_dem_host_fallback(FILE *screen, FILE *logfile, const char *reason,
                                 int &nfallbacks)
{
  if (screen) fprintf(screen,"GPU_DEM host fallback: %s\n",reason);
  if (logfile) fprintf(logfile,"GPU_DEM host fallback: %s\n",reason);
  ++nfallbacks;
}

int verify_gpu_dem_runtime(int requested_device, int precision_policy,
                           FILE *screen, FILE *logfile, char *errmsg,
                           std::size_t errmsg_size)
{
#ifdef LIGGGHTS_GPU_DEM_RUNTIME
  try {
    GPU_DEM::LiggghtsGpuRuntimeProbe probe;
    GPU_DEM::probe_liggghts_gpu_runtime(requested_device,precision_policy,probe);
    if (screen)
      fprintf(screen,"GPU_DEM CUDA runtime verified: device=%d/%d name=\"%s\"\n",
              probe.device_id,probe.device_count,probe.device_name);
    if (logfile)
      fprintf(logfile,"GPU_DEM CUDA runtime verified: device=%d/%d name=\"%s\"\n",
              probe.device_id,probe.device_count,probe.device_name);
    return 1;
  } catch (const std::exception &ex) {
    snprintf(errmsg,errmsg_size,"GPU_DEM CUDA runtime initialization failed: %s",ex.what());
    return 0;
  }
#else
  snprintf(errmsg,errmsg_size,
           "GPU_DEM CUDA runtime is not linked into this executable; rebuild the GPU target");
  return 0;
#endif
}

} // namespace

/* ---------------------------------------------------------------------- */

Update::Update(LAMMPS *lmp) : Pointers(lmp)
{
  char *str;

  ntimestep = 0;
  ntimestep_reset_since_last_run = false;
  timestep_set = false;
  atime = 0.0;
  atimestep = 0;
  first_update = 0;

  whichflag = 0;
  firststep = laststep = 0;
  beginstep = endstep = 0;
  setupflag = 0;
  multireplica = 0;

  restrict_output = 0;

  eflag_global = vflag_global = -1;

  unit_style = NULL;
  set_units("lj");

  integrate_style = NULL;
  integrate = NULL;
  minimize_style = NULL;
  minimize = NULL;
  gpu_dem_mode = 0;
  gpu_dem_mode_style = NULL;
  gpu_dem_precision = 1;
  gpu_dem_precision_style = NULL;
  gpu_dem_neighbor = 0;
  gpu_dem_neighbor_style = NULL;
  gpu_dem_device = -1;
  gpu_dem_device_style = NULL;
  gpu_dem_streams = 1;
  gpu_dem_streams_style = NULL;
  set_gpu_dem_mode(0,"off");
  set_gpu_dem_precision(1,"mixed");
  set_gpu_dem_neighbor(0,"auto");
  set_gpu_dem_device(-1,"auto");
  set_gpu_dem_streams(1,"on");

  if (lmp->cuda) {
    str = (char *) "verlet/cuda";
    create_integrate(1,&str,NULL);
  } else {
    str = (char *) "verlet";
    create_integrate(1,&str,NULL);
  }

  str = (char *) "cg";
  create_minimize(1,&str);

  force_dt_reset_ = false;
}

/* ---------------------------------------------------------------------- */

Update::~Update()
{
  delete [] unit_style;

  delete [] integrate_style;
  delete integrate;

  delete [] minimize_style;
  delete minimize;

  delete [] gpu_dem_mode_style;
  delete [] gpu_dem_precision_style;
  delete [] gpu_dem_neighbor_style;
  delete [] gpu_dem_device_style;
  delete [] gpu_dem_streams_style;
}

/* ---------------------------------------------------------------------- */

void Update::init()
{
  // if USER-CUDA mode is enabled:
  // integrate/minimize style must be CUDA variant

  if (whichflag == 1 && lmp->cuda)
    if (strstr(integrate_style,"cuda") == NULL)
      error->all(FLERR,"USER-CUDA mode requires CUDA variant of run style");
  if (whichflag == 2 && lmp->cuda)
    if (strstr(minimize_style,"cuda") == NULL)
      error->all(FLERR,"USER-CUDA mode requires CUDA variant of min style");

  // init the appropriate integrate and/or minimize class
  // if neither (e.g. from write_restart) then just return

  if (whichflag == 0) return;
  if (whichflag == 1) integrate->init();
  else if (whichflag == 2) minimize->init();

  // only set first_update if a run or minimize is being performed

  first_update = 1;

  ntimestep_reset_since_last_run = false;
}

/* ---------------------------------------------------------------------- */

void Update::set_gpu_dem_mode(int mode, const char *style)
{
  gpu_dem_mode = mode;
  delete [] gpu_dem_mode_style;
  gpu_dem_mode_style = new char[strlen(style) + 1];
  strcpy(gpu_dem_mode_style,style);
}

/* ---------------------------------------------------------------------- */

void Update::set_gpu_dem_precision(int precision, const char *style)
{
  gpu_dem_precision = precision;
  delete [] gpu_dem_precision_style;
  gpu_dem_precision_style = new char[strlen(style) + 1];
  strcpy(gpu_dem_precision_style,style);
}

/* ---------------------------------------------------------------------- */

void Update::set_gpu_dem_neighbor(int neighbor_policy, const char *style)
{
  gpu_dem_neighbor = neighbor_policy;
  delete [] gpu_dem_neighbor_style;
  gpu_dem_neighbor_style = new char[strlen(style) + 1];
  strcpy(gpu_dem_neighbor_style,style);
}

/* ---------------------------------------------------------------------- */

void Update::set_gpu_dem_device(int device, const char *style)
{
  gpu_dem_device = device;
  delete [] gpu_dem_device_style;
  gpu_dem_device_style = new char[strlen(style) + 1];
  strcpy(gpu_dem_device_style,style);
}

/* ---------------------------------------------------------------------- */

void Update::set_gpu_dem_streams(int streams, const char *style)
{
  gpu_dem_streams = streams;
  delete [] gpu_dem_streams_style;
  gpu_dem_streams_style = new char[strlen(style) + 1];
  strcpy(gpu_dem_streams_style,style);
}

/* ---------------------------------------------------------------------- */

void Update::check_gpu_dem_run_support()
{
  if (gpu_dem_mode == 0) return;

  int mpi_size = 1;
  MPI_Comm_size(world,&mpi_size);

  if (screen) {
    fprintf(screen,
            "GPU_DEM feature probe: mode=%s precision=%s neighbor=%s device=%s streams=%s atom_style=%s pair_style=%s "
            "normal_model=%s tangential_model=%s cohesion_model=%s "
            "rolling_model=%s mpi_size=%d\n",
            safe_text(gpu_dem_mode_style), safe_text(gpu_dem_precision_style),
            safe_text(gpu_dem_neighbor_style), safe_text(gpu_dem_device_style),
            safe_text(gpu_dem_streams_style),
            safe_text(atom->atom_style),
            safe_text(force->pair_style),
            safe_string_text(force->custom_contact_models.custom_normal_model),
            safe_string_text(force->custom_contact_models.custom_tangential_model),
            safe_string_text(force->custom_contact_models.custom_cohesion_model),
            safe_string_text(force->custom_contact_models.custom_rolling_model),
            mpi_size);
  }
  if (logfile) {
    fprintf(logfile,
            "GPU_DEM feature probe: mode=%s precision=%s neighbor=%s device=%s streams=%s atom_style=%s pair_style=%s "
            "normal_model=%s tangential_model=%s cohesion_model=%s "
            "rolling_model=%s mpi_size=%d\n",
            safe_text(gpu_dem_mode_style), safe_text(gpu_dem_precision_style),
            safe_text(gpu_dem_neighbor_style), safe_text(gpu_dem_device_style),
            safe_text(gpu_dem_streams_style),
            safe_text(atom->atom_style),
            safe_text(force->pair_style),
            safe_string_text(force->custom_contact_models.custom_normal_model),
            safe_string_text(force->custom_contact_models.custom_tangential_model),
            safe_string_text(force->custom_contact_models.custom_cohesion_model),
            safe_string_text(force->custom_contact_models.custom_rolling_model),
            mpi_size);
  }

  int nblockers = 0;
  int nfallbacks = 0;

  if (!atom_style_is_gpu_dem_sphere(atom->atom_style))
    write_gpu_dem_blocker(screen,logfile,
                          "only atom_style sphere/granular is supported",
                          nblockers);

  if (!fix_style_contains(modify,"nve/sphere"))
    write_gpu_dem_host_fallback(screen,logfile,
                                "non-nve/sphere integration remains on the host",
                                nfallbacks);

  if (mpi_size != 1)
    write_gpu_dem_blocker(screen,logfile,"multi-rank MPI GPU execution is not implemented",nblockers);

  if (!force->pair_style)
    write_gpu_dem_blocker(screen,logfile,"no pair_style is defined",nblockers);
  else if (!contains_text(force->pair_style,"gran") &&
           !contains_text(force->pair_style,"hooke"))
    write_gpu_dem_host_fallback(screen,logfile,
                                "non-granular pair style remains on the host",
                                nfallbacks);

  if (model_is_enabled(force->custom_contact_models.custom_normal_model) &&
      force->custom_contact_models.custom_normal_model.find("hooke") == std::string::npos &&
      force->custom_contact_models.custom_normal_model.find("hertz") == std::string::npos)
    write_gpu_dem_host_fallback(screen,logfile,
                                "unsupported normal contact remains on the host",
                                nfallbacks);

  if (model_is_enabled(force->custom_contact_models.custom_cohesion_model))
    write_gpu_dem_host_fallback(screen,logfile,"cohesion remains on the host",nfallbacks);

  if (fix_style_contains(modify,"multisphere"))
    write_gpu_dem_host_fallback(screen,logfile,"multisphere support remains on the host",nfallbacks);

  if (fix_style_contains(modify,"heat"))
    write_gpu_dem_host_fallback(screen,logfile,"thermal contact remains on the host",nfallbacks);

  if (nblockers == 0) {
    if (nfallbacks == 0) {
      char runtime_error[512];
      if (!verify_gpu_dem_runtime(gpu_dem_device,gpu_dem_precision,screen,logfile,
                                  runtime_error,sizeof(runtime_error)))
        error->all(FLERR,runtime_error);
      write_gpu_dem_line(screen,logfile,"GPU_DEM feature probe accepted this run for GPU-native execution");
      if (gpu_dem_mode == 2)
        error->all(FLERR,
                   "gpu_mode strict verified CUDA but production Verlet GPU timestep is not wired; refusing CPU timestep fallback");
      return;
    }

    if (gpu_dem_mode == 2)
      error->all(FLERR,
                 "gpu_mode strict rejected this run because host fallbacks would be required; see GPU_DEM host fallback messages above");

    write_gpu_dem_line(screen,logfile,
                       "GPU_DEM auto mode: falling back to CPU for this run");
    return;
  }

  if (gpu_dem_mode == 2)
    error->all(FLERR,"gpu_mode strict rejected this run; see GPU_DEM unsupported messages above");

  write_gpu_dem_line(screen,logfile,
                     "GPU_DEM auto mode: falling back to CPU for this run");
}

/* ---------------------------------------------------------------------- */

void Update::set_units(const char *style)
{
  // physical constants from:
  // http://physics.nist.gov/cuu/Constants/Table/allascii.txt
  // using thermochemical calorie = 4.184 J

  if (strcmp(style,"lj") == 0) {
    force->boltz = 1.0;
    force->hplanck = 0.18292026;  // using LJ parameters for argon
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 1.0;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0;
    force->femtosecond = 1.0;
    force->qelectron = 1.0;

    dt = 0.005;
    neighbor->skin = 0.3;

  } else if (strcmp(style,"real") == 0) {
    force->boltz = 0.0019872067;
    force->hplanck = 95.306976368;
    force->mvv2e = 48.88821291 * 48.88821291;
    force->ftm2v = 1.0 / 48.88821291 / 48.88821291;
    force->mv2d = 1.0 / 0.602214179;
    force->nktv2p = 68568.415;
    force->qqr2e = 332.06371;
    force->qe2f = 23.060549;
    force->vxmu2f = 1.4393264316e4;
    force->xxt2kmu = 0.1;
    force->e_mass = 1.0/1836.1527556560675;
    force->hhmrr2e = 0.0957018663603261;
    force->mvh2r = 1.5339009481951;
    force->angstrom = 1.0;
    force->femtosecond = 1.0;
    force->qelectron = 1.0;

    dt = 1.0;
    neighbor->skin = 2.0;

  } else if (strcmp(style,"metal") == 0) {
    force->boltz = 8.617343e-5;
    force->hplanck = 4.135667403e-3;
    force->mvv2e = 1.0364269e-4;
    force->ftm2v = 1.0 / 1.0364269e-4;
    force->mv2d = 1.0 / 0.602214179;
    force->nktv2p = 1.6021765e6;
    force->qqr2e = 14.399645;
    force->qe2f = 1.0;
    force->vxmu2f = 0.6241509647;
    force->xxt2kmu = 1.0e-4;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0;
    force->femtosecond = 1.0e-3;
    force->qelectron = 1.0;

    dt = 0.001;
    neighbor->skin = 2.0;

  } else if (strcmp(style,"si") == 0) {
    force->boltz = 1.3806504e-23;
    force->hplanck = 6.62606896e-34;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 8.9876e9;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0e-10;
    force->femtosecond = 1.0e-15;
    force->qelectron = 1.6021765e-19;

    dt = 1.0e-8;
    neighbor->skin = 0.001;

  } else if (strcmp(style,"cgs") == 0) {
    force->boltz = 1.3806504e-16;
    force->hplanck = 6.62606896e-27;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 1.0;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0e-8;
    force->femtosecond = 1.0e-15;
    force->qelectron = 4.8032044e-10;

    dt = 1.0e-8;
    neighbor->skin = 0.1;

  } else if (strcmp(style,"electron") == 0) {
    force->boltz = 3.16681534e-6;
    force->hplanck = 0.1519829846;
    force->mvv2e = 1.06657236;
    force->ftm2v = 0.937582899;
    force->mv2d = 1.0;
    force->nktv2p = 2.94210108e13;
    force->qqr2e = 1.0;
    force->qe2f = 1.94469051e-10;
    force->vxmu2f = 3.39893149e1;
    force->xxt2kmu = 3.13796367e-2;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.88972612;
    force->femtosecond = 0.0241888428;
    force->qelectron = 1.0;

    dt = 0.001;
    neighbor->skin = 2.0;

  } else if (strcmp(style,"micro") == 0) {
    force->boltz = 1.3806504e-8;
    force->hplanck = 6.62606896e-13;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 8.9876e30;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0e-4;
    force->femtosecond = 1.0e-9;
    force->qelectron = 1.6021765e-19;

    dt = 2.0;
    neighbor->skin = 0.1;

  } else if (strcmp(style,"nano") == 0) {
    force->boltz = 0.013806503;
    force->hplanck = 6.62606896e-4;
    force->mvv2e = 1.0;
    force->ftm2v = 1.0;
    force->mv2d = 1.0;
    force->nktv2p = 1.0;
    force->qqr2e = 8.9876e39;
    force->qe2f = 1.0;
    force->vxmu2f = 1.0;
    force->xxt2kmu = 1.0;
    force->e_mass = 0.0;    // not yet set
    force->hhmrr2e = 0.0;
    force->mvh2r = 0.0;
    force->angstrom = 1.0e-1;
    force->femtosecond = 1.0e-6;
    force->qelectron = 1.6021765e-19;

    dt = 0.00045;
    neighbor->skin = 0.1;

  } else error->all(FLERR,"Illegal units command");

  delete [] unit_style;
  int n = strlen(style) + 1;
  unit_style = new char[n];
  strcpy(unit_style,style);
}

/* ---------------------------------------------------------------------- */

void Update::create_integrate(int narg, char **arg, char *suffix)
{
  if (narg < 1) error->all(FLERR,"Illegal run_style command");

  delete [] integrate_style;
  delete integrate;

  int sflag;
  new_integrate(arg[0],narg-1,&arg[1],suffix,sflag);

  if (sflag) {
    char estyle[256];
    sprintf(estyle,"%s/%s",arg[0],suffix);
    int n = strlen(estyle) + 1;
    integrate_style = new char[n];
    strcpy(integrate_style,estyle);
  } else {
    int n = strlen(arg[0]) + 1;
    integrate_style = new char[n];
    strcpy(integrate_style,arg[0]);
  }
}

/* ----------------------------------------------------------------------
   create the Integrate style, first with suffix appended
------------------------------------------------------------------------- */

void Update::new_integrate(char *style, int narg, char **arg,
                           char *suffix, int &sflag)
{
  int success = 0;

  if (suffix && lmp->suffix_enable) {
    sflag = 1;
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);
    success = 1;

    if (0) return;

#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) \
    else if (strcmp(estyle,#key) == 0) integrate = new Class(lmp,narg,arg);
#include "style_integrate.h"
#undef IntegrateStyle
#undef INTEGRATE_CLASS

    else success = 0;
  }

  if (!success) {
    sflag = 0;

    if (0) return;

#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) \
    else if (strcmp(style,#key) == 0) integrate = new Class(lmp,narg,arg);
#include "style_integrate.h"
#undef IntegrateStyle
#undef INTEGRATE_CLASS

    else error->all(FLERR,"Illegal integrate style");
  }
}

/* ---------------------------------------------------------------------- */

void Update::create_minimize(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal min_style command");

  delete [] minimize_style;
  delete minimize;

  if (0) return;      // dummy line to enable else-if macro expansion

#define MINIMIZE_CLASS
#define MinimizeStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) minimize = new Class(lmp);
#include "style_minimize.h"
#undef MINIMIZE_CLASS

  else error->all(FLERR,"Illegal min_style command");

  int n = strlen(arg[0]) + 1;
  minimize_style = new char[n];
  strcpy(minimize_style,arg[0]);
}

/* ----------------------------------------------------------------------
   reset timestep as called from input script
------------------------------------------------------------------------- */

void Update::reset_timestep(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal reset_timestep command");
  bigint newstep = ATOBIGINT(arg[0]);
  reset_timestep(newstep);
}

/* ----------------------------------------------------------------------
   reset timestep
   set atimestep to new timestep, so future update_time() calls will be correct
   trigger reset of timestep for output and for fixes that require it
   do not allow any timestep-dependent fixes to be defined
   reset eflag/vflag global so nothing will think eng/virial are current
   reset invoked flags of computes,
     so nothing will think they are current between runs
   clear timestep list of computes that store future invocation times
   called from rerun command and input script (indirectly)
------------------------------------------------------------------------- */

void Update::reset_timestep(bigint newstep)
{
  
  ntimestep_reset_since_last_run = true;
  bigint oldtimestep = ntimestep;

  ntimestep = newstep;
  if (ntimestep < 0) error->all(FLERR,"Timestep must be >= 0");
  if (ntimestep > MAXBIGINT) error->all(FLERR,"Too big a timestep");

  atime += (ntimestep - atimestep) * dt;
  if (atime < 0)
      atime = 0;
  atimestep = ntimestep;

  output->reset_timestep(ntimestep);

  for (int i = 0; i < modify->nfix; i++) {
    if (modify->fix[i]->time_depend && !force_dt_reset_) 
      error->all(FLERR,
                 "Cannot reset timestep with a time-dependent fix defined");
    modify->fix[i]->reset_timestep(ntimestep,oldtimestep);
  }

  eflag_global = vflag_global = -1;

  for (int i = 0; i < modify->ncompute; i++) {
    modify->compute[i]->invoked_scalar = -1;
    modify->compute[i]->invoked_vector = -1;
    modify->compute[i]->invoked_array = -1;
    modify->compute[i]->invoked_peratom = -1;
    modify->compute[i]->invoked_local = -1;
  }

  for (int i = 0; i < modify->ncompute; i++)
    if (modify->compute[i]->timeflag) modify->compute[i]->clearstep();

  // NOTE: 7Jun12, adding rerun command, don't think this is required

  //for (int i = 0; i < domain->nregion; i++)
  //  if (domain->regions[i]->dynamic_check())
  //    error->all(FLERR,"Cannot reset timestep with a dynamic region defined");
}

/* ----------------------------------------------------------------------
   update elapsed simulation time
   called at end of runs or when timestep size changes
------------------------------------------------------------------------- */

void Update::update_time()
{
  atime += (ntimestep-atimestep) * dt;
  atimestep = ntimestep;
}

/* ----------------------------------------------------------------------
   memory usage of update and integrate/minimize
------------------------------------------------------------------------- */

bigint Update::memory_usage()
{
  bigint bytes = 0;
  if (whichflag == 1) bytes += integrate->memory_usage();
  else if (whichflag == 2) bytes += minimize->memory_usage();
  return bytes;
}
