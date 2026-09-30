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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include "pointers.h"
#include "lammps.h"
#include "fix_property_global.h"
#include <map>
#include <set>
#include <string>
#include "properties.h"
#include "error.h"
#include "modify.h"
#include "property_registry.h"
#include "update.h"
#include <cstring>

using namespace std;
using namespace LAMMPS_NS;

PropertyRegistry::PropertyRegistry(LAMMPS* lmp) : Pointers(lmp), properties(lmp),
  refreshing_(false), n_refreshed_(0), cached_max_type_(0)
{
}

PropertyRegistry::~PropertyRegistry()
{
  init();
}

int PropertyRegistry::max_type()
{
  // During refresh() the property sizes are fixed (the refresh copies into
  // the storage built at init), so reuse the value from the last build and
  // avoid the O(nlocal) scan + 2 MPI_Allreduce per creator call.
  if(refreshing_ && cached_max_type_ > 0) return cached_max_type_;
  cached_max_type_ = properties.max_type();
  return cached_max_type_;
}

double PropertyRegistry::min_radius()
{
  return properties.min_radius();
}

double PropertyRegistry::max_radius()
{
  return properties.max_radius();
}

LAMMPS * PropertyRegistry::getLAMMPS()
{
  return lmp;
}

FixPropertyGlobal* PropertyRegistry::getGlobalProperty(const char *varname, const char *style, const char *svmstyle, int len1, int len2, const char *caller)
{
  FixPropertyGlobal *fix = static_cast<FixPropertyGlobal*>(modify->find_fix_property(varname, style, svmstyle, len1, len2, caller));

  // record the dependency for every property currently under construction
  // (a derived property depends on everything its base properties read)
  if(fix && !build_stack_.empty() && 0 == strcmp(fix->style,"property/global"))
    for(size_t i = 0; i < build_stack_.size(); i++)
      add_dependency(build_stack_[i], fix);

  return fix;
}

/* ----------------------------------------------------------------------
   dependency tracking and in-place refresh (C-01/F-01)
------------------------------------------------------------------------- */

void PropertyRegistry::add_dependency(Entry &e, FixPropertyGlobal *fix)
{
  for(size_t i = 0; i < e.deps.size(); i++)
    if(e.deps[i].fix == fix) return;
  Dependency d;
  d.fix = fix;
  d.version = fix->version();
  e.deps.push_back(d);
}

void PropertyRegistry::inherit_dependencies(int kind, const std::string &varname)
{
  if(build_stack_.empty()) return;
  std::map<std::string,size_t>::iterator it = entry_index_[kind].find(varname);
  if(it == entry_index_[kind].end()) return;
  const std::vector<Dependency> &deps = entries_[it->second].deps;
  for(size_t i = 0; i < build_stack_.size(); i++)
    for(size_t j = 0; j < deps.size(); j++)
      add_dependency(build_stack_[i], deps[j].fix);
}

void PropertyRegistry::begin_build(int kind, const std::string &varname, const char *caller)
{
  Entry e;
  e.kind = kind;
  e.name = varname;
  e.caller = caller ? caller : "";
  build_stack_.push_back(e);
}

void PropertyRegistry::end_build()
{
  Entry e = build_stack_.back();
  build_stack_.pop_back();
  // snapshot the versions after construction (creators may have triggered
  // the first evaluation of v_ values)
  for(size_t i = 0; i < e.deps.size(); i++)
    e.deps[i].version = e.deps[i].fix->version();
  entry_index_[e.kind][e.name] = entries_.size();
  entries_.push_back(e);
}

ScalarProperty * PropertyRegistry::create_scalar(const std::string &varname, const char *caller)
{
  begin_build(KIND_SCALAR, varname, caller);
  ScalarProperty *p = (*scalar_creators[varname])(*this, caller, use_sanity_checks[varname]);
  scalars[varname] = p;
  end_build();
  return p;
}

VectorProperty * PropertyRegistry::create_vector(const std::string &varname, const char *caller)
{
  begin_build(KIND_VECTOR, varname, caller);
  VectorProperty *p = (*vector_creators[varname])(*this, caller, use_sanity_checks[varname]);
  vectors[varname] = p;
  end_build();
  return p;
}

MatrixProperty * PropertyRegistry::create_matrix(const std::string &varname, const char *caller)
{
  begin_build(KIND_MATRIX, varname, caller);
  MatrixProperty *p = (*matrix_creators[varname])(*this, caller, use_sanity_checks[varname]);
  matrices[varname] = p;
  end_build();
  return p;
}

void PropertyRegistry::rebuild(size_t k)
{
  // re-run the creator and copy the result into the existing object, so
  // every pointer handed out by connect()/get*Property() stays valid.
  // Copies are taken because a creator may append to entries_.
  const int kind = entries_[k].kind;
  const std::string name = entries_[k].name;
  const std::string caller = entries_[k].caller;

  if(KIND_SCALAR == kind)
  {
    ScalarProperty *target = scalars[name];
    ScalarProperty *fresh = (*scalar_creators[name])(*this, caller.c_str(), use_sanity_checks[name]);
    if(fresh && fresh != target)
    {
      target->data = fresh->data;
      delete fresh;
    }
    target->updateAll();   // scalar listeners hold copies of the value
  }
  else if(KIND_VECTOR == kind)
  {
    VectorProperty *target = vectors[name];
    VectorProperty *fresh = (*vector_creators[name])(*this, caller.c_str(), use_sanity_checks[name]);
    if(fresh && fresh != target)
    {
      if(fresh->cols != target->cols)
        error->all(FLERR,"internal error: property size changed during refresh");
      for(int c = 0; c < target->cols; c++) target->data[c] = fresh->data[c];
      delete fresh;
    }
  }
  else
  {
    MatrixProperty *target = matrices[name];
    MatrixProperty *fresh = (*matrix_creators[name])(*this, caller.c_str(), use_sanity_checks[name]);
    if(fresh && fresh != target)
    {
      if(fresh->rows != target->rows || fresh->cols != target->cols)
        error->all(FLERR,"internal error: property size changed during refresh");
      for(int r = 0; r < target->rows; r++)
        for(int c = 0; c < target->cols; c++)
          target->data[r][c] = fresh->data[r][c];
      delete fresh;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PropertyRegistry::refresh()
{
  // never while properties are being built (creators evaluate v_ values),
  // never recursively, and only while a run/minimize is set up: between runs
  // listeners may belong to deleted models, and Force::init() rebuilds anyway
  if(refreshing_ || !build_stack_.empty() || entries_.empty()) return;
  if(0 == update->whichflag) return;

  refreshing_ = true;

  // fixes that still exist (a dependency may have been unfixed)
  std::set<FixPropertyGlobal*> alive;
  for(int i = 0; i < modify->nfix; i++)
    if(modify->fix[i] && 0 == strcmp(modify->fix[i]->style,"property/global"))
      alive.insert(static_cast<FixPropertyGlobal*>(modify->fix[i]));

  // entries_ is in order of completed construction, so base properties
  // (e.g. youngsModulus, coefficientRestitution) are refreshed before the
  // derived ones built from them (Yeff, coeffRestLog, betaeff, ...)
  const size_t nentries = entries_.size();
  for(size_t k = 0; k < nentries; k++)
  {
    bool stale = false, valid = true;
    {
      const Entry &e = entries_[k];
      for(size_t i = 0; i < e.deps.size(); i++)
      {
        if(alive.find(e.deps[i].fix) == alive.end()) { valid = false; break; }
        if(e.deps[i].fix->version() != e.deps[i].version) stale = true;
      }
    }
    if(!valid || !stale) continue;

    rebuild(k);
    n_refreshed_++;
    Entry &e = entries_[k];
    for(size_t i = 0; i < e.deps.size(); i++)
      e.deps[i].version = e.deps[i].fix->version();
  }

  refreshing_ = false;
}

ScalarProperty * PropertyRegistry::getScalarProperty(string varname,const char *caller)
{
  if(scalars.find(varname) == scalars.end()) {
    if(scalar_creators.find(varname) != scalar_creators.end()) {
      create_scalar(varname, caller);
    } else {
      error->message(FLERR, "unknown scalar property");
    }
  }
  else
    inherit_dependencies(KIND_SCALAR, varname);
  return scalars[varname];
}

VectorProperty * PropertyRegistry::getVectorProperty(string varname,const char *caller)
{
  if(vectors.find(varname) == vectors.end()) {
    if(vector_creators.find(varname) != vector_creators.end()) {
      create_vector(varname, caller);
    } else {
      error->message(FLERR, "unknown vector property");
    }
  }
  else
    inherit_dependencies(KIND_VECTOR, varname);
  return vectors[varname];
}

MatrixProperty * PropertyRegistry::getMatrixProperty(string varname,const char *caller)
{
  if(matrices.find(varname) == matrices.end()) {
    if(matrix_creators.find(varname) != matrix_creators.end()) {
      create_matrix(varname, caller);
    } else {
      error->message(FLERR, "unknown matrix property");
    }
  }
  else
    inherit_dependencies(KIND_MATRIX, varname);
  return matrices[varname];
}

void PropertyRegistry::registerProperty(string varname, ScalarPropertyCreator creator, bool sanity_checks)
{
  if(scalar_creators.find(varname) == scalar_creators.end()) {
    scalar_creators[varname] = creator;
    use_sanity_checks[varname] = sanity_checks;
  } else if(scalar_creators[varname] != creator) {
    error->message(FLERR, "property with the same name, but different implementation registered");
  }
}

void PropertyRegistry::registerProperty(string varname, VectorPropertyCreator creator, bool sanity_checks)
{
  if(vector_creators.find(varname) == vector_creators.end()) {
    vector_creators[varname] = creator;
    use_sanity_checks[varname] = sanity_checks;
  } else if(vector_creators[varname] != creator) {
    error->message(FLERR, "property with the same name, but different implementation registered");
  }
}

void PropertyRegistry::registerProperty(string varname, MatrixPropertyCreator creator, bool sanity_checks)
{
  if(matrix_creators.find(varname) == matrix_creators.end()) {
    matrix_creators[varname] = creator;
    use_sanity_checks[varname] = sanity_checks;
  } else if(matrix_creators[varname] != creator) {
    error->message(FLERR, "property with the same name, but different implementation registered");
  }
}

void PropertyRegistry::connect(string varname, double ** & variable, const char *caller)
{
  if(matrices.find(varname) == matrices.end()) {
    if(matrix_creators.find(varname) != matrix_creators.end()) {
      create_matrix(varname, caller);
    } else {
      // ERROR unknown property
      error->message(FLERR, "unknown matrix property");
    }
  }
  matrices[varname]->connect(variable);
}

void PropertyRegistry::connect(string varname, double * & variable, const char *caller)
{
  if(vectors.find(varname) == vectors.end()) {
    if(vector_creators.find(varname) != vector_creators.end()) {
      create_vector(varname, caller);
    } else {
      // ERROR unknown property
      error->message(FLERR, "unknown vector property");
    }
  }
  vectors[varname]->connect(variable);
}

void PropertyRegistry::connect(string varname, double & variable, const char *caller)
{
  if(scalars.find(varname) == scalars.end()) {
    if(scalar_creators.find(varname) != scalar_creators.end()) {
      create_scalar(varname, caller);
    } else {
      // ERROR unknown property
      error->message(FLERR, "unknown scalar property");
    }
  }
  scalars[varname]->connect(variable);
}

void PropertyRegistry::init()
{
  for(std::map<string,ScalarProperty*>::iterator it = scalars.begin(); it != scalars.end(); ++it) {
      delete it->second;
  }
  for(std::map<string,VectorProperty*>::iterator it = vectors.begin(); it != vectors.end(); ++it) {
      delete it->second;
  }
  for(std::map<string,MatrixProperty*>::iterator it = matrices.begin(); it != matrices.end(); ++it) {
      delete it->second;
  }
  scalars.clear();
  vectors.clear();
  matrices.clear();

  entries_.clear();
  build_stack_.clear();
  for(int k = 0; k < 3; k++) entry_index_[k].clear();
}

void PropertyRegistry::print_all(FILE * out)
{
  for(std::map<string,ScalarProperty*>::iterator it = scalars.begin(); it != scalars.end(); ++it) {
      fprintf(out, " %s = ", it->first.c_str());
      it->second->print_value(out);
      fprintf(out, "\n");
  }
  for(std::map<string,VectorProperty*>::iterator it = vectors.begin(); it != vectors.end(); ++it) {
      fprintf(out, " %s = ", it->first.c_str());
      it->second->print_value(out);
      fprintf(out, "\n");
  }
  for(std::map<string,MatrixProperty*>::iterator it = matrices.begin(); it != matrices.end(); ++it) {
      fprintf(out, " %s = ", it->first.c_str());
      it->second->print_value(out);
      fprintf(out, "\n");
  }
}
