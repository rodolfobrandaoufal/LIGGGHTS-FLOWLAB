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

    Richard Berger (JKU Linz)

    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>
#include "lmptype.h"
#include "force.h"
#include <string>
#include <map>
#include <iostream>
#include "contact_interface.h"
#include "contact_model_constants.h"
#include <sstream>

namespace LIGGGHTS {
using namespace LAMMPS_NS;

namespace Utils {

  inline int64_t generate_gran_hashcode(int model, int tangential, int cohesion, int rolling, int surface)
  {
    return (((int64_t)model)           ) |
           (((int64_t)tangential) <<  6) |
           (((int64_t)cohesion)   << 12) |
           (((int64_t)rolling)    << 18) |
           (((int64_t)surface)    << 24) ;
  }

  /* decode the 6-bit fields of a granular contact-model hash (see above)
     index: 0 normal, 1 tangential, 2 cohesion, 3 rolling, 4 surface */
  inline int gran_hashcode_field(int64_t hash, int index)
  {
    return static_cast<int>((hash >> (6*index)) & ((1<<6) - 1));
  }

  /* name of one sub-model of a contact-model hash.
     as_identifier = true : the C++ identifier used in GRAN_MODEL(...) lines,
                            e.g. HERTZ, TANGENTIAL_HISTORY, COHESION_OFF
     as_identifier = false: the keyword used in the input script,
                            e.g. hertz, history, off                     */
  inline const char * gran_submodel_name(int64_t hash, int index, bool as_identifier)
  {
    const int id = gran_hashcode_field(hash, index);
    switch(index)
    {
    case 0:
      if(id == ContactModels::NORMAL_OFF) return as_identifier ? "NORMAL_OFF" : "off";
      #define NORMAL_MODEL(identifier,str,constant) \
      if(id == constant) return as_identifier ? #identifier : #str;
      #include "style_normal_model.h"
      #undef NORMAL_MODEL
      break;
    case 1:
      if(id == ContactModels::TANGENTIAL_OFF) return as_identifier ? "TANGENTIAL_OFF" : "off";
      #define TANGENTIAL_MODEL(identifier,str,constant) \
      if(id == constant) return as_identifier ? #identifier : #str;
      #include "style_tangential_model.h"
      #undef TANGENTIAL_MODEL
      break;
    case 2:
      if(id == ContactModels::COHESION_OFF) return as_identifier ? "COHESION_OFF" : "off";
      #define COHESION_MODEL(identifier,str,constant) \
      if(id == constant) return as_identifier ? #identifier : #str;
      #include "style_cohesion_model.h"
      #undef COHESION_MODEL
      break;
    case 3:
      if(id == ContactModels::ROLLING_OFF) return as_identifier ? "ROLLING_OFF" : "off";
      #define ROLLING_MODEL(identifier,str,constant) \
      if(id == constant) return as_identifier ? #identifier : #str;
      #include "style_rolling_model.h"
      #undef ROLLING_MODEL
      break;
    case 4:
      #define SURFACE_MODEL(identifier,str,constant) \
      if(id == constant) return as_identifier ? #identifier : #str;
      #include "style_surface_model.h"
      #undef SURFACE_MODEL
      break;
    }
    return "UNKNOWN";
  }

  /* "model hertz tangential history cohesion off rolling_friction epsd2 surface default" */
  inline std::string gran_hashcode_to_keywords(int64_t hash)
  {
    std::string s("model ");
    s += gran_submodel_name(hash, 0, false);
    s += " tangential ";       s += gran_submodel_name(hash, 1, false);
    s += " cohesion ";         s += gran_submodel_name(hash, 2, false);
    s += " rolling_friction "; s += gran_submodel_name(hash, 3, false);
    s += " surface ";          s += gran_submodel_name(hash, 4, false);
    return s;
  }

  /* "GRAN_MODEL(HERTZ, TANGENTIAL_HISTORY, COHESION_OFF, ROLLING_EPSD2, SURFACE_DEFAULT)" */
  inline std::string gran_hashcode_to_whitelist_entry(int64_t hash)
  {
    std::string s("GRAN_MODEL(");
    for(int i = 0; i < 5; i++)
    {
      s += gran_submodel_name(hash, i, true);
      s += (i < 4) ? ", " : ")";
    }
    return s;
  }

  /* remedy text shared by the fallback warning and the whitelist errors */
  inline std::string gran_whitelist_remedy(int64_t hash)
  {
    std::string s("To compile this combination statically, add the line\n    ");
    s += gran_hashcode_to_whitelist_entry(hash);
    s += "\n  to src/contact_model_whitelist.txt (tracked, used by Make and by CMake by default)"
         " or to src/style_contact_model_user.whitelist (local, Make route),"
         " or run src/genAutoExamplesWhitelist.sh on your input script,"
         " or configure CMake with -DLIGGGHTS_CONTACT_WHITELIST=<file> / the ENABLE_MODEL_* options;"
         " then recompile.";
    return s;
  }

  inline std::string int_to_string(int a)
  {
    // see https://www.cfdem.com/forums/error-non-const-lvalue-reference-type-basicostringstream-cannot-bind-temporary-type
    // return static_cast< std::ostringstream & >(( std::ostringstream() << std::dec << a ) ).str();
    std::ostringstream ss;
    ss << std::dec << a;
    return ss.str();
  }

  inline std::string double_to_string(double dbl)
  {
    std::ostringstream strs;
    strs << dbl;
    std::string str = strs.str();
    return str;
  }

  template <typename T>
  inline T* ptr_reduce(T** &t)
  { return &(t[0][0]); }

  template <typename T>
  inline T* ptr_reduce(T* &t)
  { return t; }

  template <typename T>
  inline T* ptr_reduce(T &t)
  { return &t; }

  template<typename Interface>
  class AbstractFactory {
    typedef typename Interface::ParentType ParentType;
    typedef Interface * (*Creator)(class LAMMPS * lmp, ParentType* parent, int64_t hash);
    typedef int64_t (*VariantSelector)(int & argc, char ** & argv, Custom_contact_models ccm);
    // finding P0-18: the key carries the full 64-bit contact-model hash
    // (generate_gran_hashcode / GranStyle::HASHCODE are int64_t); an int key
    // would truncate it, and two hashes differing only in the upper 32 bits
    // would collide
    typedef std::pair<std::string, int64_t> StyleKey;
    typedef std::map<StyleKey, Creator> StyleTable;
    typedef std::map<std::string, VariantSelector> VariantSelectorTable;
    StyleTable styleTable;
    VariantSelectorTable variantSelectorTable;
    AbstractFactory(const AbstractFactory &){}
    AbstractFactory& operator=(const AbstractFactory&){}

  protected:
    AbstractFactory() {}

  public:
    Interface * create(const std::string & name, int64_t variant, class LAMMPS * lmp, ParentType* parent) {
      StyleKey key(name, variant);
      if(styleTable.find(key) != styleTable.end()) {
        return styleTable[key](lmp, parent, variant);
      }
      // not in the static whitelist: use the runtime-composed fallback
      // ContactModel<GranStyle<> > (registered under the all-OFF hash in
      // granular_styles.h unless LIGGGHTS_NO_CONTACT_MODEL_FALLBACK is set)
      int64_t default_variant = generate_gran_hashcode(ContactModels::NORMAL_OFF, ContactModels::TANGENTIAL_OFF, ContactModels::COHESION_OFF, ContactModels::ROLLING_OFF, 0);
      StyleKey default_key(name, default_variant);
      if(styleTable.find(default_key) != styleTable.end()) {
        return styleTable[default_key](lmp, parent, variant);
      }
      return NULL;
    }

    // true if 'variant' is compiled as a static (whitelisted) combination
    bool hasStaticVariant(const std::string & name, int64_t variant) {
      StyleKey key(name, variant);
      return styleTable.find(key) != styleTable.end();
    }

    int64_t selectVariant(const std::string & name, int & argc, char ** & argv,Custom_contact_models ccm) {
      if(variantSelectorTable.find(name) != variantSelectorTable.end()) {
        return variantSelectorTable[name](argc, argv,ccm);
      }
      return 0;
    }

    void addStyle(const std::string & name, int64_t variant, Creator create) {
      StyleKey key(name, variant);
      if(styleTable.find(key) != styleTable.end()){
        std::cerr << "WARNING! Style collision detected! Duplicate entry (" << key.first << ", " << key.second << ") in style table." << std::endl;
      }
      styleTable[key] = create;
    }

    void addVariantSelector(const std::string & name, VariantSelector selector) {
      if(variantSelectorTable.find(name) != variantSelectorTable.end()){
        std::cerr << "WARNING! VariantSelector collision detected! Duplicate entry '" << name << "' in variant selector table." << std::endl;
      }
      variantSelectorTable[name] = selector;
    }
  };
}

}

#endif // UTILS_H
