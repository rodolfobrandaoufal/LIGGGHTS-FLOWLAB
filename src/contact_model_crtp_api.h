/* ----------------------------------------------------------------------
   Compile-time contact model plugin API prototype.

   The design keeps runtime selection outside the neighbor/contact loop and
   uses CRTP inside the loop, so model calls are statically bound and can be
   inlined by the compiler.
------------------------------------------------------------------------- */

#ifndef LMP_CONTACT_MODEL_CRTP_API_H
#define LMP_CONTACT_MODEL_CRTP_API_H

#include <cmath>

namespace LIGGGHTS {
namespace ContactModels {
namespace CRTP {

struct ContactInput {
  int i;
  int j;
  int itype;
  int jtype;
  double deltan;
  double radi;
  double radj;
  double meff;
  double vn;
  double vtr1;
  double vtr2;
  double vtr3;
  double en[3];
};

struct ContactForce {
  double fx_i;
  double fy_i;
  double fz_i;
  double fx_j;
  double fy_j;
  double fz_j;
  double tx_i;
  double ty_i;
  double tz_i;
  double tx_j;
  double ty_j;
  double tz_j;
  double Fn;
  double kt;
  double gamman;
  double gammat;

  void clear()
  {
    fx_i = fy_i = fz_i = 0.0;
    fx_j = fy_j = fz_j = 0.0;
    tx_i = ty_i = tz_i = 0.0;
    tx_j = ty_j = tz_j = 0.0;
    Fn = kt = gamman = gammat = 0.0;
  }
};

template <class Derived>
class ContactModelBase {
 public:
  inline void evaluate(const ContactInput &in, ContactForce &out) const
  {
    static_cast<const Derived *>(this)->evaluate_impl(in, out);
  }
};

class HertzMindlin : public ContactModelBase<HertzMindlin> {
 public:
  HertzMindlin(const double **youngs_modulus_eff,
               const double **shear_modulus_eff,
               const double **coeff_restitution_log,
               double sqrt_five_over_six)
      : Yeff_(youngs_modulus_eff),
        Geff_(shear_modulus_eff),
        beta_(coeff_restitution_log),
        sqrt_five_over_six_(sqrt_five_over_six)
  {}

  inline void evaluate_impl(const ContactInput &in, ContactForce &out) const
  {
    const double reff = (in.radi * in.radj) / (in.radi + in.radj);
    const double sqrt_delta_reff = std::sqrt(reff * in.deltan);
    const double Sn = 2.0 * Yeff_[in.itype][in.jtype] * sqrt_delta_reff;
    const double St = 8.0 * Geff_[in.itype][in.jtype] * sqrt_delta_reff;

    out.kt = St;
    out.gamman = -2.0 * sqrt_five_over_six_ *
                 beta_[in.itype][in.jtype] * std::sqrt(Sn * in.meff);
    out.gammat = -2.0 * sqrt_five_over_six_ *
                 beta_[in.itype][in.jtype] * std::sqrt(St * in.meff);

    const double Fn_contact = (4.0 / 3.0) *
                              Yeff_[in.itype][in.jtype] *
                              sqrt_delta_reff * in.deltan;
    const double Fn_damping = -out.gamman * in.vn;
    out.Fn = Fn_contact + Fn_damping;

    out.fx_i += out.Fn * in.en[0];
    out.fy_i += out.Fn * in.en[1];
    out.fz_i += out.Fn * in.en[2];
    out.fx_j -= out.Fn * in.en[0];
    out.fy_j -= out.Fn * in.en[1];
    out.fz_j -= out.Fn * in.en[2];
  }

 private:
  const double **Yeff_;
  const double **Geff_;
  const double **beta_;
  double sqrt_five_over_six_;
};

class NoCohesion : public ContactModelBase<NoCohesion> {
 public:
  inline void evaluate_impl(const ContactInput &, ContactForce &) const {}
};

class LinearTangentialDamping : public ContactModelBase<LinearTangentialDamping> {
 public:
  inline void evaluate_impl(const ContactInput &in, ContactForce &out) const
  {
    const double Ft1 = -out.gammat * in.vtr1;
    const double Ft2 = -out.gammat * in.vtr2;
    const double Ft3 = -out.gammat * in.vtr3;
    out.fx_i += Ft1;
    out.fy_i += Ft2;
    out.fz_i += Ft3;
    out.fx_j -= Ft1;
    out.fy_j -= Ft2;
    out.fz_j -= Ft3;
  }
};

template <class NormalModel, class CohesionModel, class TangentialModel>
class ContactPipeline {
 public:
  ContactPipeline(const NormalModel &normal,
                  const CohesionModel &cohesion,
                  const TangentialModel &tangential)
      : normal_(normal), cohesion_(cohesion), tangential_(tangential)
  {}

  inline void evaluate(const ContactInput &in, ContactForce &out) const
  {
    out.clear();
    normal_.evaluate(in, out);
    cohesion_.evaluate(in, out);
    tangential_.evaluate(in, out);
  }

 private:
  NormalModel normal_;
  CohesionModel cohesion_;
  TangentialModel tangential_;
};

struct AtomSoAView {
  double **x;
  double **v;
  double **f;
  double *radius;
  double *rmass;
  int *type;
  int *mask;
  int nlocal;
};

struct HalfNeighborListView {
  int inum;
  int *ilist;
  int *numneigh;
  int **firstneigh;
};

template <class Pipeline>
class PairGranCRTPKernel {
 public:
  explicit PairGranCRTPKernel(const Pipeline &pipeline) : pipeline_(pipeline) {}

  inline void compute(AtomSoAView &atom,
                      const HalfNeighborListView &list,
                      int groupbit,
                      double cutsq) const
  {
    ContactInput in;
    ContactForce out;

    for (int ii = 0; ii < list.inum; ++ii) {
      const int i = list.ilist[ii];
      if (!(atom.mask[i] & groupbit)) continue;

      const double xi = atom.x[i][0];
      const double yi = atom.x[i][1];
      const double zi = atom.x[i][2];
      const int itype = atom.type[i];
      int * const neighs = list.firstneigh[i];
      const int n = list.numneigh[i];

      for (int jj = 0; jj < n; ++jj) {
        const int j = neighs[jj];
        const double dx = xi - atom.x[j][0];
        const double dy = yi - atom.x[j][1];
        const double dz = zi - atom.x[j][2];
        const double rsq = dx*dx + dy*dy + dz*dz;
        if (rsq >= cutsq) continue;

        const double r = std::sqrt(rsq);
        const double radsum = atom.radius[i] + atom.radius[j];
        const double deltan = radsum - r;
        if (deltan <= 0.0) continue;

        const double rinv = 1.0 / r;
        in.i = i;
        in.j = j;
        in.itype = itype;
        in.jtype = atom.type[j];
        in.radi = atom.radius[i];
        in.radj = atom.radius[j];
        in.deltan = deltan;
        in.meff = (atom.rmass[i] * atom.rmass[j]) /
                  (atom.rmass[i] + atom.rmass[j]);
        in.en[0] = dx * rinv;
        in.en[1] = dy * rinv;
        in.en[2] = dz * rinv;

        const double vr1 = atom.v[i][0] - atom.v[j][0];
        const double vr2 = atom.v[i][1] - atom.v[j][1];
        const double vr3 = atom.v[i][2] - atom.v[j][2];
        in.vn = vr1*in.en[0] + vr2*in.en[1] + vr3*in.en[2];
        in.vtr1 = vr1 - in.vn*in.en[0];
        in.vtr2 = vr2 - in.vn*in.en[1];
        in.vtr3 = vr3 - in.vn*in.en[2];

        pipeline_.evaluate(in, out);

        atom.f[i][0] += out.fx_i;
        atom.f[i][1] += out.fy_i;
        atom.f[i][2] += out.fz_i;
        if (j < atom.nlocal) {
          atom.f[j][0] += out.fx_j;
          atom.f[j][1] += out.fy_j;
          atom.f[j][2] += out.fz_j;
        }
      }
    }
  }

 private:
  Pipeline pipeline_;
};

enum ContactPluginId {
  CONTACT_PLUGIN_HERTZ_MINDLIN = 0
};

template <ContactPluginId Id>
struct ContactPluginFactory;

template <>
struct ContactPluginFactory<CONTACT_PLUGIN_HERTZ_MINDLIN> {
  typedef ContactPipeline<HertzMindlin,
                          NoCohesion,
                          LinearTangentialDamping> Pipeline;

  static Pipeline create(const double **Yeff,
                         const double **Geff,
                         const double **beta,
                         double sqrt_five_over_six)
  {
    return Pipeline(HertzMindlin(Yeff, Geff, beta, sqrt_five_over_six),
                    NoCohesion(),
                    LinearTangentialDamping());
  }
};

template <ContactPluginId Id>
inline void compute_with_plugin(AtomSoAView &atom,
                                const HalfNeighborListView &list,
                                int groupbit,
                                double cutsq,
                                const double **Yeff,
                                const double **Geff,
                                const double **beta,
                                double sqrt_five_over_six)
{
  typename ContactPluginFactory<Id>::Pipeline pipeline =
      ContactPluginFactory<Id>::create(Yeff, Geff, beta, sqrt_five_over_six);
  PairGranCRTPKernel<typename ContactPluginFactory<Id>::Pipeline> kernel(pipeline);
  kernel.compute(atom, list, groupbit, cutsq);
}

} // namespace CRTP
} // namespace ContactModels
} // namespace LIGGGHTS

#endif
