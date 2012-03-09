// MOSFET LEVEL 2 MODEL
//
//          Drain 2
//                  o
//                  |
//                  |
//              |---+
//              |
// Gate 1 o-----|------o 4 Bulk
//              |
//              |---+
//                  |
//                  |
//                  o
//         Source 3
//
//
//	  Author:
//		Aaron Walker

#ifndef Mosn2_h
#define Mosn2_h 1
#define NMOS 1
#define PMOS -1

#include "../../../../network/ADInterface.h"
#include "../../../../analysis/TimeDomainSV.h"

class Mosn2 : public ADInterface
{
public:

  Mosn2(const string& iname);

  ~Mosn2() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Some constants
  double pi, epsilon0, epsilons, leff, vt, vtnom;
  double Kox;
  // Model threshold voltage factors
  double eg, ni, vfb, vbi;
  double tdev;

  // Intermediate variables
  double sqrphi, gammasqr, xlambda, factor, eta, gamma_cmp;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double type, vt0, kp, gamma, phi, lambda, rd, rs;
  double is, pb,rsh;
  double js, tox, ld, u0, fc, nsub, tpg;
  double nss, tnom, kf, af, t, l, w, alpha;
  double delta, uexp, ucrit, vmax, xj, neff;
  double nfs;

  // Computed model variables
  double nssm, cox, egfet1, fermig, beta, fermis;
  double wkfng, wkfngs, coxwl, xd;

  // Temperature modification factors
  double ratio, fact1, ratio4, kt, egfet, arg, pbfact;
  double phio, sqrtphi, tref, arg1, pbfact1, fact2, ktnom;
  double pbo, tpb, sbiarg, stphi3;

  // Temperature modified parameters
  double tkp, tu0, tphi, tvbi, tvt0;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
