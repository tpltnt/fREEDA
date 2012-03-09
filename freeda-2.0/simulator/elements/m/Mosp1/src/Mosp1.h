// MOSFET LEVEL 1 MODEL
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

#ifndef Mosp1_h
#define Mosp1_h 1
#define NMOS 1
#define PMOS -1

#include "../../../../network/ADInterface.h"
#include "../../../../analysis/TimeDomainSV.h"

class Mosp1 : public ADInterface
{
	public:

  Mosp1(const string& iname);

  ~Mosp1() {}

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
  double sqrphi, gammasqr;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double type, vt0, kp, gamma, phi, lambda, rd, rs;
  double cbd, cbs, is, pb, cgdo, cgso, cgbo, rsh, cj, mj;
  double cjsw, mjsw, js, tox, ld, u0, fc, nsub, tpg;
  double nss, tnom, kf, af, t, l, w, alpha;

  // Computed model variables
  double nssm, cox, egfet1, fermig, beta, fermis;
  double wkfng, wkfngs, coxwl;

  // Temperature modification factors
  double ratio, fact1, ratio4, kt, egfet, arg, pbfact;
  double phio, sqrtphi, tref, arg1, pbfact1, fact2, ktnom;

  // Temperature modified parameters
  double tkp, tu0, tphi, tvbi, tvt0;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
