// Mosnbsim3 MOSFET MODEL
//
//                Drain 1
//                  o
//                  |
//                  |
//              |---+
//              |
// Gate 2 o-----|------o 4 Bulk
//              |
//              |---+
//                  |
//                  |
//                  o
//               Source 3
//
//
//	  Author: Ramya Mohan

#ifndef Mosnbsim3_h
#define Mosnbsim3_h 1

#include "../../../../network/ADInterface.h"
#include "../../../../analysis/TimeDomainSV.h"

class Mosnbsim3:public ADInterface
{
	public:

  Mosnbsim3(const string& iname);

  ~Mosnbsim3() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  // Generic state variable evaluation routine
  virtual void eval(AD * x, AD * effort, AD * flow);

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double l, w, tox, toxm, cdsc, cdscb, cdscd, cit, nfactor, xj, vsat, at, a0;
  double ags, a1, a2, keta, nsub, nch, ngate, vbm, xt1, kt1, kt1l, kt2;
  double k3, k3b, w0, nlx, dvt0, dvt1, dvt2, dvt0w, dvt1w, dvt2w, drout, dsub;
  double ua, ub, uc, u0, voff, tnom, elm, delta, rdsw, prwg, prwb;
  double prt, eta0, etab, pclm, pdibl1, pdibl2, pdiblb, pscbe1;
  double pscbe2, pvag, vfb, acde, moin, noff, voffcv, lint, ll, llc, lln;
  double lw, lwc, lwn, lwl, lwlc, wr, wint, dwg, dwb, wl, wlc, wln, ww, wwc;
  double wwn, wwl, wwlc, b0, b1, clc, cle, alpha0, alpha1, beta0, ute, k1, k2;
  double temp, ua1, ub1, uc1;

  bool pmos;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
