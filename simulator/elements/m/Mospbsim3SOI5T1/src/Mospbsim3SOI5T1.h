// Mospbsim3SOI5T1 MOSFET MODEL
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

#ifndef Mospbsim3SOI5T1_h
#define Mospbsim3SOI5T1_h 1

#include "../../../../network/ADInterface.h"
#include "../../../../analysis/TimeDomainSV.h"

class Mospbsim3SOI5T1:public ADInterface
{
	public:

  Mospbsim3SOI5T1(const string& iname);

  ~Mospbsim3SOI5T1() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  // Generic state variable evaluation routine
  virtual void eval(AD * x, AD * eval, AD * flow);

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables

  double l, w, tsi, tbox, tox, toxqm, xj, nch, nsub, ngate, vth0, k1, k1w1, k1w2;
  double k2, k3, k3b, kb1, w0, nlx, dvt0, dvt1, dvt2, dvt0w, dvt1w, dvt2w, u0;
  double ua, ub, uc, vsat, a0, ags, b0, b1, keta, ketas, a1, a2, rdsw, prwb, prwg;
  double wr, nfactor, wint, lint, dwg, dwb, dwbc, voff, eta0, etab, dsub, cit;
  double cdsc, cdscb, cdscd, pclm, pdibl1, pdibl2, pvag, delta, alpha0, beta0, beta1, beta2;
  double vdsatii0, cgeo, cjswg, pbswg, mjswg, tt, csdesw, cgsl, cgdl, ckappa, clc, cle;
  double dlc, dlcb, dlbg, dwc, delvt, fbody, moin, tnom, ute, kt1, kt1l, kt2;
  double ua1, ub1, uc1, at, prt, vbm, xt1, pdiblb, ll, llc, lln, lw, lwc, lwn;
  double lwl, lwlc, wl, wlc, wln, ww, wwc, wwn, wwl, wwlc, temp, acde;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
