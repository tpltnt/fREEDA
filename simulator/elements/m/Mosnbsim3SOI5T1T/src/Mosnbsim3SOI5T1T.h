// Bsim3nsoitest MOSFET MODEL
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

#ifndef Mosnbsim3SOI5T1T_h
#define Mosnbsim3SOI5T1T_h 1

#include "../../../../network/ADInterface.h"
#include "../../../../analysis/TimeDomainSV.h"

class Mosnbsim3SOI5T1T:public ADInterface
{
public:

  Mosnbsim3SOI5T1T(const string& iname);

  ~Mosnbsim3SOI5T1T() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

private:

  // Generic state variable evaluation routine
  virtual void eval(AD * x, AD * effort, AD * flow);

  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list);
  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double l, w, dtoxcv, llc, lwc, lwlc, wlc, wwc, wwlc, tsi, tox, toxref, tbox, tnom, rbody;
  double rbsh, rsh, rhalo, wint, lint, wth0, ll, wl, lln, wln, lw, ww, lwn, wwn, lwl, wwl, ln;
  double xpart, xj, k1b, k2b, dk2b, vbsa, aigc, bigc, cigc, aigsd, bigsd, cigsd, nigc;
  double pioxedge, pigcd,  k1w1, k1w2,  k3, k3b, kb1, w0, nlx,  nsub, ngate;
  double dvt0, dvt1, dvt2, dvt0w, dvt1w, dvt2w, eta0, etab, dsub, voff, nfactor, cdsc, cdscb, cdscd;
  double cit, u0, prwg, prwb, wr, rdsw, a0, ags, a1, a2, b0, b1, vsat, keta, ketas;
  double dwg, dwb, dwbc, pclm, pdibl1, pdibl2, pdiblb, drout, pvag, delta, alpha0, beta0;
  double beta1, beta2, fbjtii, vdsatii0, tii, lii, esatii, sii0, sii1, sii2, siid, agidl;
  double bgidl, ngidl, ebg, vgb1, vgb2, voxh, deltavox, ntox, ntun, ndiode, nrecf0, nrecr0;
  double isbjt, isdif, isrec, istun, vrec0, vtun0, nbjt, lbjt0, vabjt, aely,  vevb, vecb;
  double cjswg, mjswg, pbswg, tt, ldif0, cgeo, cgso, cgdo, dlc, dwc, dlcb, dlbg, fbody, clc;
  double cle, cf, csdmin, asd, csdesw,  vsdth, delvt, acde, moin, ckappa, cgdl, cgsl;
  double ndif, rth0, cth0, tpbswg, tcjswg, kt1, kt1l, kt2, ute, ua1, ub1, uc1, prt, at, ntrecf;
  double ntrecr, xbjt, xdif, xrec, xtun, dlcig, nbc, nseg, pdbcp, psbcp, toxqm, type, toxm;
  double xt1, dvbd0, dvbd1, temp,npeak,capMod, vbm,nofffd, vofffd,moinFD,shmod ;
  double ahli,ua,ub,uc,vsdfb,k2,k1,vth0,nch;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
