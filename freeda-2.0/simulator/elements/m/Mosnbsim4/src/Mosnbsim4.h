// BSIM4 MOSFET MODEL
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
//	  Author: Nikhil Kriplani

#ifndef Mosnbsim4_h
#define Mosnbsim4_h 1

#include "../../../../network/ADInterface.h"
#include "../../../../analysis/TimeDomainSV.h"

class Mosnbsim4:public ADInterface
{
public:

  Mosnbsim4(const string& iname);

  ~Mosnbsim4() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double TOXE, TOXP, EPSROX, VFB ,VTH0, NGATE, XL, XW, NF, W, L;
  double DWG, DWB, WINT, WL, WLN, WW, WWN, WWL, LINT, LL, LLN, LW;
  double LWN, LWL, K1, K2, LPEB, LPE0, K3, K3B, W0, DVT0W, DVT0;
  double DVT1W, DVT1, DSUB, ETA0, ETAB, TOXM, T, NDEP, PHIN, VBM, NSUB;
  double DVT2W, NSD, DVT2, MINV, NFACTOR, CDSC, CDSCD, CDSCB, CIT, KETA;
  double B0, B1, A0, AGS, XJ, U0, UA, UB, UC, EU, DELTA, PDITS, FPROUT, PDITSL;
  double PDITSD, PSCBE2, PSCBE1, PDIBLCB, PVAG, PDIBL1, PDIBL2, DROUT, PCLM;
  double A1, A2, RDSWMIN, RDSW, PRWG, PRWB, WR, WLC, WWC;
  double DWJ, CLC, CLE, NOFF, VOFFCV, CF, CKAPPAD, CKAPPAS;
  double LLC, LWC, LWLC, WWLC, CGB0, VOFF, VOFFL, POXEDGE, TOXREF, NTOX, DLCIG;
  double AIGSD, BIGSD, CIGSD, MOIN, VSAT, AIGC, BIGC, CIGC, NIGC;
  double PIGCD, DVTP0, DVTP1, PRT, AT, XT, ALPHA0, ALPHA1, BETA0;
  double AGIDL, BGIDL, CGIDL, EGIDL, ACDE, DLC, DWC;
  double AIGBACC, BIGBACC, CIGBACC, NIGBACC, AIGBINV, BIGBINV, CIGBINV;
  double EIGBINV, NIGBINV, KT1l, KT1, KT2;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
