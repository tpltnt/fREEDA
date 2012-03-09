//          Drain 2
//                  o
//                  |
//                  |
//              |---+
//              |
// Gate 1 o-----|
//              |
//              |---+
//                  |
//                  |
//                  o
//         Source 3
//
//
//	  Author:
//          Shuping Zhang
//
// This is a Triquint (TOM) Intrinsic model


#ifndef MesfetTom_h
#define MesfetTom_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class MesfetTom : public ADInterface
{
	public:

  MesfetTom(const string& iname);

  ~MesfetTom() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

	// Some constants
  double VT, Veff, A1, Vnew, F1, F2, F3,  KTanh, Vp, Ids0, Vt;
  double delta, tn, Beta, Nn, is1, Vbi, Vt0, k1, k2, k3, k4, k5, k6;
  double Ebarr, EbarrN;
  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double vgsi, vdsi, vgdi, vmax, k, area, vds0;
  double beta, vt0, gama, q, delt, alfa, t, cgs0, cgd0, vbi, is, n, ib0, nr, vbd;
  double tj, t1, tnom, tbet, xti, a0, a1, a2, a3, fcc, avt0, bvt0, tm, tme, eg, m;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
