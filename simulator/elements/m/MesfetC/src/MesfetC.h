// Curtice cubic MESFET model (from Microwave Harmonica Manual)
//
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
//		  Carlos E. Christoffersen
//                Hector Gutierrez

#ifndef MesfetC_h
#define MesfetC_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class MesfetC : public ADInterface
{
	public:

  MesfetC(const string& iname);

  ~MesfetC() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  virtual void eval(AD * x,	AD * effort, AD * flow);

  // Some constants
  double k2, k3;
  double delta_T, tn, Vt, k1, k4, k5, k6, Vt0, Beta, Ebarr, EbarrN, Nn;
  double Is, Vbi;
  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double a0, a1, a2, a3, beta, vds0, gama, vt0;
  double cgs0, cgd0, is, n, ib0, nr;
  double t, vbi, fcc, vbd, area;
  double tnom, avt0, bvt0, tbet, tm, tme, eg, m, xti, tj;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
