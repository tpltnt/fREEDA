// Curtice cubic MESFET model (from Microwave Harmonica Manual)
// with power disipation and thermal effects.
//
//          Drain 2
//                  o
//                  |
//                  |
//              |---+
//              |     -----o 4 Power out
// Gate 1 o-----|
//              |     -----o 5 Thermal reference
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

#ifndef MesfetCT_h
#define MesfetCT_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class MesfetCT : public ADInterface
{
	public:

  MesfetCT(const string& iname);

  ~MesfetCT() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // Get a vector with the indexes of the local reference nodes.
  // Get terminal pointers in term_list
  // ordered by local reference node:
  //
  // t0 t1 t2 t3
  //
  //     ^     ^
  //     LRN1  LRN2
  //
  // local_ref_vec contains: {1, 3}
  //
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Some constants
  double k2, k3, bm1;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double a0, a1, a2, a3, beta, vds0, gama, vt0;
  double cgs0, cgd0, is, n, ib0, nr;
  double t, vbi, fcc, vbd, area;
  double tnom, avt0, bvt0, tbet, tm, tme, eg, m, xti;
  bool kirchhoff;
  double b, ts;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
