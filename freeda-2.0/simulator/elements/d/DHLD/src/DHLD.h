// This is a Double Heterojunction laser Diode model (The R. S. Tucker model)
//
//               Laser Diode
//
//              ++++++++++++++
//              +            +
//          o---------       -----o
//              +    |       +
//              +  -----     +
//              +  \   /     +
//              +   \ /      +
//              +   ---      +
//              +    |       +
//          o---------       -----o
//       "gnd"  +            +   "oref"
//              +            +
//              ++++++++++++++
//
//
//
// by Houssam Kanj

#ifndef DHLD_h
#define DHLD_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class DHLD : public ADInterface
{
	public:

  DHLD(const string& iname);

  ~DHLD() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  virtual void init() throw(string&);
	virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
	TerminalVector& term_list);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Some constants
  double Vpara;  // like V1 in V. Rizzoli
  double vt, q;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double rs, re, i01, i02, b, tns, c0, vd, d, a, rp, cp, sc, beta;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
