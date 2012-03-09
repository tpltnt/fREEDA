//		  Spice diode model (with charge conservation model)
//
//	      |\ |
//   0 o------| >|------o 1
//	      |/ |
// anode                 cathode
//
//
//	  Author:
//		Carlos E. Christoffersen
//

#ifndef DiodeSP_h
#define DiodeSP_h 1

#include "../../../../network/CircuitManager.h"
#include "../../../../network/ADInterface.h"

class DiodeSP : public ADInterface
{
	public:

  DiodeSP(const string& iname);

  ~DiodeSP() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Some constants
  double v1;
  //  double k1, k2, k3, k4, k5, k6;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double is, n, ibv, bv, fc, cj0, vj, m, tt, area, rs;
  bool charge;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
