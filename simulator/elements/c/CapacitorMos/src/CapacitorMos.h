// MOS capacitor model (2-terminal device) with tunneling
// Essentially behaves like a non-linear resistor as of now
//               |  |
//      0 o------|  |------o 1
//               |  |
// anode                 cathode
//
//
// Authors:
// Krishnanshu Dandu, Yawei Jin

#ifndef CapacitorMos_h
#define CapacitorMos_h 1

#include "../../../../network/CircuitManager.h"
#include "../../../../network/ADInterface.h"

class CapacitorMos : public ADInterface
{
	public:

  CapacitorMos(const string& iname);

  ~CapacitorMos() {}

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
  // double k1, k2, k3, k4, k5, k6;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double w, l, tox;
  double t, npoly, nsub;
  double moxecb, moxevb, moxhvb;
  double phibecb, phibevb, phibhvb;
  double phib0ecb, phib0evb, phib0hvb;
  double vfb, epsrox, s;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
