//	Nonlinear Thermal Boundary Conditions
//
//	      +--+
//   0 o------|  |------o 1
//	      +--+
// Discretized           Reference
// Surface               Temperature
//
//

#ifndef ThermalBlockBC1_h
#define ThermalBlockBC1_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class ThermalBlockBC1 : public ADInterface
{
	public:

  ThermalBlockBC1(const string& iname);

  ~ThermalBlockBC1() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Some constants
  static const double sigma;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double epsilon, hn, hf, T0;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
