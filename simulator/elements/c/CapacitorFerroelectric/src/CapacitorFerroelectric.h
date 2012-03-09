
// Ferroelectric capacitor model (2-terminal device)
//               |  |
//      0 o------|  |------o 1
//               |  |
// anode                 cathode
//
//
//	  Author:
//		Vrinda Haridasan

#ifndef CapacitorFerroelectric_h
#define CapacitorFerroelectric_h 1

#include "../../../../network/CircuitManager.h"
#include "../../../../network/ADInterface.h"

class CapacitorFerroelectric : public ADInterface
{
	public:

  // Constructor
  CapacitorFerroelectric(const string& iname);

  // Destructor
  ~CapacitorFerroelectric() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  // Generic state variable evaluation routine
  virtual void eval(AD * x, AD * effort, AD * flow);

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double epsid;
  double epsb0;
  double a;
  double d;
  double t;
  double k;
  double alpha3;
  double T0;
  double T;
  double beta;
  double p;

  // Init routine variables
 double ci, cmax, cmaxinv, cnst, cf, alpha1, cbmax, b, epsb;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
