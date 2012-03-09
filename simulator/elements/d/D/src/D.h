//		  Contains routines for analysis of compact software's
//		  microwave Diode.
//
//	      |\ |
//   0 o------| >|------o 1
//	      |/ |
// anode                 cathode
//
//
//	  Author:
//		Carlos E. Christoffersen
//              Mete Ozkar
//

#ifndef D_h
#define D_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class D : public ADInterface
{
	public:

  D(const string& iname);

  ~D() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // State variable transient analysis
  //virtual void svTran(TimeDomainSV *tdsv);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Some constants
  double v1;
  double k1, k2, k3, k4, k5, k6;

  // Internal state variables (charge)
  //  double qold, qnew;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double js;
  double alfa;
  double jb;
  double vb;
  double e;
  double ct0;
  double fi;
  double gama;
  double cd0;
  double afac;
  double r0;
  double t;
  double area;
  double imax;
  double eg;
  double m;
  double aro;
  double bro;
  double afag;
  double xti;
  //  double ind;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
