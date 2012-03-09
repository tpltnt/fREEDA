// The microwave diode.
//
//	      |\ |
//   0 o------| >|------o 1
//	      |/ |
// anode                 cathode
//
//
// Authors:
// Carlos E. Christoffersen and Mete Ozkar

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

#ifndef Diode_h
#define Diode_h 1

class Diode : public ADInterface
{
public:

  Diode(const string& iname);

  ~Diode() {}

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
  double k1, k2, k3, k4, k5, k6;

  // Internal state variables (charge)
  //  double qold, qnew;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double js, alfa, jb, vb, e, ct0, fi, gama, cd0, afac;
  double r0, t, area, imax, eg, m, aro, bro, afag, xti;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
