//		  Microwave diode model with charge conservation model
//                (new adolc interface)
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

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

#ifndef DiodeCompact_h
#define DiodeCompact_h 1

class DiodeCompact : public ADInterface
{
	public:

  DiodeCompact(const string& iname);

  ~DiodeCompact() {}

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
  double js;
  double alfa;
  double jb;
  double vb;
  double e;
  double ct0;
  double fi;
  double gama;
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

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
