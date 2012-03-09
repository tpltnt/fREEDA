//		  Contains routines for analysis of a tunnel diode
//
//	      |\ |
//   0 o------| >|------o 1
//	      |/ |
// anode                cathode
//
//
//	  Author:
//		Stephen Bowyer
//              Jennifer Huckaby
//

#ifndef DiodeTun_h
#define DiodeTun_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class DiodeTun : public ADInterface
{
	public:

  DiodeTun(const string& iname);

  ~DiodeTun() {}

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
  double k, q;
  double vt, vx, vth, vxth;
  double k1, k2, k4, k5, k6;
  double at1, ax1, ath1, at2, ax2, ath2, b, ct, cx, cth;

  // Internal state variables (charge)
  //  double qold, qnew;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double js;
  double ct0;
  double fi;
  double gama;
  double cd0;
  double afac;
  double r0;
  double t;
  double area;
  double jv;
  double jp;
  double vv;
  double vpk;
  double a2;
  double mt;
  double mx;
  double mth;
  double temper;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
