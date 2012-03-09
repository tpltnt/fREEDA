// ThermalInt is a two-port element used at interfaces with different
// Kirchhoff transformation parameters.
//
//	          +-----+
//             o--+     +--o
//   Medium 1     |     |     Medium 2
//	          +--+--+
//                   |
//                   o

#ifndef ThermalInt_h
#define ThermalInt_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class ThermalInt : public ADInterface
{
	public:

  ThermalInt(const string& iname);

  ~ThermalInt() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  double bm_1, bm_2;
  // Element information
  static ItemInfo einfo;
  // Number of parameters of this element
  static const unsigned n_par;
  // Parameter information
  static ParmInfo pinfo[];
  // Parameter variables
  double ts, b1, b2;
};

#endif
