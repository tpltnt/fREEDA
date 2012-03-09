// General Purpose Operational Amplifier
//
//
//                  Vdd 0
//                      o
//                      |
//                    | |
//           V+ 1 o---|V+
//                    |   \___ 3 Vout
//                    |   /
//           V- 2 o---|V-/
//                    |/|
//                      |
//                      o
//                  GND 5 o--/\/\/--o 4 Temp
//
//
//	  Author:
//		Gregory Parsons
//

#ifndef OpAmpT_h
#define OpAmpT_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

// The OpAmpT Class
class OpAmpT : public ADInterface
{
public:

  OpAmpT(const string& iname);

  ~OpAmpT() {}

  static const char* getNetlistName()
  {
     return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  private:

  // Evaluation function
  virtual void eval(AD * x, AD * effort, AD * flow);

  // Parameter Variables
  double gain, rin, rout, psr, cin, cout, pdr;
  double rint, routt, psrt, nomt;
  // Model Variables

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
