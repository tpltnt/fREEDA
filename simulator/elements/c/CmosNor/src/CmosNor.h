// CMOS NOR Gate
//
//
//                                       Vdd 1
//                                          o
//                                          |
//                                      |---+
//                                      |
//            |------------------------O|
//            |                         |
//            |                         |---+
//            |                             |
//            |                             |
//            |                         |---+
//            |                         |
//            |                   |----O|
//            |                   |     |
//            |                   |     |---+
//            |                   |         |
//            |        Input2 3 o-|         |
//            |     *-------------|---------*---o 4 Output
//            |     |             |         |
//            | |---+             |     |---+
//            | |                 |     |
// Input1 2 o-|-|                 |-----|
//              |                       |
//              |---+                   |---+
//                  |                       |
//                  |                       |
//                  ----------------------*--
//                                        |
//                                        |
//                                        o
//                                     GND 5
//
//
//	  Author:
//		Mazen M Kharbutli
//

#ifndef CmosNor_h
#define CmosNor_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

// The CmosNor Class
class CmosNor : public ADInterface
{
	public:

  CmosNor(const string& iname);

  ~CmosNor() {}

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
  double vtn, vtp, un, up, en, ep, tox, ln, wn, lp, wp, td;

  // Model Variables
  double betan, betap;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
