// CMOS NAND Gate
//
//
//                                           Vdd 1
//                                              o
//                                              |
//                      *-----------------------*
//                      |                       |
//                  |---+                   |---+
//                  |                       |
// Input1 2 o--|---O|    Input2 3 o--|-----O|
//             |    |                |      |
//             |    |---+            |      |---+
//             |        |            |          |
//             |        |            |          |
//             |        -------------|----------*--o 4 Output
//             |                     |          |
//             |                     |          |
//             |                     |          |
//             |                     |      |---+
//             |                     |      |
//             |                     |------|
//             |                            |
//             |                            |---+
//             |                                |
//             |                                |
//             |                            |---+
//             |                            |
//             |----------------------------|
//                                          |
//                                          |---+
//                                              |
//                                              |
//                                              o
//                                           GND 5
//
//
//	  Author:
//		Mazen M Kharbutli
//

#ifndef CmosNand_h
#define CmosNand_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

// The CmosNand Class
class CmosNand : public ADInterface
{
	public:

  CmosNand(const string& iname);

  ~CmosNand() {}

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
