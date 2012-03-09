// CMOS NOR Gate with Junction temperature
//
//
//                                                         Vdd 1
//                                                             o
//                                                              |
//                                                        |---+
//                                                        |
//                   |-------------------------O|
//                   |                                     |
//                   |                                     |---+
//                   |                                           |                             Thermal nodes
//                   |                                           |
//                   |                                     |---+                            Junction Temp
//                   |                                     |                                            6
//                   |                            |----O|                                           o
//                   |                            |        |                                        |
//                   |                            |        |----+                                   |
//                   |         Input2 3 o--|               |                                        |
//                   |                            |               |                                 |
//                   |            *----------|-----------*---o 4 Output                             |
//                   |             |              |               |                                 |
//                   |     |-----+             |         |---+                                      o
//                   |     |                      |         |                                       7
// Input1 2 o-|----|                      |------|                                             Thermal Ref
//                         |                               |
//                         |---+                         |---+
//                               |                               |
//                               |                               |
//                               ------------------*-----
//                                                        |
//                                                        |
//                                                        o
//                                                    GND 5
//
//
//
//	  Authors:
//		 Tony Mulder, Travis Lentz
//

#ifndef CmosNorT_h
#define CmosNorT_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

// The CmosNorT Class
class CmosNorT : public ADInterface
{
  public:

  CmosNorT(const string& iname);

  ~CmosNorT() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list);

  private:

  // Evaluation function
  virtual void eval(AD * x, AD * effort, AD * flow);

  // Parameter Variables
  double vtn, vtp, un, up, en, ep, tox, ln, wn, lp, wp, td, ts, tnom, zt;
  double c1, c2, c3, c4, c5, freq, lk;

  // Model Variables
  double betan, betap, power;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  bool thermal;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
