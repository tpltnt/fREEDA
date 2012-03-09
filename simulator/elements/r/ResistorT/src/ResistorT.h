// This is an electro-thermal resistor model
//
//           ++++++++++++++
//           +   res      +
//      o--- ----/\/\/\------o
//           +            +
//           +  *      *  +
//           +  |      |  +
//           ++ | ++++ | ++
//              |      |
//              o      o
//             tref   tout
//
// by  Houssam S. Kanj

#ifndef ResistorT_h
#define ResistorT_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"

class ResistorT : public ADInterface
{
	public:

  ResistorT(const string& iname);

  ~ResistorT() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // Get a vector with the indexes of the local reference nodes.
  // Get terminal pointers in term_list
  // ordered by local reference node:
  //
  // t0 t1 t2 t3
  //
  //     ^     ^
  //     LRN1  LRN2
  //
  // local_ref_vec contains: {1, 3}
  //
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list);

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double r0, rsh, l, w, narrow, defw, t, tnom, tc1, tc2;
  bool pdr;

  // Parameter information
  static ParmInfo pinfo[];

  // The calculated resistor value
	double r;
};

#endif
