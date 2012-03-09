// This is an electro-thermal interconect model
//
//           +++++++++++++++++++
//           + InterconnectRTSH+
//           +  +----+         +
//      o----+--+    +    +----+---o
//           +       +----+    +
//           +                 +
//           +  *           *  +
//           +  |           |  +
//           ++ | +++++++++ | ++
//              |           |
//              o           o
//             tref       tout
//
// by  Kai Li, Theodore R Harris

#ifndef InterconnectRTSH_h
#define InterconnectRTSH_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"

class InterconnectRTSH : public ADInterface
{
	public:

  InterconnectRTSH(const string& iname);

  ~InterconnectRTSH() {}

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
  double l, w, tm, rho;
  string metal;
  double t, tnom, tc;
  bool pdr;
  //double idf;

  // Parameter information
  static ParmInfo pinfo[];

  // The calculated resistance
	double r;
};

#endif
