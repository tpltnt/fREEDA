// CPW model. It uses TLinp4 to model the transmission line.
// It is considered as a 4 terminal device with 2 local reference nodes.
//
//               +--------+
//       1 o-----+        +------o 3
//               |        |
//               |        |
//       2 o-----+        +------o 4
//         |     +--------+      |
//         ^   Local references  ^
//
// by Carlos E. Christoffersen using equations from the Collin book.
//
// This is an approximate model. We need to set up (later) some
// numerical integration algorithm to improve accuracy and generality.

#ifndef CPW_h
#define CPW_h 1

#include "../../../t/Tlinp4/src/Tlinp4.h"

class CPW : public Element
{
	public:

  CPW(const string& iname);

  ~CPW() {}

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
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
	TerminalVector& term_list);

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);

	private:

  // Underlying transmission line
  Tlinp4 *tline;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double w;
  double s;
  double length;
  double er;
  double t;
  double tand;
  double sigma;
  int nsect;
  double fopt;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif

