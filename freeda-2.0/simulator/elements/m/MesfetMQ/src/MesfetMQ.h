//		  Contains routines for analysis of compact software's
//                Materka-Kacprzak MESFET model
//
//          Drain 2
//                  o
//                  |
//                  |
//              |---+
//              |
// Gate 1 o-----|       ( -----o Bulk 4    in a future release )
//              |
//              |---+
//                  |
//                  |
//                  o
//         Source 3
//
//
//	  Author:
//		  Senthil Velu
//

#ifndef MesfetMQ_h
#define MesfetMQ_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class MesfetMQ : public ADInterface
{
	public:

  MesfetMQ(const string& iname);

  ~MesfetMQ() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Some constants
  double k01, k02, k03;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double idss, vp0, gama, e, ke, sl, kg, ss, t, ig0, afag,
	ib0, afab, vbc, r10, kr, c10, k1, c1s, cf0, kf, area;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
