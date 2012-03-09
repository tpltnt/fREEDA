// This is an capacitor model
//
//                cap
//
//                 ||
//         o--+----||----+--o
//            |    ||    |
//            |          |
//            +--\/\/\/--+
//               int_g
//
// by Sonali Luniya

#ifndef CapacitorLC_h
#define CapacitorLC_h 1

#include "../../../../network/CircuitManager.h"
#include "../../../../network/ADInterface.h"

class CapacitorLC : public ADInterface
{
	public:

  CapacitorLC(const string& iname);

  ~CapacitorLC() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  double EPSILON_0;
  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double L, W, D, DELTA, GAMMA, DTIME, VC, EPSILON_PL;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif

