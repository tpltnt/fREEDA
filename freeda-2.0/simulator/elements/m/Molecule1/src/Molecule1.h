// Molecule1.h

#ifndef Molecule1_h
#define Molecule1_h 1

#include "../../../../network/CircuitManager.h"
#include "../../../../network/ADInterface.h"

class Molecule1 : public ADInterface
{
	public:
  Molecule1(const string& iname);
  ~Molecule1() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:
  virtual void eval(AD * x, AD * effort, AD * flow);
  // Element information
  static ItemInfo einfo;
  // Number of parameters of this element
  static const unsigned n_par;
  // Parameter variables
  double e0, g1, g2, ef, etha, m, u0, t, accur, area;
  // Parameter information
  static ParmInfo pinfo[];
};

#endif

