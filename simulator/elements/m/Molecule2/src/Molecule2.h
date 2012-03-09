// Molecule2.h

#ifndef Molecule2_h
#define Molecule2_h 1

#include "../../../../network/CircuitManager.h"
#include "../../../../network/ADInterface.h"

class Molecule2 : public ADInterface
{
	public:
  Molecule2(const string& iname);
  ~Molecule2() {}
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
  double e0, g1, g2, ef, etha, u0, t, accur, npoints, width, area;
  // Parameter information
  static ParmInfo pinfo[];
};

#endif

