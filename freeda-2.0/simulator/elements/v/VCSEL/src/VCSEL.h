// This is a Comprehensive VCSEL Model

#ifndef VCSEL_h
#define VCSEL_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class VCSEL : public ADInterface
{
	public:

  VCSEL(const string& iname);

  ~VCSEL() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  virtual void init() throw(string&);

  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
	double etai, beta, tn, k, g0, n0, tp, a0, a1, a2, a3, a4, rho, n, lambda0;
	double rth, tth, t0;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
