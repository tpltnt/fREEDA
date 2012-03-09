// This is an mutual inductor model
//
//                 ind
//
//          o-----CCCCC-----o
//          o-----CCCCC-----o
//
//                 ind
//
// by Wei Zheng
// added optional 3rd inductance: NMK (May '05)

#ifndef K_h
#define K_h 1

#include "../../../l/L/src/L.h"

class K : public Element
{
	public:

  K(const string& iname);

  ~K() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);

  // State variable transient analysis
  //virtual void svTran(TimeDomainSV *tdsv);
  //virtual void deriv_svTran(TimeDomainSV *tdsv);

	private:

  unsigned my_row1, my_row2;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;
  static const double factor;

  // Parameter variables
  double coupling, M;
  string l1, l2;
  bool time_d;

  L* ind1;
  L* ind2;
  double inductance1, inductance2;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif

