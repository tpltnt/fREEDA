// Open circuit element
//
//
//   0 o-             -o 1
//
//
//
//

#ifndef Open_h
#define Open_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "../../../../analysis/FreqDomainSV.h"
#include "../../../../network/ADInterface.h"

class Open : public ADInterface
{
public:

  Open(const string& iname);

  ~Open() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  virtual void eval(AD * x, AD * effort, AD * flow);
  // State variable HB.
  virtual void svHB(FreqDomainSV* fdsv);
  virtual void deriv_svHB(FreqDomainSV* fdsv);

  // State variable transient analysis
  virtual void svTran(TimeDomainSV* tdsv);
  virtual void deriv_svTran(TimeDomainSV* tdsv);

private:

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter information
  // static ParmInfo pinfo[];

};

#endif
