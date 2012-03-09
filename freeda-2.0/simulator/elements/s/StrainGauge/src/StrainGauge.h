// Contains routines for analysis of a piezoelectric strain // gauge
//                       	+---------+
//                      	|         |
//             0  o------+         +-------o 1
//                       |         |
//     Mechanical        +----+----+         Electrical
//                            |
//                o-----------+------------o 2
// Author:
// Jonathan Cantor and Richard McMunn
//

#ifndef StrainGauge_h
#define StrainGauge_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class StrainGauge : public ADInterface
{
public:

  StrainGauge(const string& iname);

  ~StrainGauge() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  virtual void init() throw(string&);

private:

  virtual void eval(AD * x, AD * effort, AD * flow);
  static ItemInfo einfo;
  static const unsigned n_par;
  double h1,h2,h3,C,e,eps;
  static ParmInfo pinfo[];
};

#endif
