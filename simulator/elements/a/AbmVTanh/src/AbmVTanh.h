// This may look like C code, but it is really -*- C++ -*-
//
// This element is aided in multisine simulation
//
// the current Iout = L*tanh(vin[idx]*g/L)/Rout
//
// Otherwise, the current is adjusted by a polynomial factor.
//
//
//           1(IN)  o------o     ------o  2(OUT)
//           		   +    |      +
//           		 Vin   ( )     Vout
//           		   -    |      -
//           	    o------------------o
//                        3(Common)
//
//  Author:
//  Minsheng Li

#ifndef AbmVTanh_h
#define AbmVTanh_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class AbmVTanh : public ADInterface
{
public:

  AbmVTanh(const string& iname);

  ~AbmVTanh() {}

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
  double g;
  double l;
  double rout;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
