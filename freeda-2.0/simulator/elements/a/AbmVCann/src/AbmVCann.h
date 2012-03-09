// This may look like C code, but it is really -*- C++ -*-
//
// This element is aided in multisine simulation
//
// the current Iout = g*vin[idx]/(1+(g/L*abs(vin[idx]))^s)^(1/s)/Rout
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

#ifndef AbmVCann_h
#define AbmVCann_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class AbmVCann : public ADInterface
{
public:

  AbmVCann(const string& iname);

  ~AbmVCann() {}

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
  double s;
  double rout;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
