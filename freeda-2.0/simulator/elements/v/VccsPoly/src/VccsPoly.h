//      This is a Voltage controlled Current source where:
//                 Ig =  k * polynomial(Vin)
//
//
//                 0                            1
//                  o--------+       +---------o
//                           |       |    |
//           	    +  	       \       |    \         +
//           	  Vin	         / Ri   <Ig>  / Ro    Vout
//           	    -	         \       |    \         -
//                           |       |    |
//                           |       |    |
//           	     o--------+        +---------o
//                2                             3
//
// Author: Jim Hall
//------------------------------------------------------------

#ifndef VccsPoly_h
#define VccsPoly_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class VccsPoly : public ADInterface
{
	public:

  VccsPoly(const string& iname) ;

  ~VccsPoly() {};

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // This is a multi-reference element.

  virtual void eval(AD * x, AD * effort, AD * flow) ;
  virtual void init() throw(string&) ;

	private:

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double k;
  double ri;
  double ro;
  DenseDoubleVector poly_coef;
  double td;
  double gamma;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
