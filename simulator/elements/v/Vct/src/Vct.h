//  Voltage-to-Current Transducer
//  This is a current limited Voltage-to-Current Transducer
//  This is a special element created to aid in modeling the ERA6 amplifier
//
// If the polynomial coefficients are not specified,
// the current Iout = kf*tanh(beta*(vin[idx] +gamma)) + kf
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
//  Mark Summers
//  C++ version and refinement by Carlos E. Christoffersen

#ifndef Vct_h
#define Vct_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class Vct : public ADInterface
{
	public:

  Vct(const string& iname);

  ~Vct() {}

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
  double gamma;
  double beta;
  double kf;
  DenseDoubleVector poly_coef;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
