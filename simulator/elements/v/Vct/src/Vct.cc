#include "Vct.h"

// Static members
const unsigned Vct::n_par = 4;

// Element information
ItemInfo Vct::einfo =
{
  "vct",
  "Voltage-to-Current Transducer",
  "Mark Summers, Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:source",
  "2000_07_20"
};

// Parameter information
ParmInfo Vct::pinfo[] =
{
  {"gamma", "Element parameter", TR_DOUBLE, true},
  {"beta", "Element parameter", TR_DOUBLE, true},
  {"kf", "Element parameter", TR_DOUBLE, true},
  {"poly", "Element parameter", TR_DOUBLE_VECTOR, false}
};


Vct::Vct(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &gamma;
  paramvalue[1] = &beta;
  paramvalue[2] = &kf;
  paramvalue[3] = &poly_coef;

  // Set the number of terminals
  setNumTerms(3);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(2);
}

void Vct::init() throw(string&)
{
  DenseIntVector var(2);
  var[1] = 1;
  DenseIntVector novar;
  DenseDoubleVector nodelay;
  initializeAD(var, novar, novar, novar, nodelay);
}

void Vct::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: vin
  // x[1]: vout
  //  Assign voltages directly.
  effort[0] = x[0];
  effort[1] = x[1];

  if (!isSet(3))
	{
    // The polynomial coefficients are not specified.
    // Use small signal model.
    // Set input current to zero.
    flow[0] = zero;

    // Equation for output current
    flow[1] = kf*tanh(beta*(x[0] +gamma)) + kf;
  }
  else
	{
    // Use polynomial
    AD i_sum_out;
    // Set input current to zero.
    flow[0] = zero;

    // Calculate p0 + vin * (p1 + vin * (p2 + vin * (p3)))
    i_sum_out = poly_coef[poly_coef.length()-1];
    for (unsigned i = poly_coef.length()-1; i > 0; i--)
		{
      i_sum_out *= x[0];
      i_sum_out += poly_coef[i-1];
    }
    // Adjust Current for dependency on output voltage
    flow[1] = i_sum_out * kf * tanh(beta * x[1] - gamma);
  }
}

