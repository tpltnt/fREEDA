#include "ResistorPhyPoly.h"

// Static members
const unsigned ResistorPhyPoly::n_par = 15;

// Element information
ItemInfo ResistorPhyPoly::einfo =
{
  "resistorphypoly",
  "Physical resistor model",
  "NCSU ECE718 student",
  DEFAULT_ADDRESS"category:lumped,electrothermal,semiconductor",
  "2003_05_15"
};

// Parameter information
ParmInfo ResistorPhyPoly::pinfo[] =
{
  {"r", "Resistance (ohms)", TR_DOUBLE, false},
  // Could not eliminate parser errors when coefficients were contained in a vector
  {"coeff0", "Constant term of conductance polynomial", TR_DOUBLE, false},
  {"coeff1", "First order coefficient of conductance polynomial", TR_DOUBLE, false},
  {"coeff2", "Second order coefficient of conductance polynomial", TR_DOUBLE, false},
  {"coeff3", "Third order coefficient of conductance polynomial", TR_DOUBLE, false},
  {"coeff4", "Fourth order coefficient of conductance polynomial", TR_DOUBLE, false},
  {"coeff5", "Fifth order coefficient of conductance polynomial", TR_DOUBLE, false},
  {"polyarg", "Polynomial model argument type.  Possible values are true (sum) or false (diff).", TR_BOOLEAN, false},
  {"tc1", "Linear temperature coefficient of resistor (1/C)", TR_DOUBLE, false},
  {"tc2", "Quadratic temperature coefficient of resistor (1/C^2)", TR_DOUBLE, false},
  {"tc1c", "Linear temperature coefficient of linear capacitor (1/C)", TR_DOUBLE, false},
  {"tc2c", "Quadratic temperature coefficient of linear capacitor (1/C^2)", TR_DOUBLE, false},
  {"tnom", "Parameter measurement temperature (K)", TR_DOUBLE, false},
  {"tdev", "Device operating temperature (K)", TR_DOUBLE, false},
  {"c", "Linear capacitance (F)", TR_DOUBLE, false}

};

ResistorPhyPoly::ResistorPhyPoly(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
	// Set default parameter values
  paramvalue[0]  = &(r = 1.0E+9);
  paramvalue[1]  = &(coeff0 = 1.0);
  paramvalue[2]  = &(coeff1 = 0.0);
  paramvalue[3]  = &(coeff2 = 0.0);
  paramvalue[4]  = &(coeff3 = 0.0);
  paramvalue[5]  = &(coeff4 = 0.0);
  paramvalue[6]  = &(coeff5 = 0.0);
  paramvalue[7]  = &(polyarg = true);
  paramvalue[8]  = &(tc1 = 0.0);
  paramvalue[9]  = &(tc2 = 0.0);
  paramvalue[10]  = &(tc1c = 0.0);
  paramvalue[11]  = &(tc2c = 0.0);
  paramvalue[12] = &(tnom = 300.0);
  paramvalue[13] = &(tdev = 300.0);
  paramvalue[14] = &(c = 0.0);

  // Set the number of terminals
  setNumTerms(3);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(2);
}


void ResistorPhyPoly::init() throw(string&)
{
  DenseIntVector var(2);
  var[0] = 0;
  var[1] = 1;
  DenseIntVector dvar(2);
  dvar[0] = 0;
  dvar[1] = 1;
  initializeAD(var, dvar);
}

void ResistorPhyPoly::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: v1 = vterminal0-vterminal2
  // x[1]: v2 = vterminal1-vterminal2
  // x[2]: dv1/dt
  // x[3]: dv2/dt
  // effort[0]: v1 , flow[0]: current flowing into terminal0
  // effort[1]: v2 , flow[1]: current flowing into terminal1

  // Assign known output voltages
  effort[0] = x[0];
  effort[1] = x[1];

  // Declaration of double variables used in this routine
  AD vctrl, rofv, ir, icap1, icap2, condpoly;

  // Calculate temperature adjustment for resistor and capacitors
  double rtemp = 1.0 + tc1*(tdev-tnom) + tc2*(tdev-tnom)*(tdev-tnom);
  double ctemp = 1.0 + tc1c*(tdev-tnom) + tc2c*(tdev-tnom)*(tdev-tnom);

  // Calculate conductance polynomial and then resistance.
  // Coefficients should be chosen carefully to avoid negative
  // resistance for all anticipated values of vctrl.
  if(polyarg)
	{
    vctrl = (x[0] + x[1]) / 2.0;
    condpoly = coeff0;
    condpoly += coeff1*vctrl;
    condpoly += coeff2*vctrl*vctrl;
    condpoly += coeff3*vctrl*vctrl*vctrl;
    condpoly += coeff4*vctrl*vctrl*vctrl*vctrl;
    condpoly += coeff5*vctrl*vctrl*vctrl*vctrl*vctrl;
    rofv = r*rtemp/condpoly;
  }
  else
	{
    vctrl = x[0] - x[1];
    condpoly = coeff0;
    condpoly += (1.0/2.0)*coeff1*vctrl;
    condpoly += (1.0/3.0)*coeff2*vctrl*vctrl;
    condpoly += (1.0/4.0)*coeff3*vctrl*vctrl*vctrl;
    condpoly += (1.0/5.0)*coeff4*vctrl*vctrl*vctrl*vctrl;
    condpoly += (1.0/6.0)*coeff5*vctrl*vctrl*vctrl*vctrl*vctrl;
    rofv = r*rtemp/condpoly;
  }

  // Calculate resistor current
  ir = ( x[0] - x[1] ) / rofv;

  // calculate current into cap between terminals 0 and 2, i = c*dv/dt
  icap1 = c*ctemp*x[2];

  // calculate current into cap between terminals 1 and 2, i = c*dv/dt
  icap2 = c*ctemp*x[3];

  // Calculate the output currents
  flow[0]  = ir + icap1;
  flow[1]  = -1.0*ir +icap2;
}

