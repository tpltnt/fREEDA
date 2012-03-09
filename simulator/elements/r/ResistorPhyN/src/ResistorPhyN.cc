#include "ResistorPhyN.h"

// Static members
const unsigned ResistorPhyN::n_par = 23;

// Element information
ItemInfo ResistorPhyN::einfo =
{
  "resistorphyn",
  "Physical resistor model",
  "NCSU ECE718 student",
  DEFAULT_ADDRESS"category:lumped,electrothermal,semiconductor",
  "2003_05_15"
};

// Parameter information
ParmInfo ResistorPhyN::pinfo[] =
{
	// Resistance parameters
  {"r", "Resistance (ohms)", TR_DOUBLE, false},
  // Could not eliminate parser errors when coefficients were contained in a vector
  {"coeff0", "Constant term of conductance polynomial", TR_DOUBLE, false},
  {"coeff1", "First order coefficient of conductance polynomial", TR_DOUBLE, false},
  {"coeff2", "Second order coefficient of conductance polynomial", TR_DOUBLE, false},
  {"coeff3", "Third order coefficient of conductance polynomial", TR_DOUBLE, false},
  {"coeff4", "Fourth order coefficient of conductance polynomial", TR_DOUBLE, false},
  {"coeff5", "Fifth order coefficient of conductance polynomial", TR_DOUBLE, false},
  {"polyarg", "Polynomial model argument type.  Possible values are true (sum) or false (diff).", TR_BOOLEAN, false},
	// Temperature Effects Parameters
  {"tc1", "Linear temperature coefficient of resistor (1/C)", TR_DOUBLE, false},
  {"tc2", "Quadratic temperature coefficient of resistor (1/C^2)", TR_DOUBLE, false},
  {"tnom", "Parameter measurement temperature (K)", TR_DOUBLE, false},
  {"tdev", "Device operating temperature (K)", TR_DOUBLE, false},
	// Junction Diode Model Parameters
  {"is", "Saturation current (A)", TR_DOUBLE, false},
  {"n", "Emission coefficient", TR_DOUBLE, false},
  {"ibv", "Current magnitude at the reverse breakdown voltage (A)", TR_DOUBLE, false},
  {"bv", "Junction reverse breakdown voltage (V)", TR_DOUBLE, false},
  {"fc", "Coefficient for forward-bias depletion capacitance", TR_DOUBLE, false},
  {"cj0", "Zero-bias junction capacitance (F)", TR_DOUBLE, false},
  {"vj", "Junction built-in potential (V)", TR_DOUBLE, false},
  {"m", "Junction grading coefficient", TR_DOUBLE, false},
  {"tt", "Transit time (s)", TR_DOUBLE, false},
  {"area", "Diode area multiplier", TR_DOUBLE, false},
  {"rs", "Diode series resistance (ohms)", TR_DOUBLE, false}
};


ResistorPhyN::ResistorPhyN(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
	// Set default parameter values
	// r should be infinity, making it 1Tohm
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
  paramvalue[10] = &(tnom = 300.0);
  paramvalue[11] = &(tdev = 300.0);
	// Using default values from SPDiode.cc
  paramvalue[12] = &(is = 1E-14);
  paramvalue[13] = &(n = 1.0);
  paramvalue[14] = &(ibv = 1E-10);
  paramvalue[15] = &(bv = 0.0);
  paramvalue[16] = &(fc = 0.5);
  paramvalue[17] = &(cj0 = 0.0);
  paramvalue[18] = &(vj = 1.0);
  paramvalue[19] = &(m = 0.5);
  paramvalue[20] = &(tt = 0.0);
  paramvalue[21] = &(area = 1.0);
  paramvalue[22] = &(rs = 0.0);

  // Set the number of terminals
  setNumTerms(3);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(2);
}


void ResistorPhyN::init() throw(string&)
{
  // Looks like a lot of temperature and constant
  // calculations that only need to be done once
  DenseIntVector var(2);
  var[0] = 0;
  var[1] = 1;
  DenseIntVector dvar(2);
  dvar[0] = 0;
  dvar[1] = 1;
  initializeAD(var, dvar);
}

void ResistorPhyN::eval(AD * x, AD * effort, AD * flow)
{
  // Model is parameterized for convergence (necessary because of diodes), the state variables are no longer voltages.
  //
  // Diode voltages and currents are calculated from x[0] (diode1) and x[1] (diode 2). Output voltages effort[0]
  // and effort[1] are then calculated from their respective diode voltages (including drop across diode series r).
  // Using the output voltages, the conductance polynomial is evaluated, then voltage-dependent resistance and
  // current through that resistor are calculated.
  //
  // x[0]: State variable used to calculate diode1 voltage and current
  // x[1]: State variable used to calculate diode2 voltage and current
  // x[2]: dx[0]/dt
  // x[3]: dx[1]/dt
  // effort[0]: vterminal0-vterminal2 , flow[0]: current flowing into terminal0
  // effort[1]: vterminal1-vterminal2 , flow[1]: current flowing into terminal1

  // Declaration of double and AD variables used in this routine
  double alpha = eCharge / n / kBoltzman / tdev;
  double v1 = log(5e8 / alpha) / alpha; // normal is .5e9
  double k3 = exp(alpha * v1);
  AD vd1, vd2, id1, id2, dvd1_dx, dvd2_dx, cd1, cd2;
  AD vctrl, rofv, ir, icap1, icap2, condpoly;

  // Calculate temperature adjustment for resistor
  double rtemp = 1.0 + tc1*(tdev-tnom) + tc2*(tdev-tnom)*(tdev-tnom);

  // Calculate diode1 voltage
  if (v1 > x[0])
    vd1 = x[0] + 0.0;
  else
    vd1 = v1 + log(1.0 + alpha*(x[0]-v1))/alpha;

  // Calculate diode2 voltage
  if (v1 > x[1])
    vd2 = x[1] + 0.0;
  else
    vd2 = v1 + log(1.0 + alpha*(x[1]-v1))/alpha;

  // Calculate diode1 current
  if (v1 > x[0])
    id1 = is * (exp(alpha * x[0]) - 1.0);
  else
    id1 = is * k3 * (1.0 + alpha * (x[0] - v1)) - is;
  // subtract the breakdown current
  id1 -= ibv * exp(-alpha * (vd1 + bv));

  // Calculate diode2 current
  if (v1 > x[1])
    id2 = is * (exp(alpha * x[1]) - 1.0);
  else
    id2 = is * k3 * (1.0 + alpha * (x[1] - v1)) - is;
  // subtract the breakdown current
  id2 -= ibv * exp(-alpha * (vd2 + bv));

  // Calculate diode1 junction capacitance
  if (v1 > x[0])
    dvd1_dx = 1.0;
  else
    dvd1_dx = 1.0 / (1.0 + alpha*(x[0]-v1));
  if (isSet(&cj0))
  {
    if (fc * vj > vd1)
      cd1 = cj0 * pow(1.0 - vd1 / vj, -m);
    else
      cd1 = cj0 * pow(1.0 - fc, - m - 1.0) * ((1.0 - fc * (1.0 + m)) + m * vd1 / vj);
  }
  else
    cd1 = zero;
  if (isSet(&tt))
    cd1 += alpha * tt * id1;

  // Calculate diode2 junction capacitance
  if (v1 > x[1])
    dvd2_dx = 1.0;
  else
    dvd2_dx = 1.0 / (1.0 + alpha*(x[1]-v1));
  if (isSet(&cj0))
  {
    if (fc * vj > vd2)
      cd2 = cj0 * pow(1.0 - vd2 / vj, -m);
    else
      cd2 = cj0 * pow(1.0 - fc, - m - 1.0) * ((1.0 - fc * (1.0 + m)) + m * vd2 / vj);
  }
  else
    cd2 = zero;
  if (isSet(&tt))
    cd2 += alpha * tt * id2;

  // Add junction capacitor currents into diode currents, ic = c*(dv/dx)*(dx/dt) = c*dv/dt
  id1 += cd1 * dvd1_dx * x[2];
  id2 += cd2 * dvd2_dx * x[3];

  // Scale diode currents by area multiplier
  id1 *= area;
  id2 *= area;

  // Calculate output voltages
  effort[0] = vd1 + id1*rs;
  effort[1] = vd2 + id2*rs;


  // Calculate conductance polynomial and then resistance
  if(polyarg)
	{
    vctrl = (effort[0] + effort[1]) / 2.0;
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
    // take absolute value of control voltage, assuming symmetry -- may take this out
    if (effort[0] > effort[1])
      vctrl = effort[0] - effort[1];
    else
      vctrl = -1.0*(effort[0] - effort[1]);
    condpoly = coeff0;
    condpoly += (1.0/2.0)*coeff1*vctrl;
    condpoly += (1.0/3.0)*coeff2*vctrl*vctrl;
    condpoly += (1.0/4.0)*coeff3*vctrl*vctrl*vctrl;
    condpoly += (1.0/5.0)*coeff4*vctrl*vctrl*vctrl*vctrl;
    condpoly += (1.0/6.0)*coeff5*vctrl*vctrl*vctrl*vctrl*vctrl;
    rofv = r*rtemp/condpoly;
  }

  // Calculate resistor current
  ir = ( effort[0] - effort[1] ) / rofv;

  // Calculate the output currents
  flow[0]  = ir + id1;
  flow[1]  = -1.0*ir + id2;
}

