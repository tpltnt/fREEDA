// General Purpose Operational Amplifier
#include "OpAmpT.h"

// Define the number of parameters
const unsigned OpAmpT :: n_par = 11;

// Element Information
ItemInfo OpAmpT::einfo =
{
  "opampt",
  "LM741 OpAmp (Electro-thermal)",
  "Gregory S Parsons",
  DEFAULT_ADDRESS"category:behavioral",
  "2004_05_10"
};

// Parameter Information
ParmInfo OpAmpT::pinfo[] =
{
  {"gain","opamp Gain (A/V)", TR_DOUBLE, false},
  {"rin","Linear Input Resistance (ohm)", TR_DOUBLE, false},
  {"rint","Input Thermal Effect (ohm/degree)", TR_DOUBLE, false},
  {"rout","Linear Output Resistance (ohm)", TR_DOUBLE, false},
  {"routt","Output Thermal Effect (ohm/degree)", TR_DOUBLE, false},
  {"psr","PS Resistance Vdd to gnd (ohm)", TR_DOUBLE, false},
  {"psrt","PS Resistance Thermal Effect (ohm/degree)", TR_DOUBLE, false},
  {"cin","Input Capacitance (F)", TR_DOUBLE, false},
  {"cout","Output Capacitance (F)", TR_DOUBLE, false},
  {"nomt","Nominal Temperature (K)", TR_DOUBLE, false},
  {"pdr","Thermal Dependent Device (bool)", TR_BOOLEAN, false}
};

// Constructor
OpAmpT::OpAmpT(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(gain = 1);
  paramvalue[1] = &(rin = 2000000);
  paramvalue[2] = &(rint = 12500);
  paramvalue[3] = &(rout = 80);
  paramvalue[4] = &(routt = -0.4);
  paramvalue[5] = &(psr = 11000);
  paramvalue[6] = &(psrt = -27);
  paramvalue[7] = &(cin = 0.00000005);
  paramvalue[8] = &(cout = 0.000005);
  paramvalue[9] = &(nomt = 25);
  paramvalue[10] = &(pdr = false);
}

// Initialization Function
void OpAmpT::init() throw(string&)
{
  if(pdr)
  {
    setNumTerms(6);
    setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
    setNumberOfStates(5);
    DenseIntVector var(5);
    var[1] = 1;
    var[2] = 2;
    var[3] = 3;
    var[4] = 4;
    initializeAD(var, var);
  }
  else
  {
    setNumTerms(5);
    setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
    setNumberOfStates(4);
    DenseIntVector var(4);
    var[1] = 1;
    var[2] = 2;
    var[3] = 3;
    initializeAD(var, var);
  }
}

// Evaluate Function
void OpAmpT::eval(AD * x, AD * effort, AD * flow)
{
  //x[0]  : Vdd   x[5]  : dVdd
  //x[1]  : V+    x[6]  : dV+
  //x[2]  : V-    x[7]  : dV-
  //x[3]  : Vout  x[8]  : dVout
  //x[4]  : Temp  x[9]  : dTemp
  //effort[0] : Vdd   flow[0] : Idd
  //effort[1] : V+    flow[1] : I+
  //effort[2] : V-    flow[2] : I-
  //effort[3] : Vout  flow[3] : Iout
  //effort[4] : Temp  flow[4] : Temp_Flux

  AD Iout, Iin1, Iin2, Idd, IoutA, IoutB, IoutC, Tout;

  if(pdr)
  {
    // if Thermal
    // Region A: V- < V+ - Vdiff
    IoutA = (x[3] - x[0]) / (rout + (routt * (x[4]-nomt)));

    // Region C: V+ > V- + Vdiff
    IoutC = x[3] / (rout + (routt * (x[4]-nomt)));
  }
  else
  {
    // if NOT Thermal
    // Region A: V- < V+ - Vdiff
    IoutA = (x[3] - x[0]) / rout;
    // Region C: V+ > V- + Vdiff
    IoutC = x[3] / rout;
  }

  // Region B: strobe
  IoutB = (x[2] - x[1]) * gain;

  // Find maximum of Region A and Region B
  if (IoutB > IoutA)
    Iout = IoutB;
  else
    Iout = IoutA;

  // Find Minimum of Prev-Result and Region C
  if (IoutC > Iout)
    Iout = Iout;
  else
    Iout = IoutC;

  // Determine Current Consumed
  if (Iout > 0.0)
    Idd = 0.0;
  else
    Idd = -Iout;
  // result positive if sink, neg if source

  if(pdr)
  {
    // If Thermal
    if (Iout < 0.0)
      Tout = Iout*x[3];
    else
      Tout = -Iout * (x[0] - x[3]);

    if (Tout < 0.0)
      Tout = -Tout; // Sanity check

    // Generate and Write Temperate Elements
    Idd += x[0] / (psr + (psrt * (x[4]-nomt)));
    // Include PSR Current on T
    Tout += x[0] / (psr + (psrt * (x[4]-nomt)));
    effort[4] = x[4];
    flow[4] = -Tout;
  }
  else
  {
    Idd += x[0] / psr; // include PSR Current
  }

  effort[0] = x[0];
  effort[1] = x[1];
  effort[2] = x[2];
  effort[3] = x[3];

  flow[0] = Idd;
  flow[1] = (x[6] * cin); // input capacitance
  flow[2] = (x[7] * cin);
  flow[3] = Iout + (x[8] * cout); // Output Capacitance
}
