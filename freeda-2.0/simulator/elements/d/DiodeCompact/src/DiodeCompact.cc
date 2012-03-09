#include "DiodeCompact.h"

// Static members
const unsigned DiodeCompact::n_par = 18;

// Element information
ItemInfo DiodeCompact::einfo =
{
  "diodecompact",
  "Spice diode model (conserves charge)",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:diode",
  "2000_05_15"
};

// Parameter information
ParmInfo DiodeCompact::pinfo[] =
{
  {"js", "Saturation current (A)", TR_DOUBLE, false},
  {"alfa", "Slope factor of conduction current (1/V)", TR_DOUBLE, false},
  {"jb", "Breakdown saturation current (A)", TR_DOUBLE, false},
  {"vb", "Breakdown voltage (V)", TR_DOUBLE, false},
  {"e", "Power-law parameter of breakdown current", TR_DOUBLE, false},
  {"ct0", "Zero-bias depletion capacitance (F)", TR_DOUBLE, false},
  {"fi", "Built-in barrier potential (V)", TR_DOUBLE, false},
  {"gama", "Capacitance power-law parameter", TR_DOUBLE, false},
  {"r0", "Bias-dependent part of series resistance in forward-bias (Ohms)",
	TR_DOUBLE, false},
  {"t", "Intrinsic time constant of depletion layer (s)", TR_DOUBLE, false},
  {"area", "Area multiplier", TR_DOUBLE, false},
  {"imax", "Maximum forward and reverse current (A)", TR_DOUBLE, false},
  {"eg", "Barrier height at 0 K (eV)", TR_DOUBLE, false},
  {"m", "Grading coefficient", TR_DOUBLE, false},
  {"aro", "r0 linear temperature coefficient (1/K)", TR_DOUBLE, false},
  {"bro", "r0 quadratic temperature coefficient (1/K^2)", TR_DOUBLE, false},
  {"afag", "Temperature-related coefficient", TR_DOUBLE, false},
  {"xti", "Js temperature exponent", TR_DOUBLE, false}
};


DiodeCompact::DiodeCompact(const string& iname)
: ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(js = 1e-16);
  paramvalue[1] = &(alfa = 38.696);
  paramvalue[2] = &(jb = 1e-5);
  paramvalue[3] = &(vb = -1e20);
  paramvalue[4] = &(e = 10.);
  paramvalue[5] = &(ct0 = zero);
  paramvalue[6] = &(fi = .8);
  paramvalue[7] = &(gama = .5);
  paramvalue[8] = &(r0 = zero);
  paramvalue[9] = &(t = zero);
  paramvalue[10] = &(area = one);
  paramvalue[11] = &(imax = zero);
  paramvalue[12] = &(eg = .8);
  paramvalue[13] = &(m = .5);
  paramvalue[14] = &(aro = zero);
  paramvalue[15] = &(bro = zero);
  paramvalue[16] = &(afag = one);
  paramvalue[17] = &(xti = 2.);

  // Set the number of terminals
  setNumTerms(2);
  // Set number of states
  setNumberOfStates(1);
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void DiodeCompact::init() throw(string&)
{
  DenseIntVector var(1);
  initializeAD(var, var);
}

void DiodeCompact::eval(AD * x, AD * effort, AD * flow)
{
  double v1 = log(5e8 / alfa) / alfa; // normal is .5e9
  double k1 = ct0 / pow(.2, gama);
  double k2 = fi * .8;
  double k3 = exp(alfa * v1);

  // x[0]: state variable
  // x[1]: time derivative of x[0]

  AD vj, qj, rs, itmp, cj, idynamic;

  if (v1 > x[0])
    vj = x[0] + zero;
  else
    vj = v1 + log(one + alfa*(x[0]-v1))/alfa;

  // Calculate charge across the diode
  if (isSet(&ct0))
  {
    if (vj > k2)
      qj = ct0 * fi * (one - pow(.2, one - gama)) / (one - gama)
           + k1 * fi * (vj / fi - .8);
    else
      qj = ct0 * fi * (one - pow(one - vj / fi, one - gama)) / (one - gama);
  }

  // Now use the state variable again to calculate the total
	// current.  This way, we save some exp() calls.  The total
	// current is the current through the capacitor plus the ideal
	// diode current.
	if (v1 > x[0])
	  itmp = js * (exp(alfa * x[0]) - one);
	else
	  itmp = js * k3 * (one + alfa * (x[0] - v1)) - js;
  flow[0] = itmp;

  // subtract the breakdown current
  if (vj - vb - one > zero)
    itmp = zero;
  else
    itmp = jb * pow(one + vb - vj, e);
  flow[0] -= itmp;

  // Calculate capacitance as usual
  if (isSet(&ct0))
  {
    if (vj > k2)
      cj = k1;
    else
      cj = ct0 / pow(one - vj/fi, gama);

    if (t/cj > r0)
      rs = zero;
    else
      rs = r0 - t/cj;
  }
  else
    rs = r0;

  // Add dynamic current contribution, dqj/dt
  // dqj/dt = dqj/dx * dx/dt
  idynamic = qj.fastAccessDx(0) * x[1];
  flow[0] += idynamic;

  effort[0] = vj + flow[0] * rs;
  // scale the current according to area.
  flow[0] *= area;
}
