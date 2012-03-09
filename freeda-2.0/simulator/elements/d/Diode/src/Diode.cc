#include "Diode.h"

// Static members
const unsigned Diode::n_par = 20;

// Element information
ItemInfo Diode::einfo =
{
  "diode",
  "Microwave diode",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:diode",
  "2000_06_15"
};

// Parameter information
ParmInfo Diode::pinfo[] =
{
  {"js", "Saturation current (A)", TR_DOUBLE, false},
  {"alfa", "Slope factor of conduction current (1/V)", TR_DOUBLE, false},
  {"jb", "Breakdown saturation current (A)", TR_DOUBLE, false},
  {"vb", "Breakdown voltage (V)", TR_DOUBLE, false},
  {"e", "Power-law parameter of breakdown current", TR_DOUBLE, false},
  {"ct0", "Zero-bias depletion capacitance (F)", TR_DOUBLE, false},
  {"fi", "Built-in barrier potential (V)", TR_DOUBLE, false},
  {"gama", "Capacitance power-law parameter", TR_DOUBLE, false},
  {"cd0", "Zero-bias diffusion capacitance (F)", TR_DOUBLE, false},
  {"afac", "Slope factor of diffusion capacitance (1/V)", TR_DOUBLE, false},
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


Diode::Diode(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
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
  paramvalue[8] = &(cd0 = zero);
  paramvalue[9] = &(afac = 38.696);
  paramvalue[10] = &(r0 = zero);
  paramvalue[11] = &(t = zero);
  paramvalue[12] = &(area = one);
  paramvalue[13] = &(imax = zero);
  paramvalue[14] = &(eg = .8);
  paramvalue[15] = &(m = .5);
  paramvalue[16] = &(aro = zero);
  paramvalue[17] = &(bro = zero);
  paramvalue[18] = &(afag = one);
  paramvalue[19] = &(xti = 2.);

  // Set the number of terminals
  setNumTerms(2);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(1);
}

void Diode::init() throw(string&)
{
  v1 = log(5e8 / alfa) / alfa; // normal is .5e9
  k1 = ct0 / pow(.2, gama);
  k2 = fi * .8;
  k3 = exp(alfa * v1);
  k4 = - ct0 * fi / (one - gama);
  k5 = k4 * (.2 * k1 / ct0 - one);
  k6 = cd0 / afac;

  // initialize automatic differentiation
  DenseIntVector var(1);
  initializeAD(var, var);
}

void Diode::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: state variable
  // x[1]: time derivative of x[0]

  AD vj, dvj_dx, cj, rs, itmp;

  if (v1 > x[0])
  {
    vj = x[0] + zero;
    dvj_dx = one;
  }
  else
  {
    vj = v1 + log(one + alfa*(x[0] - v1))/alfa;
    dvj_dx = one / (one + alfa*(x[0] - v1));
  }

  // Calculate the junction capacitance (experimental)
  // Use a modified function with continous first and second deriv.
  cj = zero;
  if (isSet(&ct0))
  {
    AD exp1 = exp(10. * (vj-k2));
    const double k14 = ct0 * gama / fi;
    const double k15 = k4 * (gama - one) / fi;
    if (vj < zero)
      cj = ct0 / pow(one - vj/fi, gama);
    else
      cj = (ct0 + vj * (k14 + k15 * vj)) / (one + exp1) +
           k1 * exp1 / (one + exp1);
  }

  // add depletion capacitance
  if (isSet(&cd0))
    cj += cd0 * exp(afac * vj);

  // Now calculate the current through the capacitor.
  // Using the chain rule:
  // dq/dt = cj(vj) * dvj/dx * dx/dt
  // x[1] is dx/dt
  flow[0] = cj * dvj_dx * x[1];

  // Now use the state variable again to calculate the total
  // current.  This way, we save some exp() calls.  The total
  // current is the current through the capacitor plus the ideal
  // diode current.
  if (v1 > x[0])
    itmp = js * (exp(alfa * x[0]) - one);
  else
    itmp = js * k3 * (one + alfa * (x[0] - v1)) - js;
  flow[0] += itmp;

  // subtract the breakdown current
  if (vj - vb > one)
    itmp = zero;
  else
    itmp = jb * pow(one + vb - vj, e);
  flow[0] -= itmp;

  // Calculate Rs
  if (cj != zero)
  {
    if (t/cj > r0)
      rs = zero;
    else
      rs = r0 - t/cj;
  }
  else
    rs = r0;

  effort[0] = vj + flow[0] * rs;

  // scale the current according to area. All the calculations were made
  // for a unit area diode.
  flow[0] *= area;
}
