#include "DiodeTun.h"
#include <cstdio>

// Static members
const unsigned DiodeTun::n_par = 18;

// Element information
ItemInfo DiodeTun::einfo =
{
  "diodetun",
  "Tunnel Diode",
  "Stephen Bowyer and Jennifer Huckaby",
  DEFAULT_ADDRESS"category:diode",
  "2003_05_15"
};

// Parameter information
ParmInfo DiodeTun::pinfo[] =
{
  {"js", "Saturation current (A)", TR_DOUBLE, false},
  {"ct0", "Zero-bias depletion capacitance (F)", TR_DOUBLE, false},
  {"fi", "Built-in barrier potential (V)", TR_DOUBLE, false},
  {"gama", "Capacitance power-law parameter", TR_DOUBLE, false},
  {"cd0", "Zero-bias diffusion capacitance (F)", TR_DOUBLE, false},
  {"afac", "Slope factor of diffusion capacitance (1/V)", TR_DOUBLE, false},
  {"r0", "Bias-dependent part of series resistance in forward-bias (Ohms)", TR_DOUBLE, false},
  {"t", "Intrinsic time constant of depletion layer (s)", TR_DOUBLE, false},
  {"area", "Area multiplier", TR_DOUBLE, false},
  {"jv", "Valley current at the valley voltage (A)", TR_DOUBLE, false},
  {"jp", "Peak current at the peak voltage (A)", TR_DOUBLE, false},
  {"vv", "Valley volatge at the valley current (V)", TR_DOUBLE, false},
  {"vpk", "Peak voltage at the peak current (V)", TR_DOUBLE, false},
  {"a2", "Prefactor in the excess current exponent", TR_DOUBLE, false},
  {"mt", "Slope factor of tunnel current (1/V)", TR_DOUBLE, false},
  {"mx", "Slope factor of excess current (1/V)", TR_DOUBLE, false},
  {"mth", "Slope factor of thermal current (1/V)", TR_DOUBLE, false},
  {"temper", "Temperature (K)", TR_DOUBLE, false}
};


DiodeTun::DiodeTun(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(js = 1e-16);
  paramvalue[1] = &(ct0 = zero);
  paramvalue[2] = &(fi = .8);
  paramvalue[3] = &(gama = .5);
  paramvalue[4] = &(cd0 = zero);
  paramvalue[5] = &(afac = 38.696);
  paramvalue[6] = &(r0 = zero);
  paramvalue[7] = &(t = zero);
  paramvalue[8] = &(area = one);
  paramvalue[9] = &(jv = 1e-04);
  paramvalue[10] = &(jp = 1e-03);
  paramvalue[11] = &(vv = 0.5);
  paramvalue[12] = &(vpk = 0.1);
  paramvalue[13] = &(a2 = 30);
  paramvalue[14] = &(mt = -1.);
  paramvalue[15] = &(mx = 1.);
  paramvalue[16] = &(mth = 1.);
  paramvalue[17] = &(temper = 300);

  // Set the number of terminals
  setNumTerms(2);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(1);
}

void DiodeTun::init() throw(string&)
{
	//Define Electronic Charge (C)
  q = 1.60218e-19;

	//Define Boltzmann's Constant (J/K)
  k = 1.3807e-23;

	// Defining constants based on the parameters
  k1 = ct0 / pow(.2, gama);
  k2 = fi * .8;
  k4 = - ct0 * fi / (one - gama);
  k5 = k4 * (.2 * k1 / ct0 - one);
  k6 = cd0 / afac;

  vt = vpk - vpk * log((-mt * vpk) / jp);
  vx = vv + (one / a2) * log(mx / (jv * a2));
  vth = ((k * temper) / q) * log((mth * k * temper) / (js * q));

	//Defining constant based on constraints based on various
	//possible relationships between vt, vx, and vth
  if(vx > vth)
	{
    vxth = vth;
    at2 = jp * exp(one - (vth / vpk));
    ax2 = jv * exp(a2 * (vth - vv));
    ath2 = js * exp((q * vth) / (k * temper));
    b = q / (k * temper);
    ct = -(k * temper) / (q * vpk);
    cx = (a2 * k * temper) / q;
    cth = one;
	}
  else
	{
    vxth = vx;
    at2 = jp * exp(one - (vx / vpk));
    ax2 = jv * exp(a2 * (vx - vv));
    ath2 = js * exp((q * vx) / (k * temper));
    b = a2;
    ct = -one / (a2 * vpk);
    cx = one;
    cth = q / (k * temper * a2);
	}

  if(vt > vxth)
	{
    vt = vxth;
	}

	//Expressions for constants for the three currents
  at1 = jp * exp(one - (vt / vpk));
  ax1 = jv * exp(a2 * (vt - vv));
  ath1 = js * exp((q * vt) / (k * temper));

  DenseIntVector var(1);
  initializeAD(var, var);
}

void DiodeTun::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: state variable
  // x[1]: time derivative of x[0]
  AD vj, dvj_dx, cj, rs, itmp, it, ix, ith;

  //Calculate the voltage of the tunnel diode
  if (vxth > x[0])
    vj = x[0] + zero;
  else
    vj = vxth + (one / b) * log(one + b*(x[0] - vxth));

  if (vt > x[0])
    vj -= vpk * log(one - ((x[0] - vt) / vpk));

  //Calculate the derivative of the voltage
  if (vxth > x[0])
    dvj_dx = one;
  else
    dvj_dx = one / (one + b*(x[0]-vxth));

  if (vt > x[0])
    dvj_dx = one / (one - (one/vpk) * (x[0] - vt));

  // Calculate the junction capacitance (experimental)
  // Use a modified function with continous first and second deriv.
  cj = zero;
  if (isSet(&ct0))
	{
		AD exp1 = exp(10. * (vj-k2));
		const double k14 = ct0 * gama / fi;
		const double k15 = k4 * (gama - one) / fi;
    if (vj < 0.0)
      cj = ct0 / pow(one - vj/fi, gama);
    else
      cj = (ct0 + vj * (k14 + k15 * vj)) / (one + exp1) + k1 * exp1 / (one + exp1);
	}

  // add depletion capacitance
  if (isSet(&cd0))
    cj += cd0 * exp(afac * vj);

  // Now calculate the current through the capacitor.
  // Using the chain rule:
  //
  // dq/dt = cj(vj) * dvj/dx * dx/dt
  //
  // x[1] is dx/dt
  flow[0] = cj * dvj_dx * x[1];

  // Now use the state variable again to calculate the total
  // current.  This way, we save some exp() calls.  The total
  // current is the current through the capacitor plus the tunnel
  // current plus the excess current plus the thermal current of
  // the diode.

  // Tunnel current
  if (vxth > x[0])
    it = jp * (exp(one - (x[0] / vpk)));
  else
    it = at2 * pow(one + b * (x[0] - vxth), ct);

  if (vt > x[0])
    it = at1 * (one - (one / vpk) * (x[0] - vt));

  itmp = it * (vj / vpk);
  //Add to the total current
  flow[0] += itmp;

  // Excess Current
  if (vxth > x[0])
    ix = jv * exp(a2 * (x[0] - vv));
  else
    ix = ax2 * pow(one + b * (x[0] - vxth), cx);

  if (vt > x[0])
    ix = ax1 * pow(one - (one / vpk) * (x[0] - vt), -a2 * vpk);

  // Add to the total current
  flow[0] += ix;

  // Thermal Current
  if (vxth > x[0])
    ith = js * exp(q * x[0] / (k * temper)) - js;
  else
    ith = ath2 * pow(one + b * (x[0] - vxth), cth) - js;

  if (vt > x[0])
    ith = ath1 * pow(one - (one / vpk) * (x[0] - vt), ((-q * vpk) / (k * temper))) - js;

  //Add to the total current
  flow[0] += ith;

	/***This resistance equation is based on the normal microwave diode.***/
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

  //Total voltage is diode voltage plus
  //the voltage across series resistor
  effort[0] = vj + flow[0]*rs;

  // scale the current according to area. All the calculations were made
  // for a unit area diode.
  flow[0] *= area;
}
