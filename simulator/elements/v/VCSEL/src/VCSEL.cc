#include "VCSEL.h"

// Static members
const unsigned VCSEL::n_par = 18;

// Element information
ItemInfo VCSEL::einfo =
{
  "vcsel",
  "The VCSEL Model",
  "Houssam Kanj",
  DEFAULT_ADDRESS"category:laser diode,diode>laser"
};

// Parameter information
ParmInfo VCSEL::pinfo[] =
{
  {"etai", "Injection Efficiency", TR_DOUBLE, false},
  {"beta", "Spontaneous Emission Coupling Coefficient", TR_DOUBLE, false},
  {"tn", "Carrier Reconbination Lifetime (nano sec)", TR_DOUBLE, false},
  {"k", "Scalin factor accounting for the output coupling efficiency (Watts)", TR_DOUBLE, false},
  {"g0", "Gain Coefficient (1/sec)", TR_DOUBLE, false},
  {"n0", "Carrier Transparency number", TR_DOUBLE, false},
  {"tp", "Photon Lifetime (sec)", TR_DOUBLE, false},
  {"a0", "1st temperature coefficient of the offset current (A)", TR_DOUBLE, false},
  {"a1", "2nd temperature coefficient of the offset current (A/K)", TR_DOUBLE, false},
  {"a2", "3rd temperature coefficient of the offset current (A/K^2)", TR_DOUBLE, false},
  {"a3", "4nd temperature coefficient of the offset current (A/K^3)", TR_DOUBLE, false},
  {"a4", "5th temperature coefficient of the offset current (A/K^4)", TR_DOUBLE, false},
  {"rho", "Refractive index change", TR_DOUBLE, false},
  {"n", "Refractive Index", TR_DOUBLE, false},
  {"lambda0", "Wavelength(meters)", TR_DOUBLE, false},
  {"rth", "Thermal Impedence (C/mW)", TR_DOUBLE, false},
  {"tth", "Thermal time constant (sec)", TR_DOUBLE, false},
  {"t0", "Ambient Temperature (C)", TR_DOUBLE, false}
};

VCSEL::VCSEL(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(etai = 1);
  paramvalue[1] = &(beta = 1e-6);
  paramvalue[2] = &(tn = 5e-9);
  paramvalue[3] = &(k = 2.6e-8);
  paramvalue[4] = &(g0 = 1.6e4);
  paramvalue[5] = &(n0 = 1.94e7);
  paramvalue[6] = &(tp = 2.28e-12);
  paramvalue[7] = &(a0 = 1.246e-3);
  paramvalue[8] = &(a1 = -2.545e-5);
  paramvalue[9] = &(a2 = 2.908e-7);
  paramvalue[10] = &(a3 = -2.531e-10);
  paramvalue[11] = &(a4 = 1.022e-12);
  paramvalue[12] = &(rho = 2.4e-9);
  paramvalue[13] = &(n = 3.5);
  paramvalue[14] = &(lambda0 = 863e-9);
  paramvalue[15] = &(rth = 2.6e3);
  paramvalue[16] = &(tth = 1e-6);
  paramvalue[17] = &(t0 = 20);

  // Set the number of terminals
  setNumTerms(12);

  // Set flags
  setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(6);
}

void VCSEL::init() throw(string&)
{
  DenseIntVector var(6);
  DenseIntVector novar;
  DenseDoubleVector nodelay;
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  var[3] = 3;
  var[4] = 4;
  var[5] = 5;
  initializeAD(var, var, novar, novar, nodelay);
}

void VCSEL::getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1)); // Local reference terminal
  local_ref_vec.push_back(1); // Local reference index

  term_list.push_back(getTerminal(2));
  term_list.push_back(getTerminal(3)); // Local reference terminal
  local_ref_vec.push_back(3); // Local reference index

  term_list.push_back(getTerminal(4));
  term_list.push_back(getTerminal(5)); // Local reference terminal
  local_ref_vec.push_back(5); // Local reference index

  term_list.push_back(getTerminal(6));
  term_list.push_back(getTerminal(7)); // Local reference terminal
  local_ref_vec.push_back(7); // Local reference index

  term_list.push_back(getTerminal(8));
  term_list.push_back(getTerminal(9)); // Local reference terminal
  local_ref_vec.push_back(9); // Local reference index

  term_list.push_back(getTerminal(10));
  term_list.push_back(getTerminal(11)); // Local reference terminal
  local_ref_vec.push_back(11); // Local reference index
}

void VCSEL::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: terminal current, I
  // x[1]: carrier density, vn
  // x[2]: temperature, td
  // x[3]: photon density, vm
  // x[4]: output power, p0
  // x[5]: lambda, wavelenth
  // x[6]: dI/dt
  // x[7]: dvn/dt
  // x[8]: dtd/dt
  // x[9]: dvm/dt
  // x[10]: dp0/dt
  // x[11]: dlambda/dt

  AD delta = 1e-10;
  AD zn = 1e7;
  AD q = 1.6e-19;
  AD epsilon = 3.4e-23;

  effort[0] = 1.721 + 275*x[0] - 2.439e4*x[0]*x[0] + 1.338e6*x[0]*x[0]*x[0] - 4.154e7*pow(x[0], 4)
	+ 6.683e8*pow(x[0], 5) - 4.296e9*pow(x[0], 6);

  flow[0] = x[0];

	AD Ioff = a0+a1*x[2]+a2*x[2]*x[2]+a3*x[2]*x[2]*x[2]+a4*x[2]*x[2]*x[2]*x[2];

	AD N = zn*x[1] + n0;
	AD dN_dt = zn*x[7] + 0;

	AD S = (x[3]+0.001)/k;
	AD dS_dt = (1/k)*x[9];

	flow[1] = (flow[0]-Ioff - q*g0*(N-n0)*S/(1+epsilon*S) -q*N/tn -q*dN_dt);
	effort[1] = N/zn;

  //this implementation is to pull the thermal code into the model
  //this implementation is in degree celcius
	AD cth = tth/rth;
	AD Ith = (t0 + (flow[0]*effort[0] - (x[3]+delta)*(x[3]+delta))*rth);
	flow[2] = (Ith -x[2] -tth*x[8])/rth;
	effort[2] = x[2];

  effort[3] = x[3];

	flow[3] = tp*k*(-S/tp + beta*N/tn + g0*(N-n0)*S/(1+epsilon*S)- dS_dt );

  flow[4] = x[4];
  effort[4] = (k+0.25e-8)*S;

  flow[5] = x[5];
  effort[5] = lambda0*(1 - (rho/n)*(N-n0));
}

