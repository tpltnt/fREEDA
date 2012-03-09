#include "DHLD.h"

// Static members
const unsigned DHLD::n_par = 14;

// Element information
ItemInfo DHLD::einfo =
{
  "dhld",
  "Double Heterojunction Laser Diode",
  "Houssam Kanj",
  DEFAULT_ADDRESS"category:laser diode,diode>laser"
};

// Parameter information
ParmInfo DHLD::pinfo[] =
{
  {"rs", "Series Resistance value (Ohms)", TR_DOUBLE, false},
  {"re", "Nonlinear Series Resistance value (Ohms)", TR_DOUBLE, false},
  {"i01", "Equivalent Diode1 Saturation Current (Ampere)", TR_DOUBLE, false},
  {"i02", "Equivalent Diode2 Saturation Current (Ampere)", TR_DOUBLE, false},
  {"b", "Current Controlled Current Source gain (1/Ampere)", TR_DOUBLE, false},
  {"tns", "Equivqlent Reconbination Lifetime (Sec)", TR_DOUBLE, false},
  {"c0", "Zero-Bias Capacitance (Farad)", TR_DOUBLE, false},
  {"vd", "Diode Built-in potential (Volt)", TR_DOUBLE, false},
  {"d", "Constant relate the radiative recombination current per unit volume to the optical gain (meter^6/Ampere/Volt)", TR_DOUBLE, false},
  {"a", "nonradiative recombination lifetime/(nonradiative recombination lifetime + low level injection spontaneous recombination lifetime) = tn/(ts+tn)", TR_DOUBLE, false},
  {"rp", " --- (Ohms)", TR_DOUBLE, false},
  {"cp", " --- (Farad)", TR_DOUBLE, false},
  {"sc", "Photon density normalisation constant (meter^-3)", TR_DOUBLE, false},
  {"beta", "fraction of spontaneous emission coupled into the lasing mode", TR_DOUBLE, false}

};

DHLD::DHLD(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(rs = 2);
  paramvalue[1] = &(re = 0.468);
  paramvalue[2] = &(i01 = 2.54e-25);
  paramvalue[3] = &(i02 = 18.13e-3);
  paramvalue[4] = &(b = 6.92);
  paramvalue[5] = &(tns = 2.25e-9);
  paramvalue[6] = &(c0 = 10e-12);
  paramvalue[7] = &(vd = 1.65);
  paramvalue[8] = &(d = 1.79e-29);
  paramvalue[9] = &(a = 0.125);
  paramvalue[10] = &(rp = 29.4);
  paramvalue[11] = &(cp = 0.102e-12);
  paramvalue[12] = &(sc = 1e21);
  paramvalue[13] = &(beta = 1e-3);

  // Set the number of terminals
  setNumTerms(4);

  // Set flags
  setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(2);
}

void DHLD::init() throw(string&)
{
  vt = 25.6802271e-3; q = 1.6022e-19;
  Vpara = log(vt/i01)*vt;

  DenseIntVector var(2);
  DenseIntVector novar;
  DenseDoubleVector nodelay;
  var[0] = 0;
  var[1] = 1;
  initializeAD(var, var, novar, novar, nodelay);
}

void DHLD::eval(AD * x, AD * effort, AD * flow)
{
  // x[0] state variable to compute v1 & i1
  // x[1]=Sn
  // x[2]=dx[0]/dt  time derivative of state variable
  // x[3]=dSn/dt

  AD v1, dv1_dx, i1, dv1_dt;

  if (Vpara > x[0])
  {
    v1 = x[0] + zero;
    dv1_dx = one;
    i1 = i01*(exp(x[0]/vt) - one);
  }
  else
  {
    v1 = Vpara + vt*log(one + (1/vt)*(x[0]-Vpara));
    dv1_dx = one/(one + (1/vt)*(x[0]-Vpara));
    i1 = i01*exp(Vpara/vt)*(one + (1/vt)*(x[0] - Vpara)) - i01;
  }
	dv1_dt = dv1_dx * x[2];

  AD v2 = vt*log(i1/i02 +1);
  AD vre = i1 * re;
  AD vj = vre + v1 + v2;
  AD ie = a*i1 + b*i1*i1;
  AD G = d*pow((ie/(cp/q/sc) - 2e13),2);

	AD x1, dx1_dt;

	x1 = pow((x[1] - 2), 2) - .01; //best parametrisation
	dx1_dt = (2*pow((x[1] - 2),1))*x[3];//best parametrization

  flow[0] = i1 + b*i1*i1 + tns*((i1+i01)/vt)*dv1_dt +
	c0*pow((1-vj/vd), -0.5)*( dv1_dt + (i01/i02)*exp((v1-v2)/vt)*dv1_dt
	+ re*((i1+i01)/vt)*dv1_dt) + G*x1;

  effort[0] = rs*flow[0] + vj;

  flow[1] = -G*x1 - beta*ie + x1/rp + cp*dx1_dt;

  effort[1] = x1;
}

void DHLD::getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list)
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
}
