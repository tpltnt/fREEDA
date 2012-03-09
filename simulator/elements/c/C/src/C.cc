#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "C.h"

// Static members
const unsigned C::n_par = 3;

// Element information
ItemInfo C::einfo =
{
  "c",
  "Linear capacitor",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:lumped",
  "2000_07_20"
};

// Parameter information
ParmInfo C::pinfo[] =
{
  {"c", "Capacitance value (F)", TR_DOUBLE, true},
  {"int_g", "Internal conductance value (S).", TR_DOUBLE, false},
  {"time_d", "Flag, if true, calculate in the time domain.", TR_BOOLEAN, false}
};


C::C(const string& iname) :
Element(&einfo, pinfo, n_par, iname)
{
  // Value of cap is required
  paramvalue[0] = &cap;
  // Default value of internal resistor. (large)
  paramvalue[1] = &(int_g = 1e-8);
  paramvalue[2] = &(time_d = false);

  // Set the number of terminals
  setNumTerms(2);
  // Set number of states
  setNumberOfStates(1);
  // Set flags
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);
}

void C::init() throw(string&)
{
  if (time_d)
    setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void C::fillMNAM(FreqMNAM* mnam)
{
  // Calculate admittance at current frequency
  double_complex y(int_g, twopi * mnam->getFreq() * cap);

  // Ask my terminals the row numbers
  mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), y);
}

void C::fillMNAM(TimeMNAM* mnam)
{
  // Ask my terminals the row numbers
  if (int_g)
    mnam->
	setMAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), int_g);
  mnam->setMpAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), cap);
}

void C::svTran(TimeDomainSV* tdsv)
{
  // The state variable is the voltage
  tdsv->u(0) = tdsv->getX(0);
  // i = C dv/dt
  tdsv->i(0) = cap * tdsv->getdX_dt(0) + tdsv->u(0) * int_g;
}

void C::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = one;
  tdsv->getJi()(0,0) = cap * tdsv->getdx_dtFactor() + int_g;
}

