#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "Vtwotone.h"


// Static members
const unsigned Vtwotone::n_par = 7;


// Element information
ItemInfo Vtwotone::einfo =
{
  "vtwotone",
  "General two-tone sinusoidal voltage source",
  "Aravind K. Mikkilineni",
  DEFAULT_ADDRESS"category:source",
  "2008_10_28"
};


// Parameter information
ParmInfo Vtwotone::pinfo[] = {
  {"vac1", "AC voltage peak amplitude (V) for f1", TR_DOUBLE, false},
  {"f1", "AC frequency 1 (Hz)", TR_DOUBLE, false},
  {"phase1", "Source phase (degrees) for f1", TR_DOUBLE, false},
  {"vac2", "AC voltage peak amplitude (V) for f2", TR_DOUBLE, false},
  {"f2", "AC frequency 2 (Hz)", TR_DOUBLE, false},
  {"phase2", "Source phase (degrees) for f2", TR_DOUBLE, false},
  {"delay", "Delay to AC start (s)", TR_DOUBLE, false},
};


Vtwotone::Vtwotone(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(vac1 = zero);
  paramvalue[1] = &(f1 = zero);
  paramvalue[2] = &(phase1 = zero);
  paramvalue[3] = &(vac2 = zero);
  paramvalue[4] = &(f2 = zero);
  paramvalue[5] = &(phase2 = zero);
  paramvalue[6] = &(delay = zero);

  // Set the number of terminals
  setNumTerms(2);
  // Set flags
  setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN | SOURCE);
  // Set number of states
  setNumberOfStates(1);

  my_row = 0;
}


void Vtwotone::init() throw(string&)
{
  omega1 = twopi * f1;
  omega2 = twopi * f2;
  phase_rad1 = deg2rad * phase1;
  phase_rad2 = deg2rad * phase2;
}


unsigned Vtwotone::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Keep the equation number assigned to this element
  my_row = eqn_number;
  // Add one extra RC
  return 1;
}


void Vtwotone::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row);
  first_eqn = my_row;
  n_rows = 1;
}


void Vtwotone::fillMNAM(FreqMNAM* mnam)
{
  assert(my_row);
  // Ask my terminals the row numbers
  mnam->setOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);

  return;
}


void Vtwotone::fillMNAM(TimeMNAM* mnam)
{
  assert(my_row);
  // Ask my terminals the row numbers
  mnam->setMOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);
  return;
}


void Vtwotone::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime();
  double e = 0.0;
    e += vac1 * cos(omega1 * (ctime - delay) + phase_rad1);
    e += vac2 * cos(omega2 * (ctime - delay) + phase_rad2);
  mnam->setSource(my_row, e);
  return;
}


void Vtwotone::svTran(TimeDomainSV* tdsv)
{
  // Calculate voltage
  double& e = tdsv->u(0);

  const double& ctime = tdsv->getCurrentTime();
  e = 0.0;
  e += vac1 * cos(omega1 * (ctime - delay) + phase_rad1);
  e += vac2 * cos(omega2 * (ctime - delay) + phase_rad2);

  // Scale state variable for numerical stability
  tdsv->i(0) = tdsv->getX(0) * 1e-2;
}


void Vtwotone::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = zero;
  tdsv->getJi()(0,0) = 1e-2;
}

