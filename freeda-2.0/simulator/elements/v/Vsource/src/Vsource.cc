#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "Vsource.h"


// Static members
const unsigned Vsource::n_par = 7;


// Element information
ItemInfo Vsource::einfo =
{
  "vsource",
  "General DC and sinusoidal voltage source",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:source",
  "2000_07_20"
};


// Parameter information
ParmInfo Vsource::pinfo[] = 
{
  {"vdc", "DC voltage (V)", TR_DOUBLE, false},
  {"vac", "AC voltage peak amplitude (V)", TR_DOUBLE, false},
  {"f", "AC frequency (Hz)", TR_DOUBLE, false},
  {"phase", "Source phase (degrees)", TR_DOUBLE, false},
  {"delay", "Delay to AC start (s)", TR_DOUBLE, false},
  {"tr", "Rising time for DC component (s)", TR_DOUBLE, false},
  {"periods", "num of periods for voltage to ramp linearly", TR_DOUBLE, false} 
};


Vsource::Vsource(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(vdc = zero);
  paramvalue[1] = &(vac = zero);
  paramvalue[2] = &(frequency = zero);
  paramvalue[3] = &(phase = zero);
  paramvalue[4] = &(delay = zero);
  paramvalue[5] = &(tr = zero);
  paramvalue[6] = &(periods = 0.0);

  // Set the number of terminals
  setNumTerms(2);
  // Set flags
  setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN | SOURCE);
  // Set number of states
  setNumberOfStates(1);

  my_row = 0;
}


void Vsource::init() throw(string&)
{
  omega = twopi * frequency;
  phase_rad = deg2rad * phase;
  ramp_slope = frequency / periods;
}


unsigned Vsource::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Keep the equation number assigned to this element
  my_row = eqn_number;
  // Add one extra RC
  return 1;
}


void Vsource::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row);
  first_eqn = my_row;
  n_rows = 1;
}


void Vsource::fillMNAM(FreqMNAM* mnam)
{
  assert(my_row);
  // Ask my terminals the row numbers
  mnam->setOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);

  const double& freq(mnam->getFreq());
  // Set source vector
  if (freq == zero)
  {
    mnam->setSource(my_row, vdc);
  }
  else if ((freq && abs(freq - frequency)/freq < 1e-14 ) || !frequency)
  {
    if (phase == zero)
      mnam->setSource(my_row, vac);
    else
    {
      double_complex vs(vac * cos(phase_rad),	vac * sin(phase_rad));
      mnam->setSource(my_row, vs);
    }
  }
  // Since the source vector is assumed to be initialized to zero,
  // we do not need to fill it if freq != frequency.
  return;
}


void Vsource::fillMNAM(TimeMNAM* mnam)
{
  assert(my_row);
  // Ask my terminals the row numbers
  mnam->setMOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);
  return;
}


void Vsource::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime();
  double e = (tr && ctime < tr) ? vdc * ctime / tr : vdc;
  if ((ctime > delay) && (ctime < periods/frequency))
    e += vac * cos(omega * (ctime - delay) + phase_rad) * (ramp_slope * ctime);
  else if (ctime > delay)
    e += vac * cos(omega * (ctime - delay) + phase_rad);
  mnam->setSource(my_row, e);
  return;
}


void Vsource::svTran(TimeDomainSV* tdsv)
{
  // Calculate voltage
  double& e = tdsv->u(0);

  e = (tr) ? zero : vdc;
  if (!tdsv->DC()) 
  {
    const double& ctime = tdsv->getCurrentTime();
    e = (tr && ctime < tr) ? vdc * ctime / tr : vdc;
    if ((ctime > delay) && (ctime < periods/frequency))
      e += vac * cos(omega * (ctime - delay) + phase_rad) * (ramp_slope * ctime);
    else if (ctime > delay)
      e += vac * cos(omega * (ctime - delay) + phase_rad);
  }

  // Scale state variable for numerical stability
  tdsv->i(0) = tdsv->getX(0) * 1e-2;
}


void Vsource::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = zero;
  tdsv->getJi()(0,0) = 1e-2;
}

