#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "VAM.h"

// Static members
const unsigned VAM::n_par = 5 ;

// Element information
ItemInfo VAM::einfo =
{
  "vam",
  "Amplitude Modulated voltage source",
  "Satish V. Uppathil",
  DEFAULT_ADDRESS"category:source",
  "2001_06_20"
};

// Parameter information
ParmInfo VAM::pinfo[] =
{
  {"oc", "Offset Constant (Dimensionless)", TR_DOUBLE, false},
  {"sa", "Signal amplitude (V)", TR_DOUBLE, false},
  {"fcarrier", "Carrier frequency (Hz)", TR_DOUBLE, false},
  {"fmod", "modulation frequency(Hz)", TR_DOUBLE, false},
  {"td", "Time Delay (seconds)", TR_DOUBLE, false}
} ;

VAM::VAM(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(oc = zero) ;
  paramvalue[1] = &(sa = zero) ;
  paramvalue[2] = &(fcarrier = zero) ;
  paramvalue[3] = &(fmod = zero) ;
  paramvalue[4] = &(td = zero) ;

  // Set the number of terminals
  setNumTerms(2) ;

  // Set flags
  setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN | SOURCE);
  // Set number of states
  setNumberOfStates(1) ;
  my_row = 0 ;
}

unsigned VAM::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Keep the equation number assigned to this element
  my_row = eqn_number ;
  // Add one extra RC
  return 1 ;
}

void VAM::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row) ;
  first_eqn = my_row ;
  n_rows = 1 ;
}

void VAM::fillMNAM(FreqMNAM* mnam)
{
  assert(my_row) ;
  // Ask my terminals the row numbers
  mnam->setOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row) ;
  mnam->setSource(my_row, sa) ;
  return ;
}


void VAM::fillMNAM(TimeMNAM* mnam)
{
  assert(my_row) ;
  // Ask my terminals the row numbers
  mnam->setMOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row) ;
  return ;
}

void VAM::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime() ;
  double e = sa * (oc + sin(2 * pi * fmod * (ctime - td))) * sin(2 * pi * fcarrier * (ctime - td)) ;
  mnam->setSource(my_row, e) ;
  return ;
}

void VAM::svTran(TimeDomainSV* tdsv)
{
  // Calculate voltage
  double& e = tdsv->u(0) ;

  if (!tdsv->DC())
	{
    const double& ctime = tdsv->getCurrentTime() ;
    e = sa * (oc + sin(2 * pi * fmod * (ctime - td))) * sin(2 * pi * fcarrier * (ctime - td)) ;
  }

  // Scale state variable for numerical stability
  tdsv->i(0) = tdsv->getX(0) * 1e-2 ;
}

void VAM::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = zero ;
  tdsv->getJi()(0,0) = 1e-2 ;
}

