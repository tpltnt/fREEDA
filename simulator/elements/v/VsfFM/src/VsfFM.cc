#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "VsfFM.h"
#include "../../../../analysis/TimeDomainSV.h"

// Static members
const unsigned VsfFM::n_par = 5 ;

// Element information
ItemInfo VsfFM::einfo =
{
  "vsffm",
  "Single Frequency FM voltage source",
  "Satish V. Uppathil",
  DEFAULT_ADDRESS"category:source",
  "2001_06_20"
} ;

// Parameter information
ParmInfo VsfFM::pinfo[] =
{
  {"vo", "Output voltage Offset (V)", TR_DOUBLE, false},
  {"va", "Output voltage amplitude (V)", TR_DOUBLE, false},
  {"fcarrier", "Carrier frequency (Hz)", TR_DOUBLE, false},
  {"mdi", "modulation index (No Dimension)", TR_DOUBLE, false},
  {"fsignal", "Signal frequency (Hz)", TR_DOUBLE, false}
} ;

VsfFM::VsfFM(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(vo) ;
  paramvalue[1] = &(va) ;
  paramvalue[2] = &(fcarrier = zero) ;
  paramvalue[3] = &(mdi = zero) ;
  paramvalue[4] = &(fsignal = zero) ;

  // Set the number of terminals
  setNumTerms(2) ;

  // Set flags
  setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN | SOURCE) ;
  // Set number of states
  setNumberOfStates(1) ;

  my_row = 0 ;
}

unsigned VsfFM::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Keep the equation number assigned to this element
  my_row = eqn_number ;
  // Add one extra RC
  return 1 ;
}

void VsfFM::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row) ;
  first_eqn = my_row ;
  n_rows = 1 ;
}

void VsfFM::fillMNAM(FreqMNAM* mnam)
{
  assert(my_row) ;
  // Ask my terminals the row numbers
  mnam->setOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row) ;

	mnam->setSource(my_row, vo) ;
  return ;
}

void VsfFM::fillMNAM(TimeMNAM* mnam)
{
  assert(my_row) ;
  // Ask my terminals the row numbers
  mnam->setMOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row) ;
  return ;
}

void VsfFM::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime() ;

  double e = vo + va * sin(2 * pi * fcarrier * ctime + mdi * sin(2 * pi * fsignal * ctime)) ;

  mnam->setSource(my_row, e) ;
  return ;
}

void VsfFM::svTran(TimeDomainSV* tdsv)
{
  // Calculate voltage
  double& e = tdsv->u(0) ;


  if (!tdsv->DC()) {
    const double& ctime = tdsv->getCurrentTime() ;

		e = vo + va * sin(2 * pi * fcarrier * ctime + mdi * sin(2 * pi * fsignal * ctime)) ;
  }

  // Scale state variable for numerical stability
  tdsv->i(0) = tdsv->getX(0) * 1e-2 ;
}

void VsfFM::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = zero ;
  tdsv->getJi()(0,0) = 1e-2 ;
}

