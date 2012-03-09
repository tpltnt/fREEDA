#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "IsfFM.h"

// Static members
const unsigned IsfFM::n_par = 5 ;

// Element information
ItemInfo IsfFM::einfo = 
{
  "isffm",
  "Single Frequency FM current source",
  "Satish V. Uppathil",
  DEFAULT_ADDRESS"source"
} ;

// Parameter information
ParmInfo IsfFM::pinfo[] = 
{
  {"io", "Offset current (V)", TR_DOUBLE, false},
  {"ia", "RMS current amplitude (V)", TR_DOUBLE, false},
  {"fcarrier", "AC frequency (Hz)", TR_DOUBLE, false},
  {"mdi", "modulation index (No Dimension)", TR_DOUBLE, false},
  {"fsignal", "Signal frequency (Hz)", TR_DOUBLE, false}
  } ;

IsfFM::IsfFM(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(io) ;
  paramvalue[1] = &(ia) ;
  paramvalue[2] = &(fcarrier = zero) ;
  paramvalue[3] = &(mdi = zero) ;
  paramvalue[4] = &(fsignal = zero) ;

  // Set the number of terminals
  setNumTerms(2) ;

  // Set flags
  setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN | SOURCE) ;
  // Set number of states
  setNumberOfStates(1) ;

}

void IsfFM::fillMNAM(FreqMNAM* mnam)
{

  const double& freq(mnam->getFreq()) ;
  // Set source vector

  if (freq == zero) {
    mnam->addToSource(getTerminal(0)->getRC(),
		      getTerminal(1)->getRC(),
		      io) ;
    mnam->addToSource(getTerminal(0)->getRC(),
			getTerminal(1)->getRC(),
			ia) ;}

  // Since the source vector is assumed to be initialized to zero,
  // we do not need to fill it if freq != frequency.
  return ;
}
  
void IsfFM::fillMNAM(TimeMNAM* mnam)
{
  // Nothing to be done.
  return ;
}

void IsfFM::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime() ;
  
  double is = io + ia * sin(2 * pi * fcarrier * ctime + mdi * sin(2 * pi * fsignal * ctime)) ;
  
  mnam->addToSource(getTerminal(0)->getRC(),
		    getTerminal(1)->getRC(),
		    is) ;
  return ;
}

void IsfFM::svTran(TimeDomainSV* tdsv)
{
  // Calculate current
  double& is = tdsv->i(0) ;

  is = io ;
  if (!tdsv->DC()) {
    const double& ctime = tdsv->getCurrentTime() ;

    is = io + ia * sin(2 * pi * fcarrier * ctime + mdi * sin(2 * pi * fsignal * ctime)) ;
  }
  
  // State variable is voltage
  tdsv->u(0) = tdsv->getX(0) ;
}

void IsfFM::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = one ;
  tdsv->getJi()(0,0) = zero ;
}
