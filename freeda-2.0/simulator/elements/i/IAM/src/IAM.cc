#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "IAM.h"

// Static members
const unsigned IAM::n_par = 5;

// Element information
ItemInfo IAM::einfo = 
{
  "iam",
  "Amplitude Modulated current source",
  "Satish V. Uppathil",
  DEFAULT_ADDRESS"source",
  "2001_06_20"
};

// Parameter information
ParmInfo IAM::pinfo[] = 
{
  {"oc", "Offset Constant (Dimensionless)", TR_DOUBLE, false},
  {"sa", "Signal amplitude (V)", TR_DOUBLE, false},
  {"fcarrier", "Carrier frequency (Hz)", TR_DOUBLE, false},
  {"fmod", "modulation frequency(Hz)", TR_DOUBLE, false},
  {"td", "Time Delay (seconds)", TR_DOUBLE, false}
};

IAM::IAM(const string& iname) : Element(&einfo, pinfo, n_par, iname)
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
  setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN | SOURCE) ;
  // Set number of states
  setNumberOfStates(1) ;

}

const char* IAM::getNetlistName()
{
  return einfo.name;
}


//void IAM::fillMNAM(FreqMNAM* mnam)
//{
//    // This is not correct. The spectrum of the signal must be calculated.
//    // Set source vector
  
//    mnam->addToSource(getTerminal(0)->getRC(),
//  		    getTerminal(1)->getRC(),
//  		    sa) ;
//    mnam->addToSource(getTerminal(0)->getRC(),
//  		    getTerminal(1)->getRC(),
//  		    sa) ;

//    // Since the source vector is assumed to be initialized to zero,
//    return ;
//}
  
void IAM::fillMNAM(TimeMNAM* mnam)
{
  // Nothing to be done.
  return ;
}

void IAM::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime() ;
  
  double is = sa * (oc + sin(2 * pi * fmod * (ctime - td))) 
    * sin(2 * pi * fcarrier * (ctime - td)) ;
  
  mnam->addToSource(getTerminal(0)->getRC(),
		    getTerminal(1)->getRC(),
		    is) ;
  return ;
}

void IAM::svTran(TimeDomainSV* tdsv)
{
  // Calculate current
  double& is = tdsv->i(0) ;

  is = sa ;
  if (!tdsv->DC()) {
    const double& ctime = tdsv->getCurrentTime() ;

    is = sa * (oc + sin(2 * pi * fmod * (ctime - td))) 
      * sin(2 * pi * fcarrier * (ctime - td)) ;
  }
  
  // State variable is voltage
  tdsv->u(0) = tdsv->getX(0) ;
}

void IAM::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = one ;
  tdsv->getJi()(0,0) = zero ;
}
