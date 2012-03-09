#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "Iexp.h"

// Static members
const unsigned Iexp::n_par = 6 ;

// Element information
ItemInfo Iexp::einfo = 
{
  "iexp",
  "Amplitude Modulated current source",
  "Satish V. Uppathil",
  DEFAULT_ADDRESS"source"
};

// Parameter information
ParmInfo Iexp::pinfo[] = 
{
  {"v1", "Inital voltage (V)", TR_DOUBLE, false},
  {"v2", "Final voltage (V)", TR_DOUBLE, false},
  {"tdr", "Rise Time delay (Sec)", TR_DOUBLE, false},
  {"tdf", "Fall Time delay (Sec)", TR_DOUBLE, false},
  {"tcr", "Rise Time Constant (Sec)", TR_DOUBLE, false},
  {"tcf", "Fall Time Constant (Sec)", TR_DOUBLE, false}
  } ;

Iexp::Iexp(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(v1) ;
  paramvalue[1] = &(v2) ;
  paramvalue[2] = &(tdr = zero) ;
  paramvalue[3] = &(tdf = zero) ;
  paramvalue[4] = &(tcr = zero) ;
  paramvalue[5] = &(tcf = zero) ;
  // Set the number of terminals
  setNumTerms(2) ;

  // Set flags
  setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN | SOURCE) ;
  // Set number of states
  setNumberOfStates(1) ;

}

void Iexp::fillMNAM(FreqMNAM* mnam)
{

  // Set source vector

    mnam->addToSource(getTerminal(0)->getRC(),
		      getTerminal(1)->getRC(),
		      v1) ;
    mnam->addToSource(getTerminal(0)->getRC(),
			getTerminal(1)->getRC(),
			v1) ;

  // Since the source vector is assumed to be initialized to zero,
  return ;
}
  
void Iexp::fillMNAM(TimeMNAM* mnam)
{
  // Nothing to be done.
  return ;
}

void Iexp::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime() ;
  
  double is = v1 ;
  
  if ((ctime > tdr) && (ctime < tdf))
    is = v1 + (v2 - v1) * (1 - exp(-(ctime - tdr) / tcr)) ;
    
  if (ctime > tdf)
    is = v1 + (v2 - v1) * (1 - exp(-(tdf - tdr) / tcr)) * exp((ctime - tdf) / tcf) ;
  
  mnam->addToSource(getTerminal(0)->getRC(),
		    getTerminal(1)->getRC(),
		    is) ;
  return ;
}

void Iexp::svTran(TimeDomainSV* tdsv)
{
  // Calculate current
  double& is = tdsv->i(0) ;

  is = v1 ;
  if (!tdsv->DC()) {
    const double& ctime = tdsv->getCurrentTime() ;

    if ((ctime > tdr) && (ctime < tdf))
    is = v1 + (v2 - v1) * (1 - exp(-(ctime - tdr) / tcr)) ;
    
  if (ctime > tdf)
    is = v1 + (v2 - v1) * (1 - exp(-(tdf - tdr) / tcr)) * exp((ctime - tdf) / tcf) ;
  }
  
  // State variable is voltage
  tdsv->u(0) = tdsv->getX(0) ;
}

void Iexp::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = one ;
  tdsv->getJi()(0,0) = zero ;
}
