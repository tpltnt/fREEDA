#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Vexp.h"
#include "../../../../analysis/TimeDomainSV.h"

// Static members
const unsigned Vexp::n_par = 6 ;

// Element information
ItemInfo Vexp::einfo =
{
  "vexp",
  "Exponential voltage source",
  "Satish V. Uppathil",
  DEFAULT_ADDRESS"category:source"
} ;

// Parameter information
ParmInfo Vexp::pinfo[] = {
  {"v1", "Inital voltage (V)", TR_DOUBLE, false},
  {"v2", "Final voltage (V)", TR_DOUBLE, false},
  {"tdr", "Rise Time delay (Sec)", TR_DOUBLE, false},
  {"tdf", "Fall Time delay (Sec)", TR_DOUBLE, false},
  {"tcr", "Rise Time Constant (Sec)", TR_DOUBLE, false},
  {"tcf", "Fall Time Constant (Sec)", TR_DOUBLE, false}
} ;

Vexp::Vexp(const string& iname) : Element(&einfo, pinfo, n_par, iname)
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

  my_row = 0 ;
}

unsigned Vexp::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Keep the equation number assigned to this element
  my_row = eqn_number ;

  // Add one extra RC
  return 1 ;
}

void Vexp::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row) ;
  first_eqn = my_row ;
  n_rows = 1 ;
}

void Vexp::fillMNAM(TimeMNAM* mnam)
{
  assert(my_row) ;

  // Ask my terminals the row numbers
  mnam->setMOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row) ;
  return ;
}

void Vexp::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime() ;
  double e = v1 ;
  if ((ctime > tdr) && (ctime < tdf))
    e = v1 + (v2 - v1) * (1 - exp(-(ctime - tdr) / tcr)) ;
  if (ctime > tdf)
    e = v1 + (v2 - v1) * (1 - exp(-(tdf - tdr) / tcr)) * exp((ctime - tdf) / tcf) ;
  mnam->setSource(my_row, e) ;
  return ;
}

void Vexp::svTran(TimeDomainSV* tdsv)
{
  // Calculate voltage
  double& e = tdsv->u(0) ;
  e = v1 ;
  if (!tdsv->DC()) {
    const double& ctime = tdsv->getCurrentTime() ;
		if ((ctime > tdr) && (ctime < tdf))
			e = v1 + (v2 - v1) * (1 - exp(-(ctime - tdr) / tcr)) ;
		if (ctime > tdf)
			e = v1 + (v2 - v1) * (1 - exp(-(tdf - tdr) / tcr)) * exp((ctime - tdf) / tcf) ;
  }

  // Scale state variable for numerical stability
  tdsv->i(0) = tdsv->getX(0) * 1e-2 ;
}

void Vexp::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = zero ;
  tdsv->getJi()(0,0) = 1e-2 ;
}

