#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Resistor.h"

// Static members
const unsigned Resistor::n_par = 1;

// Element information
ItemInfo Resistor::einfo = 
{
  "resistor",
  "Resistor",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:lumped",
  "2000_07_15"
};

// Parameter information
ParmInfo Resistor::pinfo[] = 
{
  {"r", "Resistance value (Ohms)", TR_DOUBLE, true}
};

Resistor::Resistor(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Value of r is required
  paramvalue[0] = &res;
	
  // Set the number of terminals
  setNumTerms(2);
	
  // Set flags
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);
}

void Resistor::fillMNAM(FreqMNAM* mnam)
{
  // Ask my terminals the row numbers
  mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
	one/res);
}

void Resistor::fillMNAM(TimeMNAM* mnam)
{
  // Ask my terminals the row numbers
  mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
	one/res);
}

