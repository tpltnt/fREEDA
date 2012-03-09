#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Vccs.h"

// Static members
const unsigned Vccs::n_par = 3;

// Element information
ItemInfo Vccs::einfo =
{
  "vccs",
  "Voltage-controlled current source",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"elements/Vccs.h.html",
  "2000_07_20"
};

// Parameter information
ParmInfo Vccs::pinfo[] =
{
  {"g", "Transconductance (Siemens)", TR_DOUBLE, true},
  {"ri", "Input resistance value (Ohms)", TR_DOUBLE, false},
  {"ro", "Output resistance value (Ohms)", TR_DOUBLE, false}
};


Vccs::Vccs(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Value of r is required
  paramvalue[0] = &g;
  paramvalue[1] = &(ri = zero);
  paramvalue[2] = &(ro = zero);

  // Set the number of terminals
  setNumTerms(4);

  // Set flags
  setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN);
}

void Vccs::getLocalRefIdx(UnsignedVector& local_ref_vec,
TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  // in the netlist, the output terminals are specified first
  // followed by the input terminals
  // so the input terminals are specified in the 2nd and 3rd 
  // position in the declaration of this element, and the output
  // terminals are specified in the 0th and 1st position
  term_list.push_back(getTerminal(2));
  term_list.push_back(getTerminal(3)); // Local reference terminal
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1)); // Local reference terminal

  local_ref_vec.push_back(1); // Local reference index
  local_ref_vec.push_back(3); // Local reference index
}

void Vccs::fillMNAM(FreqMNAM* mnam)
{
  if (ri)
    mnam->setAdmittance(getTerminal(2)->getRC(), getTerminal(3)->getRC(),
	one/ri);

  if (ro)
    mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
	one/ro);

  mnam->setQuad(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
	getTerminal(2)->getRC(), getTerminal(3)->getRC(),
	g);
}

void Vccs::fillMNAM(TimeMNAM* mnam)
{
  if (ri)
    mnam->setMAdmittance(getTerminal(2)->getRC(), getTerminal(3)->getRC(),
	one/ri);

  if (ro)
    mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
	one/ro);

  mnam->setMQuad(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
	getTerminal(2)->getRC(), getTerminal(3)->getRC(),
	g);
}

