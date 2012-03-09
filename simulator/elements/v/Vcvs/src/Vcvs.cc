#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Vcvs.h"

// Static members
const unsigned Vcvs::n_par = 3;

// Element information
ItemInfo Vcvs::einfo =
{
  "vcvs",
  "Voltage-controlled voltage source",
  "Satish V. Uppathil, Justin Lowry, Nikhil Kriplani, Chris Saunders",
  DEFAULT_ADDRESS"elements/Vcvs.h.html",
  "2001_06_20"
} ;

// Parameter information
ParmInfo Vcvs::pinfo[] =
{
  {"k", "gain", TR_DOUBLE, false},
  {"ri", "Input resistance value (Ohms)", TR_DOUBLE, false},
  {"ro", "Output resistance value (Ohms)", TR_DOUBLE, false}
};

Vcvs::Vcvs(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Value of r is required
  paramvalue[0] = &(k=1) ;
  paramvalue[1] = &(ri = zero) ;
  paramvalue[2] = &(ro = zero) ;
  
  setNumTerms(4);

  // Set flags
  setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN) ;

  my_row = 0;
}

unsigned Vcvs::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Keep the equation number assigned to this element
  my_row = eqn_number;
  // Add one extra RC
  return 1;
}

void Vcvs::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row);
  first_eqn = my_row;
  n_rows = 1;
}

void Vcvs::getLocalRefIdx(UnsignedVector& local_ref_vec,
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


void Vcvs::fillMNAM(FreqMNAM* mnam)
{
  if (ri)
    mnam->setAdmittance(getTerminal(2)->getRC(), getTerminal(3)->getRC(),
	one/ri) ;

  if (ro)
    mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
	one/ro) ;

  mnam->setOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row) ;
  mnam->setElement(my_row, getTerminal(2)->getRC(), -k) ;
  mnam->setElement(my_row, getTerminal(3)->getRC(), k) ;
}

void Vcvs::fillMNAM(TimeMNAM* mnam)
{
  if (ri)
    mnam->setMAdmittance(getTerminal(2)->getRC(), getTerminal(3)->getRC(),
	one/ri) ;

  if (ro)
    mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
	one/ro) ;

  //assert(my_row) ;
  mnam->setMOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row) ;
  mnam->setMElement(my_row, getTerminal(2)->getRC(), -k) ;
  mnam->setMElement(my_row, getTerminal(3)->getRC(), k) ;
}
