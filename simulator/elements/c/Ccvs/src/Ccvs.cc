#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Ccvs.h"

// Static members
const unsigned Ccvs::n_par = 3;

// Element information
ItemInfo Ccvs::einfo =
{
  "ccvs",
  "Current-controlled voltage source",
  "Andrew Mugisha",
  DEFAULT_ADDRESS"elements/Ccvs.h.html",
  "2008_08_29"
} ;

// Parameter information
ParmInfo Ccvs::pinfo[] =
{
  {"r", "Transresistance", TR_DOUBLE, false},
  {"ri", "Input resistance value (Ohms)", TR_DOUBLE, false},
  {"ro", "Output resistance value (Ohms)", TR_DOUBLE, false}
};

Ccvs::Ccvs(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Value of r is required
  paramvalue[0] = &(r=1) ;
  paramvalue[1] = &(ri = zero) ;
  paramvalue[2] = &(ro = zero) ;
  
  setNumTerms(4);

  // Set flags
  setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN) ;
  
  my_row = 0;
}

unsigned Ccvs::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Keep the equation number assigned to this element
  my_row = eqn_number;
  // Add one extra RC
  return 2;

}

void Ccvs::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row);
  first_eqn = my_row;
  n_rows = 2;
}

void Ccvs::getLocalRefIdx(UnsignedVector& local_ref_vec,
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


void Ccvs::fillMNAM(FreqMNAM* mnam)
{
  if (ri)
    mnam->setAdmittance(getTerminal(2)->getRC(), getTerminal(3)->getRC(),
	one/ri) ;

  if (ro)
    mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
	one/ro) ;

  mnam->setOnes(getTerminal(2)->getRC(), getTerminal(3)->getRC(), my_row) ;
  mnam->setOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row+1) ;
  mnam->setElement(my_row+1, my_row, -r) ;

}

void Ccvs::fillMNAM(TimeMNAM* mnam)
{
  if (ri)
    mnam->setMAdmittance(getTerminal(2)->getRC(), getTerminal(3)->getRC(),
	one/ri) ;

  if (ro)
    mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
	one/ro) ;

  //assert(my_row) ;
  mnam->setMOnes(getTerminal(2)->getRC(), getTerminal(3)->getRC(), my_row) ;
  mnam->setMOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row+1) ;
  mnam->setMElement(my_row+1, my_row, -r) ;

}
