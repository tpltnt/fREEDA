#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Isolator.h"

// Static members
const unsigned Isolator::n_par = 2;

// Element information
ItemInfo Isolator::einfo =
{
  "isolator",
  "Isolator",
  "Don Widdows - Isac Lima - Daryl Lindsey",
  DEFAULT_ADDRESS"elements/Isolator.h.html",
  "2003_05_15"
};


// Parameter information
ParmInfo Isolator::pinfo[] =
{
  {"r", "Resistance looking into both ports of isolator (Ohms)", TR_DOUBLE, false},
  {"g", "Conductance looking into both ports of isolator (S)", TR_DOUBLE, false}
};


Isolator::Isolator(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  paramvalue[0] = &(r=0);
  paramvalue[1] = &(g=0);

  //future variations may require inductance, since isolators contain
  //ferrites

  // Set the number of terminals
  setNumTerms(4);

	// Set flags
  setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN); //must put multi-ref
	// for local ref term's!
  //Linear means it's going to check for fillMNAM functions.
  //fREEDA is written to have certain conventions about functions.
  //It has to know which functions to call from within the main body
  //  of the code when certain flags are set.
}

void Isolator::init() throw(string&)
{
  if(r) //(order of precedence is r, then g)
		{g=1/r;}
  else if (!g)
		{g=.02;} //if no resistance or conductance specified,
	//make conductance = .02 S
	//this is set up so that the value of the conductance will default to .02S
	// (resistance defaults to 50 ohms) if neither parameter is specified.
}


/*  For future implementation of non-ideal isolator
//unsigned Isolator::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
//{
	//   for (int i=0; i < 2; i++)
	//    my_row[i] = eqn_number + i;
	//
	//   // Add 2 extra RCs
	//    return 2;
	//
//}

//void Isolator::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
//{
	//   assert(my_row[0]);
	//   first_eqn = my_row[0];
	//   n_rows = 2;
	//   assert(my_row[0]);
	//
//}
*/


void Isolator::getLocalRefIdx(UnsignedVector& local_ref_vec,
TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1));
  term_list.push_back(getTerminal(2)); //local reference terminal
  term_list.push_back(getTerminal(3)); // local reference terminal

  local_ref_vec.push_back(2); // Local reference index
  local_ref_vec.push_back(3); // Local reference index
}


void Isolator::fillMNAM(FreqMNAM* mnam)
{
  //initialize admittance values
  //(values may be complex in a future non-ideal implementation)

  double y=g;

  // Ask my terminals the row numbers
	const unsigned term1=getTerminal(0)->getRC();
	const unsigned term2=getTerminal(1)->getRC();
	const unsigned term3=getTerminal(2)->getRC();
	const unsigned term4=getTerminal(3)->getRC();

  //2-Port (4-terminal) isolator model
  mnam->setElement(term1,term1,y);
  mnam->setElement(term1,term2,0);
  mnam->setElement(term1,term3,-y);
  mnam->setElement(term1,term4,0);
  mnam->setElement(term2,term1,-2*y);
  mnam->setElement(term2,term2,y);
  mnam->setElement(term2,term3,2*y);
  mnam->setElement(term2,term4,-y);
  mnam->setElement(term3,term1,-y);
  mnam->setElement(term3,term2,0);
  mnam->setElement(term3,term3,y);
  mnam->setElement(term3,term4,0);
  mnam->setElement(term4,term1,2*y);
  mnam->setElement(term4,term2,-y);
  mnam->setElement(term4,term3,-2*y);
  mnam->setElement(term4,term4,y);
}


void Isolator::fillMNAM(TimeMNAM* mnam)
{
  //got to get the terminal numbers
  const unsigned term1=getTerminal(0)->getRC();
  const unsigned term2=getTerminal(1)->getRC();
  const unsigned term3=getTerminal(2)->getRC();
  const unsigned term4=getTerminal(3)->getRC();

  //for now, the isolator is implemented as a purely
  //resistive (or conductive) element

  //2-Port (4-terminal) isolator Model
  mnam->setMElement(term1,term1,g);
  mnam->setMElement(term1,term2,0);
  mnam->setMElement(term1,term3,-g);
  mnam->setMElement(term1,term4,0);
  mnam->setMElement(term2,term1,-2*g);
  mnam->setMElement(term2,term2,g);
  mnam->setMElement(term2,term3,2*g);
  mnam->setMElement(term2,term4,-g);
  mnam->setMElement(term3,term1,-g);
  mnam->setMElement(term3,term2,0);
  mnam->setMElement(term3,term3,g);
  mnam->setMElement(term3,term4,0);
  mnam->setMElement(term4,term1,2*g);
  mnam->setMElement(term4,term2,-g);
  mnam->setMElement(term4,term3,-2*g);
  mnam->setMElement(term4,term4,g);
}

