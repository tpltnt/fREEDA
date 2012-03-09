#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Transformer.h"

// Static members
const unsigned Transformer::n_par = 3; 

// Element information
ItemInfo Transformer::einfo =
{
	"transformer",
	"Transformer 2-port",
	"Olga Andreescu",
	"category:passive, 2-port transformer",
	"2008_April"
};

// Parameter information
ParmInfo Transformer::pinfo[] =
{
	{"n1", "Number of turns on primary winding", TR_DOUBLE, false},
	{"n2", "Number of turns on secondary winding", TR_DOUBLE, false},
	{"R", "Leakage resistance", TR_DOUBLE, false}
};

Transformer::Transformer(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
	// Values not required: 
	paramvalue[0] = &(n1=1.);
	paramvalue[1] = &(n2=1.);
	paramvalue[2] = &(R=1.e10);

	// Set the number of terminals; this model is 4 terminals
	setNumTerms(4);

	// Set flags for multi reference local terms
	setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN);
}

void Transformer::init() throw(string&)
{				
	return;
}

unsigned Transformer::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
   	my_row = eqn_number;
	return 1;	
}

void Transformer::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
	assert(my_row);
	first_eqn = my_row;
	n_rows = 1;			//add one extra row to matrix
}

void Transformer::getLocalRefIdx(UnsignedVector& local_ref_vec,
TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));	// Input Terminal np
  term_list.push_back(getTerminal(1));	// local reference terminal
  term_list.push_back(getTerminal(2)); 	// Output Terminal ns
  term_list.push_back(getTerminal(3)); 	// Local reference terminal 

  local_ref_vec.push_back(1); // Local reference index
  local_ref_vec.push_back(3); // Local reference index
}

void Transformer::fillMNAM(FreqMNAM* mnam)
{
  // Ask my terminals the row/column numbers 				
	const unsigned term1=getTerminal(0)->getRC();
	const unsigned term2=getTerminal(1)->getRC();
	const unsigned term3=getTerminal(2)->getRC();
	const unsigned term4=getTerminal(3)->getRC();

      // There is no need to set the conductances. 
	// Scaling provides a better condition number
	// for analysis run.
      double turns_ratio = n2/n1;
	double_complex cxg(1./R,0);
		
	assert(my_row);
	//set the extra current vector
	mnam->setElement(term1, my_row, cxg);
	mnam->setElement(term2, my_row, -cxg);
	mnam->setElement(term3, my_row, -cxg/turns_ratio);
	mnam->setElement(term4, my_row, cxg/turns_ratio);

	//set the extra voltage vector
	mnam->setElement(my_row, term1 , turns_ratio*cxg);
	mnam->setElement(my_row, term2 , -turns_ratio*cxg);
	mnam->setElement(my_row, term3 , -cxg);
	mnam->setElement(my_row, term4 , cxg);
}

void Transformer::fillMNAM(TimeMNAM* mnam)
{
  	// Ask my terminals the row/column numbers				
  	const unsigned term1=getTerminal(0)->getRC();
 	const unsigned term2=getTerminal(1)->getRC();
  	const unsigned term3=getTerminal(2)->getRC();
  	const unsigned term4=getTerminal(3)->getRC();
  
      // The linear state variable is the current through the 
	// primary coil. "g" will scale this current, so that the enteries
	// in MNAM are typical of those from other elements in the circuit.
      double turns_ratio = n2/n1;
	double g = 1./R;
	
	assert(my_row);
	//set the primary coil conductance
	mnam->setMElement(term1, term1, g);
	mnam->setMElement(term1, term2, -g);
	mnam->setMElement(term2, term1, -g);
	mnam->setMElement(term2, term2, g);

	//set the secondary coil conductance
	mnam->setMElement(term3, term3, g);
	mnam->setMElement(term3, term4, -g);
	mnam->setMElement(term4, term3, -g);
	mnam->setMElement(term4, term4, g);

	//set the extra current column
	mnam->setMElement(term1, my_row, g );
	mnam->setMElement(term2, my_row, -g);
	mnam->setMElement(term3, my_row, -g/turns_ratio);
	mnam->setMElement(term4, my_row, g/turns_ratio);
	
	//set the extra voltage row
	mnam->setMElement(my_row, term1, turns_ratio*g);	
	mnam->setMElement(my_row, term2, -turns_ratio*g);
	mnam->setMElement(my_row, term3, -g);
	mnam->setMElement(my_row, term4, g);
}

