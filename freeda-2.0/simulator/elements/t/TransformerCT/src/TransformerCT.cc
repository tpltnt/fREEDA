#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "TransformerCT.h"

// Static members
const unsigned TransformerCT::n_par = 6; 

// Element information
ItemInfo TransformerCT::einfo =
{
	"transformerct",
	"Center Tap Transformer",
	"Olga Andreescu",
	"category:passive, 5-terminals center tap transformer",
	"2008_April"
};

// Parameter information
ParmInfo TransformerCT::pinfo[] =
{
	{"n1", "Number of turns on primary winding", TR_DOUBLE, false},
	{"n2", "Number of turns on 1st secondary winding", TR_DOUBLE, false},
	{"n3", "Number of turns on 2nd secondary winding", TR_DOUBLE, false},
	{"k12", "coupling coefficient between coils 1&2", TR_DOUBLE, false},
	{"k13", "coupling coefficient between coils 1&3", TR_DOUBLE, false},
	{"R", "Leakage resistance", TR_DOUBLE, false}
};

TransformerCT::TransformerCT(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
	// Value not required, default assignment
	paramvalue[0] = &(n1=1.);
	paramvalue[1] = &(n2=1.);
	paramvalue[2] = &(n3=1.);
	paramvalue[3] = &(k12=1.);
	paramvalue[4] = &(k13=1.);
	paramvalue[5] = &(R=1.e10);

	// Set the number of terminals 
	setNumTerms(5);
	
	// Set flags for multi reference local terms
	setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN);
}

void TransformerCT::init() throw(string&)
{	//add 2 extra rows&columns to MNA matrix			
	my_row.resize(2);
}

unsigned TransformerCT::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
	for (int i=0; i < 2; i++)
   	my_row[i] = eqn_number + i;
	return 2;	// Add 2 extra RCs 
}

void TransformerCT::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
	assert(my_row[0]);
	first_eqn = my_row[0];
	n_rows = 2;
}

void TransformerCT::getLocalRefIdx(UnsignedVector& local_ref_vec,
TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));	// Input Terminal
  term_list.push_back(getTerminal(1));	// local reference terminal 
  term_list.push_back(getTerminal(2));	// Output1 Terminal
  term_list.push_back(getTerminal(3)); 	// Local reference terminal 
  term_list.push_back(getTerminal(4)); 	// Output2 Terminal

  local_ref_vec.push_back(1); // Local reference index
  local_ref_vec.push_back(3); // Local reference index
}

void TransformerCT::fillMNAM(FreqMNAM* mnam)
{
  // Ask my terminals the row/column numbers 				
	const unsigned term1=getTerminal(0)->getRC();
	const unsigned term2=getTerminal(1)->getRC();
	const unsigned term3=getTerminal(2)->getRC();
	const unsigned term4=getTerminal(3)->getRC();
	const unsigned term5=getTerminal(4)->getRC();

	double_complex cxg(1./R,0);
	double t1=(n2/n1)*k12;
	double t2=(n3/n1)*k13;

	// There is no need to set the conductances. 
	// Scaling provides a better condition number
	// for analysis run.
	
	assert(my_row[0]);
	//set the extra current vector									
	mnam->setElement(term1, my_row[0], cxg);
	mnam->setElement(term2, my_row[0], -cxg);
	mnam->setElement(term3, my_row[0], -cxg/t1);
	mnam->setElement(term4, my_row[0], cxg/t1);
	mnam->setElement(term5, my_row[0], 0);

	assert(my_row[1]);
	//set the extra current vector
	mnam->setElement(term1, my_row[1], cxg);
	mnam->setElement(term2, my_row[1], -cxg);
	mnam->setElement(term3, my_row[1], 0);
	mnam->setElement(term4, my_row[1], -cxg/t2);
	mnam->setElement(term5, my_row[1], cxg/t2);

	//set the extra voltage vectors			
	mnam->setElement(my_row[0], term1, cxg*t1);
	mnam->setElement(my_row[0], term2, -cxg*t1);
	mnam->setElement(my_row[0], term3, -cxg);
	mnam->setElement(my_row[0], term4, cxg);
	mnam->setElement(my_row[0], term5, 0);
	mnam->setElement(my_row[1], term1, cxg*t2);
	mnam->setElement(my_row[1], term2, -cxg*t2);
	mnam->setElement(my_row[1], term3, 0);
	mnam->setElement(my_row[1], term4, -cxg);
	mnam->setElement(my_row[1], term5, cxg);
}

void TransformerCT::fillMNAM(TimeMNAM* mnam)
{
  	// Ask my terminals the row numbers				
  	const unsigned term1=getTerminal(0)->getRC();
 	const unsigned term2=getTerminal(1)->getRC();
  	const unsigned term3=getTerminal(2)->getRC();
  	const unsigned term4=getTerminal(3)->getRC();
  	const unsigned term5=getTerminal(4)->getRC();

	// The linear state variable is the current through the 
	// primary coil. "g" will scale this current, so that the enteries
	// in MNAM are typical of those from other elements in the circuit.
	double g = 1/R;

	// t=turn ratio times coupling coefficient
	double t1=(n2/n1)*k12;
	double t2=(n3/n1)*k13;

	//set the primary coil conductance
	mnam->setMElement(term1, term1, g);
	mnam->setMElement(term1, term2, -g);
	mnam->setMElement(term2, term1, -g);
	mnam->setMElement(term2, term2, g);

	//set the 1st secondary coil conductance
	mnam->setMElement(term3, term3, g);
	mnam->setMElement(term3, term4, -g);
	mnam->setMElement(term4, term3, -g);
	mnam->setMElement(term4, term4, g);

	//set the 2nd secondary coil conductance
	mnam->setMElement(term4, term4, g);
	mnam->setMElement(term4, term5, -g);
	mnam->setMElement(term5, term4, -g);
	mnam->setMElement(term5, term5, g);

	assert(my_row[0]);
  	//set the extra current vector					
	mnam->setMElement(term1, my_row[0], g);
	mnam->setMElement(term2, my_row[0], -g);
	mnam->setMElement(term3, my_row[0], -g/t1);
	mnam->setMElement(term4, my_row[0], g/t1);
	mnam->setMElement(term5, my_row[0], 0);

	assert(my_row[1]);
	//set the extra current vector
	mnam->setMElement(term1, my_row[1], g);
	mnam->setMElement(term2, my_row[1], -g);
	mnam->setMElement(term3, my_row[1], 0);
	mnam->setMElement(term4, my_row[1], -g/t2);
	mnam->setMElement(term5, my_row[1], g/t2);

	//set the extra voltage vectors			
	mnam->setMElement(my_row[0], term1, g*t1);
	mnam->setMElement(my_row[0], term2, -g*t1);
	mnam->setMElement(my_row[0], term3, -g);
	mnam->setMElement(my_row[0], term4, g);
	mnam->setMElement(my_row[0], term5, 0);
	mnam->setMElement(my_row[1], term1, g*t2);
	mnam->setMElement(my_row[1], term2, -g*t2);
	mnam->setMElement(my_row[1], term3, 0);
	mnam->setMElement(my_row[1], term4, -g);
	mnam->setMElement(my_row[1], term5, g);
}
