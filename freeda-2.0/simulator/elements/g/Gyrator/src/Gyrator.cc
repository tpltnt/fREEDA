#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Gyrator.h"

// Static members
const unsigned Gyrator::n_par = 7;

// Element information
ItemInfo Gyrator::einfo =
{
  "gyrator",
  "Gyrator",
  "Don Widdows - Isac Lima - Daryl Lindsey",
  DEFAULT_ADDRESS"category:behavioral,functional",
  "2003_05_15"
};

// Parameter information
ParmInfo Gyrator::pinfo[] =
{
  {"r1", "Resistance looking into Port 1 (Ohms)", TR_DOUBLE, false},
  {"r2", "Resistance looking into Port 2 (Ohms)", TR_DOUBLE, false},
  {"l1", "Inductance looking into Port 1 (H)", TR_DOUBLE, false},
  {"l2", "Inductance looking into Port 2 (H)", TR_DOUBLE, false},
  {"x1", "Reactance looking into Port 1 (Ohms)", TR_DOUBLE, false},
  {"x2", "Reactance looking into Port 2 (Ohms)", TR_DOUBLE, false},
  {"f", "Center Frequency (Hz)", TR_DOUBLE, false},
};

Gyrator::Gyrator(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  paramvalue[0] = &(r1=0);
  paramvalue[1] = &(r2=0);
  paramvalue[2] = &(l1=0);
  paramvalue[3] = &(l2=0);
  paramvalue[4] = &(x1=0);
  paramvalue[5] = &(x2=0);
  paramvalue[6] = &(f=1e9);

	setNumTerms(4);
	// Set flags
  setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN); //must put multi-ref
	// for local ref term's!
  //Linear means it's going to check for fillMNAM functions.
  //fREEDA is written to have certain conventions about functions.
  //It has to know which functions to call from within the main body
  //  of the code when certain flags are set.
}

void Gyrator::init() throw(string&)
{
  if(r1 && !r2)
		r2=r1;
  else if (r2 && !r1)
		r1=r2;

  if(l1 && !l2)
		l2=l1;
  else if(l2 && !l1)
		l1=l2;

	if(x1 && !x2)
		x2=x1;
	else if(x2 && !x1)
		x1=x2;

	parameter_type=0;

  //----------------------
  if(l1)
		parameter_type=1;
	//define inductance  //Priority: Ind, Reactance
  else if(x1)  //if frequency is not specified, it defaults to 1 GHz.
		parameter_type=2;
	//define reactance
  else if(r1)
		parameter_type=1;
  else
	{
		r1=50;
		r2=50;
		//l1=8e-6;  //default condition if nothing on command line
		//l2=8e-6;  //for now, default to all real condition
		parameter_type=1; //define inductance
	}

  //-----------------------
	if(parameter_type==2) //reactance
	{
		l1=x1/(twopi*f);
		l2=x2/(twopi*f);
		parameter_type=1;
	} //it's now an inductor problem


  //END OF THE VARIABLE CONDITIONS
  //-------------------------------------
  //setNumberOfStates(2);

  my_row.resize(2);
}

unsigned Gyrator::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
	for (int i=0; i < 2; i++)
    my_row[i] = eqn_number + i;

	// Add 2 extra RCs
	return 2;
}

void Gyrator::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
	assert(my_row[0]);
	first_eqn = my_row[0];
	n_rows = 2;

	assert(my_row[0]);
}

void Gyrator::getLocalRefIdx(UnsignedVector& local_ref_vec,
TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1));
  term_list.push_back(getTerminal(2)); //local reference terminal
  term_list.push_back(getTerminal(3)); // Local reference terminal

  local_ref_vec.push_back(2); // Local reference index
  local_ref_vec.push_back(3); // Local reference index
}

void Gyrator::fillMNAM(FreqMNAM* mnam)
{
  // Ask my terminals the row numbers
	const unsigned term1=getTerminal(0)->getRC();
	const unsigned term2=getTerminal(1)->getRC();
	const unsigned term3=getTerminal(2)->getRC();
	const unsigned term4=getTerminal(3)->getRC();

	const double& freq(mnam->getFreq());

	//WE WILL FILL THE MATRIX USING IMPEDANCE STAMPS
	//(LOOKING INTO EACH SIDE OF THE GYRATOR, THERE IS A
	//  RESISTANCE AND AN INDUCTANCE WHICH FORM AN IMPEDANCE)
	assert(my_row[0]);

	//set the current vectors
	mnam->setElement(term1, my_row[0], 1);
	mnam->setElement(term3, my_row[0], -1);
	mnam->setElement(term2, my_row[1], 1);
	mnam->setElement(term4, my_row[1], -1);

	//set the voltage vectors
	mnam->setElement(my_row[0], term2 , 1);
	mnam->setElement(my_row[0], term4 , -1);
	assert(my_row[1]);
	mnam->setElement(my_row[1], term1 , 1);
	mnam->setElement(my_row[1], term3 , -1);

	//set the impedances
	double_complex z2(r2, twopi*freq*l2);
	mnam->setElement(my_row[0], my_row[0], z2);
	double_complex z1(r1, twopi*freq*l1);
	mnam->setElement(my_row[1], my_row[1], -z1);
}

void Gyrator::fillMNAM(TimeMNAM* mnam)
{
	//  //got to get the terminal numbers
  const unsigned term1=getTerminal(0)->getRC();
  const unsigned term2=getTerminal(1)->getRC();
  const unsigned term3=getTerminal(2)->getRC();
  const unsigned term4=getTerminal(3)->getRC();

  //WE WILL FILL THE MATRIX USING IMPEDANCE STAMPS
  //(LOOKING INTO EACH SIDE OF THE GYRATOR, THERE IS A
  //  RESISTANCE AND AN INDUCTANCE WHICH FORM AN IMPEDANCE)

	//set the current vectors
	mnam->setMElement(term1, my_row[0], 1);
	mnam->setMElement(term3, my_row[0], -1);
	mnam->setMElement(term2, my_row[1], 1);
	mnam->setMElement(term4, my_row[1], -1);

	//set the voltage vectors and impedances
	assert(my_row[0]);
	mnam->setMElement(my_row[0], term2 , 1);
	mnam->setMElement(my_row[0], term4 , -1);
	mnam->setMElement(my_row[0], my_row[0], r2);
	mnam->setMpElement(my_row[0], my_row[0], l2);

	assert(my_row[1]);
	mnam->setMElement(my_row[1], term1 , 1);
	mnam->setMElement(my_row[1], term3 , -1);
	mnam->setMElement(my_row[1], my_row[1], -r1);
	mnam->setMpElement(my_row[1], my_row[1], -l1);
}

