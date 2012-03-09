#include "../../../../network/ElementManager.h"
#include "../../../../network/Circuit.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "K.h"

// Static members
const unsigned K::n_par = 4;
const double K::factor = 1e-2;

// Element information
ItemInfo K::einfo =
{
  "k",
  "Transformer",
  "Wei Zheng",
  "elements/K.h.html",
  "2003_05_15"
};

// Parameter information
ParmInfo K::pinfo[] =
{
  {"l1", "Inductor 1", TR_STRING, true},
  {"l2", "Inductor 2", TR_STRING, true},
  {"coupling", "Inductance value (H)", TR_DOUBLE, true},
  {"time_d", "Flag, if true, calculate in the time domain.", TR_BOOLEAN, false}
};


K::K(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Value of ind is required
  paramvalue[0] = &l1;
  paramvalue[1] = &l2;
  paramvalue[2] = &coupling;
  paramvalue[3] = &(time_d = false);

  // Set the number of terminals
  setNumTerms(0);
  // Set number of states
  // setNumberOfStates(0);
  // Set flags
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);

  my_row1 = 0;
  my_row2 = 0;
}

void K::init() throw(string&)
{
  if (time_d)
    setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN);

  ind1 = (L*)(getCircuit()->getElement(l1));
  ind2 = (L*)(getCircuit()->getElement(l2));
  if (ind1==NULL || ind2==NULL)
  {
		// inductor not defined, terminate program
  	throw(string("Inductor not defined: ") + *(string*)(paramvalue[0]) + string(" or ") + *(string*)(paramvalue[1]));
  }
  else
	{
		inductance1 = ind1->getInd();
		inductance2 = ind2->getInd();
		M = coupling*sqrt(inductance1*inductance2);
	}
}

void K::fillMNAM(FreqMNAM* mnam)
{
  const double& freq(mnam->getFreq());
  double_complex zl(0, -twopi * freq * M);
  mnam->setElement(ind1->getMyrow(), ind2->getMyrow(), zl);
  mnam->setElement(ind2->getMyrow(), ind1->getMyrow(), zl);
}

void K::fillMNAM(TimeMNAM* mnam)
{
  //assert(my_row);
  mnam->setMpElement(ind1->getMyrow(), ind2->getMyrow(), -M);
  mnam->setMpElement(ind2->getMyrow(), ind1->getMyrow(), -M);
}

