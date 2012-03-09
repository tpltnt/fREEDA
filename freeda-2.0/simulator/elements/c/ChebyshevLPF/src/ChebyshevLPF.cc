#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "ChebyshevLPF.h"
#include <cstdio>
#include <math.h>

// Static members
const unsigned ChebyshevLPF::n_par = 5;

// Element information
ItemInfo ChebyshevLPF::einfo =
{
  "chebyshevlpf",
  "Lumped Lowpass Chebyshev Filter",
  "Michael Bucher",
  "category:filter",
  "2008_04_21"
};

// Parameter information
ParmInfo ChebyshevLPF::pinfo[] =
{
  {"n", "Filter order", TR_INT, false},
  {"f0", "Filter corner frequency", TR_DOUBLE, true},
  {"z0", "Characteristic impedance (ohms)", TR_DOUBLE, true},
  {"type", "Type of filter, lowpass(0) or highpass(1)", TR_INT, true},
  {"ripple", "Filter ripple (in dB)", TR_DOUBLE, false}
};


ChebyshevLPF::ChebyshevLPF(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  //
  paramvalue[0] = &(n = 3);
  paramvalue[1] = &(f0);
  paramvalue[2] = &(z0);
  paramvalue[3] = &(type);
  paramvalue[4] = &(ripple = 1);


  // Set the number of terminals
  setNumTerms(3);

  setFlags(LINEAR | TR_FREQ_DOMAIN);
}

void ChebyshevLPF::init() throw(string&)
{
  // Clear flags so this element is not called to fill the MNAM
  setFlags(MULTI_REF);

  // Start out with a parallel capacitor
  int parallel = 1;
  double beta = log(1/tanh(ripple/17.37));
  olda = 0;
  oldg=0;
  oldb=0;
  double R = 1.0;
  // Expand nsect sections
  Circuit* cir = getCircuit();
  unsigned term_id1 = getTerminal(0)->getId();
  unsigned tref_id = getTerminal(1)->getId();
  for (int i=0; i<n; i++)
    {
      char loopit[10];
      sprintf(loopit, "%d", i);
      alpha=sin((2*(i+1)-1)*pi/(2*n));
      gamma = sinh(beta/(2*n));
      b=gamma*gamma + sin((i+1)*pi/n)*sin((i+1)*pi/n);
      if(i==0) {
	g = 2*alpha/gamma;
      }
      else {
	g = 4*olda*alpha/(oldb*oldg);
      }
      if(type == 0) 
	{
	  if(parallel == 1)
	    {
	      C = g/(2*pi*f0*z0);
	      // Add parallel branch (capacitor with nopropagate flag)
	      unsigned newelem_id =
		cir->addElement("capacitor", getInstanceName() + ":capacitor:" + loopit, true);
	      // Connect to previous terminal and reference
	      cir->connect(newelem_id, term_id1);
	      cir->connect(newelem_id, tref_id);
	      // Get capacitor pointer.
	      Element* elem = cir->getElement(newelem_id);
	      // Set parallel parameters
	      elem->setParam("c", &C, TR_DOUBLE);


	      // Init the capacitor
	      elem->init();
	      parallel = 0;
	    }
	  else {
	    L = g*z0/(2*pi*f0);
	    // Add serial branch (inductor with nopropagate flag)
	    unsigned newelem_id =
	      cir->addElement("inductor", getInstanceName() + ":inductor:" + loopit, true);
	    // Connect to previous terminal
	    cir->connect(newelem_id, term_id1);
	    // create new terminal (with nocheck) unless this is the last one
	    unsigned term_id2;
	    if (i != n-1 && i != n-2)
	      {
		term_id2 = cir->addTerminal(getInstanceName() + ":" + loopit, true);
	      }
	    else
	      {
		term_id2 = getTerminal(2)->getId();
	      }
	    // Connect to term_id2
	    cir->connect(newelem_id, term_id2);
	    // Get inductor pointer.
	    Element* elem = cir->getElement(newelem_id);
	    // Set serial parameters
	    elem->setParam("l", &L, TR_DOUBLE);

	    // Init the inductor
	    elem->init();
	    parallel = 1;
	    // Prepare for next loop
	    term_id1 = term_id2;
	  }
	}
      else // Highpass filter, so use parallel inductors and series caps
	{
	  if(parallel == 1)
	    {
	      L = (1/g)*z0/(2*pi*f0);
	      // Add parallel branch (inductor with nopropagate flag)
	      unsigned newelem_id =
		cir->addElement("inductor", getInstanceName() + ":inductor:" + loopit, true);
	      // Connect to previous terminal and reference
	      cir->connect(newelem_id, term_id1);
	      cir->connect(newelem_id, tref_id);
	      // Get inductor pointer.
	      Element* elem = cir->getElement(newelem_id);
	      // Set parallel parameters
	      elem->setParam("l", &L, TR_DOUBLE);
	      // Init the inductor
	      elem->init();
	      parallel = 0;
	    }
	  else {
	    C = (1/g)/(2*pi*f0*z0);
	    // Add serial branch (capacitor with nopropagate flag)
	    unsigned newelem_id =
	      cir->addElement("capacitor", getInstanceName() + ":capacitor:" + loopit, true);
	    // Connect to previous terminal
	    cir->connect(newelem_id, term_id1);
	    // create new terminal (with nocheck) unless this is the last one
	    unsigned term_id2;
	    if (i != n-1 && i != n-2)
	      term_id2 = cir->addTerminal(getInstanceName() + ":" + loopit, true);
	    else
	      term_id2 = getTerminal(2)->getId();
	    // Connect to term_id2
	    cir->connect(newelem_id, term_id2);
	    // Get capacitor pointer.
	    Element* elem = cir->getElement(newelem_id);
	    // Set serial parameters
	    elem->setParam("c", &C, TR_DOUBLE);
	    // Init the capacitor
	    elem->init();
	    parallel = 1;
	    // Prepare for next loop
	    term_id1 = term_id2;
	  }
	}
      olda=alpha;
      oldb=b;
      oldg=g;
    }
  if(n%2==0) 
    {
      double R = tanh(beta/4)*tanh(beta/4);
      char last[10];
      sprintf(last, "%d", n);
      // Add parallel resistor
      unsigned newelem_id =
	cir->addElement("resistor", getInstanceName() + ":resistor:" + last, true);
      // Connect to previous terminal
      cir->connect(newelem_id, term_id1);
      // Connect to tref_id
      cir->connect(newelem_id, tref_id);
      // Get resistor pointer.
      Element* elem = cir->getElement(newelem_id);
      // Set serial parameters
      elem->setParam("r", &R, TR_DOUBLE);
      // Init the resistor
      elem->init();
    }
}

void ChebyshevLPF::getLocalRefIdx(UnsignedVector& local_ref_vec,
			    TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1)); // Local reference terminal
  term_list.push_back(getTerminal(2));
  //  term_list.push_back(getTerminal(3)); // Local reference terminal

  local_ref_vec.push_back(1); // Local reference index
  //  local_ref_vec.push_back(3); // Local reference index
}


void ChebyshevLPF::fillMNAM(FreqMNAM* mnam)
{
  // If nsect is set, there is nothing to do
  return;
}

void ChebyshevLPF::fillMNAM(TimeMNAM* mnam)
{
  // The work is done by the expanded circuit
  return;
}

