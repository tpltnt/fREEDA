#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "ChebyshevHPF.h"
#include <cstdio>
#include <math.h>

extern "C"
{
#include "../../../../inout/report.h"
}

// Static members
const unsigned ChebyshevHPF::n_par = 7;

// Element information
ItemInfo ChebyshevHPF::einfo =
{
  "chebyshevhpf",
  "Lumped Highpass Chebyshev Filter",
  "Michael Bucher, Michael Steer",
  "category:filter:lumped",
  "2008_05_11"
};

// Parameter information
ParmInfo ChebyshevHPF::pinfo[] =
{
  {"n", "Filter order (must be odd)", TR_INT, false},
  {"f0", "Filter corner frequency", TR_DOUBLE, true},
  {"z0", "Characteristic impedance (ohms)", TR_DOUBLE, true},
  {"ripple", "Filter ripple (in dB)", TR_DOUBLE, false},
  {"ql", "Q of inductors (at f0)", TR_DOUBLE, false},
  {"qc", "Q of capacitors (at f0)", TR_DOUBLE, false},
  {"type", "Filter Type  (Type 1 or Type 2)", TR_INT, false}
};


ChebyshevHPF::ChebyshevHPF(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  //
  paramvalue[0] = &(n = 3);
  paramvalue[1] = &(f0 = 1);
  paramvalue[2] = &(z0 = 50);
  paramvalue[3] = &(ripple = 1);
  paramvalue[4] = &(ql = 10000);
  paramvalue[5] = &(qc = 10000);
  paramvalue[6] = &(type = 2);

  // Set the number of terminals
  setNumTerms(3);

  setFlags(LINEAR | TR_FREQ_DOMAIN);
}

void ChebyshevHPF::init() throw(string&)
{
  double olda,oldb,oldg,alpha,gamma,b,g,R,C,L;

  // Clear flags so this element is not called to fill the MNAM
  setFlags(MULTI_REF);

  // Start out with a parallel capacitor
  int parallel = 1;
  double beta = log(1/tanh(ripple/17.37));
  olda = 0;
  oldg=0;
  oldb=0;

  // Do some parameter testing. This would be best done when parsing so
  // the error is identified in the netlist.

  if(n%2==0) 
   report(FATAL,"Error in use of ChebyshevHPF element, order (n) must be odd.");
  if(type!=2) 
    report(FATAL,"Error in use of ChebyshevHPF element, Only Type II filters.");
  

  // Expand nsect sections
  if(type == 2) // Now for Filter Type II
    {
    Circuit* cir = getCircuit();
    unsigned term_id1 = getTerminal(0)->getId();
    unsigned tref_id = getTerminal(1)->getId();

    for (int i=0; i<n; i++)
      {
      char loopit[10];
      sprintf(loopit, "%d", i);
  
      // Get filter coefficient
      alpha=sin((2*(i+1)-1)*pi/(2*n));
      gamma = sinh(beta/(2*n));
      b=gamma*gamma + sin((i+1)*pi/n)*sin((i+1)*pi/n);
      if(i==0)
	g = 2*alpha/gamma;
      else
	g = 4*olda*alpha/(oldb*oldg);
  
      // Highpass filter, so use parallel inductors and series caps
      if(parallel == 1)
        {
	L = (1/g)*z0/(2*pi*f0);
	// Add parallel branch (inductor with nopropagate flag)
	unsigned newelem_id = cir->addElement("l",
                     getInstanceName() + ":inductor:" + loopit, true);
	// Connect to previous terminal and reference
	cir->connect(newelem_id, term_id1);
	cir->connect(newelem_id, tref_id);
	// Get inductor pointer.
	Element* elem = cir->getElement(newelem_id);
	// Set parallel parameters
	elem->setParam("l", &L, TR_DOUBLE);
        // If the inductor has a finite Q then add internal resistor
        if(ql > 0)
          {
          // Get series resistance
          // R = (2*pi*f0*L)/ql;
	  R = g * z0 / ql;
	  elem->setParam("int_res", &R, TR_DOUBLE);
          }
  
	// Init the inductor
	elem->init();
	parallel = 0;
	}
      else
        {
	C = (1/g)/(2*pi*f0*z0);
	// Add serial branch (capacitor with nopropagate flag)
	unsigned newelem_id = cir->addElement("c",
              getInstanceName() + ":capacitor:" + loopit, true);
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
        // If the capacitor has a finite Q then add shunt resistor
        if(qc > 0)
          {
          // Get shunt resistance
          // R = qc/(2*pi*f0*C);
	  R = qc * z0/g;
	  unsigned newelem_id = cir->addElement("r",
                  getInstanceName() + ":R:" + loopit, true);
	  // Connect to previous terminal and reference
	  cir->connect(newelem_id, term_id1);
	  cir->connect(newelem_id, term_id2);
	  // Get resistor pointer.
	  Element* elem = cir->getElement(newelem_id);
	  // Set parameter
	  elem->setParam("r", &R, TR_DOUBLE);
	  // Init the resistor
	  elem->init();
          }
  
	parallel = 1;
	// Prepare for next loop
	term_id1 = term_id2;
	}

      olda=alpha;
      oldb=b;
      oldg=g;
      }
    }

}

void ChebyshevHPF::getLocalRefIdx(UnsignedVector& local_ref_vec,
			    TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1)); // Local reference terminal
  term_list.push_back(getTerminal(2));

  local_ref_vec.push_back(1); // Local reference index
}

void ChebyshevHPF::fillMNAM(FreqMNAM* mnam)
{
  // If nsect is set, there is nothing to do
  return;
}

void ChebyshevHPF::fillMNAM(TimeMNAM* mnam)
{
  // The work is done by the expanded circuit
  return;
}

