#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "ChebyshevBSF.h"
#include <cstdio>

//---------------------------------------------------------------------
// This is an analog nth-order type-2 chebyshev bandstop filter model
//
//         --------------
//         |   order,   |
//  n1 o---|   f0, bw,  |---o n3
//         | z0, ripple |
//         |            |
//         --------------
//                |
//                |
//                o n2
//
// Michael Steer September 2008
// Adapted from ChebyshevBSF
//
// This model designs a type-2 Cauer topology chebyshev bandstop
// filter, using parallel inductors and capacitors and series
// inductors and capacitors in a ladder structure.  The filter can be
// of any order up to a maximum order of 100.  The filter has a center
// frequency, f0, a bandwidth, bw, and a passband ripple, ripple. The
// filter is designed so that the input and output impedance, z0, of
// the filter is the same.  
//
// The native fREEDA netlist format is
// Chebyshevbsf:1 n1 n2 n3 order=? f0=? bw=? z0=? ripple=?
// where:
// Chebyshevbsf indicates the model to use
// Chebyshevbsf:1 is the name of this instance
// n1 is the input terminal (Can be a string or an integer)
// n2 is the reference terminal (Can be a string or an integer)
// n3 is the output terminal (Can be a string or an integer)
// order is the filter order
// f0 is the center frequency of the filter in Hz
// bw is the bandwidth of the filter in Hz
// z0 is the input/output impedance in Ohms
// ripple is the passband ripple in decibels
//----------------------------------------------------------------------

// Static members
// Set the number of parameters which is 6 for this element
const unsigned ChebyshevBSF::n_par = 6;

// Element information
ItemInfo ChebyshevBSF::einfo =
{
  "chebyshevbsf",
  "Chebyshev BandStop Filter (BSF) (lumped element)",
  "Michael Steer, Shawn Evans",
  "category:filter:lumped",
  "2008_08_14"
};

// Parameter information
// The size of ParmInfo is automatically assigned here.
// The information here enables automatic documentation.
// So here are the elements of the "order" parameter
// "n"  the name of the first parameter. This is the name that must
//      appear in the netlist such as "n = 3"
// "Filter Order" This is the string used in documentation. It is
//      a fuller name of the element with units.
// "TR_INT" This indicates the type of the element. In this case it is
//      a int. The full list of valid input parameter types
//      is given in .../network/NetListItem.h
// "true" indicates that the parameter is required and must be specified in
//      the netlist. The alternative is "false" indicating that the parameter
//      is optional in which case a default value must be sepcified.
ParmInfo ChebyshevBSF::pinfo[] =
{
  {"n", "Filter Order", TR_INT, false},
  {"f0", "Center Frequency (Hz)", TR_DOUBLE, true},
  {"z0", "Input/Output Impedance (Ohms)", TR_DOUBLE, false},
  {"q", "Resonator Q", TR_DOUBLE, false},
  {"ripple", "Stopband Ripple (dB)", TR_DOUBLE, false},
  {"bw", "Bandwidth (Hz)", TR_DOUBLE, true}
};


ChebyshevBSF::ChebyshevBSF(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Now list all of the parameters and give them meaningful variable names
  // to use in the codea. These names are the same names used in the header
  // file.
  paramvalue[0] = &(n=11);
  paramvalue[1] = &f0;
  paramvalue[2] = &(z0=50);
  paramvalue[3] = &(q=10000);
  paramvalue[4] = &(ripple=0.1);
  paramvalue[5] = &bw;

  // Set the number of terminals which for this model is 3
  setNumTerms(3);

  // Set flags
  // The flags indicate
  // Linear = This is a linear element
  // ONE_REF = There is only one reference terminal
  // TR_FREQ_DOMAIN = Element works in both the time domain and
  //           frequency domain. (As all good elements should.)
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);
}

void ChebyshevBSF::init() throw(string&)
{
  double a[100], b[100], g[100], C[100], L[100], R;
  double wc = f0 * 2.0 * pi;      // center frequency in radians/sec
  double bw_rad = bw * 2.0 * pi;  // bandwidth in radians/sec
  assert(n>100); // maximum order is 1000
  assert((n+1)%2); // order must be odd

  // Initialize, calculate, and display all of the L and C values
  double beta = log(cosh(ripple/17.37)/sinh(ripple/17.37));
  double gamma = sinh(beta/(2.0*n));
  for (int h=1; h<(n+1); h++)
  {
    // calculate coefficients of a chebyshev lowpass prototype filter
    // normalized to a radian corner frequency of 1 radian/sec and a 
    // 1 Ohm system impedance
    a[h]=sin(((2.0*h-1)*pi)/(2.0*n));
    b[h]=gamma*gamma+sin(h*pi/n)*sin(h*pi/n);
    if (h == 1)
      g[1]=2*a[1]/gamma;
    else
      g[h]=4.0*a[h-1]*a[h]/(b[h-1]*g[h-1]);
    if (h % 2 == 1) // (odd) calculate series L and C network values
      {
      C[h]=g[h]*bw_rad/(wc*wc*z0);
      L[h]=z0/(bw_rad*g[h]);
      }
    else // (even) calculate parallel L and C network values
      {
      C[h]=1./(bw_rad*g[h]*z0);
      L[h]=(g[h]*bw_rad*z0)/(wc*wc);
      }
  }

  // Clear flags so this element is not called to fill the MNAM
  setFlags(ONE_REF);

  Circuit* cir = getCircuit();
  unsigned term_id1 = getTerminal(0)->getId();
  unsigned term_id2 = getTerminal(2)->getId();
  unsigned term_id_last = getTerminal(2)->getId();
  unsigned tref_id = getTerminal(1)->getId();
  for (int i=1; i<(n+1); i++)
  {
    char loopit[10];
    sprintf(loopit, "%d", i);

    if (i % 2 == 1) // (odd Add series L and C in shunt leg
    {  
      // Add series inductor in shunt leg
      unsigned newelem_id =
        cir->addElement("l", getInstanceName() + ":L:" + loopit, true);
      // Connect to previous terminal
      cir->connect(newelem_id, term_id1);
      // Add an internal terminal
      term_id2 = cir->addTerminal(getInstanceName() + ":" + loopit, true);
      // Connect to internal terminal
      cir->connect(newelem_id, term_id2);
      // Get inductor pointer.
      Element* elem = cir->getElement(newelem_id);
      // Set inductor value
      elem->setParam("l", &L[i], TR_DOUBLE);
      // If the resonator has a finite Q then add internal resistor
      if(q > 0)
        {
        R = sqrt(L[i]/C[i])/q;
        elem->setParam("int_res", &R, TR_DOUBLE);
        }
      // Init the inductor
      elem->init();

      // Add series capacitor in shunt leg
      newelem_id =
        cir->addElement("c", getInstanceName() + ":C:" + loopit, true);
      // Connect to previous terminal and reference
      cir->connect(newelem_id, term_id2);
      cir->connect(newelem_id, tref_id);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C[i], TR_DOUBLE);
      // Init the capacitor
      elem->init();
  
    } 
    else
    {
    // Add parallel L and C in series leg

      // Add parallel capacitor in series leg
      unsigned newelem_id =
        cir->addElement("c", getInstanceName() + ":C:" + loopit, true);
      // Connect to previous internal terminal
      cir->connect(newelem_id, term_id1);
      // determine if the capacitor should be connected to the output or not
      if (i == n || i == (n-1))
        term_id2 = term_id_last;
      else
        term_id2 = cir->addTerminal(getInstanceName() + ":x" + loopit, true);
      // connect to the next terminal
      cir->connect(newelem_id, term_id2);
      // Get capacitor pointer.
      Element* elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C[i], TR_DOUBLE);
      // Init the capacitor
      elem->init();

      // Add parallel inductor
      newelem_id =
      cir->addElement("l", getInstanceName() + ":L:" + loopit, true);
      // Connect to previous terminal
      cir->connect(newelem_id, term_id1);
      // Connect to next terminal
      cir->connect(newelem_id, term_id2);
      // Get inductor pointer.
      elem = cir->getElement(newelem_id);
      // Set inductor value
      elem->setParam("l", &L[i], TR_DOUBLE);
      // Init the inductor
      elem->init();

      // If the resonator has a finite Q then add shunt resistor
      if(q > 0)
        {
        R = q * sqrt(L[i]/C[i]);
	newelem_id =
          cir->addElement("r", getInstanceName() + ":R:" + loopit, true);
	// Connect to previous terminal and reference
	cir->connect(newelem_id, term_id2);
	cir->connect(newelem_id, tref_id);
	// Get resistor pointer.
	Element* elem = cir->getElement(newelem_id);
	// Set parameter
	elem->setParam("r", &R, TR_DOUBLE);
	// Init the resistor
	elem->init();
        }

      // Prepare for next loop
      term_id1 = term_id2;
    }

    // if the order is 1 add a tiny resistor to short input terminal to output
    // terminal if input and output terminals are the same
    if (n == 1) // Add resistor to short input terminal to output terminal
      {
      if(term_id1 != term_id_last)
        {  
        // add resistor
        unsigned newelem_id =
          cir->addElement("r", getInstanceName() + ":R:" + loopit, true);
        // Connect to input terminal
        cir->connect(newelem_id, term_id1);
        // Connect to output terminal
        cir->connect(newelem_id, term_id_last);
        // Get resistor pointer.
        Element* elem = cir->getElement(newelem_id);
        // Set resistor value to be a short
        double R=1e-12;
        elem->setParam("r", &R, TR_DOUBLE);
        // Init the resistor
        elem->init();
        }
      }
  }
}

void ChebyshevBSF::getLocalRefIdx(UnsignedVector& local_ref_vec,
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

void ChebyshevBSF::fillMNAM(FreqMNAM* mnam)
{
  // If nsect is set, there is nothing to do
  if (n)
    return;
}

void ChebyshevBSF::fillMNAM(TimeMNAM* mnam)
{
  // The work is done by the expanded circuit
  assert(n);
}

