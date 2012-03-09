#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "ButterworthBPF.h"
#include <cstdio>

//---------------------------------------------------------------------
// This is an analog nth-order type-2 butterworth bandpass filter model
//
//         --------------
//         |   order,   |
//  n1 o---|   fc, bw,  |---o n2
//         |     z0     |
//         |            |
//         --------------
//                |
//                |
//                o n3
//
// by Shawn D. Evans
//
// This model designs a type-2 Cauer topology butterworth bandpass
// filter, using parallel inductors and capacitors and series
// inductors and capacitors in a ladder structure.  The filter can be
// of any order up to a maximum order of 100.  The filter has a center
// frequency, fc, and a bandwidth, bw. The filter is designed so that
// the input and output impedance, z0, of the filter is the same.  
//
// The native fREEDA netlist format is
// Butterworthbpf:1 n1 n2 n3 order=? fc=? bw=? z0=?
// where:
// Butterworthbpf indicates the model to use
// Butterworthbpf:1 is the name of this instance
// n1 is the input terminal (Can be a string or an integer)
// n2 is the output terminal (Can be a string or an integer)
// n3 is the reference terminal (Can be a string or an integer)
// order is the filter order
// fc is the center frequency of the filter in Hz
// bw is the bandwidth of the filter in Hz
// z0 is the input/output impedance in Ohms
//----------------------------------------------------------------------

// Static members
// Set the number of parameters which is 4 for this element
const unsigned ButterworthBPF::n_par = 4;

// Element information
ItemInfo ButterworthBPF::einfo =
{
  "butterworthbpf",
  "Butterworth bandpass filter (BPF)",
  "Shawn D. Evans",
  "category:filter:bpf",
  "2008_04_21"
};

// Parameter information
// The size of pinfa[]o is automatically assigned here.
// The information here enables automatic documentation.
// So here are the elements of the "order" parameter
// "order"  the name of the first parameter. This is the name that must
//      appear in the netlist such as "order = 3"
// "Filter Order" This is the string used in documentation. It is
//      a fuller name of the element with units.
// "TR_INT" This indicates the type of the element. In this case it is
//      a int. The full list of valid input parameter types
//      is given in .../network/NetListItem.h
// "true" indicates that the parameter is required and must be specified in
//      the netlist. The alternative is "false" indicating that the parameter
//      is optional in which case a default value must be sepcified.
ParmInfo ButterworthBPF::pinfo[] =
{
  {"order", "Filter Order", TR_INT, true},
  {"fc", "Center Frequency (Hz)", TR_DOUBLE, true},
  {"bw", "Bandwidth (Hz)", TR_DOUBLE, true},
  {"z0", "Input/Output Impedance (Ohms)", TR_DOUBLE, true}
};


ButterworthBPF::ButterworthBPF(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Now list all of the parameters and give them meaningful variable names
  // to use in the codea. These names are the same names used in the header
  // file.
  paramvalue[0] = &order;
  paramvalue[1] = &fc;
  paramvalue[2] = &bw;
  paramvalue[3] = &z0;
  // If order had a default value of 3 we would have used:
  // paramvalue[0] = &(order = 3);

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

void ButterworthBPF::init() throw(string&)
{
  wc = fc * 2.0 * pi;      // center frequency in radians/sec
  bw_rad = bw * 2.0 * pi;  // bandwidth in radians/sec

  // Initialize, calculate, and display all of the L and C values
  cout << endl << endl;
  cout << "----------------------------------------" << endl;
  cout << "Butterworth Bandpass Filter Elements" << endl;
  for (int h=1; h<(order+1); h++)
  {
    // calculate coefficients of a butterworth lowpass prototype filter
    // normalized to a radian corner frequency of 1 radian/sec and a 
    // 1 Ohm system impedance
    g[h]=2*sin(((2.0*h-1)*pi)/(2.0*order));
    
    if (h % 2 == 1) // calculate parallel L and C network values
    {
      C[h]=g[h]/(bw_rad*z0);
      L[h]=(bw_rad*z0)/(g[h]*wc*wc);
      cout << "C" << h << " = " << C[h] << " F" << endl;
      cout << "L" << h << " = " << L[h] << " H" << endl;
    }
    else // calculate series L and C network values
    {
      C[h]=bw_rad/(wc*wc*g[h]*z0);
      L[h]=(g[h]*z0)/bw_rad;
      cout << "C" << h << " = " << C[h] << " F" << endl;
      cout << "L" << h << " = " << L[h] << " H" << endl;
    }
  }
  cout << "----------------------------------------" << endl;

  // Clear flags so this element is not called to fill the MNAM
  setFlags(ONE_REF);
  Circuit* cir = getCircuit();
  unsigned term_id1 = getTerminal(0)->getId();
  unsigned term_id2 = getTerminal(1)->getId();
  unsigned tref_id = getTerminal(2)->getId();
  for (int i=1; i<(order+1); i++)
  {
    char loopit[10];
    sprintf(loopit, "%d", i);
    
    if (i % 2 == 1) // Add parallel L and C
    {  
      // Add parallel inductor
      unsigned newelem_id =
      cir->addElement("inductor", getInstanceName() + ":inductor:" + loopit, true);
      // Connect to previous terminal
      cir->connect(newelem_id, term_id1);
      // Connect to tref_id
      cir->connect(newelem_id, tref_id);
      // Get inductor pointer.
      Element* elem = cir->getElement(newelem_id);
      // Set inductor value
      elem->setParam("l", &L[i], TR_DOUBLE);
      // Init the inductor
      elem->init();

      // Add parallel capacitor
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":capacitor:" + loopit, true);
      // Connect to previous terminal and reference
      cir->connect(newelem_id, term_id1);
      cir->connect(newelem_id, tref_id);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C[i], TR_DOUBLE);
      // Init the capacitor
      elem->init();
    } 
    else // Add series L and C
    {
      // Add series inductor
      unsigned newelem_id =
      cir->addElement("inductor", getInstanceName() + ":inductor:" + loopit, true);
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
      // Init the inductor
      elem->init();

      // Add series capacitor
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":capacitor:" + loopit, true);
      // Connect to previous internal terminal
      cir->connect(newelem_id, term_id2);
      // determine if the capacitor should be connected to the output or not
      if (i == order || i == (order-1))
        term_id2 = getTerminal(1)->getId();
      else
      {
        sprintf(loopit, "%d", (i+1));
        term_id2 = cir->addTerminal(getInstanceName() + ":" + loopit, true);
      }
      // connect to the next terminal
      cir->connect(newelem_id, term_id2);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C[i], TR_DOUBLE);
      // Init the capacitor
      elem->init();

      // Prepare for next loop
      term_id1 = term_id2;
    }

    if (order == 1) // Add resistor to short input terminal to output terminal
    {  
      // add resistor
      unsigned newelem_id =
      cir->addElement("resistor", getInstanceName() + ":resistor:" + loopit, true);
      // Connect to input terminal
      cir->connect(newelem_id, term_id1);
      // Connect to output terminal
      cir->connect(newelem_id, term_id2);
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

void ButterworthBPF::getLocalRefIdx(UnsignedVector& local_ref_vec,
			    TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1));
  term_list.push_back(getTerminal(2)); // Local reference terminal

  local_ref_vec.push_back(2); // Local reference index
}

void ButterworthBPF::fillMNAM(FreqMNAM* mnam)
{
  // If nsect is set, there is nothing to do
  if (order)
    return;
}

void ButterworthBPF::fillMNAM(TimeMNAM* mnam)
{
  // The work is done by the expanded circuit
  // Set order to some value
  assert(order);
}

