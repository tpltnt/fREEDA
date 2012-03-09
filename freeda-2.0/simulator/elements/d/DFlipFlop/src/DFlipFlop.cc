#include "../../../../network/ElementManager.h"
#include "../../../../network/CircuitManager.h"
#include "../../../../network/Element.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "DFlipFlop.h"
#include <cstdio>

// Set the number of parameters which is 4 for this element
const unsigned DFlipFlop::n_par = 4;

// Element information
ItemInfo DFlipFlop::einfo =
{
  "dflipflop",
  "Master Slave Negative Edge Triggered D-Flip Flop",
  "Shivam Priyadarshi",
  "category: Digital Logic",
  "2009_04_20"
};

ParmInfo DFlipFlop::pinfo[] =
{
        {"ln","Channel Width of NMOS (m)", TR_DOUBLE, false},
        {"wn","Channel Length of NMOS (m)", TR_DOUBLE, false},
        {"lp","Channel Width of PMOS (m)", TR_DOUBLE, false},
        {"wp","Channel Length of PMOS (m)", TR_DOUBLE, false},
 /* Later it need to be modified. The input parameters should come as nmodel and pmodel for nmos and pmos. char* nmodel and char *pmodel type of construct should be used*/

};


DFlipFlop::DFlipFlop(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
        // Set default parameter values
        paramvalue[0] = &(ln = 1.0e-6);
        paramvalue[1] = &(wn = 1.0e-6);
        paramvalue[2] = &(lp = 1.0e-6);
        paramvalue[3] = &(wp = 1.0e-6);

  // Set the number of terminals which for this model is 6
  setNumTerms(6);

  // Set flags
  // The flags indicate
  // Non Linear = This is a non linear element
  // ONE_REF = There is only one reference terminal
  // TR_TIME_DOMAIN = Element works in  the time domain
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void DFlipFlop::init() throw(string&)
{
   C = 1e-12;
   Rs = 1;
   // Clear flags so this element is not called to fill the MNAM
   setFlags(ONE_REF);

  Circuit* cir = getCircuit();

  unsigned term_id1 = getTerminal(0)->getId();  // Vdd
  unsigned term_id2 = getTerminal(1)->getId();  // D
  unsigned term_id3 = getTerminal(2)->getId();  // CLK
  unsigned term_id4 = getTerminal(3)->getId();  // Q
  unsigned term_id5 = getTerminal(4)->getId();  // Qb
  unsigned tref_id = getTerminal(5)->getId();   // Gnd 

     // Add Nand Gate 1 : At the Input 

     unsigned  newelem_id =
      cir->addElement("cmos2nand", getInstanceName() + ":Nand1:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id2);   // D  Connect
      cir->connect(newelem_id, term_id3);   // CLK Connect
       // Add an internal Terminal
     unsigned  term_id_i4 = cir->addTerminal(getInstanceName() + ":4", true);     
      cir->connect(newelem_id, term_id_i4);   // Output
      cir->connect(newelem_id, tref_id);    // Ground Conenct
      // Get Nand pointer.
      Element*  elem = cir->getElement(newelem_id);
       // Set Nand Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the Nand
      elem->init();
       
   // Add  Capacitance at Internal Node 4 
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C4:", true);
      // Connect to Internal Terminal 4 and reference
      cir->connect(newelem_id, term_id_i4);
      cir->connect(newelem_id, tref_id);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C, TR_DOUBLE);
      // Init the capacitor
      elem->init();

     // Add Nand Gate 3

      newelem_id =
      cir->addElement("cmos2nand", getInstanceName() + ":Nand3:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id_i4);   // Terminal 4  Connect
      cir->connect(newelem_id, term_id3);   // CLK Connect
       // Add an internal Terminal
     unsigned  term_id_i6 = cir->addTerminal(getInstanceName() + ":6", true);
      cir->connect(newelem_id, term_id_i6);   // Output
      cir->connect(newelem_id, tref_id);    // Ground Conenct
      // Get Nand pointer.
       elem = cir->getElement(newelem_id);
       // Set Nand Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the Nand
      elem->init();

   // Add  Capacitance at Internal Node 6 
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C6:", true);
      // Connect to Internal Terminal 6 and reference
      cir->connect(newelem_id, term_id_i6);
      cir->connect(newelem_id, tref_id);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C, TR_DOUBLE);
      // Init the capacitor
      elem->init();

  // Add Nand Gate 2
      newelem_id =
      cir->addElement("cmos2nand", getInstanceName() + ":Nand2:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id_i4);   // In1
      // Add Internal Terminal 5,7,8
      unsigned  term_id_i5 = cir->addTerminal(getInstanceName() + ":5", true);
      unsigned  term_id_i7 = cir->addTerminal(getInstanceName() + ":7", true);
      unsigned  term_id_i8 = cir->addTerminal(getInstanceName() + ":8", true);

      cir->connect(newelem_id, term_id_i7);   // In2
      cir->connect(newelem_id, term_id_i5);   // Output
      cir->connect(newelem_id, tref_id);    // Ground Conenct
      // Get Nand pointer.
       elem = cir->getElement(newelem_id);
       // Set Nand Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the Nand
      elem->init();

      // Add Nand Gate 4
      newelem_id =
      cir->addElement("cmos2nand", getInstanceName() + ":Nand4:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id_i5);   // In1
      cir->connect(newelem_id, term_id_i6);   // In2
      cir->connect(newelem_id, term_id_i7);   // Output
      cir->connect(newelem_id, tref_id);    // Ground Conenct
      // Get Nand pointer.
       elem = cir->getElement(newelem_id);
       // Set Nand Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the Nand
      elem->init();

// Add  Capacitance at Internal Node 5 
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C5:", true);
      // Connect to Internal Terminal 5 and reference
      cir->connect(newelem_id, term_id_i5);
      cir->connect(newelem_id, tref_id);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C, TR_DOUBLE);
      // Init the capacitor
      elem->init();

// Add  Capacitance at Internal Node 7 
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C7:", true);
      // Connect to Internal Terminal 7 and reference
      cir->connect(newelem_id, term_id_i7);
      cir->connect(newelem_id, tref_id);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C, TR_DOUBLE);
      // Init the capacitor
      elem->init();

   // Add Inverter
      newelem_id =
      cir->addElement("cmosinv", getInstanceName() + ":Inv1:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id3);   // CLK
      cir->connect(newelem_id, term_id_i8);   // CLKB
      cir->connect(newelem_id, tref_id);    // Ground Conenct
      // Get INV  pointer.
       elem = cir->getElement(newelem_id);
       // Set INV Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the INV
      elem->init();

// Add  Capacitance at Internal Node 8 
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C8:", true);
      // Connect to Internal Terminal 8 and reference
      cir->connect(newelem_id, term_id_i8);
      cir->connect(newelem_id, tref_id);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C, TR_DOUBLE);
      // Init the capacitor
      elem->init();

// Add Internal Node 9, 10, 

      unsigned  term_id_i9 = cir->addTerminal(getInstanceName() + ":9", true);
      unsigned  term_id_i10 = cir->addTerminal(getInstanceName() + ":10", true);

  // Add  Nand Gate 5

  newelem_id =
      cir->addElement("cmos2nand", getInstanceName() + ":Nand5:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id_i8);   // In1
      cir->connect(newelem_id, term_id_i5);   // In2
      cir->connect(newelem_id, term_id_i9);   // Output
      cir->connect(newelem_id, tref_id);    // Ground Conenct
      // Get Nand pointer.
       elem = cir->getElement(newelem_id);
       // Set Nand Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the Nand
      elem->init();

 // Add  Nand Gate 6
  newelem_id =
      cir->addElement("cmos2nand", getInstanceName() + ":Nand6:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id_i8);   // In1
      cir->connect(newelem_id, term_id_i7);   // In2
      cir->connect(newelem_id, term_id_i10);   // Output
      cir->connect(newelem_id, tref_id);    // Ground Conenct
      // Get Nand pointer.
       elem = cir->getElement(newelem_id);
       // Set Nand Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the Nand
      elem->init();

     //    Add Capacitance at node 9
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C9:", true);
      // Connect to internal node 9 and reference
      cir->connect(newelem_id, term_id_i9);
      cir->connect(newelem_id, tref_id);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C, TR_DOUBLE);
      // Init the capacitor
      elem->init();

        //    Add Capacitance at node 10
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C10:", true);
      // Connect to internal node 10 and reference
      cir->connect(newelem_id, term_id_i10);
      cir->connect(newelem_id, tref_id);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C, TR_DOUBLE);
      // Init the capacitor
      elem->init();

// Add  Nand Gate 7
  newelem_id =
      cir->addElement("cmos2nand", getInstanceName() + ":Nand7:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id_i9);   // In1
      cir->connect(newelem_id, term_id5);   // In2
      cir->connect(newelem_id, term_id4);   // Q
      cir->connect(newelem_id, tref_id);    // Ground Conenct
      // Get Nand pointer.
       elem = cir->getElement(newelem_id);
       // Set Nand Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the Nand
      elem->init();


// Add  Nand Gate 8
  newelem_id =
      cir->addElement("cmos2nand", getInstanceName() + ":Nand8:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id_i10);   // In1
      cir->connect(newelem_id, term_id4);   // In2
      cir->connect(newelem_id, term_id5);   // Qb
      cir->connect(newelem_id, tref_id);    // Ground Conenct
      // Get Nand pointer.
       elem = cir->getElement(newelem_id);
       // Set Nand Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the Nand
      elem->init();


    //    Add Load Capacitance to output :Q
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C11:", true);
      // Connect to Q and reference
      cir->connect(newelem_id, term_id4);
      cir->connect(newelem_id, tref_id);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C, TR_DOUBLE);
      // Init the capacitor
      elem->init();

    //    Add Load Capacitance to output :Qb
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C12:", true);
      // Connect to internal terminal and reference
      cir->connect(newelem_id, term_id5);
      cir->connect(newelem_id, tref_id);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &C, TR_DOUBLE);
      // Init the capacitor
      elem->init();

}
