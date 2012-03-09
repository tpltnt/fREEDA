#include "../../../../network/ElementManager.h"
#include "../../../../network/CircuitManager.h"
#include "../../../../network/Element.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "XOR.h"
#include <cstdio>

// Set the number of parameters which is 5 for this element
const unsigned XOR::n_par = 5;

// Element information
ItemInfo XOR::einfo =
{
  "xor",
  "XOR Gate based on NAND Gate",
  "Shivam Priyadarshi",
  "category: Digital Logic",
  "2009_04_20"
};

ParmInfo XOR::pinfo[] =
{
        {"ln","Channel Width of NMOS (m)", TR_DOUBLE, false},
        {"wn","Channel Length of NMOS (m)", TR_DOUBLE, false},
        {"lp","Channel Width of PMOS (m)", TR_DOUBLE, false},
        {"wp","Channel Length of PMOS (m)", TR_DOUBLE, false},
        {"cshunt","Intermediate Nodal Capacitance(F)", TR_DOUBLE, false},

 /* Later it need to be modified. The input parameters should come as nmodel and pmodel for nmos and pmos. char* nmodel and char *pmodel type of construct should be used*/

};


XOR::XOR(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
        // Set default parameter values
        paramvalue[0] = &(ln = 1.0e-6);
        paramvalue[1] = &(wn = 1.0e-6);
        paramvalue[2] = &(lp = 1.0e-6);
        paramvalue[3] = &(wp = 1.0e-6);
        paramvalue[4] = &(cshunt = 1.0e-12);

  // Set the number of terminals which for this model is 5
  setNumTerms(5);

  // Set flags
  // The flags indicate
  // Non Linear = This is a non linear element
  // ONE_REF = There is only one reference terminal
  // TR_TIME_DOMAIN = Element works in  the time domain
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void XOR::init() throw(string&)
{
   // Clear flags so this element is not called to fill the MNAM
   setFlags(ONE_REF);

  Circuit* cir = getCircuit();

  unsigned term_id1 = getTerminal(0)->getId();  // VDD
  unsigned term_id2 = getTerminal(1)->getId();  // IN1
  unsigned term_id3 = getTerminal(2)->getId();  // IN2
  unsigned term_id4 = getTerminal(3)->getId();  // OUT
  unsigned term_id5 = getTerminal(4)->getId();  // GND

     // Add Nand Gate 1 : At the Input
  unsigned newelem_id = cir->addElement("cmos2nand", getInstanceName() + ":Nand1:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id2);   // IN1  Connect
      cir->connect(newelem_id, term_id3);   // IN2 Connect
       // Add an internal Terminal
     unsigned  term_id_i1 = cir->addTerminal(getInstanceName() + ":5", true);
      cir->connect(newelem_id, term_id_i1);   // Output
      cir->connect(newelem_id, term_id5);    // Ground Conenct
      // Get Nand pointer.
      Element* elem = cir->getElement(newelem_id);
       // Set Nand Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the Nand
      elem->init();

   // Add  Capacitance at Internal Node 1
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C1:", true);
      // Connect to Internal Terminal 1 and reference
      cir->connect(newelem_id, term_id_i1);
      cir->connect(newelem_id, term_id5);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &cshunt, TR_DOUBLE);
      // Init the capacitor
      elem->init();
     // Add Nand Gate 2

      newelem_id =
      cir->addElement("cmos2nand", getInstanceName() + ":Nand2:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id2);   // IN1  Connect
      cir->connect(newelem_id, term_id_i1);   // IN2 Connect
       // Add an internal Terminal
     unsigned  term_id_i2 = cir->addTerminal(getInstanceName() + ":6", true);
      cir->connect(newelem_id, term_id_i2);   // Output
      cir->connect(newelem_id, term_id5);    // Ground Conenct
      // Get Nand pointer.
       elem = cir->getElement(newelem_id);
       // Set Nand Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the Nand
      elem->init();

   // Add  Capacitance at Internal Node 2
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C2:", true);
      // Connect to Internal Terminal 6 and reference
      cir->connect(newelem_id, term_id_i2);
      cir->connect(newelem_id, term_id5);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &cshunt, TR_DOUBLE);
      // Init the capacitor
      elem->init();

  // Add Nand Gate 3
      newelem_id =
      cir->addElement("cmos2nand", getInstanceName() + ":Nand3:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id3);   // In1
      cir->connect(newelem_id, term_id_i1);   // In1
      // Add Internal Terminal 3
      unsigned  term_id_i3 = cir->addTerminal(getInstanceName() + ":7", true);
      cir->connect(newelem_id, term_id_i3);   // Output
      cir->connect(newelem_id, term_id5);    // Ground Conenct
      // Get Nand pointer.
       elem = cir->getElement(newelem_id);
       // Set Nand Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the Nand
      elem->init();

// Add  Capacitance at Internal Node 3
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C3:", true);
      // Connect to Internal Terminal 5 and reference
      cir->connect(newelem_id, term_id_i3);
      cir->connect(newelem_id, term_id5);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &cshunt, TR_DOUBLE);
      // Init the capacitor
      elem->init();

  // Add  Nand Gate 4

  newelem_id =
      cir->addElement("cmos2nand", getInstanceName() + ":Nand5:", true);
      // Connect to terminals
      cir->connect(newelem_id, term_id1);   //Vdd Connect
      cir->connect(newelem_id, term_id_i2);   // In1
      cir->connect(newelem_id, term_id_i3);   // In2
      cir->connect(newelem_id, term_id4);   // Final Output
      cir->connect(newelem_id, term_id5);    // Ground Conenct
      // Get Nand pointer.
       elem = cir->getElement(newelem_id);
       // Set Nand Parameter Values
       elem->setParam("ln",&ln, TR_DOUBLE);
       elem->setParam("wn",&wn, TR_DOUBLE);
       elem->setParam("lp",&lp, TR_DOUBLE);
       elem->setParam("wp",&wp, TR_DOUBLE);
      // Init the Nand
      elem->init();

     //    Add Capacitance at node output
      newelem_id =
      cir->addElement("capacitor", getInstanceName() + ":C4:", true);
      // Connect to internal node 9 and reference
      cir->connect(newelem_id, term_id4);
      cir->connect(newelem_id, term_id5);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set capacitor value
      elem->setParam("c", &cshunt, TR_DOUBLE);
      // Init the capacitor
     elem->init();
}
