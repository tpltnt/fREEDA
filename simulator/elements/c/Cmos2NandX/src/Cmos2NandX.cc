#include "../../../../network/ElementManager.h"
#include "../../../../network/CircuitManager.h"
#include "../../../../network/Element.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Cmos2NandX.h"
#include <cstdio>

// Set the number of parameters which is 4 for this element
const unsigned Cmos2NandX::n_par = 4;

// Element information
ItemInfo Cmos2NandX::einfo =
{
  "cmos2nandx",
  "Cmos Nand Gate",
  "Shivam Priyadarshi",
  "category: Logic",
  "2009_04_17"
};

ParmInfo Cmos2NandX::pinfo[] =
{
  {"ln","Channel Width of NMOS (m)", TR_DOUBLE, false},
  {"wn","Channel Length of NMOS (m)", TR_DOUBLE, false},
  {"lp","Channel Width of PMOS (m)", TR_DOUBLE, false},
  {"wp","Channel Length of PMOS (m)", TR_DOUBLE, false},
};


Cmos2NandX::Cmos2NandX(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(ln = 1.0e-6);
  paramvalue[1] = &(wn = 1.0e-6);
  paramvalue[2] = &(lp = 1.0e-6);
  paramvalue[3] = &(wp = 1.0e-6);

  // Set the number of terminals which for this model is 6
  setNumTerms(5);

  // Set flags
  // The flags indicate
  // Non Linear = This is a non linear element
  // ONE_REF = There is only one reference terminal
  // TR_TIME_DOMAIN = Element works in  the time domain
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void Cmos2NandX::init() throw(string&)
{
  C = 1e-12;
  Rs = 10;
  double uo_n = 400;
  double uo_p = 200;

  // Clear flags so this element is not called to fill the MNAM
  setFlags(ONE_REF);

  Circuit* cir = getCircuit();
  unsigned term_id1 = getTerminal(1)->getId();  // First Input
  unsigned term_id2 = getTerminal(2)->getId();  // Second Input
  unsigned term_id3 = getTerminal(3)->getId();  // Output
  unsigned term_id4 = getTerminal(0)->getId();  // Vdd
  unsigned tref_id = getTerminal(4)->getId();   // Gnd

  // Add Source Resistance 1
  unsigned newelem_id =
    cir->addElement("resistor", getInstanceName() + ":resistors1:", true);
  // Connect to input terminal
  cir->connect(newelem_id, term_id1);
  // Add an internal Terminal
  unsigned term_id_i1 = cir->addTerminal(getInstanceName() + ":1", true);
  // Connect to internal terminal
  cir->connect(newelem_id, term_id_i1);
  // Get resistor pointer.
  Element* elem = cir->getElement(newelem_id);
  // Set resistor value to be a short
  elem->setParam("r", &Rs, TR_DOUBLE);
  // Init the resistor
  elem->init();

  // Add Source Resistance 2
  newelem_id =
    cir->addElement("resistor", getInstanceName() + ":resistors2:", true);
  // Connect to input terminal
  cir->connect(newelem_id, term_id2);
  // Add an internal Terminal
  unsigned  term_id_i2 = cir->addTerminal(getInstanceName() + ":2", true);
  // Connect to internal terminal
  cir->connect(newelem_id, term_id_i2);
  // Get resistor pointer.
  elem = cir->getElement(newelem_id);
  // Set resistor value to be a short
  elem->setParam("r", &Rs, TR_DOUBLE);
  // Init the resistor
  elem->init();


  // Add Source Capacitance 1
  newelem_id =
    cir->addElement("capacitor", getInstanceName() + ":capacitors1:", true);
  // Connect to Input terminal1 and reference
  cir->connect(newelem_id, term_id_i1);
  cir->connect(newelem_id, tref_id);
  // Get capacitor pointer.
  elem = cir->getElement(newelem_id);
  // Set capacitor value
  elem->setParam("c", &C, TR_DOUBLE);
  // Init the capacitor
  elem->init();

  // Add Source Capacitance 2
  newelem_id =
    cir->addElement("capacitor", getInstanceName() + ":capacitors2:", true);
  // Connect to Input terminal2 and reference
  cir->connect(newelem_id, term_id_i2);
  cir->connect(newelem_id, tref_id);
  // Get capacitor pointer.
  elem = cir->getElement(newelem_id);
  // Set capacitor value
  elem->setParam("c", &C, TR_DOUBLE);
  // Init the capacitor
  elem->init();

  // Add EKV Pmos MP1
  newelem_id =
    cir->addElement("mospekv", getInstanceName() + "MP1", true);

  cir->connect(newelem_id,term_id3); // Connect to D input
  cir->connect(newelem_id,term_id_i1); // Connect to G
  cir->connect(newelem_id,term_id4); // Connect to S = VDD
  cir->connect(newelem_id,term_id4); // Connect To Bulk = VDD

  elem = cir->getElement(newelem_id);
  // Set parameter of pmos
  elem->setParam("lp",&lp, TR_DOUBLE);
  elem->setParam("wp",&wp, TR_DOUBLE);
  // elem->setParam("uo",&uo_p, TR_DOUBLE);
  elem->init();

  // Add EKV Pmos MP2
  newelem_id =
    cir->addElement("mospekv", getInstanceName() + "MP2", true);

  cir->connect(newelem_id,term_id3); // Connect to D input
  cir->connect(newelem_id,term_id_i2); // Connect to G
  cir->connect(newelem_id,term_id4); // Connect to S = VDD
  cir->connect(newelem_id,term_id4); // Connect To Bulk = VDD

  elem = cir->getElement(newelem_id);
  // Set parameter of pmos
  elem->setParam("lp",&lp, TR_DOUBLE);
  elem->setParam("wp",&wp, TR_DOUBLE);
  //   elem->setParam("uo",&uo_p, TR_DOUBLE);
  elem->init();

  //Add EKV nmos MN1
  newelem_id =
    cir->addElement("mosnekv", getInstanceName() + "MN1", true);

  cir->connect(newelem_id,term_id3); // Connect to D input
  cir->connect(newelem_id,term_id_i1); // Connect to G

  // Add an internal Terminal
  unsigned  term_id_i3 = cir->addTerminal(getInstanceName() + ":3", true);

  cir->connect(newelem_id,term_id_i3); // Connect to S = Internal terminal
  cir->connect(newelem_id,term_id_i3); // Connect To Bulk = Internal terminal Ignore Body Effect

  elem = cir->getElement(newelem_id);
  // Set parameter of nmos
  elem->setParam("ln",&ln, TR_DOUBLE);
  elem->setParam("wn",&wn, TR_DOUBLE);
  //     elem->setParam("uo",&uo_n, TR_DOUBLE);
  elem->init();

  //    Add capacitance to internal node of nand
  newelem_id =
    cir->addElement("capacitor", getInstanceName() + ":capacitors3:", true);
  // Connect to internal terminal and reference
  cir->connect(newelem_id, term_id_i3);
  cir->connect(newelem_id, tref_id);
  // Get capacitor pointer.
  elem = cir->getElement(newelem_id);
  // Set capacitor value
  elem->setParam("c", &C, TR_DOUBLE);
  // Init the capacitor
  elem->init();


  //Add EKV nmos MN2
  newelem_id =
    cir->addElement("mosnekv", getInstanceName() + "MN2", true);

  cir->connect(newelem_id,term_id_i3); // Connect to D input
  cir->connect(newelem_id,term_id_i2); // Connect to G
  cir->connect(newelem_id,tref_id); // Connect to S = Gnd
  cir->connect(newelem_id,tref_id); // Connect To Bulk = Gnd

  elem = cir->getElement(newelem_id);
  // Set parameter of nmos
  elem->setParam("ln",&ln, TR_DOUBLE);
  elem->setParam("wn",&wn, TR_DOUBLE);
  //       elem->setParam("uo",&uo_n, TR_DOUBLE);
  elem->init();

  //    Add Load Capacitance to output
  newelem_id =
    cir->addElement("capacitor", getInstanceName() + ":capacitors4:", true);
  // Connect to Output and reference
  cir->connect(newelem_id, term_id3);
  cir->connect(newelem_id, tref_id);
  // Get capacitor pointer.
  elem = cir->getElement(newelem_id);
  // Set capacitor value
  elem->setParam("c", &C, TR_DOUBLE);
  // Init the capacitor
  elem->init();
}

