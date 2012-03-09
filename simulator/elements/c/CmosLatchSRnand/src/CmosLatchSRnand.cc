#include "../../../../network/ElementManager.h"
#include "../../../../network/CircuitManager.h"
#include "../../../../network/Element.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "CmosLatchSRnand.h"
#include <cstdio>

// Set the number of parameters which is 4 for this element
const unsigned CmosLatchSRnand::n_par = 4;

// Element information
ItemInfo CmosLatchSRnand::einfo =
{
  "cmoslatchsrnand",
  "Cmos Latch Based on Nand gate",
  "Shivam Priyadarshi",
  "category: Digital Logic",
  "2009_04_18"
};

ParmInfo CmosLatchSRnand::pinfo[] =
{
  {"ln","Channel Width of NMOS (m)", TR_DOUBLE, false},
  {"wn","Channel Length of NMOS (m)", TR_DOUBLE, false},
  {"lp","Channel Width of PMOS (m)", TR_DOUBLE, false},
  {"wp","Channel Length of PMOS (m)", TR_DOUBLE, false}
  // Later it need to be modified. The input parameters should come as nmodel and pmodel for nmos
  // and pmos. char* nmodel and char *pmodel type of construct should be used
};


CmosLatchSRnand::CmosLatchSRnand(const string& iname) : Element(&einfo, pinfo, n_par, iname)
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

void CmosLatchSRnand::init() throw(string&)
{
  C = 1e-12;
  Rs = 10;
  // Clear flags so this element is not called to fill the MNAM
  setFlags(ONE_REF);

  Circuit* cir = getCircuit();

  unsigned term_id1 = getTerminal(0)->getId();  // Vdd
  unsigned term_id2 = getTerminal(1)->getId();  // Sb
  unsigned term_id3 = getTerminal(2)->getId();  // Rb
  unsigned term_id4 = getTerminal(3)->getId();  // Q
  unsigned term_id5 = getTerminal(4)->getId();  // Qb
  unsigned tref_id = getTerminal(5)->getId();   // Gnd

  // Add Source Resistance 1
  unsigned newelem_id =
    cir->addElement("resistor", getInstanceName() + ":resistors1:", true);
  // Connect to input terminal
  cir->connect(newelem_id, term_id2);
  // Add an internal Terminal
  unsigned term_id_i2 = cir->addTerminal(getInstanceName() + ":1", true);
  // Connect to internal terminal
  cir->connect(newelem_id, term_id_i2);
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
  cir->connect(newelem_id, term_id3);
  // Add an internal Terminal
  unsigned  term_id_i3 = cir->addTerminal(getInstanceName() + ":2", true);
  // Connect to internal terminal
  cir->connect(newelem_id, term_id_i3);
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
  cir->connect(newelem_id, term_id_i2);
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
  cir->connect(newelem_id, term_id_i3);
  cir->connect(newelem_id, tref_id);
  // Get capacitor pointer.
  elem = cir->getElement(newelem_id);
  // Set capacitor value
  elem->setParam("c", &C, TR_DOUBLE);
  // Init the capacitor
  elem->init();


  // Add Nand Gate 1 order of Nand Terminals are : Vdd In1 In2 Out Gnd
  newelem_id =
    cir->addElement("cmos2nand", getInstanceName() + ":Nand1:", true);
  // Connect to terminals
  cir->connect(newelem_id, term_id1);   //Vdd Connect
  cir->connect(newelem_id, term_id_i2);   // Sb Connect
  cir->connect(newelem_id, term_id5);   // Qb Connect
  cir->connect(newelem_id, term_id4);   // Q Connect
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

  // Add Another Nand Gate

  newelem_id =
    cir->addElement("cmos2nand", getInstanceName() + ":Nand2:", true);
  // Connect to terminals
  cir->connect(newelem_id, term_id1);   //Vdd Connect
  cir->connect(newelem_id, term_id_i3);   // Rb Connect
  cir->connect(newelem_id, term_id4);   // Q Connect
  cir->connect(newelem_id, term_id5);   // Qb Connect
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

  //    Add Load Capacitance to output : Q
  newelem_id =
    cir->addElement("capacitor", getInstanceName() + ":capacitors4:", true);
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
    cir->addElement("capacitor", getInstanceName() + ":capacitors5:", true);
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

