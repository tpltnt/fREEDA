#include "../../../../network/ElementManager.h"
#include "../../../../network/CircuitManager.h"
#include "../../../../network/Element.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "CmosLatchSRnandX.h"
#include <cstdio>

// Set the number of parameters which is 4 for this element
const unsigned CmosLatchSRnandX::n_par = 4;

// Element information
ItemInfo CmosLatchSRnandX::einfo =
{
  "cmoslatchsrnandx",
  "Latch Based on Nand gate",
  "Shivam Priyadarshi",
  "category: Digital Logic",
  "2009_04_18"
};

ParmInfo CmosLatchSRnandX::pinfo[] =
{
  {"ln","Channel Width of NMOS (m)", TR_DOUBLE, false},
  {"wn","Channel Length of NMOS (m)", TR_DOUBLE, false},
  {"lp","Channel Width of PMOS (m)", TR_DOUBLE, false},
  {"wp","Channel Length of PMOS (m)", TR_DOUBLE, false},
};


CmosLatchSRnandX::CmosLatchSRnandX(const string& iname) : Element(&einfo, pinfo, n_par, iname)
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

void CmosLatchSRnandX::init() throw(string&)
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

  // Add Nand Gate 1 order of Nand Terminals are : Vdd In1 In2 Out Gnd
  unsigned newelem_id =
    cir->addElement("cmos2nandx", getInstanceName() + ":Nand1:", true);
  // Connect to terminals
  cir->connect(newelem_id, term_id1);   //Vdd Connect
  cir->connect(newelem_id, term_id2);   // Sb Connect
  cir->connect(newelem_id, term_id5);   // Qb Connect
  cir->connect(newelem_id, term_id4);   // Q Connect
  cir->connect(newelem_id, tref_id);    // Ground Conenct
  // Get Nand pointer.
  Element* elem = cir->getElement(newelem_id);
  // Set Nand Parameter Values
  elem->setParam("ln",&ln, TR_DOUBLE);
  elem->setParam("wn",&wn, TR_DOUBLE);
  elem->setParam("lp",&lp, TR_DOUBLE);
  elem->setParam("wp",&wp, TR_DOUBLE);
  // Init the Nand
  elem->init();

  // Add Another Nand Gate

  newelem_id =
    cir->addElement("cmos2nandx", getInstanceName() + ":Nand2:", true);
  // Connect to terminals
  cir->connect(newelem_id, term_id1);   //Vdd Connect
  cir->connect(newelem_id, term_id3);   // Rb Connect
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
}

