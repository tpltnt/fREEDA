// George Feller
#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Lspiral.h"
#include <cstdio>

// Static members
const unsigned Lspiral::n_par = 15;

// element type indicator - used by addLRC routine
const int LTYPE=0; // inductor
const int RTYPE=1; // resistor
const int CTYPE=2; // capacitor

// Element information
ItemInfo Lspiral::einfo =
{
  "lspiral",
  "planar spiral inductor for RF IC",
  "George Feller",
  DEFAULT_ADDRESS"category:inductor",
  "2008_04_07"
};

// Parameter information
ParmInfo Lspiral::pinfo[] =
{
  {"rdc", "", TR_DOUBLE, true},
  {"ldc", "", TR_DOUBLE, true},
  {"rs1", "", TR_DOUBLE, true},
  {"ls1", "", TR_DOUBLE, true},
  {"ms1", "", TR_DOUBLE, true},
  {"cw", "", TR_DOUBLE, true},
  {"cox1", "", TR_DOUBLE, true},
  {"csub11", "", TR_DOUBLE, true},
  {"rsub11", "", TR_DOUBLE, true},
  {"cox2", "", TR_DOUBLE, true},
  {"csub21", "", TR_DOUBLE, true},
  {"rsub21", "", TR_DOUBLE, true},
  {"cox3", "", TR_DOUBLE, true},
  {"csub31", "", TR_DOUBLE, true},
  {"rsub31", "", TR_DOUBLE, true},
 };

// yet another data blob
struct grfBLOB
{
  const char *param_name;
  const char *element_name;
  const char *type_suffix;
};

// yet another array of blobs
grfBLOB EBLOBS[]=
{
  {"l","inductor",":inductor:"}, 
  {"r","resistor",":resistor:"}, 
  {"c","capacitor",":capacitor:"} 
};

Lspiral::Lspiral(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  paramvalue[0] = &rdc;
  paramvalue[1] = &ldc;
  paramvalue[2] = &rs1;
  paramvalue[3] = &ls1;
  paramvalue[4] = &ms1;
  paramvalue[5] = &cw;
  paramvalue[6] = &cox1;
  paramvalue[7] = &csub11;
  paramvalue[8] = &rsub11;
  paramvalue[9] = &cox2;
  paramvalue[10] = &csub21;
  paramvalue[11] = &rsub21;
  paramvalue[12] = &cox3;
  paramvalue[13] = &csub31;
  paramvalue[14] = &rsub31;

  // Set the number of terminals
  setNumTerms(3);

  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);
}

void Lspiral::init() throw(string&)
{
  unsigned et1=getTerminal(0)->getId(); // external term
  unsigned et2=getTerminal(1)->getId();
  unsigned rt1=getTerminal(2)->getId(); // reference term
 
  Circuit *theCircuit = getCircuit();
  assert(theCircuit);

  // The following is an example of poor coding technique
  unsigned it1 = theCircuit->addTerminal(getInstanceName()+":1",true);
  unsigned it2 = theCircuit->addTerminal(getInstanceName()+":2",true);
  unsigned it3 = theCircuit->addTerminal(getInstanceName()+":3",true);
  unsigned it4 = theCircuit->addTerminal(getInstanceName()+":4",true);
  unsigned it5 = theCircuit->addTerminal(getInstanceName()+":5",true);
  unsigned it6 = theCircuit->addTerminal(getInstanceName()+":6",true);
  unsigned it7 = theCircuit->addTerminal(getInstanceName()+":7",true);
  unsigned it8 = theCircuit->addTerminal(getInstanceName()+":8",true);

  sid = 0; // no constituent elements added yet

  // interwinding capacitance
  addLRC(CTYPE,&cw,et1,et2); // C1 

  // frequency dependent inductor left pi-network	
  Element* l1 = addLRC(LTYPE,&ldc,et1,it1); // L1
  Element* l2 = addLRC(LTYPE,&ls1,it1,it2); // L2
  addLRC(RTYPE,&rs1,it2,it3); // R2
  addLRC(RTYPE,&rdc,it1,it3); // R1

  // inductive coupling
  coupleInductors(l1,l2);

  // frequency dependent inductor right pi-network
  l1 = addLRC(LTYPE,&ldc,it3,it4); // L3
  l2 = addLRC(LTYPE,&ls1,it4,it5); // L4
  addLRC(RTYPE,&rs1,it5,et2); // R4
  addLRC(RTYPE,&rdc,it4,et2); // R3

  // inductive coupling
  coupleInductors(l1,l2);

  // substrate parasitics
  addLRC(CTYPE,&cox1,et1,it6);
  addLRC(CTYPE,&csub11,it6,rt1);
  addLRC(RTYPE,&rsub11,it6,rt1);	

  addLRC(CTYPE,&cox2,it3,it7);
  addLRC(CTYPE,&csub21,it7,rt1);
  addLRC(RTYPE,&rsub21,it7,rt1);	

  addLRC(CTYPE,&cox3,et2,it8);
  addLRC(CTYPE,&csub31,it8,rt1);
  addLRC(RTYPE,&rsub31,it8,rt1);	
}

void Lspiral::fillMNAM(FreqMNAM* mnam)
{
  // Internal elements handle this
}

void Lspiral::fillMNAM(TimeMNAM* mnam)
{
  // Internal elements handle this
}


// Helper routine for adding inductors, resistors and capacitors
Element* Lspiral::addLRC(int type, double* value,unsigned int t1,unsigned int t2) throw(string&)
{  
   char csid[10]={0};
   sprintf(csid,"%d",sid);

   Circuit *theCircuit = getCircuit();
   assert(theCircuit);
   
   // The items under consideration are all two terminal elements
   unsigned int eid = theCircuit->addElement(EBLOBS[type].element_name,
	getInstanceName() + EBLOBS[type].type_suffix + csid,
	true);
#ifdef DEBUG
   fprintf("Adding element: %s %s",EBLOBS[type].element_name,std::string(getInstanceName() + EBLOBS[type].type_suffix+csid).c_str());
#endif

   theCircuit->connect(eid,t1);
   theCircuit->connect(eid,t2);
   
   Element *theElement = theCircuit->getElement(eid);
   assert(theElement);
	
   // All parameters are of the same type and ought to be set	
   theElement->setParam(EBLOBS[type].param_name,value,TR_DOUBLE);
   theElement->init();
   
   sid++; // next element gets a new id.
   return theElement;
}

// Helper routine for coupling inductors
void Lspiral::coupleInductors(Element* l1, Element* l2) throw(string&)
{ 
   char csid[10]={0};
   sprintf(csid,"%d",sid);
 
   assert(l1);
   assert(l2);

   Circuit *theCircuit = getCircuit();
   assert(theCircuit);
   
   // The items under consideration are all two terminal elements
   unsigned int eid = theCircuit->addElement("k",
	getInstanceName() + ":k:" + csid,true);
	   
   Element *theElement = theCircuit->getElement(eid);
   assert(theElement);

   std::string l1name(l1->getInstanceName());
   std::string l2name(l2->getInstanceName());
	
   // All parameters are of the same type and ought to be set	
   theElement->setParam("coupling",&ms1,TR_DOUBLE);
   theElement->setParam("l1",&(l1name),TR_STRING);
   theElement->setParam("l2",&(l2name),TR_STRING);
   theElement->init();
   
   sid++; // next element gets a new id.
}

