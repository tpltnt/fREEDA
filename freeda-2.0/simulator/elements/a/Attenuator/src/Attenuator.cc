#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "math.h"
#include "Attenuator.h"

extern "C"
{
#include "../../../../inout/parserGlobals.h"
}

extern int thisTermCount;

// Static members
const unsigned Attenuator::n_par = 2;

// Element information
ItemInfo Attenuator::einfo = 
{
  "attenuatorgen",
  "Resistive Attenuator, 3 or 4 terminal",
  "Shravanthi Ippagunta, Michael Steer",
  "category:passive",
  "2008_04_29"
};

// Parameter information
ParmInfo Attenuator::pinfo[] = 
{
  {"zo", "System characteristic impedances (Ohms)", TR_DOUBLE, true},
  {"alpha", "Attenuation (dB)", TR_DOUBLE, true}
};

Attenuator::Attenuator(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Value of Z0 and alpha are required
  paramvalue[0] = &(zo);
  paramvalue[1] = &(alpha);

  // Check to see if this is really a three terminal element by checking
  // If the second and fourth terminals are the same
  if(thisTermCount == 4)
  {
    if(thisTermList[1].val == thisTermList[3].val)
      numberOfTerminals = 3;
    else
      numberOfTerminals = 4;
  }
  else
    ErrMsg("Incorrect number of terminals.");

  if(thisTermList[0].val == thisTermList[2].val)
    ErrMsg("Incorrect terminals.");

  // Four terminals must specified, even if the second and
  //  fourth are the same. 
  setNumTerms(4);

  // Set flags
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN); 
};

//Resistor and Attenuation Calculations 
void Attenuator::init() throw(string&)
{ 
  double y = (alpha/20);
  double w = (alpha/10);
  double x = 10;

  if(numberOfTerminals==3)
  {
    r1 = zo*((pow(x,y)+1)/(pow(x,y)-1));
    r2 = zo*((pow(x,y)+1)/(pow(x,y)-1));
    r3 = zo*((pow(x,w)-1)/(2*pow(x,y))); 

  }
  else if(numberOfTerminals==4)
  {
    r1 = zo*((pow(x,y)+1)/(pow(x,y)-1));
    r2 = zo*((pow(x,y)+1)/(pow(x,y)-1));
    r3 = zo*((pow(x,w)-1)/(4*pow(x,y))); 
    r4 = zo*((pow(x,w)-1)/(4*pow(x,y))); 

  } 
};

void Attenuator::fillMNAM(FreqMNAM* mnam)
{
  // Ask my terminals the row numbers
  if(numberOfTerminals==3)
  {
    mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(2)->getRC(), 1/r3);
    mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), 1/r1);
    mnam->setAdmittance(getTerminal(2)->getRC(), getTerminal(1)->getRC(), 1/r2); 
  }
  else if(numberOfTerminals==4)
  {
    mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(2)->getRC(), 1/r3);
    mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), 1/r1);
    mnam->setAdmittance(getTerminal(2)->getRC(), getTerminal(3)->getRC(), 1/r2);
    mnam->setAdmittance(getTerminal(1)->getRC(), getTerminal(3)->getRC(), 1/r4);  
  }

}

void Attenuator::fillMNAM(TimeMNAM* mnam)
{
  // Ask my terminals the row numbers
  if(numberOfTerminals==3)
  {
    mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(2)->getRC(), 1/r3);
    mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), 1/r1);
    mnam->setMAdmittance(getTerminal(2)->getRC(), getTerminal(1)->getRC(), 1/r2); 
  }
  else if(numberOfTerminals==4)
  {
    mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(2)->getRC(), 1/r3);
    mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), 1/r1);
    mnam->setMAdmittance(getTerminal(2)->getRC(), getTerminal(3)->getRC(), 1/r2);
    mnam->setMAdmittance(getTerminal(1)->getRC(), getTerminal(3)->getRC(), 1/r4); 
  }
}

