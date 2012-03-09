#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Circulator.h"

// Static members
const unsigned Circulator::n_par = 2;

// Element information
ItemInfo Circulator::einfo =
{
  "circulator",
  "Circulator",
  "Don Widdows - Isac Lima - Daryl Lindsey",
  DEFAULT_ADDRESS"category:behavioral,functional",
  "2003_05_15"
};

// Parameter information
ParmInfo Circulator::pinfo[] =
{
  {"r", "Resistance looking into each circulator port (Ohms)", TR_DOUBLE, false},
  {"g", "Conductance looking into each circulator port (S)", TR_DOUBLE, false}
};

Circulator::Circulator(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Value either of r or g may be input.
  paramvalue[0] = &(r=0);
  paramvalue[1] = &(g=0);

  // Set the number of terminals
  setNumTerms(6);

  // Set flags
  setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN); //must put multi-ref
                                                 // for local ref term's!
  //Linear means it's going to check for fillMNAM functions.
  //fREEDA is written to have certain conventions about functions.
  //It has to know which functions to call from within the main body
  //  of the code when certain flags are set.
}

void Circulator::init() throw(string&)
{
  if(r)
   g=1/r;
  else if (!g)
     g=.02;

  //If no resistance or conductance specified, make conductance = .02 S.
  //This is set up so that the value of the conductance will default to .02S
  // (resistance defaults to 50 ohms) if neither parameter is specified.
}

/*These functions may be useful in a later implementation of a non-ideal
circulator
//unsigned Circulator::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
//{
//   for (int i=0; i < 2; i++)
//    my_row[i] = eqn_number + i;
//
//   // Add 2 extra RCs
//    return 2;
//
//}

//void Circulator::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
//{
//   assert(my_row[0]);
//   first_eqn = my_row[0];
//   n_rows = 2;
//   assert(my_row[0]);
//
//}
*/


void Circulator::getLocalRefIdx(UnsignedVector& local_ref_vec,
			    TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1));
  term_list.push_back(getTerminal(2)); //local reference terminal
  term_list.push_back(getTerminal(3)); //local reference terminal
  term_list.push_back(getTerminal(4));
  term_list.push_back(getTerminal(5)); //local reference terminal

  local_ref_vec.push_back(2); // Local reference index
  local_ref_vec.push_back(3); // Local reference index
  local_ref_vec.push_back(5); // Local reference index
}


void Circulator::fillMNAM(FreqMNAM* mnam)
{
  //initialize complex admittance values
  //this may be useful in a later non-ideal implementation:
  //const double& freq(mnam->getFreq());
  //double_complex z1(r1,twopi*freq*l1);
  //double_complex z2(r2,twopi*freq*l2);
  //double_complex y1=1/z1;
  //double_complex y2=1/z2;

  double y1=g;
  double y2=g;
  //During the process of formulating the admittance matrix,
  //we used the admittance stamps of 3 ideal gyrators.  Adding
  //the stamps together, we produced a fairly dense admittance matrix.
  //Typically, a circulator will have a single admittance value.
  //The reason we have y1 and y2 in the matrix instead of just y is that
  //the circulator is gyrator-implemented, and each gyrator has a
  //certain conductance looking into each port.  We are leaving the matrix as it
  //is--with both values of admittance--so that this fact will be better
  //understood.

  // Ask my terminals the row numbers
   const unsigned term1=getTerminal(0)->getRC();
   const unsigned term2=getTerminal(1)->getRC();
   const unsigned term3=getTerminal(2)->getRC();
   const unsigned term4=getTerminal(3)->getRC();
   const unsigned term5=getTerminal(4)->getRC();
   const unsigned term6=getTerminal(5)->getRC();

  //3-Port (6-terminal) Circulator Model
  //It is essentially comprised of 3 gyrators
  mnam->setElement(term1,term1,0);
  mnam->setElement(term1,term2,y1);
  mnam->setElement(term1,term3,-y2);
  mnam->setElement(term1,term4,0);
  mnam->setElement(term1,term5,-y1);
  mnam->setElement(term1,term6,y2);
  mnam->setElement(term2,term1,-y2);
  mnam->setElement(term2,term2,0);
  mnam->setElement(term2,term3,y1);
  mnam->setElement(term2,term4,y2);
  mnam->setElement(term2,term5,0);
  mnam->setElement(term2,term6,-y1);
  mnam->setElement(term3,term1,y1);
  mnam->setElement(term3,term2,-y2);
  mnam->setElement(term3,term3,0);
  mnam->setElement(term3,term4,-y1);
  mnam->setElement(term3,term5,y2);
  mnam->setElement(term3,term6,0);
  mnam->setElement(term4,term1,0);
  mnam->setElement(term4,term2,-y1);
  mnam->setElement(term4,term3,y2);
  mnam->setElement(term4,term4,0);
  mnam->setElement(term4,term5,y1);
  mnam->setElement(term4,term6,-y2);
  mnam->setElement(term5,term1,y2);
  mnam->setElement(term5,term2,0);
  mnam->setElement(term5,term3,-y1);
  mnam->setElement(term5,term4,-y2);
  mnam->setElement(term5,term5,0);
  mnam->setElement(term5,term6,y1);
  mnam->setElement(term6,term1,-y1);
  mnam->setElement(term6,term2,y2);
  mnam->setElement(term6,term3,0);
  mnam->setElement(term6,term4,y1);
  mnam->setElement(term6,term5,-y2);
  mnam->setElement(term6,term6,0);
}

void Circulator::fillMNAM(TimeMNAM* mnam)
{
  double g1=g;
  double g2=g;
  //see the comments on y1 and y2 in the FreqMNAM matrix above.
  //The same idea applies here.
  //have to get the terminal numbers
  const unsigned term1=getTerminal(0)->getRC();
  const unsigned term2=getTerminal(1)->getRC();
  const unsigned term3=getTerminal(2)->getRC();
  const unsigned term4=getTerminal(3)->getRC();
  const unsigned term5=getTerminal(4)->getRC();
  const unsigned term6=getTerminal(5)->getRC();

  //3-Port (6-terminal) Circulator Model
  //It is essentially comprised of 3 gyrators
  mnam->setMElement(term1,term1,0);
  mnam->setMElement(term1,term2,g1);
  mnam->setMElement(term1,term3,-g2);
  mnam->setMElement(term1,term4,0);
  mnam->setMElement(term1,term5,-g1);
  mnam->setMElement(term1,term6,g2);
  mnam->setMElement(term2,term1,-g2);
  mnam->setMElement(term2,term2,0);
  mnam->setMElement(term2,term3,g1);
  mnam->setMElement(term2,term4,g2);
  mnam->setMElement(term2,term5,0);
  mnam->setMElement(term2,term6,-g1);
  mnam->setMElement(term3,term1,g1);
  mnam->setMElement(term3,term2,-g2);
  mnam->setMElement(term3,term3,0);
  mnam->setMElement(term3,term4,-g1);
  mnam->setMElement(term3,term5,g2);
  mnam->setMElement(term3,term6,0);
  mnam->setMElement(term4,term1,0);
  mnam->setMElement(term4,term2,-g1);
  mnam->setMElement(term4,term3,g2);
  mnam->setMElement(term4,term4,0);
  mnam->setMElement(term4,term5,g1);
  mnam->setMElement(term4,term6,-g2);
  mnam->setMElement(term5,term1,g2);
  mnam->setMElement(term5,term2,0);
  mnam->setMElement(term5,term3,-g1);
  mnam->setMElement(term5,term4,-g2);
  mnam->setMElement(term5,term5,0);
  mnam->setMElement(term5,term6,g1);
  mnam->setMElement(term6,term1,-g1);
  mnam->setMElement(term6,term2,g2);
  mnam->setMElement(term6,term3,0);
  mnam->setMElement(term6,term4,g1);
  mnam->setMElement(term6,term5,-g2);
  mnam->setMElement(term6,term6,0);
}

