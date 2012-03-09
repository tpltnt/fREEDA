#include "Element.h"

Element::Element(ItemInfo* einfo, ParmInfo* param_desc, const int& numparms,
const string& iname)
: GraphNode(iname), NetListItem(einfo, param_desc, numparms)
{
  my_circuit = NULL;

  // Set flags to zero by default
  myflags = 0;

  // Set graph processing flag to default value
  black = false;
  nopropagate = false;
  // Terminal initialization
  // 0 means unknown
  numterms = 0;

  // Default number of states is zero
  my_nstates = 0;

  // Set elemdata to null until the number of terminals are set.
  elemdata = NULL;
}

Element::~Element()
{
  // Disconnect element from terminals
  for (unsigned i=0; i < getTermCount(); i++)
    getTerminal(i)->detachElement(this);

  // This memory is allocated when numterms is set.
  if (numterms)
    delete elemdata;
}

void Element::setNumTerms(const unsigned& numterms)
{
  this->numterms = numterms;
  // Allocate output instance
  elemdata = new ElementData(numterms);
}

void Element::connect(Terminal* terminal)
{
  link(terminal);
  terminal->attachElement(this);
}

void Element::init() throw(string&)
{
  // By this time, the number of terminals must be already set.
  assert(numterms);
  // do nothing by default
  // We can use this routine to set flags once the element is built
  // If something goes wrong, an exception can be thrown.
}


void Element::check() throw(string&)
{
  if (!checkConnect())
    throw(getInstanceName() + ": Wrong number of connections.");
  try
	{
    checkParams();
  }
  catch(string& msg)
	{
    throw(getInstanceName() + ": " + msg);
  }
}

bool Element::checkConnect() const
{
  if (numterms)
  {
    return (numterms == getGraphNodeCount());
  }
  else
    // We can not say anything yet... (numterms can be set in init()).
	return true;
}


bool Element::satisfies(const ElemFlag& mask) const
{
  if ( (myflags & mask) == mask )
    return true;
  else
    return false;
}

// Get the number of secondary state variables
unsigned Element::getNumberOfSecStates() const
{
  // Return zero by default
  return 0;
}


// Frequency MNAM methods
unsigned Element::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Default is not to add extra rows and columns.
  return 0;
}

void Element::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  // Default is not to add extra rows and columns.
  first_eqn = n_rows = 0;
}

void Element::fillMNAM(FreqMNAM* mnam)
{
  // It is not desirable this function being called.
  // Use flags to select the proper elements
  assert(false);
}


void Element::fillMNAM(TimeMNAM* mnam)
{
  // It is not desirable this function being called.
  // Use flags to select the proper elements
  assert(false);
}

void Element::fillSourceV(TimeMNAM* mnam)
{
  // It is not desirable this function being called.
  // Use flags to select the proper elements
  assert(false);
}

void Element::setLastResult(DenseDoubleVector& res, const double& time)
{
  // It is not desirable this function being called.
  // Use flags to select the proper elements
  assert(false);
}

void Element::getTimeInfo(double & step, double & stop)
{
  // It is not desirable this function being called.
  // Use flags to select the proper elements
  assert(false);
}

void Element::getLocalRefIdx(UnsignedVector& local_ref_vec,
TerminalVector& term_list)
{
  // The element must have only one local reference
  assert(satisfies(ONE_REF));

  // Clean vectors first
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());
  // Assume the last terminal is the reference
  local_ref_vec.push_back(numterms-1);
  // Fill terminal vector with terminals in the netlist order
  for (unsigned i = 0; i < numterms; i++)
    term_list.push_back(getTerminal(i));
}

void Element::setGroup(unsigned ref_id)
{
  if (nopropagate || black || satisfies(MULTI_REF))
    return;

  black = true;
  // Otherwise, propagate external local reference terminal
  for (unsigned i=0; i < numterms; i++)
    getTerminal(i)->setGroup(ref_id);
  // After processing return black to normal state.
  black = false;
}

// ---------------------------------------
// Analysis Methods
// ---------------------------------------
// State-variable related methods
void Element::setNumberOfStates(const unsigned& n)
{
  my_nstates = n;
}

void Element::svHB(FreqDomainSV* fdsv)
{
  // control should never reach this point
  assert(false);
}

void Element::deriv_svHB(FreqDomainSV* fdsv)
{
  // control should never reach this point
  assert(false);
}

void Element::svTran(TimeDomainSV* tdsv)
{
  // control should never reach this point
  assert(false);
}

void Element::deriv_svTran(TimeDomainSV* tdsv)
{
  // control should never reach this point
  assert(false);
}

void Element::svWav(WaveletDomainSV* wdsv)
{
  // control should never reach this point
  assert(false);
}

void Element::deriv_svWav(WaveletDomainSV* wdsv)
{
  // control should never reach this point
  assert(false);
}

DenseDoubleVector Element::getDelayVec()
{
  cout << "\n*** Warning: Should not call getDelayVec() on base class 'Element'  ***\n";
  return 0;
}

