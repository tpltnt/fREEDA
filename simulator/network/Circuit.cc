#include "Circuit.h"
#include "../elements/x/Xsubckt/src/Xsubckt.h"

ElementManager* Circuit::the_EM = ElementManager::getElementManager();

Circuit::Circuit(const string& iname, bool subckt = false):Instanciable(iname)
{
  last_element = NULL;

  subcircuit = subckt;

	// reserve space in connections vector
  if (subcircuit)
	  connect_term.reserve(64);

  current_element_it = elem_set.begin();
  current_terminal_it = term_set.begin();

  // reserve space for subckt instances
  subckt_list.reserve(64);

  flattened = false;
  used_by = 0;

  // Set default mask to be zero (i.e.: take all elements)
  elem_mask = 0;
}

Circuit::~Circuit()
{
  // Free all the elements in this circuit
  ElementMap::iterator i = elem_set.begin();
  while (i != elem_set.end()) 
	{
    delete ((*i).second);
    ++i;
  }
  
	// Free all the terminals in this circuit
  TerminalMap::iterator j = term_set.begin();
  while (j != term_set.end()) 
	{
    delete ((*j).second);
    ++j;
  }
}

unsigned Circuit::addElement(const string& elem_type, const string& iname,
bool nopropagate)
throw(string&)
{
  // Check if the element is already included in the circuit.
  Element* n_elem = getElement(iname);
  if (!n_elem) 
	{
    n_elem = the_EM->createElement(elem_type, iname);
    elem_set[n_elem->getId()] = n_elem;
    n_elem->setCircuit(this);
    if (nopropagate)
      n_elem->noPropagate();
    last_element = n_elem;
    return n_elem->getId();
  }
  else
    throw ("Element " + iname + " duplicated in circuit " +
	this->getInstanceName());
}


void Circuit::removeElement(const unsigned& id)
{
  if (last_element && last_element->getId() == id)
    last_element = NULL;

  ElementMap::iterator emi = elem_set.find(id);
  // Be sure that the element is in the circuit
  assert(emi != elem_set.end());

  // Get pointer to element
  Element* elem = (*emi).second;
  // Erase the element from set
  elem_set.erase(emi);
  // Delete element
  delete elem;
}

void Circuit::connect(const unsigned& elem_id, const unsigned& term_id)
{
  Terminal* terminal = getTerminal(term_id);
  assert(terminal);

  if (!last_element || elem_id != last_element->getId()) 
	{
    last_element = getElement(elem_id);
    assert(last_element);
  }

  // Tell the element to connect itself to the terminal.
  last_element->connect(terminal);
}

Element* Circuit::getElement(const string& elem_name)
{
  // For speed check first with last_element
  if (last_element && elem_name == last_element->getInstanceName())
    return last_element;

  for (ElementMap::iterator j = elem_set.begin(); j != elem_set.end() ; j++)
    if ((*j).second->getInstanceName() == elem_name)
      return (*j).second;

  return NULL;
}


Element* Circuit::getElement(const unsigned& id)
{
  // For speed check first with last_element
  if (last_element && id == last_element->getId())
    return last_element;

  ElementMap::iterator j = elem_set.find(id);
  if (j != elem_set.end())
    return((*j).second);
  else
    return NULL;
}


void Circuit::setFirstElement(ElemFlag mask = 0)
{
  current_element_it = elem_set.begin();
  elem_mask = mask;
}

Element* Circuit::nextElement()
{
  // Travel the list until we find the end or an element with the
  // desired flags.
  while(current_element_it != elem_set.end() &&
		! ((*current_element_it).second)->satisfies(elem_mask))
	current_element_it++;

  // If the list is not at the end, then return the element found.
  if (current_element_it != elem_set.end())
    return (*current_element_it++).second;
  else
    return NULL;
}


unsigned Circuit::addTerminal(const string& term_name, bool nocheck)
{
  // Check if the terminal is already created in the circuit.
  Terminal* n_term = getTerminal(term_name);
  if (!n_term) 
	{
    n_term = new Terminal(term_name);
    term_set[n_term->getId()] = n_term;
  }
  if (nocheck)
    n_term->noCheck();
  return n_term->getId();
}

void Circuit::removeTerminal(const unsigned& id)
{
  TerminalMap::iterator emi = term_set.find(id);
  // Be sure that the terminal is in the circuit
  assert(emi != term_set.end());

  // Get pointer to terminal
  Terminal* term = (*emi).second;

  // Assert the terminal is not connected to elements
  assert(term->getNumberOfConnect() == 0);

  // Erase the terminal from set
  term_set.erase(emi);
  // Delete terminal
  delete term;
}

void Circuit::setRefTerm(unsigned& id)
{
  Terminal* term = getTerminal(id);
  assert(term);
  if (!term->isRef()) 
	{
    term->setRef();
    ref_terms.push_back(term);
  }
}

Terminal* Circuit::getTerminal(const string& name)
{
  for (TerminalMap::iterator j = term_set.begin(); j != term_set.end() ; j++)
    if ((*j).second->getInstanceName() == name)
      return (*j).second;

  return NULL;
}

Terminal* Circuit::getTerminal(const unsigned& id)
{
  TerminalMap::iterator j = term_set.find(id);
  if (j != term_set.end())
    return((*j).second);
  else
    return NULL;
}

Terminal* Circuit::nextTerminal()
{
  if(current_terminal_it != term_set.end())
    return (*current_terminal_it++).second;
  else
    return NULL;
}


void Circuit::check() throw(string&)
{
  // check all elements
  ElementMap::iterator j = elem_set.begin();
  while (j != elem_set.end()) 
	{
    ((*j).second)->check();
    ++j;
  }
  // Check connectivity of terminals
	checkTerms();
}

void Circuit::init() throw(string&)
{
  // init all elements currently in the circuit.
  // Some of the elements may expand into more elements.
  ElementVector elem_vec;
  unsigned nelem(getNumberOfElements());
  elem_vec.reserve(nelem);
  ElementMap::iterator j = elem_set.begin();
  unsigned i(0);

  while (j != elem_set.end()) 
	{
    elem_vec[i] = ((*j).second);
    ++j;
    i++;
  }

  for (i=0; i < nelem; i++)
    elem_vec[i]->init();
}

void Circuit::checkReferences() throw(string&)
{
  // This method only work if the circuit has been flattened
  assert(flattened);

  unsigned nref = ref_terms.size();

  if (!nref)
    throw(string("Circuit ") + getInstanceName() +
	" has no reference terminal.");

  for (unsigned i=0; i<nref; i++)
    // Propagate each reference
	ref_terms[i]->propagate();

  // Check groups in multi reference elements.
  // Loop through all selected elements in circuit.
  setFirstElement(MULTI_REF);
  Element* elem = nextElement();
  while(elem) 
	{
    UnsignedVector local_ref_vec;
    TerminalVector term_list;
    // Get local groups information
    elem->getLocalRefIdx(local_ref_vec, term_list);
    unsigned ngroups = local_ref_vec.size();
    // jbase is the first terminal index in each group
    unsigned jbase = 0;
    for (unsigned i=0; i<ngroups; i++) 
		{
      unsigned current_group = 0;
      for (unsigned j = jbase; j <= local_ref_vec[i]; j++) 
			{
				current_group = term_list[j]->getGroup();
				if (current_group)
					break;
      }
      // Now propagate the group reference, if any found
      if (current_group)
				for (unsigned j = jbase; j <= local_ref_vec[i]; j++)
					term_list[j]->setGroup(current_group);

      // update jbase
      jbase = local_ref_vec[i] + 1;
    }
    // get next element pointer
    elem = nextElement();
  }

  // Check that there are no floating groups of terminals.
  TerminalMap::iterator ti = term_set.begin();
  // Check each terminal
  while (ti != term_set.end()) 
	{
    if (!((*ti).second)->getGroup())
      throw("Terminal \"" + ((*ti).second)->getInstanceName() +
		"\" belongs to a group without a reference terminal.");
    ++ti;
  }
}

void Circuit::checkTerms() throw(string&)
{
  TerminalMap::iterator i = term_set.begin();
  // Check each terminal
  while (i != term_set.end()) 
	{
    if (!((*i).second)->checkConnect())
      throw("Terminal \"" + ((*i).second)->getInstanceName() +
		"\" connected to less than two elements.");
    ++i;
  }
}

void Circuit::expand()
{
  if (!flattened) 
	{
    // Tell all subcircuits to expand themselves
    unsigned sc = subckt_list.size();
    for (unsigned i=0; i < sc; i++) 
		{
      ((Xsubckt*)(subckt_list[i]))->expandToCircuit(this);
      // Remove subckt instance
      removeElement(subckt_list[i]->getId());
    }
    // Leave subckt_list as an empty vector
    subckt_list = ElementVector();
    flattened = true;
  }
}


void Circuit::delUser()
{
  used_by--;
  assert(used_by >= 0);
}

void Circuit::setConnection(const unsigned& id)
{
  assert(subcircuit);
  Terminal* terminal = getTerminal(id);
  // Set the terminal as extern.
  terminal->setExtern();
  // add the new terminal to the connections vector
  connect_term.push_back(terminal);
}


Terminal* Circuit::getConnect(const unsigned& index)
{
  assert(index < connect_term.size());
  return connect_term[index];
}

unsigned Circuit::getConnectIndex(const unsigned& id) const
{
  unsigned tc = getNumberOfConnections();
  for (unsigned i=0; i < tc; i++)
    if (connect_term[i]->getId() == id)
      return i;

  // if control reaches here, that means that the terminal is not external
  assert(false);
}

