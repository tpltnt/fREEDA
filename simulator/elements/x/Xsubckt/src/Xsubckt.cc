#include <cassert>
#include "../../../../network/CircuitManager.h"
#include "Xsubckt.h"

// Static members
const unsigned Xsubckt::n_par = 0;

// Element information
ItemInfo Xsubckt::einfo =
{
  "xsubckt",
  "Subcircuit instance",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:subcircuit",
  "2000_07_20"
};

Xsubckt::Xsubckt(const string& iname) : Element(&einfo, NULL, n_par, iname)
{
	circuit_def = NULL;
  // Set flags for the subcircuit element
  // We do not know anything by now.
  setFlags(0);
}

Xsubckt::~Xsubckt()
{
  if (circuit_def)
	{
    circuit_def->delUser();
    // By now, call the circuit manager to erase the definition
    // (if unused).
    if (!circuit_def->getNumberOfUsers())
      the_CM->deleteCircuit(circuit_def->getInstanceName());
  }
}

void Xsubckt::init() throw(string&)
{
  // At this point everything should be Ok so,
  // set the number of terminals.
  setNumTerms(getTermCount());
}


bool Xsubckt::checkConnect() const
{
  assert(circuit_def);
  // We check if the the circuit definition has the same number of
  // terminals
  return (getTermCount() == circuit_def->getNumberOfConnections());
}

void Xsubckt::attachDefinition(Circuit *circuit)
{
  assert(!circuit_def);
  circuit_def = circuit;
  circuit_def->addUser();
}


void Xsubckt::expandToCircuit(Circuit* target_c)
{
  assert(target_c);
  // Firtst of all expand circuit_def
  circuit_def->expand();

  // Set mask to zero (want all elements)
  Terminal* term;
  circuit_def->setFirstElement(0);
  Element* elem = circuit_def->nextElement();
  unsigned newelem_id;
  Element* newelem;
  while(elem)
	{
    // Do not forget to copy the element!!!
    newelem_id = target_c->addElement(elem->getName(),
		getInstanceName() + ":" +
		elem->getInstanceName());
    // Get element pointer.
    newelem = target_c->getElement(newelem_id);
    // Copy element parameters.
    newelem->getParamsFrom(elem);
    // Look at the terminals of the element
    unsigned tc = elem->getNumTerms();
    for (unsigned i=0; i < tc; i++)
		{
      term = elem->getTerminal(i);
      assert(term);
      if (term->isExternal())
			{
				// search for external terminal index.
				unsigned j = circuit_def->getConnectIndex(term->getId());
				// Now connect to existing terminal.
				target_c->connect(newelem_id, getTerminal(j)->getId());
				// Set reference
				if (term->isRef())
					getTerminal(j)->setRef();
      }
      else
			{
				// The name of the new terminal has the subcircuit
				// name preceding it.
				unsigned term_id = target_c->addTerminal(getInstanceName() + ":" +
				term->getInstanceName());
				target_c->connect(newelem_id, term_id);
				// Set reference
				if (term->isRef())
					target_c->getTerminal(term_id)->setRef();
      }
    }
    // Init the element just added
    //    newelem->init();

    // get next element pointer
    elem = circuit_def->nextElement();
  }
}

