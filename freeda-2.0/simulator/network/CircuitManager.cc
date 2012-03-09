#include "CircuitManager.h"

// Create the circuit manager
CircuitManager *the_CM = CircuitManager::getCircuitManager();

CircuitManager* CircuitManager::cm = NULL;

CircuitManager::CircuitManager()
{
  // Reserve space in the vectors
  circuit_list.reserve(64);
  circuit_stack.reserve(64);
  current_circuit = NULL;
}

CircuitManager::~CircuitManager()
{
  // Free all circuits in the vector
  for(unsigned i=0; i<circuit_list.size(); i++)
    delete circuit_list[i];
}


// It can exist only one circuit manager
CircuitManager* CircuitManager::getCircuitManager()
{
  if (!cm)
    cm = new CircuitManager;

  return cm;
}


Circuit* CircuitManager::getCircuit(const string& iname)
{
  // Find if the circuit is already created.
  for(unsigned i=0; i<circuit_list.size(); i++)
	if (circuit_list[i]->getInstanceName() == iname)
	{
		current_circuit = circuit_list[i];
		return current_circuit;
	}
  // Create a new circuit
  current_circuit = createCircuit(iname);
  return current_circuit;
}

Circuit* CircuitManager::getCircuit(unsigned& i)
{
  assert(i<circuit_list.size());
  return circuit_list[i];
}

Circuit* CircuitManager::getSubCircuit(const string& iname)
{
  // Find if the circuit is already created.
  for(unsigned i=0; i<circuit_list.size(); i++)
    if (circuit_list[i]->getInstanceName() == iname)
      return(circuit_list[i]);

  // Create a new circuit
  return createCircuit(iname, true);
}

Circuit* CircuitManager::createCircuit(const string& iname, bool subckt)
{
  // create the new circuit
  Circuit* current = new Circuit(iname, subckt);

  // insert the circuit in the vector
  circuit_list.push_back(current);

  return current;
}

void CircuitManager::deleteCircuit(const string& iname)
{
  // Find if the circuit is already created.
  unsigned i;
  for(i=0; i<circuit_list.size(); i++)
	{
    if (circuit_list[i]->getInstanceName() == iname)
      break;
  }
  assert(circuit_list[i]->getInstanceName() == iname);
  delete circuit_list[i];
  circuit_list.erase(circuit_list.begin() + i);
}


void CircuitManager::check() throw(string&)
{
  // Set up all circuits.
  for(unsigned i=0; i<circuit_list.size(); i++)
	{
    try
		{
      circuit_list[i]->check();
    }
    catch(string& errmsg)
		{
      throw("In circuit \"" + circuit_list[i]->getInstanceName() +
	    "\":\n" + errmsg);
    }
  }
}

void CircuitManager::init() throw(string&)
{
  // init all circuits.
  for(unsigned i=0; i<circuit_list.size(); i++)
	{
    try
		{
      circuit_list[i]->init();
    }
    catch(string& errmsg)
		{
      throw("In circuit \"" + circuit_list[i]->getInstanceName() +
	    "\":\n" + errmsg);
    }
  }
}


// The following are stack functions.
// They are used for subcircuit processing.
unsigned CircuitManager::push(const string& name)
{
  Circuit* cir = NULL;

  // Find if the circuit is already created.
  for(unsigned i=0; i<circuit_list.size(); i++)
	if (circuit_list[i]->getInstanceName() == name)
	{
		cir = circuit_list[i];
		break;
	}

  if (!cir)
    // Create a new subcircuit (flag == true)
	cir = createCircuit(name, true);

  // put the circuit at the end of the vector
  circuit_stack.push_back(cir);

  assert(circuit_stack.size() > 0);
  return circuit_stack.size();
}


unsigned CircuitManager::pop()
{
  if(circuit_stack.size() == 0)
    return 0;

  circuit_stack.erase(circuit_stack.end()-1);
  return circuit_stack.size();
}


Circuit* CircuitManager::getCurrent()
{
  if(circuit_stack.size() == 0)
    return current_circuit;

  return circuit_stack.back();
}

