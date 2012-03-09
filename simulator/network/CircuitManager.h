// Keep track of all the circuits created
// by Carlos E. Christoffersen

#ifndef CircuitManager_h
#define CircuitManager_h 1

#include "Circuit.h"
#include <cassert>

typedef std::vector<Circuit*> CircuitVector;

class CircuitManager
{
	public:
  
  ~CircuitManager();
	
  // Only one element manager can exist.
  static CircuitManager* getCircuitManager();
  
  // Get a Circuit by name. Create if it does not exist.
  Circuit* getCircuit(const string& iname);
	
  // Get a subcircuit by name. Create if it does not exist.
  Circuit* getSubCircuit(const string& iname);
	
  // Remove one circuit
  void deleteCircuit(const string& iname);
	
  // Setup circuits and check for errors
  void check() throw(string&);
	
  // Init circuits
  void init() throw(string&);
	
  // Get total number of circuits and subcircuits
  inline unsigned getNumberOfCircuits() const
	{
		return circuit_list.size();
	}
	
  // Get circuit by index
  Circuit* getCircuit(unsigned& i);
	
  // The following are stack functions.
  // They are used for subcircuit processing.
	
  // push one circuit in the stack
  unsigned push(const string& name);
	
  // pull one circuit from the stack
  unsigned pop();
	
  // get the circuit at the top of the stack.
  Circuit* getCurrent();
	
	private:
  
  CircuitManager();
	
  // Create a new Circuit (or subcircuit) object.
  Circuit* createCircuit(const string& iname, bool subckt = false);
	
  static CircuitManager* cm;
	
  // Pointer to hold the current circuit not in the stack
  // i.e. not subcircuit.
  Circuit* current_circuit;
	
  // The keys are the instance names.
  CircuitVector circuit_list;
	
  // This vector is to implement the stack
  CircuitVector circuit_stack;
};

// Global pointer to the CM instance.
extern CircuitManager* the_CM;

#endif

