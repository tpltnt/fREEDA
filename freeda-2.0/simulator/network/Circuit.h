// Each circuit in the program (including subcircuits) is contained
// by an object of type Circuit.
// by Carlos E. Christoffersen

#ifndef Circuit_h
#define Circuit_h 1

#include "ElementManager.h"
#include <cassert>
#include <map>

// This class is required by the container "map".
struct ltunsigned
{
  bool operator()(const unsigned& id1, const unsigned& id2) const
  {
    return id1 < id2;
  }
};

typedef std::map<unsigned, Element*, ltunsigned> ElementMap;
typedef std::map<unsigned, Terminal*, ltunsigned> TerminalMap;

class Circuit : public Instanciable
{
	public:
  Circuit(const string& iname, bool subckt);
  ~Circuit();

  // Elements ---------------
  // Add one element to the circuit if it is not
  // already included. Returns element id.
  unsigned addElement(const string& elem_type, const string& iname,
	bool nopropagate = false) throw(string&);

  // remove element from the circuit
  void removeElement(const unsigned& id);

  // Connect one element to one terminal.
  void connect(const unsigned& elem_id, const unsigned& term_id);

  // Get one element by Id
  Element* getElement(const unsigned& id);

  // Get one element by name.
  Element* getElement(const string& elem_name);

  // Set the next element to the begining and set mask to
  // select elements.
  void setFirstElement(ElemFlag mask);

  // Get next element
  Element* nextElement();

  // Get the current number of elements in the circuit
  inline unsigned getNumberOfElements()
	{
    return elem_set.size();
  }


  // Terminals -------------
  // Create a new circuit terminal if it does not exist.
  // Returns terminal id.
  unsigned addTerminal(const string& term_name, bool nocheck = false);

  // Remove a terminal not connected to any element by Id
  void removeTerminal(const unsigned& id);

  // Set one terminal as a local reference.
  void setRefTerm(unsigned& id);

  // Get a terminal by name
  Terminal* getTerminal(const string& name);

  // Get a terminal by Id
  Terminal* getTerminal(const unsigned& id);

  // Set the next terminal to the begining
  inline void setFirstTerminal()
  {
    current_terminal_it = term_set.begin();
  }

  // Get next terminal
  Terminal* nextTerminal();

  // get the number of terminals of the circuit
  inline unsigned getNumberOfTerminals()
  {
    return term_set.size();
  }

  // General -------------
  // Initialize all the elements in circuit
  void init() throw(string&);

  // Set-up all circuit elements after the circuit is complete.
  void check() throw(string&);

  // Check consistency of local (global) reference terminals.
  void checkReferences() throw(string&);

  //--------------------------------------
  // Methods to treat subcircuit instances
  //--------------------------------------
  // Add one subcircuit to the instance list.
  // It would be nice if this could be done automagically,
  // but that would complicate element handling and it is
  // not such a big deal, I think.
  void addXinstance(Element* sub_instance)
  {
    subckt_list.push_back(sub_instance);
  }

  // Expand all subcircuits
  void expand();

  //--------------------
  // Subcircuit methods
  //--------------------

  // Returns true is this is a subcircuit
  inline bool isSubcircuit() const
  {
    return subcircuit;
  }

  // Add one Xsubckt user to this definition
  inline void addUser()
  {
    used_by++;
  }

  // Remove one user
  void delUser();

  // Ask the number of Xsubckt instances that use this definition
  inline int getNumberOfUsers() const
  {
    return used_by;
  }

  // Set the connection terminals by name
  void setConnection(const unsigned& id);

  // Get the number of subcircuit connections
  inline unsigned getNumberOfConnections() const
  {
    return connect_term.size();
  }

  // Get connection terminal from index number
  Terminal* getConnect(const unsigned& index);

  // Get the external terminal index given a pointer to it.
  unsigned getConnectIndex(const unsigned& id) const;

	private:
  // Vector to hold the local reference terminals of the circuit.
  // The circuit may also have only one global reference.
  TerminalVector ref_terms;

  // Check the connectivity of the circuit
  void checkTerms() throw(string&);

  // Current mask to get elements
  ElemFlag elem_mask;

  static ElementManager * the_EM;

  // keep a pointer to the element used in the last operation
  Element* last_element;

  // Set to hold the terminals.
  // Use a set to find by name quickly
  TerminalMap term_set;

  // Set to hold the elements.
  ElementMap elem_set;

  // iterator to hold the current element position (for method nextElement).
  ElementMap::iterator current_element_it;

  // iterator to hold the current terminal position (for method nextTerminal).
  TerminalMap::iterator current_terminal_it;

  //--------------------------------------
  // Subcircuit-instance related variables
  //--------------------------------------

  // vector to hold subcircuit instances
  ElementVector subckt_list;

  // Flag to mark when the subcircuits are expanded.
  bool flattened;

  //-----------------------------
  // Subcircuit-related variables
  //-----------------------------

  // Vector to hold the connection terminals for subcircuits.
  // It is empty for normal circuits.
  TerminalVector connect_term;

  // Flag to indicate if this is a subcircuit (true), and therefore
  // has external connections.
  bool subcircuit;

  // Number of Xsubckt that are using this definition
  int used_by;
};

#endif

