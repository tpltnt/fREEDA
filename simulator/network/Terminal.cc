#include "GraphNode.h"
#include "Element.h"

Terminal::Terminal(const string& iname) : GraphNode(iname) 
{
  external = false;
  myrow = 0;
  my_group = 0;
  termdata = new TerminalData();
  //  x = y = z = 0.;
}

Terminal::~Terminal()
{
  delete termdata;
}

bool Terminal::checkConnect() const 
{
  unsigned count = getGraphNodeCount();
  if (external)
    return (count != 0);
  else
    return (count > 1);
}

void Terminal::setGroup(unsigned ref_id) throw(string&)
{
  // If the terminal has been visited or is in process, just return.
  if (my_group == ref_id)
    return;
  else 
	{ 
    // If the group has been set that implies an error in the network
    // topology.
    if (my_group)
      throw (getInstanceName() + ": conflict between reference terminals.");
  }
	
  // So this is the first time this terminal is called
  my_group = ref_id;
  // Now call all the elements in the adjacency list
  unsigned nelem = getGraphNodeCount();
  Element* elem;
  for (unsigned i = 0; i < nelem; i++) 
	{
    elem = (Element*)(getGNode(i));
    elem->setGroup(ref_id);
  }
}

void Terminal::propagate()
{
  // This method is used to start the propagation from the reference 
  // terminals.
  assert(isRef());
	
  // Now call all the elements in the adjacency list
  unsigned nelem = getGraphNodeCount();
  Element* elem;
  for (unsigned i = 0; i < nelem; i++) 
	{
    elem = (Element*)(getGNode(i));
    elem->setGroup(my_group);
  }
}

