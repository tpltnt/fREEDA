// This is a base class for elements and terminals, which are stored 
// in a graph.
// by Carlos E. Christoffersen

#ifndef GraphNode_h
#define GraphNode_h 1

#include <cassert>
#include "Instanciable.h"

class GraphNode : public Instanciable
{
	public:
  inline unsigned getId() const
	{
		return id;
	}

	protected:
  GraphNode(const string& iname);
  
  ~GraphNode();
	
  // Link with other node (add the other to the list)
  inline void link(GraphNode* gn)
	{
		node_list.push_back(gn);
	}
	
  // Remove the other node from the list
  void unlink(GraphNode* gn);
	
  // Return number of connections
  inline unsigned getGraphNodeCount() const
	{
		return node_list.size();
	}
	
  // Get terminal by index
  GraphNode* getGNode(const unsigned& i);
	
	private:
  // List of nodes connected
  std::vector<GraphNode*> node_list;
  
  // counter for nextGNode()
  unsigned node_idx;
	
  // Id of the node.
  unsigned id;
  // Next Id number
  static unsigned next_id;

};

#endif

