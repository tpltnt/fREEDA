#include "GraphNode.h"

unsigned GraphNode::next_id = 1;

GraphNode::GraphNode(const string& iname)
: Instanciable(iname) 
{
  // reserve space for 6 connections (enough for most nodes).
  node_list.reserve(6);
  node_idx = 0;
	
  // set the id number
  id = next_id++;
}

GraphNode::~GraphNode()
{ }

void GraphNode::unlink(GraphNode* gn)
{
  unsigned nn = node_list.size();
  for (unsigned i=0; i < nn; i++)
	{
    if (node_list[i] == gn) 
	  {
		  // Remove element from vector
		  // For speed, swap element to be removed with last element
		  // and remove last.
		  node_list[i] = node_list[nn-1];
		  node_list.pop_back();
		  return;
	  }
	}
  assert(false);
}


GraphNode* GraphNode::getGNode(const unsigned& i)
{
  assert(i < node_list.size());
  return node_list[i];
}

