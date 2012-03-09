#ifndef Terminal_h
#define Terminal_h 1

#include "TerminalData.h"

class Terminal : public GraphNode
{
	public:
	
  Terminal(const string& iname);
  ~Terminal();
  
  // attach element to terminal
  inline void attachElement(GraphNode* element) 
	{
		link(element);
	}
	
  // detach element from terminal
  inline void detachElement(GraphNode* element) 
	{
		unlink(element);
	}
	
  // check if the terminal is connected to more that one element
  bool checkConnect() const;
	
  // Return the number of elements connected to terminal
  inline unsigned getNumberOfConnect() const
	{
		return getGraphNodeCount();
	}
	
  // Make the terminal external
  inline void setExtern()
	{
		external = true;
	}
	
  // Returns true if the terminal is external
  inline bool isExternal() const
	{
		return external;
	}
	
  // Set row/column number for NAM
  inline void setRC(unsigned n)
	{
		myrow = n;
	}
	
  // Get the row/column number for NAM
  inline unsigned getRC() const
	{
		return myrow;
	}
	
  // Get terminal data instance
  inline TerminalData* getTermData()
	{
		return termdata;
	}
	
  // Set the group of the terminal and propagate.
  void setGroup(unsigned ref_id) throw(string&);
	
  // Returns the group of the terminal.
  inline unsigned getGroup() const
  {
    return my_group;
  }
	
  // Init the propagation if the terminal is a reference.
  void propagate();
	
  // Do not check if terminal is floating
  // Use this to avoid checking of internal nodes
  // (created by the expansion of an element model).
  inline void noCheck()
  {
    // getGroup() does not returns 0
    my_group = getId() + 1;
  }
	
  // Set this instance as a local reference terminal.
  inline void setRef()
  {
    my_group = getId();
  }
  
  inline bool isRef() const
  {
    return (my_group == getId());
  }

  // begin set and get XYZ modifications
  inline void setX(double X)
  {
    xpos = X;
  }

  inline double getX() const
  {
    return xpos;
  }

  inline void setY(double Y)
  {
    ypos = Y;
  }

  inline double getY() const
  {
    return ypos;
  }

  inline void setZ(double Z)
  {
    zpos = Z;
  }

  inline double getZ() const
  {
    return zpos;
  }
	
	private:
	
  // Local group where the terminal belongs
  unsigned my_group;
	
  // Output variables
  TerminalData* termdata;
	
  // row/column number in NAM
  unsigned myrow;
	
  // Flag to mark if the terminal is 
  // external => it is allowed to have only one connection
  bool external;

  // positional variables
  double xpos, ypos, zpos;
	
};

// Define the terminal vector type
typedef std::vector<Terminal*> TerminalVector;

#endif

