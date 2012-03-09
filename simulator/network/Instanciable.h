// Instaciable is a simple base class for objects that can be instanciated,
// like elements, terminals, circuits.
// by Carlos E. Christoffersen

#ifndef Instanciable_h
#define Instanciable_h 1

#include "NetListItem.h"

class Instanciable
{
	public:
  Instanciable(const string& iname)
	{
		instance_name = iname;
	}
	
  ~Instanciable() {}
	
  // Get the instance name
  inline string getInstanceName() const 
	{
		return instance_name;
	}  
	
	private:
  // Name of this instance.
  string instance_name; 
};

#endif

