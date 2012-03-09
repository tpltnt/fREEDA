// The function of this class has been reduced to only create elements
// The management is made inside each circuit instance.
// by Carlos E. Christoffersen

#ifndef ElementManager_h
#define ElementManager_h 1

#include "Element.h"
#include <fstream>
#include <cassert>

class ElementManager
{
	public:
  ~ElementManager();
	
  // It can exist only one element manager
  static ElementManager* getElementManager();
  
  // Create a new element object
  Element* createElement(const string& elem_type, const string& iname)
	throw(string&);
	
  // Print element catalog to a file (rudimentary)
  void printCatalog() const;
  void printCatalogElement(char* elementName, int launchBrowserFlag) const;
	
	private:
  ElementManager();
  static ElementManager* em;
};

#endif

