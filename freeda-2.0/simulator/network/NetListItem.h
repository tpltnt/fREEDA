// NetListItem is the base class for everything that can be included 
// in a netlist.
// by Carlos E. Christoffersen

#ifndef NetListItem_h
#define NetListItem_h 1

#include <sstream>
#include <cassert>
#include <cstdio>
#include "../containers.h"

// This is the definition of the different types for parameters
// TR_INT, integer, valid netlist parameter
// TR_LONG, long integer
// TR_FLOAT, float
// TR_DOUBLE, double, valid netlist parameter
// TR_CHAR, single character, valid netlist parameter
// TR_STRING, string, valid netlist parameter
// TR_COMPLEX, complex number, valid netlist parameter
// TR_BOOLEAN, boolean integer, valid netlist parameter
// TR_INT_VECTOR, integer vector, valid netlist parameter
// TR_DOUBLE_VECTOR, double vector, valid netlist parameter
// TR_INT_MATRIX, integer matrix, valid netlist parameter
// TR_DOUBLE_MATRIX double matrix, valid netlist parameter
enum ParamType {TR_INT, TR_LONG, TR_FLOAT, TR_DOUBLE, TR_CHAR, TR_STRING,
	TR_COMPLEX, TR_BOOLEAN, TR_INT_VECTOR, TR_DOUBLE_VECTOR, TR_INT_MATRIX, 
  TR_DOUBLE_MATRIX};

#define DEFAULT_ADDRESS "http://www.freeda.org/"

// Structure to hold static information about the item.
struct ItemInfo 
{
  const char* name; // Name of netlist item
  const char* description; // Expanded name of netlist item
  const char* author; // Acknowledge the author
  const char* documentation; // documentation string
  const char* version; // include version number yyyy_mm_dd
};

// Structure to hold parameter information
struct ParmInfo 
{
  const char*      name;          // parameter name
  const char*      comment;       // Expanded comment about parameter
  ParamType  type;          // parameter data type
  bool       required;      // is this parameter required? (YES = true)
};


class NetListItem 
{
	public:
  
  //---------------------------------------
  // Retrieve Information
  //---------------------------------------
  inline string getName() const 
	{
    return string(i_info->name);
  }  
  
	inline string getDescription() const 
	{
    return string(i_info->description);
  } 
  
	inline string getAuthor() const 
	{
    return string(i_info->author);
  }
	
  inline string getDocumentation() const 
	{
    return string(i_info->documentation);
  }
	
	
  // ---------------------------------------
  // Parameter stuff
  // ---------------------------------------
	
  // returns the number of parameters
  inline unsigned getNumberOfParams() const 
	{
		return numparms;
	}
  
  // returns name and type of a parameter given the parameter index
  void getParamSpec(const unsigned& index, ParamType& p_type, 
	const char*& param_name) const;
	
  // returns a parameter type given the name of it
  ParamType askParamType(const string& parname) const throw(string);
  
  // Set parameters from other element
  void getParamsFrom(NetListItem* src_item);
	
  // Set the value of a parameter by index
  // The flag "copy" is to avoid special treatment with
  // strings and vectors. It is used by getParamsFrom().
  void setParam(const unsigned& index, const void* parvalue,
	const ParamType& p_type);
	
  // set parameter value by name
  void setParam(const string& parname, const void* parvalue, 
	const ParamType& type);
	
  // Return true if a parameters was set given its address
  bool isSet(const void* par_address) const;
	
  // Return true if a parameters was set given its index
  bool isSet(const unsigned& index) const;
	
  // Check element parameters. 
  // Returns true if parameters are Ok.
  void checkParams() const throw(string);
	
  // Returns a string with parameter description
  void getParamDesc(const unsigned& index,
	string& parname,
	string& comment,
	string& type,
	string& dflt_val,
	string& required) const;
	
	protected:
  // Only derived classes need this constructor
  NetListItem(ItemInfo* i_info, ParmInfo* param_desc, const int& numparms);
	
  // Destructor
  ~NetListItem();
	
  //---------------------------------------
  // Parameters
  //---------------------------------------
  
  // Later this should be made private and add more checking.
  // Parameter value vector
  void ** paramvalue;
	
	private:
  //---------------------------------------
  // Information
  //---------------------------------------
  ItemInfo* i_info;
	
  //---------------------------------------
  // Parameters
  //---------------------------------------
	
  // Number of parameters
  unsigned numparms;
  // Parameter information vector
  ParmInfo* param_desc;
  // Vector to check which parameters are set
  bool* param_set;
};

#endif

