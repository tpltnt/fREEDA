#include "NetListItem.h"

NetListItem::NetListItem(ItemInfo* i_info,
ParmInfo* param_desc, const int& numparms)
{
  this->param_desc = NULL;
  this->paramvalue = NULL;
  this->i_info = i_info;
  this->numparms = numparms;
  if (numparms > 0)
	{
    paramvalue = new void*[numparms];
    this->param_desc = param_desc;
    if (numparms)
		{
      // allocate memory for param_check
      param_set = new bool[numparms];
      // Set flags to false
      for (int i=0; i<numparms; i++)
				param_set[i] = false;
    }
  }
}

NetListItem::~NetListItem()
{
  if (numparms)
	{
    // Free strings and vectors
    // -------------------------------------
    // Not needed with current types
    // -------------------------------------
		//      for (unsigned i = 0; i < numparms; i++) {
			//        switch (param_desc[i].type) {
				//        case CHARARRAY:
				//  	delete [] (char*)(paramvalue[i]);
				//  	break;
				//        default:
				//  	;	// do nothing
			//        }
		//    }
    // Free vector with addreses
    delete [] paramvalue;
    // free vector with flags
    delete [] param_set;
  }
}

void NetListItem::getParamSpec(const unsigned& index, ParamType& p_type,
const char*& param_name) const
{
  assert(numparms);
  p_type = param_desc[index].type;
  param_name = param_desc[index].name;
}


ParamType NetListItem::askParamType(const string& parname)
const throw(string)
{
  // check the number of parameters.
  if (!numparms)
    throw("Item " + getName() + " has no parameters.");

  for (unsigned i = 0; i < numparms; i++)
    if (string (param_desc[i].name) == parname)
			return param_desc[i].type;

  // if control reaches here, the parameter does not exists
  throw(getName() + ": unknown parameter.");
}

void NetListItem::getParamsFrom(NetListItem* src_item)
{
  // Both items must be of the same type
  assert(getName() == src_item->getName());

  // Set all parameters. Use flag to allow direct copy
  // of vectors and strings.
  for (unsigned i=0; i < numparms; i++)
	{
    setParam(i, src_item->paramvalue[i], param_desc[i].type);
    param_set[i] = src_item->param_set[i];
  }
}


// Here is where the dirty work is done.
// Be careful.
void NetListItem::setParam(const unsigned& index, const void* parvalue,
const ParamType& p_type)
{
  assert(numparms);
  bool flag = true;
  // By now, the parser need to be sure about the type.
  assert(param_desc[index].type == p_type);
  // Set the parameter value
  switch (p_type)
	{
		case TR_INT:
    *(int*)(paramvalue[index]) = *(int*)(parvalue);
    break;
		case TR_LONG:
    *(long*)(paramvalue[index]) = *(long*)(parvalue);
    break;
		case TR_FLOAT:
    *(float*)(paramvalue[index]) = *(float*)(parvalue);
    break;
		case TR_DOUBLE:
    *(double*)(paramvalue[index]) = *(double*)(parvalue);
    break;
		case TR_CHAR:
    *(char*)(paramvalue[index]) = *(char*)(parvalue);
    break;
		case TR_STRING:
    *(string*)(paramvalue[index]) = *(string*)(parvalue);
    break;
		case TR_COMPLEX:
    *(double_complex*)(paramvalue[index]) = *(double_complex*)(parvalue);
    break;
		case TR_BOOLEAN:
    *(bool*)(paramvalue[index]) = *(bool*)(parvalue);
    break;
		case TR_INT_VECTOR:
    if (((DenseIntVector*)(parvalue))->length())
      (*(DenseIntVector*)(paramvalue[index])) = (*(DenseIntVector*)(parvalue));
    break;
		case TR_DOUBLE_VECTOR:
    if (((DenseDoubleVector*)(parvalue))->length())
      (*(DenseDoubleVector*)(paramvalue[index])) = (*(DenseDoubleVector*)(parvalue));
    break;
		case TR_INT_MATRIX:
    if (((IntDenseMatrix*)(parvalue))->numRows())
      (*(IntDenseMatrix*)(paramvalue[index])) = (*(IntDenseMatrix*)(parvalue));
    break;
		case TR_DOUBLE_MATRIX:
    if (((DoubleDenseMatrix*)(parvalue))->numRows())
      (*(DoubleDenseMatrix*)(paramvalue[index])) = (*(DoubleDenseMatrix*)(parvalue));
    break;
		default:
    flag = false;
  }
  assert(flag == true);
  // Mark the parameter as set
  param_set[index] = true;
}

void NetListItem::getParamDesc(const unsigned& index,
string& parname,
string& comment,
string& type,
string& dflt_val,
string& required) const
{
  // First, check that the index is valid
  assert(index < numparms);

  parname = param_desc[index].name;
  comment = param_desc[index].comment;
  required = (param_desc[index].required) ? "yes" : "no";

  // Create an output string
  char s[200];

  // No choice but to switch for each case
  switch (param_desc[index].type)
	{
		case TR_INT:
    type = "INTEGER";
    if (!param_desc[index].required)
		{
      sprintf(s, "%d", *(int*)(paramvalue[index]) );
      dflt_val = s;
    }
    else
      dflt_val = "n/a";
    break;
		case TR_LONG:
    type = "LONG INTEGER";
    if (!param_desc[index].required)
		{
      sprintf(s, "%ld", *(long*)(paramvalue[index]) );
      dflt_val = s;
    }
    else
      dflt_val = "n/a";
    break;
		case TR_FLOAT:
    type = "FLOAT";
    if (!param_desc[index].required)
		{
      sprintf(s, "%g", *(float*)(paramvalue[index]) );
      dflt_val = s;
    }
    else
      dflt_val = "n/a";
    break;
		case TR_DOUBLE:
    type = "DOUBLE";
    if (!param_desc[index].required)
		{
      sprintf(s, "%g", *(double*)(paramvalue[index]) );
      dflt_val = s;
    }
    else
      dflt_val = "n/a";
    break;
		case TR_CHAR:
    type = "CHAR";
    if (!param_desc[index].required)
		{
      sprintf(s, "%c", *(char*)(paramvalue[index]) );
      dflt_val = s;
    }
    else
      dflt_val = "n/a";
    break;
		case TR_STRING:
    type = "STRING";
    if (!param_desc[index].required)
		{
      sprintf(s, "%s", ((string*)(paramvalue[index]))->c_str() );
      dflt_val = s;
    }
    else
      dflt_val = "n/a";
    break;
		case TR_COMPLEX:
    type = "DOUBLE COMPLEX";
    if (!param_desc[index].required)
		{
      sprintf(s, "%g + j %g",
			((double_complex*)(paramvalue[index]))->real(),
			((double_complex*)(paramvalue[index]))->imag());
      dflt_val = s;
    }
    else
      dflt_val = "n/a";
    break;
		case TR_BOOLEAN:
    type = "BOOLEAN";
    if (!param_desc[index].required)
		{
      if (*(bool*)(paramvalue[index]))
				sprintf(s, "true");
      else
				sprintf(s, "false");
      dflt_val = s;
    }
    else
      dflt_val = "n/a";
    break;
		case TR_INT_VECTOR:
    type = "INTEGER VECTOR";
    if (!param_desc[index].required)
		{
      sprintf(s, "See source file.");
      dflt_val = s;
    }
    else
      dflt_val = "n/a";
    break;
		case TR_DOUBLE_VECTOR:
    type = "DOUBLE VECTOR";
    if (!param_desc[index].required)
		{
      sprintf(s, "See source file.");
      dflt_val = s;
    }
    else
      dflt_val = "n/a";
    break;
		case TR_INT_MATRIX:
    type = "INTEGER MATRIX";
    if (!param_desc[index].required)
		{
      sprintf(s, "See source file.");
      dflt_val = s;
    }
    else
      dflt_val = "n/a";
    break;
		case TR_DOUBLE_MATRIX:
    type = "DOUBLE MATRIX";
    if (!param_desc[index].required)
		{
      sprintf(s, "See source file.");
      dflt_val = s;
    }
    else
      dflt_val = "n/a";
    break;
		default:
    type = "Unknown";
    dflt_val = "n/a";
  }
}


void NetListItem::setParam(const string& parname, const void* parvalue,
const ParamType& type)
{
  // numparms must be different than zero
  assert(numparms);

  // Flag to mark if the parameter was found
  bool flag = false;

  for (unsigned i = 0; i < numparms && !flag; i++)
	{
    if (param_desc[i].name == parname)
		{
      setParam(i, parvalue, type);
      flag = true;
    }
  }
  assert(flag == true);
}

bool NetListItem::isSet(const void* par_address) const
{
  unsigned i;
  // Search parameter and check flag
  for (i = 0; i < numparms; i++)
		if (paramvalue[i] == par_address)
			break;
  assert(i < numparms);

  return param_set[i];
}

bool NetListItem::isSet(const unsigned& index) const
{
  assert(index < numparms);
  return param_set[index];
}


void NetListItem::checkParams() const throw(string)
{
  // Go through all the parameters and check flags
  for (unsigned i = 0; i < numparms; i++)
    if (!param_set[i] && param_desc[i].required)
      throw(string("Missing required parameter: ") + param_desc[i].name
	+ " - " + param_desc[i].comment);
}

