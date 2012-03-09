#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "R.h"

// The native fREEDA netlist format is
// R:1 n1 n2 r=50
// where
// R indicates the model to use
// R:1 is the name of this instance of the R model
// n1 is the first terminal (Can be a string or an integer)
// n2 is the second terminal (Can be a string or an integer)
// r is the parameter (there is only one and there is not a default value)

// Static members
// Set the number of parameters
const unsigned R::n_par = 1;

// Element information
ItemInfo R::einfo =
{
  "r",
  "Resistor",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:lumped",
  "2001_06_15"
};

// Parameter information
// The size of pinfo[] is automatically assigned here.
// The information here enables automatic documentation.
// So here are the parameters of the "r" element
// "r" the name of the zeroth parameter. This is the name that must
// appear in the netlist such as "r = 50"
// "Resistance value (Ohms)" This is the string used in documentation. It is
// a fuller name of the element with units.
// "r" the name of the zeroth parameter
// "TR_DOUBLE" This indicates the type of the element. In this case it is
// a double. The full list of valid input parameter types
// is given in .../network/NetListItem.h
// "true" indicates that the parameter is required and must be specified in
// the netlist. The alternative is "false" indicating that the parameter
// is optional in which case a default value must be sepcified.
ParmInfo R::pinfo[] =
{
  {"r", "Resistance value (Ohms)", TR_DOUBLE, true}
};


R::R(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Now list all of the parameters and give them meaningful variable names
  // to use in the codea These names are the same names used in the header file.
  // Value of r is required so there is not a default value.
  paramvalue[0] = &res;

  // Set the number of terminals
  setNumTerms(2);

  // Set flags
  // The flags indicate
  // Linear = This is a linear element
  // ONE_REF = There is only one reference terminal
  // TR_FREQ_DOMAIN = Element works in both the time domain and
  // frequency domain. (As all good elements should.)
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);
}


void R::fillMNAM(FreqMNAM* mnam)
{
  // Set the 4 element conductance stamp using a common routine.
  // There are two terminals and we need to get the row and column index
  // corresponding to each terminal.
  // one/res is the value to add to the entry in the MNA
  mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
	one/res);
}

void R::fillMNAM(TimeMNAM* mnam)
{
  // Set the 4 element conductance stamp using a common routine.
  // There are two terminals and we need to get the row and column index
  // corresponding to each terminal.
  // one/res is the value to add to the entry in the MNA
  mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
	one/res);
}

