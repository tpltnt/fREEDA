#include "ResistorT.h"

// Static members
const unsigned ResistorT::n_par = 11;

// Element information
ItemInfo ResistorT::einfo =
{
  "resistort",
  "Resistor, Electro-thermal",
  "Houssam S. Kanj",
  DEFAULT_ADDRESS"category:lumped,electrothermal",
  "2002_06_20"
};

// Parameter information
ParmInfo ResistorT::pinfo[] =
{
  {"r0", "Resistance value (Ohms)", TR_DOUBLE, false},
  {"l", "length (meters)", TR_DOUBLE, false},
  {"w", "width (meters)", TR_DOUBLE, false},
  {"t", "system Temperature (Celsius)", TR_DOUBLE, false},
  {"rsh", "sheat resistance (Ohms/sq)", TR_DOUBLE, false},
  {"defw", "default device width (meters)", TR_DOUBLE, false},
  {"narrow", "narrowing due to side etching (meters)", TR_DOUBLE, false},
  {"tnom", "initial Temperature (Celsius)", TR_DOUBLE, false},
  {"tc1", "Temperature Coefficient (1/Celsius)", TR_DOUBLE, false},
  {"tc2", "Temperature Coefficient (1/Celsius)", TR_DOUBLE, false},
  {"pdr", "Power Depenent Resistor", TR_BOOLEAN, false}
};

ResistorT::ResistorT(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set parameters
  paramvalue[0] = &(r0);
  paramvalue[1] = &(l);
  paramvalue[2] = &(w);
  paramvalue[3] = &(t);
  paramvalue[4] = &(rsh);
  paramvalue[5] = &(defw = 1e-6);
  paramvalue[6] = &(narrow = zero);
  paramvalue[7] = &(tnom = 27);
  paramvalue[8] = &(tc1= zero);
  paramvalue[9] = &(tc2= zero);
  paramvalue[10] = &(pdr=false);
}

void ResistorT::init() throw(string&)
{
  if (!isSet(&w))
    w = defw;
  if (!isSet(&t))
    t = tnom;

  if (pdr)
	{
    // Set the number of terminals
    setNumTerms(4);
    // Set flags
    setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);
    // Set number of states
    setNumberOfStates(2);

    if (!isSet(&r0))
		{
      if (!isSet(&l))
				throw("r0 or l must be specified for the resistance");
      if (!isSet(&rsh))
				throw("rsh must be specified for the resistance");
    }
    DenseIntVector var(2);
    var[1] = 1;
    initializeAD(var);
  }
  else
	{  // pdr==false
    // Set the number of terminals
    setNumTerms(2);
    // Set flags
    setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);
    // calculate the resistor value
    if (isSet(&r0))
		{
      r = r0*(1+tc1*(t-tnom)+tc2*(t-tnom)*(t-tnom));
    }
    else
		{
      if (!isSet(&l))
				throw("r0 or l must be specified for the resistance");
      if (!isSet(&rsh))
				throw("rsh must be specified for the resistance");
      r = rsh*((l-narrow)/(w-narrow))*(1+tc1*(t-tnom)+tc2*(t-tnom)*(t-tnom));
    }
  }
}

void ResistorT::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: resistor voltage
  // x[1]: deltatemp in deg. Celsius

  effort[0]=x[0];
  effort[1]=x[1]+tnom+273; //effort[1]==tp[1] in Kelvin

  AD res;
  if (isSet(&r0))
	{
    res = r0 * (one + tc1 * x[1] + tc2 * x[1] * x[1]);
  }
  else
	{
    res = rsh * ((l-narrow) / (w-narrow))
		* (one + tc1*x[1] + tc2 * x[1] * x[1]);
  }

  flow[0] = x[0] / res;
  flow[1] = - x[0] * flow[0]; //flow[1]==pp[0];
}

void ResistorT::getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1)); // Local reference terminal
  local_ref_vec.push_back(1); // Local reference index

  if (pdr)
	{
    term_list.push_back(getTerminal(2));
    term_list.push_back(getTerminal(3)); // Local reference terminal
    local_ref_vec.push_back(3); // Local reference index
  }
}

void ResistorT::fillMNAM(FreqMNAM* mnam)
{
  // Ask my terminals the row numbers
  mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), one/r);
}

void ResistorT::fillMNAM(TimeMNAM* mnam)
{
  // Ask my terminals the row numbers
  mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), one/r);
}
