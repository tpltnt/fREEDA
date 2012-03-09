#include "AbmMixer.h"

// Static members
const unsigned AbmMixer::n_par = 2;

// Element information
ItemInfo AbmMixer::einfo =
{
  "abmmixer",
  "Ideal behavioral mixer",
  "Mark Buff",
  DEFAULT_ADDRESS"category:behavioral",
  "2003_05_15"
};

// Parameter information
ParmInfo AbmMixer::pinfo[] =
{
  {"op", "Mixer operation", TR_INT, false},
  {"inputs", "Number of inputs to mixer", TR_INT,	false}
};
// operator defaults to multiply. values: 0=mul, 1=div, 2=add, 3=sub
// inputs defaults to two, no limit

AbmMixer::AbmMixer(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(op	   = 0);
  paramvalue[1] = &(inputs = 2);

  // Set the number of terminals, dependent on num of inputs
  setNumTerms(inputs + 2);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states, one less than num of terms
  setNumberOfStates(inputs + 1);
}

void AbmMixer::init() throw(string&)
{
  // Need to add conditional for more than 2 inputs and div/sub operation
  if((inputs != 2) && ((op == 1)||(op == 3)))
  {
		std::cout << std::endl << std::endl
      << "Division or subtraction not allowed with more than 2 inputs..." << std::endl;
		throw("Division or subtraction not allowed with more than 2 inputs.");
  }

  DenseIntVector var(inputs + 1);
  int count;
  for(count = 0; count <= inputs; count++)
	{
  	var[count] = count;
  }
  DenseIntVector 		novar;
  DenseDoubleVector 		nodelay;
  initializeAD(var, novar, novar, novar, nodelay);
}

void AbmMixer::eval(AD * x, AD * effort, AD * flow)
{
  // assign voltages directly
  for(int st_var = 0; st_var <= inputs; st_var++)
	{
		effort[st_var] = x[st_var];
  }

  // assign the input currents to 0
  for(int st_var = 0; st_var < inputs; st_var++)
	{
		flow[st_var] = zero;
  }

  // assign the output current
  flow[inputs] = x[0];  // will always have x[0]

  // ideal mixer
	for(int st_var = 1; st_var < inputs; st_var++)
	{
		if(op==0)                    // multiply
			flow[inputs] *= x[st_var];
		if(op==1)                    // divide
			flow[inputs] /= x[st_var];
		if(op==2)                    // add
			flow[inputs] += x[st_var];
		if(op==3)                    // subtract
			flow[inputs] -= x[st_var];
	}

  flow[inputs] = (flow[inputs]);
}
