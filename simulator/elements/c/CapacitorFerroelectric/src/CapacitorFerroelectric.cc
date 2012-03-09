#include "CapacitorFerroelectric.h"

// CapacitorFerroelectric
//
// This element models a non-linear ferroelectric capacitor.
//
// Netlist form
// capacitorferroelectric:cf1 n1 n2 <parameter list>
//
// A parallel plate physical model of the capacitor is considered
// although this concept can be further extended to model ferroelectric
// IDCs (inter-digitated capacitors) and gap capacitors.
// This model accounts for interfacial, bulk capacitance and fringing
// capacitance. It also considers thickness and temperature dependence of
// the high-permittivity ferroelectric material.
//
// This model can be used during time-domain analysis.

// Static members

// Set the number of parameters.  This number should be equal to the
// number of parameters in ParmInfo below.
const unsigned CapacitorFerroelectric::n_par = 11;

// Set the element information. First entry should be in lower case and
// defines the name of the element model. The rest is used in self
// documentation.
ItemInfo CapacitorFerroelectric::einfo =
{
  "capacitorferroelectric",
  "Ferroelectric capacitor",
  "Vrinda Haridasan",
  "category:capacitor",
  "2008_04_14"
};

// Parameter information
ParmInfo CapacitorFerroelectric::pinfo[] =
{
  {"epsid", "Interfacial capacitance density (F/m^2)", TR_DOUBLE, false},
  {"epsb0", "Bulk dielectric material zero-bias permittivity (F/m)", TR_DOUBLE, false},
  {"a", "Area of cross-section of the parallel plate capacitor (m^2)", TR_DOUBLE, true},
  {"d", "Total capacitor thickness (m)", TR_DOUBLE, true},
  {"t", "2 * (Interfacial Capacitor thickness at the plate dielectric interface) (m)", TR_DOUBLE, false},
  {"k", "Fringing capacitance constant (F)", TR_DOUBLE, false},
  {"alpha3", "Describes the non-linearity of the material in the Landau-Devonshire-Ginzburg model (m^2/C^2F)", TR_DOUBLE, false},
  {"T0", "Curie-Weiss temperature for a particular BST film thickness (deg C)", TR_DOUBLE, false},
  {"T", "Current Temperature of the sample (deg C)", TR_DOUBLE, false},
  {"beta", "Temperature Coefficient of Capacitance (TCC) at zero bias (ppm/deg C)", TR_DOUBLE, false},
  {"p", "Device periphery for fringing capacitance calculations (m)", TR_DOUBLE, true}
};

// The creator. Called when creating a new element instance.
CapacitorFerroelectric::CapacitorFerroelectric(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Parameters are stored in the paraminfo vector. Here we set default
  // values and the working name of the parameter. If all parameters
  // are not set in the netlist element call, default parameter values
  // will be used for the analysis.
  paramvalue[0] = &(epsid = 32e-3);
  paramvalue[1] = &(epsb0 = 170);
  paramvalue[2] = &a;
  paramvalue[3] = &d;
  paramvalue[4] = &(t = 5e-9);
  paramvalue[5] = &(k = 1.6e-15);
  paramvalue[6] = &(alpha3 = 3.3e-3);
  paramvalue[7] = &(T0 = -167);
  paramvalue[8] = &(T = -73);
  paramvalue[9] = &(beta = 600);
  paramvalue[10] = &p;

  // Set the number of terminals that appears in the netlist.
  setNumTerms(2);

  // Set flags to tell the analysis routines what is special about this
  // element. (Here we are setting bits.)
  // NONLINEAR = this is a non-linear element.
  // ONEREF = this element has one reference terminal.
  // TR_TIME_DOMAIN = this element is used in transient domain analysis.
  // Both NONLINEAR and TR_TIME_DOMAIN affect the type of analysis to be used.
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
  // Set number of states
  setNumberOfStates(1);

}

// The initialization routine.  Called prior to an analysis.  Here the circuit
// can be changed and flags reset. All the parameters that are unchanged for
// a single simulation run are initialized here depending on the netlist element call.
void CapacitorFerroelectric::init() throw(string&)
{
  ci = epsid * a;
  cf = k*p/d;
  epsb = epsb0 * (1 - beta/1000000 *(T - T0));
  alpha1 = 1/epsb;
  cbmax = alpha1*(d-t)/a;
  cmax = 1/(1/ci + 1/cbmax) + cf;
  cmaxinv = 1/cmax;
  cnst = alpha3*(d-t)/(a*a*a);
  DenseIntVector var(1,0);
  initializeAD(var, var);
}

// The eval() routine is called by createTape() from inside AdolcElement.cc file.
// The eval() function is the implementation of the device equations.
// It is in here that the voltages and currents are related to
// the state variables.
void CapacitorFerroelectric::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: state variable --> q
  // x[1] --> i = dq/dt: time derivative of x[0]
	AD vfnq, cfnq, dvfnq_dx;
	//vfnq = x[0] * cmaxinv;
  vfnq = x[0]  * cmaxinv + cnst*x[0]*x[0]*x[0];
 	//cfnq = 1/((cmaxinv) + cnst * 3*x[0]*x[0]);
 	//dvfnq_dx = cmaxinv + cnst*3*x[0]*x[0];
 	effort[0] = vfnq;

 	flow[0] = x[1];

 	//flow[0] = cfnq*dvfnq_dx*x[1];
}
