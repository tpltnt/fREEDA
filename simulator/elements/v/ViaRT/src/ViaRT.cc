#include "ViaRT.h"

// Static members
const unsigned ViaRT::n_par = 9;

// Element information
ItemInfo ViaRT::einfo =
{
  "viart",
  "ViaRT, Electro-thermal via",
  "Kai Li",
  "category:lumped,electrothermal",
  "2008_04_21"
};

// Parameter information
ParmInfo ViaRT::pinfo[] =
{
  {"d", "Diameter (meters)", TR_DOUBLE, false},
  {"h", "Height (meters)", TR_DOUBLE, true},
  {"tm", "Thickness (meters)", TR_DOUBLE, false},
  {"rho", "Resistivity (Ohms-meters)", TR_DOUBLE, false},
  {"metal", "Metal: Silver, Copper, Gold, Aluminum", TR_STRING, false},
  {"t", "System Temperature (Celsius)", TR_DOUBLE, false},
  {"tnom", "Initial Temperature (Celsius)", TR_DOUBLE, false},
  {"tc", "Temperature Coefficient (1/Celsius)", TR_DOUBLE, false},
  {"pdr", "Power Depenent Resistor", TR_BOOLEAN, false}
  //{"idf", "ideality factor of resistance", TR_DOUBLE, false}
};

ViaRT::ViaRT(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set parameters
  paramvalue[0] = &(d = 1e-6 );
  paramvalue[1] = &(h);
  paramvalue[2] = &(tm = 3e-7 );
  paramvalue[3] = &(rho);
  paramvalue[4] = &(metal = "copper");
  paramvalue[5] = &(t);
  paramvalue[6] = &(tnom = 20);
  paramvalue[7] = &(tc = zero);
  paramvalue[8] = &(pdr = false);
  //paramvalue[9] = &(idf = 0.8);
}

void ViaRT::init() throw(string&)
{
  if (!isSet(&t))
    t = tnom;
  // calculate area of via
  double inner = d - 2*tm;
  if ( inner < 0 )
    inner = 0;
  area = pi / 4 * ( d * d - inner * inner);
  static const int M_SIZE = 4;
  string metals[M_SIZE] = {"silver","copper","gold","aluminum"};
  double metals_rho[M_SIZE] = { 1.63e-8, 1.83e-8, 2.44e-8, 2.64e-8 };
  double metals_tc[M_SIZE] =  {0.003819, 0.004041, 0.003715, 0.004308 };
  bool found =false;
  if (!isSet(&rho) || isSet(&metal))
  {
    for (int i = 0; i < M_SIZE; i++)
    {
      if ( metals[i] == metal )
      {
        rho = metals_rho[i];
        tc = metals_tc[i];
        found = true;
      }
    }
    if (!found)
      throw(getInstanceName() + ":rho or metal must be specified for the resistivity");
  }
  if (pdr)
	{
    // Set the number of terminals
    setNumTerms(4);
    // Set flags
    setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);
    // Set number of states
    setNumberOfStates(2);
    DenseIntVector var(2,0);
    var[1] = 1;
    initializeAD(var);
  }
  else
	{
    // pdr==false
    // Set the number of terminals
    setNumTerms(2);
    // Set flags
    setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);
    // calculate the resistor value
    r = rho * h / area * ( 1+tc*(t-tnom) );
  }
}

void ViaRT::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: resistor voltage
  // x[1]: deltatemp in deg. Celsius
  effort[0] = x[0];
  effort[1] = x[1]+tnom+273; //effort[1]==tp[1] in Kelvin

  AD res;
  res = rho * h / area * (one + tc*x[1]);

  flow[0] = x[0] / res;
  flow[1] = - x[0] * flow[0]; //flow[1]==pp[0];
}

void ViaRT::getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list)
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

void ViaRT::fillMNAM(FreqMNAM* mnam)
{
  // Ask my terminals the row numbers
  mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), one/r);
}

void ViaRT::fillMNAM(TimeMNAM* mnam)
{
  // Ask my terminals the row numbers
  mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), one/r);
}
