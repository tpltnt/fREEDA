#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "Inductor.h"

// Static members
const unsigned Inductor::n_par = 3;
const double Inductor::factor = 1e-2;

// Element information
ItemInfo Inductor::einfo =
{
  "inductor",
  "Inductor",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:lumped",
  "2003_05_15"
};

// Parameter information
ParmInfo Inductor::pinfo[] =
{
  {"l", "Inductance value (H)", TR_DOUBLE, true},
  {"int_res", "Internal Resistance value (Ohms)", TR_DOUBLE, false},
  {"time_d", "Flag, if true, calculate in the time domain.", TR_BOOLEAN, false}
};


Inductor::Inductor(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Value of ind is required
  paramvalue[0] = &ind;
  // Default value of internal resistor. (small)
  paramvalue[1] = &(int_res = 1e-8);
  paramvalue[2] = &(time_d = false);

  // Set the number of terminals
  setNumTerms(2);
  // Set number of states 
  setNumberOfStates(1);
  // Set flags
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);

  my_row = 0;
}

void Inductor::init() throw(string&)
{
  if (time_d) 
    setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN);
}

unsigned Inductor::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Keep the equation number assigned to this element
  my_row = eqn_number;
  // Add one extra RC
  return 1;

}

void Inductor::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row);
  first_eqn = my_row;
  n_rows = 1;
}

void Inductor::fillMNAM(FreqMNAM* mnam)
{
  assert(my_row);
  const double& freq(mnam->getFreq());
  // Ask my terminals the row numbers
  mnam->setOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);
  double_complex zl(-int_res, -twopi * freq * ind);
  mnam->setElement(my_row, my_row, zl);
}

void Inductor::fillMNAM(TimeMNAM* mnam)
{
  assert(my_row);
  // Ask my terminals the row numbers
  mnam->setMOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);
  mnam->setMElement(my_row, my_row, -int_res);
  mnam->setMpElement(my_row, my_row, -ind);
}
  

void Inductor::svTran(TimeDomainSV* tdsv)
{
  // The state variable is the current/factor (for numerical stability).
  tdsv->i(0) = tdsv->getX(0) * factor;
  // v = L di/dt
  tdsv->u(0) = ind * tdsv->getdX_dt(0) * factor + int_res * tdsv->i(0);
}

void Inductor::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = factor * (ind * tdsv->getdx_dtFactor() + int_res);
  tdsv->getJi()(0,0) = factor;
}
