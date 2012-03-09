#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Vpulse.h"
#include "../../../../analysis/TimeDomainSV.h"

// Static members
const unsigned Vpulse::n_par = 7;

// Element information
ItemInfo Vpulse::einfo =
{
  "vpulse",
  "Pulsed voltage source (transient analysis only)",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:source",
  "2003_08_28" // Fixed td behavior, taking it out of repeating per.  FPH.
};

// Parameter information
ParmInfo Vpulse::pinfo[] =
{
  {"v1", "Initial value (V)", TR_DOUBLE, false},
  {"v2", "Pulsed value (V)", TR_DOUBLE, false},
  {"td", "Delay time (s)", TR_DOUBLE, false},
  {"tr", "Rise time (s)", TR_DOUBLE, false},
  {"tf", "Fall time (s)", TR_DOUBLE, false},
  {"pw", "Pulse width (s)", TR_DOUBLE, false},
  {"per", "Period (s)", TR_DOUBLE, false}
};

Vpulse::Vpulse(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(v1 = zero);
  paramvalue[1] = &(v2 = zero);
  paramvalue[2] = &(td = zero);
  paramvalue[3] = &(tr = zero);
  paramvalue[4] = &(tf = zero);
  paramvalue[5] = &(pw = zero);
  paramvalue[6] = &(per = zero);

  // Set the number of terminals
  setNumTerms(2);

  // Set flags
  setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN | SOURCE);
  // Set number of states
  setNumberOfStates(1);

  my_row = 0;
  int_per_time = zero;
}

unsigned Vpulse::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Keep the equation number assigned to this element
  my_row = eqn_number;
  // Add one extra RC
  return 1;
}

void Vpulse::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row);
  first_eqn = my_row;
  n_rows = 1;
}

void Vpulse::fillMNAM(FreqMNAM* mnam)
{
  assert(my_row);
  // Ask my terminals the row numbers
  mnam->setOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);

  // This element behaves as a short circuit for AC/HB analysis.
  // Since the source vector is assumed to be initialized to zero,
  // we do not need to fill it if freq != frequency.
  return;
}


void Vpulse::fillMNAM(TimeMNAM* mnam)
{
  assert(my_row);
  // Ask my terminals the row numbers
  mnam->setMOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);
  return;
}


void Vpulse::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime();
  double reltime = ctime - td - int_per_time;
  if (reltime > per || reltime <0) {
    int_per_time = int((ctime - td)/per) * per;
    reltime = ctime - td - int_per_time;
  }
  double e = zero;
  if (int_per_time == zero && ctime < td)
    e = v1;
  else {
    if (reltime < tr)
      e = v1 + (v2 - v1) * ((reltime)/tr);
    else if (reltime < tr + pw)
      e = v2;
    else if (reltime < tr + pw + tf)
      e = v2 + (v1 - v2) * (reltime - tr - pw) / tf;
    else
      e = v1;
  }
  mnam->setSource(my_row, e);
  return;
}


void Vpulse::svTran(TimeDomainSV* tdsv)
{
  // Calculate voltage
  double& e = tdsv->u(0);

  if (tdsv->DC())
    e = v1;
  else {
    const double& ctime = tdsv->getCurrentTime();
    double reltime = ctime - td - int_per_time;
    if (reltime > per || reltime <0) {
      int_per_time = int((ctime-td)/per) * per;
      reltime = ctime - td - int_per_time;
    }
    double e = zero;
    if (int_per_time == zero && ctime < td)
      e = v1;
    else {
      if (reltime < tr)
        e = v1 + (v2 - v1) * (reltime/tr);
      else if (reltime < tr + pw)
        e = v2;
      else if (reltime < tr + pw + tf)
        e = v2 + (v1 - v2) * (reltime - tr - pw) / tf;
      else
        e = v1;
    }
  }
  // Scale state variable for numerical stability
  tdsv->i(0) = tdsv->getX(0) * 1e-2;
}

void Vpulse::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = zero;
  tdsv->getJi()(0,0) = 1e-2;
}

