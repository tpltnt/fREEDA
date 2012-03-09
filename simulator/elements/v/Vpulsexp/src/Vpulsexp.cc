#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Vpulsexp.h"
#include "../../../../analysis/TimeDomainSV.h"

// Static members
const unsigned Vpulsexp::n_par = 7;

// Element information
ItemInfo Vpulsexp::einfo =
{
  "vpulsexp",
  "Pulsed voltage source with exponential taper (transient analysis only)",
  "Frank P. Hart"
  DEFAULT_ADDRESS"category:source",
  "2003_08_28"
};

// Parameter information
ParmInfo Vpulsexp::pinfo[] =
{
  {"v1", "Initial value (V)", TR_DOUBLE, false},
  {"v2", "Pulsed value (V)", TR_DOUBLE, false},
  {"td", "Delay time (s)", TR_DOUBLE, false},
  {"tr", "Rise time (s)", TR_DOUBLE, false}, // 10% to 90%
  {"tf", "Fall time (s)", TR_DOUBLE, false}, // 90% to 10%
  {"pw", "Pulse width (s)", TR_DOUBLE, false},
  {"per", "Period (s)", TR_DOUBLE, false}
};


Vpulsexp::Vpulsexp(const string& iname) : Element(&einfo, pinfo, n_par, iname)
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

unsigned Vpulsexp::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Keep the equation number assigned to this element
  my_row = eqn_number;
  // Add one extra RC
  return 1;
}

void Vpulsexp::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row);
  first_eqn = my_row;
  n_rows = 1;
}

void Vpulsexp::fillMNAM(FreqMNAM* mnam)
{
  assert(my_row);
  // Ask my terminals the row numbers
  mnam->setOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);

  // This element behaves as a short circuit for AC/HB analysis.
  // Since the source vector is assumed to be initialized to zero,
  // we do not need to fill it if freq != frequency.
  return;
}

void Vpulsexp::fillMNAM(TimeMNAM* mnam)
{
  assert(my_row);
  // Ask my terminals the row numbers
  mnam->setMOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);
  return;
}

void Vpulsexp::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime();
  // Added the stuff below.  8/27/03 FPH.
  vn = v2 - v1;	// normalized non-zero voltage.  Other, by def, is zero.
  tau_r = tr/(2.0*log(5.0));
  tau_f = tf/(2.0*log(5.0));
  pwh = pw;
  if (per <= zero) { // per not specified, so one-shot
    if (tau_r > tau_f)
      pwl = 14.0*tau_r;
    else
      pwl = 14.0*tau_f;
	}
  else
    pwl = per - tf - pw - tr; // per specified, pulse width low.

  // Added td to the expression below. FPH 8/26/03.
  double reltime = ctime - td - int_per_time;
  if (reltime > per) {
    int_per_time = int((ctime - td)/per) * per;
    reltime = ctime - td - int_per_time;
  }
  double e = zero;
  if (int_per_time == zero && ctime < td)
    e = v1;
  else {
    // function 1 - increasing with + 2nd deriv.
    if (reltime < pwl/2.0 + tr/2.0)
      e = v1 +
		0.5*vn*exp((reltime - pwl/2.0 - tr/2.0)/tau_r);

    // function 2 - increasing with - 2nd deriv.
    else if (reltime < pwl/2.0 + tr + pwh/2.0)
      e = v1 + 0.5*vn +
	  0.5*vn*(1-exp((-reltime + pwl/2.0 + tr/2.0)/tau_r));

    // function 3 - decreasing with - 2nd deriv.
    else if (reltime < pwl/2.0 + tr + pwh + tf/2.0)
      e = v1 + 0.5*vn +
	  0.5*vn*(1-exp((reltime - pwl/2.0 - tr - pwh - tf/2.0)/tau_f));

    // function 4 - decreasing with + 2nd deriv.
    else if (reltime < pwl + tr + pwh +tf)
      e = v1 +
	  0.5*vn*exp((-reltime + pwl/2.0 + tr + pwh + tf/2.0)/tau_f);

    // if get here, should be between function 4 and function 1, so v1.
    else
      e = v1;
  }
  // cout << "ctime = " << ctime << ";   e = " << e << endl;
  mnam->setSource(my_row, e);
  return;
}

void Vpulsexp::svTran(TimeDomainSV* tdsv)
{
  // Calculate voltage
  double& e = tdsv->u(0);

  if (tdsv->DC())
    e = v1;
  else {
    // Added the stuff below.  8/27/03 FPH.
    vn = v2 - v1;
    tau_r = tr/(2.0*log(5.0));
    tau_f = tf/(2.0*log(5.0));
    pwh = pw;
    if (per <= zero) { // per not specified, so one-shot
      if (tau_r > tau_f)
        pwl = 14.0*tau_r;
      else
        pwl = 14.0*tau_f;
		}
    else
      pwl = per - tf - pw - tr; // per specified, pulse width low.

    const double& ctime = tdsv->getCurrentTime();
    double reltime = ctime - td - int_per_time;
    if (reltime > per) {
      int_per_time = int((ctime-td)/per) * per;
      reltime = ctime - td - int_per_time;
    }
    double e = zero;
    if (int_per_time == zero && ctime < td)
      e = v1;
    else {
			// function 1 - increasing with + 2nd deriv.
			if (reltime < pwl/2.0 + tr/2.0)
				e = v1 +
			0.5*vn*exp((reltime - pwl/2.0 - tr/2.0)/tau_r);

			// function 2 - increasing with - 2nd deriv.
			else if (reltime < pwl/2.0 + tr + pwh/2.0)
				e = v1 + 0.5*vn +
			0.5*vn*(1-exp((-reltime + pwl/2.0 + tr/2.0)/tau_r));

			// function 3 - decreasing with - 2nd deriv.
			else if (reltime < pwl/2.0 + tr + pwh + tf/2.0)
				e = v1 + 0.5*vn +
			0.5*vn*(1-exp((reltime - pwl/2.0 - tr - pwh - tf/2.0)/tau_f));

			// function 4 - decreasing with + 2nd deriv.
			else if (reltime < pwl + tr + pwh +tf)
				e = v1 +
			0.5*vn*exp((-reltime + pwl/2.0 + tr + pwh + tf/2.0)/tau_f);

			// if get here, should be between function 4 and function 1, so v1.
			else
        e = v1;
    }
  }
  // Scale state variable for numerical stability
  tdsv->i(0) = tdsv->getX(0) * 1e-2;
}

void Vpulsexp::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = zero;
  tdsv->getJi()(0,0) = 1e-2;
}

