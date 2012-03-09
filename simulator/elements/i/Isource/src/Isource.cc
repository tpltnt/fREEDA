#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "Isource.h"

// Static members
const unsigned Isource::n_par = 7;

// Element information
ItemInfo Isource::einfo =
{
  "isource",
  "General DC and sinusoidal current source",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"source",
  "2000_07_20"
};

// Parameter information
ParmInfo Isource::pinfo[] =
{
  {"idc", "DC current (V)", TR_DOUBLE, false},
  {"iac", "AC current peak amplitude (V)", TR_DOUBLE, false},
  {"f", "Source frequency (Hz)", TR_DOUBLE, false},
  {"phase", "Source phase (degrees)", TR_DOUBLE, false},
  {"delay", "Time before which the output voltage is set to DC value (s)",
	TR_DOUBLE, false},
  {"tr", "Rising time for DC component (s)", TR_DOUBLE, false},
  {"periods", "num of periods for current to ramp linearly", TR_DOUBLE, false}
};

Isource::Isource(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(idc = zero);
  paramvalue[1] = &(iac = zero);
  paramvalue[2] = &(frequency = zero);
  paramvalue[3] = &(phase = zero);
  paramvalue[4] = &(delay = zero);
  paramvalue[5] = &(tr = zero);
  paramvalue[6] = &(periods = 0.0);

  // Set the number of terminals
  setNumTerms(2);

  // Set flags
  setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN | SOURCE);
  // Set number of states
  setNumberOfStates(1);
}

void Isource::init() throw(string&)
{
  omega = twopi * frequency;
  phase_rad = deg2rad * phase;
  ramp_slope = frequency / periods;
}

void Isource::fillMNAM(FreqMNAM* mnam)
{
  const double& freq(mnam->getFreq());
  // Set source vector
  if (freq == zero)
	{
    mnam->addToSource(getTerminal(0)->getRC(),
		getTerminal(1)->getRC(),
		idc);
  }
  else if ((freq && abs(freq - frequency)/freq < 1e-14 ) || !frequency)
	{
    if (phase == zero)
      mnam->addToSource(getTerminal(0)->getRC(),
		getTerminal(1)->getRC(),
		iac);
    else
		{
      double_complex is(iac * cos(phase_rad),
			iac * sin(phase_rad));
      mnam->addToSource(getTerminal(0)->getRC(),
			getTerminal(1)->getRC(),
			is);
    }
  }

  // Since the source vector is assumed to be initialized to zero,
  // we do not need to fill it if freq != frequency.
  return;
}

void Isource::fillMNAM(TimeMNAM* mnam)
{
  // Nothing to be done.
  return;
}

void Isource::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime();
  double is = (tr && ctime < tr) ? idc * ctime / tr : idc;
  if ((ctime > delay) && (ctime < periods/frequency))
    is += iac * cos(omega * (ctime - delay) + phase_rad) * (ramp_slope * ctime);
  else if (ctime > delay)
    is += iac * cos(omega * (ctime - delay) + phase_rad);
  mnam->addToSource(getTerminal(0)->getRC(),
	getTerminal(1)->getRC(),
	is);
  return;
}

void Isource::svTran(TimeDomainSV* tdsv)
{
  // Calculate current
  double& is = tdsv->i(0);

  is = (tr) ? zero : - idc;
  if (!tdsv->DC())
	{
    const double& ctime = tdsv->getCurrentTime();
    is = (tr && ctime < tr) ? -idc * ctime / tr : -idc;
    if ((ctime > delay) && (ctime < periods/frequency))
      is += iac * cos(omega * (ctime - delay) + phase_rad) * (ramp_slope * ctime);
    else if (ctime > delay)
      is += iac * cos(omega * (ctime - delay) + phase_rad);
  }
  // State variable is voltage
  tdsv->u(0) = tdsv->getX(0) * 1e-2;
}

void Isource::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = one;
  tdsv->getJi()(0,0) = zero;
}

