#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "Vlfmpulse.h"
#include "../../../../analysis/TimeDomainSV.h"

// Static members
const unsigned Vlfmpulse::n_par = 9;

// Element information
ItemInfo Vlfmpulse::einfo = 
{
  "vlfmpulse",
  "Linear FM Chirp pulse voltage source",
  "Frank P. Hart"
    DEFAULT_ADDRESS"elements/Vlfmpulse.html",
  "2004_05_18"
};

// Parameter information
ParmInfo Vlfmpulse::pinfo[] = 
{
  {"vo", "Offset value (V)", TR_DOUBLE, false},
  {"va", "Pulse amplitude (V)", TR_DOUBLE, false},
  {"td", "Delay time (s)", TR_DOUBLE, false},
  {"fo", "Center frequency (Hz)", TR_DOUBLE, false}, 
  {"deltaf", "Sweep frequency range (Hz)", TR_DOUBLE, false}, 
  {"chirpdir", "Chirp direction ()", TR_INT, false}, 
  {"phi", "Fixed phase (degrees)",TR_DOUBLE,false},
  {"tau", "Pulse width (s)", TR_DOUBLE, false},
  {"per", "Period (s)", TR_DOUBLE, false}
};


Vlfmpulse::Vlfmpulse(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(vo = zero);
  paramvalue[1] = &(va = zero);
  paramvalue[2] = &(td = zero);
  paramvalue[3] = &(fo = 1.0);          // Set to unity by default
  paramvalue[4] = &(deltaf = 0.5);      // Half of fo
  paramvalue[5] = &(chirpdir = 1);	// Upchirp by default
  paramvalue[6] = &(phi = zero);
  paramvalue[7] = &(tau = zero);
  paramvalue[8] = &(per = zero);

  // Set the number of terminals
  setNumTerms(2);

  // Set flags
  setFlags(LINEAR | ONE_REF | TR_TIME_DOMAIN | SOURCE);
  // Set number of states
  setNumberOfStates(1);

  my_row = 0;
  int_per_time = zero;
}

unsigned Vlfmpulse::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Keep the equation number assigned to this element
  my_row = eqn_number;
  // Add one extra RC
  return 1;
}

void Vlfmpulse::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row);
  first_eqn = my_row;
  n_rows = 1;
}

void Vlfmpulse::fillMNAM(FreqMNAM* mnam)
{
  assert(my_row);
  // Ask my terminals the row numbers
  mnam->setOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);

  // This element behaves as a short circuit for AC/HB analysis.
  // Since the source vector is assumed to be initialized to zero,
  // we do not need to fill it if freq != frequency.
  return;
}

void Vlfmpulse::fillMNAM(TimeMNAM* mnam)
{
  assert(my_row);
  // Ask my terminals the row numbers
  mnam->setMOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);
  return;
}

void Vlfmpulse::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime();

  // set reltime to modulo of per, subtracting init delay td
  // so reltime is zero at the beginning of pulse, tau at end
  double reltime = ctime - td - int_per_time;
  if (reltime > per) 
  {
    int_per_time = int((ctime - td)/per) * per;
    reltime = ctime - td - int_per_time;
  }

  double e = zero;
  double fchirp = zero;
  double freq = zero;
  double phir = pi*phi/180.0;
  if (int_per_time == zero && ctime < td)
    e = vo;
  else 
  {
    if (chirpdir == 1) // upchirp
      fchirp = -deltaf/2.0 + (deltaf*reltime/tau)/2.0;
    else 	     // downchirp - chirpdir = 0
      fchirp = deltaf/2.0 - (deltaf*reltime/tau)/2.0;
    freq = fo + fchirp;

    if (reltime <= tau)
      e = vo + va*cos(2*pi*freq*reltime+phir);
    else
      e = vo;
  }
  mnam->setSource(my_row, e);
  return;
}

void Vlfmpulse::svTran(TimeDomainSV* tdsv)
{
  // Calculate voltage
  double& e = tdsv->u(0);

  if (tdsv->DC())
    e = vo;
  else 
  {
    const double& ctime = tdsv->getCurrentTime();

    // set reltime to modulo of per, subtracting init delay td
    // so reltime is zero at the beginning of pulse, tau at end
    double reltime = ctime - td - int_per_time;
    if (reltime > per) 
    {
      int_per_time = int((ctime-td)/per) * per;
      reltime = ctime - td - int_per_time;
    }

    double e = zero;
    double fchirp = zero;
    double freq = zero;
    double phir = pi*phi/180.0;
    if (int_per_time == zero && ctime < td)
      e = vo;
    else 
    {
      if (chirpdir == 1) // upchirp
        fchirp = -deltaf/2.0 + (deltaf*reltime/tau)/2.0;
      else 	     // downchirp - chirpdir = 0
        fchirp = deltaf/2.0 - (deltaf*reltime/tau)/2.0;
      freq = fo + fchirp;

      if (reltime <= tau)
        e = vo + va*cos(2*pi*freq*reltime+phir);
      else
        e = vo;
    }
  }
  // Scale state variable for numerical stability
  tdsv->i(0) = tdsv->getX(0) * 1e-2;
}

void Vlfmpulse::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = zero;
  tdsv->getJi()(0,0) = 1e-2;
}

