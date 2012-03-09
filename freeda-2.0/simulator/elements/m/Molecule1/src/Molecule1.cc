#include "Molecule1.h"

// Static members
const unsigned Molecule1::n_par = 10;      //number of Model-Parameters

// Element information
ItemInfo Molecule1::einfo =
{
  "molecule1",
  "Two terminal molecular device simulating conduction through a single molecular energy level",
  "Nikhil Kriplani",
  DEFAULT_ADDRESS"category:molecular electronic",
  "2004_03_01"
};

// Parameter information
ParmInfo Molecule1::pinfo[] =
{
  {"e0", "Equilibrium energy level (eV)", TR_DOUBLE, false},
  {"g1", "Energy level broadening due to contact one (eV)", TR_DOUBLE, false},
  {"g2", "Energy level broadening due to contact two (eV)", TR_DOUBLE, false},
  {"ef", "Fermi reference level (eV)", TR_DOUBLE, false},
  {"etha", "Voltage division factor (unity)", TR_DOUBLE, false},
  {"m", "Proportionality constant in energy-time-uncertainty (unity)", TR_DOUBLE, false},
  {"u0", "Single electron charging energy (eV)", TR_DOUBLE, false},
  {"t", "Operation Temperature (K)", TR_DOUBLE, false},
  {"accur", "Accuracy of the predictor-corrector arrangement", TR_DOUBLE, false},
  {"area", "Area of device", TR_DOUBLE, false}
};


Molecule1::Molecule1(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(e0 = -4.5);
  paramvalue[1] = &(g1 = 0.1);
  paramvalue[2] = &(g2 = 0.1);
  paramvalue[3] = &(ef = -5.1);
  paramvalue[4] = &(etha = 0.5);
  paramvalue[5] = &(m = 1.0);
  paramvalue[6] = &(u0 = 1.3);
  paramvalue[7] = &(t = 300.0);
  paramvalue[8] = &(accur = 1e-6);
  paramvalue[9] = &(area = 1.0);

  setNumTerms(2);
  setNumberOfStates(1);
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
};

void Molecule1::init() throw(string&)
{
  DenseIntVector var(1);
  initializeAD(var, var);
};

void Molecule1::eval(AD * x, AD * effort, AD * flow)
{
  double h = 6.634e-34;

  // Number of electrons on level to make molecule neutral
  double n0 = 2.0 / (1.0 + exp(eCharge * (e0 - ef)/(kBoltzman * t)));

  AD n, iout, U_scf, dU, U_tmp, E_scf, f1, f2, mu1, mu2;

  dU = 1.0;

  mu1 = ef - etha * x[0];
  mu2 = ef + (1.0 - etha) * x[0];

  // the predictor corrector loop
  while (dU > accur)
  {
    E_scf = e0 + U_scf;
    f1 = 1.0 / (1.0 + exp(eCharge * (E_scf - mu1) / (kBoltzman * t)));
    f2 = 1.0 / (1.0 + exp(eCharge * (E_scf - mu2) / (kBoltzman * t)));
    n = 2.0 * (g1 * f1 + g2 * f2) / (g1 + g2);
    U_tmp = U_scf;
    U_scf = U_tmp + 0.05 * (u0*(n - n0) - U_tmp);
    dU = fabs(U_scf - U_tmp);
  }

  iout = (1/m) * (4.0*pi * eCharge*eCharge / h) * g1*g2 / (g1 + g2) * (f2 - f1);

  //assigning output values
  effort[0] = x[0];      // voltage across molecular device
  flow[0] = area * iout; // current through molecular device
}

