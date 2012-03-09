#include "Molecule2.h"

// Static members
const unsigned Molecule2::n_par = 11;      //number of Model-Parameters

// Element information
ItemInfo Molecule2::einfo =
{
  "molecule2",
  "Two terminal molecular device simulating conduction through a single molecular energy level",
  "Nikhil Kriplani",
  DEFAULT_ADDRESS"category:molecular electronic",
  "2004_03_09"
};

// Parameter information
ParmInfo Molecule2::pinfo[] =
{
  {"e0", "Equilibrium energy level (eV)", TR_DOUBLE, false},
  {"g1", "Energy level broadening due to contact one (eV)", TR_DOUBLE, false},
  {"g2", "Energy level broadening due to contact two (eV)", TR_DOUBLE, false},
  {"ef", "Fermi reference level (eV)", TR_DOUBLE, false},
  {"etha", "Voltage division factor (unity)", TR_DOUBLE, false},
  {"u0", "Single electron charging energy (eV)", TR_DOUBLE, false},
  {"t", "Operation Temperature (K)", TR_DOUBLE, false},
  {"accur", "Accuracy of the predictor-corrector arrangement", TR_DOUBLE, false},
  {"npoints", "Number of integration points", TR_DOUBLE, false},
  {"width", "Half the width of the integration interval", TR_DOUBLE, false},
  {"area", "Device area scaling parameter", TR_DOUBLE, false}
};


Molecule2::Molecule2(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(e0 = -4.5);
  paramvalue[1] = &(g1 = 0.2);
  paramvalue[2] = &(g2 = 0.2);
  paramvalue[3] = &(ef = -5.1);
  paramvalue[4] = &(etha = 0.5);
  paramvalue[5] = &(u0 = 1.3);
  paramvalue[6] = &(t = 300.0);
  paramvalue[7] = &(accur = 1e-4);
  paramvalue[8] = &(npoints = 400.0);
  paramvalue[9] = &(width = 3.0);
  paramvalue[10] = &(area = 1.0);

  setNumTerms(2);
  setNumberOfStates(1);
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
};

void Molecule2::init() throw(string&)
{
  DenseIntVector var(1);
  initializeAD(var, var);
};

void Molecule2::eval(AD * x, AD * effort, AD * flow)
{
  double h = 6.634e-34;
  double g = g1 + g2;
  double n0 = 2.0 / (1.0 + exp(eCharge * (e0 - ef)/(kBoltzman * t)));
  double E1, E2, dE;

  AD n, iout, U_scf, dU, U_tmp, E_scf, mu1, mu2;
  AD f1, f2, D;

  dU = 1.0; U_scf = 0.0;
  mu1 = ef - etha * x[0];
  mu2 = ef + (1.0 - etha) * x[0];

	//   // NEW TECHNIQUE (is very slow)
	//   E1 = e0 - width/2;
	//   E2 = e0 + width/2;
	//   dE = fabs(E1 - E2) / npoints;
	//
	//   // fill up the energy interval
	//   double *E;
	//   E = new double[200];
	//   for (int i = 0; i < 200; i++)
	//   {
		//     E[i] = E1;
		//     E1 += dE;
	//   }
	//
	//   AD *D, *f1, *f2;
	//   D = new AD[200];
	//   // fill up f1 and f2
	//   f1 = new AD[200];
	//   f2 = new AD[200];
	//
	//   for (int i = 0; i < 200; i++)
	//   {
		//     f1[i] = 1.0 / (1.0 + exp(eCharge * (E[i] - mu1) / (kBoltzman * t)));
		//     f2[i] = 1.0 / (1.0 + exp(eCharge * (E[i] - mu2) / (kBoltzman * t)));
	//   }
	//
	//   // the predictor corrector loop
	//   while (dU > accur)
	//   {
		//     E_scf = e0 + U_scf;
		//     // fill the lorentzian
		//     for (int i = 0; i < 200; i++)
		//     {
			//       D[i] = 1.0/(2.0*pi) * g / ((E1- E_scf)*(E1 - E_scf) + g*g/4.0);
			//       n += 2.0 * D[i] * (g1*f1[i] + g2*f2[i]) / (g1 + g2) * dE;
		//     }
		//     cout << "left for loop \n";
		//     U_tmp = U_scf;
		//     U_scf = U_tmp + 0.2 * (u0*(n - n0) - U_tmp);
		//     dU = fabs(U_scf - U_tmp);
		//     cout << "dU = " << dU << "\n";
	//   }
	//
	//   // now the current
	//   for (int i = 0; i < 200; i++)
	//     iout += D[i] * (4.0*pi * eCharge * eCharge / h) * g1*g2 / (g1 + g2) * (f2[i] - f1[i]) * dE;
	//
	//   // delete used heap memory in reverse order
	//   delete []f2;
	//   delete []f1;
	//   delete []D;


  //OLD TECHNIQUE
  // the predictor corrector loop
  while (dU > accur)
  {
    E_scf = e0 + U_scf;
    n = 0.0;
    iout = 0.0;
    // width = width of the Lorentzian pulse
    E1 = e0 - width/2;
    E2 = e0 + width/2;
    dE = fabs(E1 - E2) / npoints;
    while(E1 < E2)
    {
      f1 = 1.0 / (1.0 + exp(eCharge * (E1 - mu1) / (kBoltzman * t)));
      f2 = 1.0 / (1.0 + exp(eCharge * (E1 - mu2) / (kBoltzman * t)));
      // the lorentzian
      D = 1.0/(2.0*pi) * g / ((E1- E_scf)*(E1 - E_scf) + g*g/4.0);
      n += 2.0 * D * (g1*f1 + g2*f2) / (g1 + g2) * dE;
      iout += D * (4.0*pi * eCharge * eCharge / h) * g1*g2 / (g1 + g2) * (f2 - f1) * dE;
      E1 += dE;
    }
    U_tmp = U_scf;
    U_scf = U_tmp + 0.2 * (u0*(n - n0) - U_tmp);
    dU = fabs(U_scf - U_tmp);
  }

  //assigning output values
  effort[0] = x[0];
  flow[0] = area * iout;
}
