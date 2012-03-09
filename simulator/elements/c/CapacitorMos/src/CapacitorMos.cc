#include "CapacitorMos.h"

// Static members
const unsigned CapacitorMos::n_par = 18;

// Element information
ItemInfo CapacitorMos::einfo =
{
  "capacitormos",
  "MOS Capacitor with tunneling",
  "Krishnanshu Dandu and Yawei Jin",
  DEFAULT_ADDRESS"category:lumped,semiconductor,mosfet",
  "2003_05_15"
};

// Parameter information
ParmInfo CapacitorMos::pinfo[] =
{
  {"w", "Width of the device (in cm)", TR_DOUBLE, true},
  {"l", "Length of the gate (in cm)", TR_DOUBLE, true},
  {"tox", "Oxide thickness (cm)", TR_DOUBLE, false},
  {"t", "Temperature (Kelvin)", TR_DOUBLE, false},
  {"npoly", "Poly Doping /cm3", TR_DOUBLE, false},
  {"nsub", "Substrate Doping /cm3", TR_DOUBLE, false},
  {"moxecb", "Effective mass of electron for ECB (kg)", TR_DOUBLE, false},
  {"moxevb", "Effective mass of electron for EVB(kg)", TR_DOUBLE, false},
  {"moxhvb", "Effective mass of electron for HVB (kg)", TR_DOUBLE, false},
  {"phibecb", "Electron Tunneling Barrier Height (eV)", TR_DOUBLE, false},
  {"phibevb", "Electron Tunneling Barrier Height (eV)", TR_DOUBLE, false},
  {"phibhvb", "Electron Tunneling Barrier Height (eV)", TR_DOUBLE, false},
  {"phib0ecb", "Si/SiO2 Barrier Height for ECB (eV)", TR_DOUBLE, false},
  {"phib0evb", "Si/SiO2 Barrier Height for EVB (eV)", TR_DOUBLE, false},
  {"phib0hvb", "Si/SiO2 Barrier Height for HVB(eV)", TR_DOUBLE, false},
  {"vfb", "Flat-Band Voltage (V)", TR_DOUBLE, false},
  {"epsrox", "Relative Permittivity of Oxide", TR_DOUBLE, false},
  {"s","subthreshold swing (mv/dec)", TR_DOUBLE , false},
};


CapacitorMos::CapacitorMos(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(w = 1e-4);
  paramvalue[1] = &(l = 1e-4);
  paramvalue[2] = &(tox = 2.48e-7);
  paramvalue[3] = &(t = 300);
  paramvalue[4] = &(npoly = 9e19);
  paramvalue[5] = &(nsub = 4.7e17);
  paramvalue[6] = &(moxecb = 3.64e-31);
  paramvalue[7] = &(moxevb = 2.73e-31);
  paramvalue[8] = &(moxhvb = 2.91e-31);
  paramvalue[9] = &(phibecb = 3.1);
  paramvalue[10] = &(phibevb = 4.2);
  paramvalue[11] = &(phibhvb = 4.5);
  paramvalue[12] = &(phib0ecb = 3.1);
  paramvalue[13] = &(phib0evb = 3.1);
  paramvalue[14] = &(phib0hvb = 4.5);
  paramvalue[15] = &(vfb = -0.9);
  paramvalue[16] = &(epsrox = 3.9);
  paramvalue[17] = &(s=75e-3);

	// Set the number of terminals
	setNumTerms(3); //extra dummy terminal to handle error function for Vox
	// Set number of states
	setNumberOfStates(2);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void CapacitorMos::init() throw(string&)
{
	DenseIntVector var(2);
	var[1] = 1;
	initializeAD(var);
}

void CapacitorMos::eval(AD * x, AD * effort, AD * flow)
{
  //x[0] vg; x[1] vox

  double k = kBoltzman;
  double q = eCharge;
  double h = 6.634e-34;
  double vt = k*t/q;
  double e0 = epsilon0*1e-2;
  double esi = 11.7 * e0;
  double eox = epsrox * e0;
  double pi = 3.14;

  double Eg0 = 1.16 - (7.02e-4 * t * t / (t + 1108.0));

  double ni = 1.45e10 * (t/300.15) * sqrt(t/300.15) * exp(21.5565981 - Eg0/(2.0 * vt));

  double phi_f = vt*log(nsub/ni);
  double phi_s0 = 2* phi_f;
  double phig = Eg0;

  double ninv = 1;
  double nacc = s/vt;
  double nevb = 3;

  double gamma = tox*sqrt(2*esi*q*nsub)/eox;
  double vth = vfb + phi_s0 + gamma*sqrt(phi_s0);

  double A1ecb = (q*q)/(8*pi*h*phibecb*eox);
  double A1evb = (q*q)/(8*pi*h*phibevb*eox);
  double A1hvb = (q*q)/(8*pi*h*phibhvb*eox);

  double A2ecb = (-8*pi*sqrt(2*moxecb)*pow((phibecb*q),1.5)*tox/(100*3*h*q));
  double A2evb = (-8*pi*sqrt(2*moxevb)*pow((phibevb*q),1.5)*tox/(100*3*h*q));
  double A2hvb = (-8*pi*sqrt(2*moxhvb)*pow((phibhvb*q),1.5)*tox/(100*3*h*q));

  double alfaecb = -0.58*A2ecb/(phibecb*phibecb);
  double alfaevb = -0.58*A2evb/(phibevb*phibevb);
  double alfahvb = -0.58*A2hvb/(phibhvb*phibhvb);

  double K0ecb = exp(A2ecb*1.6/phibecb);
  double K0evb = exp(A2evb*1.6/phibevb);
  double K0hvb = exp(A2hvb*1.6/phibhvb);

  double v1 = log(1/(K0ecb*alfaecb))/alfaecb;
  double k3 = exp(alfaecb * v1);

  AD y,idecb, idevb, idhvb;

  //vg = x[0]; x as in rizzolli's equation is x[1]
  // Y is the oxide voltage calculated using parameterization
  // idecb calculates the exponential term in the current equation
  if (v1 > x[1])
  {
    y = x[1] + zero;
    idecb = K0ecb * (exp(alfaecb * x[1]));
  }
  else
  {
    y = v1 + log(one + alfaecb*(x[1]-v1))/alfaecb;
    idecb = K0ecb * k3 * (one + alfaecb * (x[1] - v1));
  }

  //calculating the expnonential terms for evb and hvb from ecb
  idevb = K0evb * pow((idecb/K0ecb),(alfaevb/alfaecb));
  idhvb = K0hvb * pow((idecb/K0ecb),(alfahvb/alfaecb));

  AD v_ge = vfb + phi_s0 + (q*esi*npoly*tox*tox)*(sqrt((1+2*eox*eox*(x[0]-vfb-phi_s0)/(q*esi*npoly*tox*tox)))-1)/(eox*eox);

  AD phi_s = (gamma*gamma/4)*(-1 + sqrt(1 + (4*(x[0] - v_ge - vfb)/(gamma*gamma))))*(-1 + sqrt(1 + (4*(x[0] - v_ge - vfb)/(gamma*gamma))));

  // this is the Vox calculated from Vg
  AD v_ox = v_ge - phi_s - vfb;

  // calculating the absolute value of Vox - required for C and to generate the error term
  AD absvox;
  if (v_ox > 0.0)
    absvox = v_ox;
  else
    absvox = -v_ox;

  AD Necbhvb = eox*(ninv*vt*log(1+exp((v_ge-vth)/(ninv*vt))) + nacc*vt*log(1+exp(-(x[0]-vfb)/(nacc*vt))))/tox;
  AD Nevb = eox*(nevb*vt*log(1+exp((absvox-phig)/(nevb*vt))))/tox;
  // evb current can exist only if |Vox| > phig
  if (absvox <= phig)
    Nevb = 0.0;

  AD Cecb = exp(20*pow(((absvox-phibecb)/phib0ecb+1),0.6)*(1-absvox/phibecb)/phibecb)*x[0]*Necbhvb/tox;
  AD Cevb = exp(20*pow(((absvox-phibevb)/phib0evb+1),1)*(1-absvox/phibevb)/phibevb)*x[0]*Nevb/tox;
  AD Chvb = exp(20*pow(((absvox-phibhvb)/phib0hvb+1),0.4)*(1-absvox/phibhvb)/phibhvb)*x[0]*Necbhvb/tox;

  AD Jgecb = A1ecb*idecb*Cecb;
  AD Jgevb = A1evb*idevb*Cevb;
  AD Jghvb = A1hvb*idhvb*Chvb;

  AD Jtotal = Jgecb + Jgevb + Jghvb;
  AD Itunnel = w*l*Jtotal;

  effort[0] = x[0]; // gate voltage
  flow[0] = Itunnel; // tunneling current
  effort[1] = 0; // dummy terminal voltage
  flow[1] = y - absvox;  //error expression
}

