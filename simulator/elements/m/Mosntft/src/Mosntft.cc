// Based on the HP a-Si TFT model
#include "Mosntft.h"

// Static members
const unsigned Mosntft::n_par = 36;
const double E0 = 8.854194e-12;
const double q = 1.6e-19;

// Element information
ItemInfo Mosntft::einfo =
{
  "mosntft",
  "HP a-Si Mosfet model level 16",
  "ECE718 Class Project",
  DEFAULT_ADDRESS"transistor>mosfet",
  "2003_05_15"
};

// Parameter information
ParmInfo Mosntft::pinfo[] =
{
  {"L", "Channel length (m)",TR_DOUBLE, false},
  {"LD", "Lateral diffusion (m)", TR_DOUBLE, false},
  {"W", "Channel width (m)", TR_DOUBLE, false},
  {"WD", "Lateral diffusion width (m)", TR_DOUBLE, false},
  {"U0", "Surface mobility (cm^2/V-s)", TR_DOUBLE, false},
  {"VT0", "Zero-bias threshold voltage (V)", TR_DOUBLE, false},
  {"PHI", "Surface inversion potential (V)", TR_DOUBLE, false},
  {"NFS", "Fast surface state density (cm^-2)", TR_DOUBLE, false},
  {"NSS", "Surface state density (cm^-2)",TR_DOUBLE,false},
  {"T1", "Film Thickness (m)", TR_DOUBLE, false},
  {"T2", "Film Thickness (m)", TR_DOUBLE, false},
  {"E1", "Dielectric constant of film1", TR_DOUBLE, false},
  {"E2", "Dieleectric constant of film2", TR_DOUBLE, false},
  {"THETA", "Mobility modulation (V^-1)", TR_DOUBLE, false},
  {"ETA", "Static feedback on threshold voltage (V^-1) ", TR_DOUBLE, false},
  {"VMAX", "Saturated velocity (m/sec)", TR_DOUBLE, false},
  {"G0", "Conductance (ohm-1)", TR_DOUBLE, false},
  {"DEFF", "Coefficient of drain leakage", TR_DOUBLE, false},
  {"VX", "potential (V)", TR_DOUBLE, false},
  {"TVST", "Voltage supply time (sec)", TR_DOUBLE, false},
  {"PSI", "temperature exponential part", TR_DOUBLE, false},
  {"GAMMA", "first order temperature gradient", TR_DOUBLE, false},
  {"VTIME", "voltage stress (sec)", TR_DOUBLE, false},
  {"TREF", "Nominal temperature (K)", TR_DOUBLE, false},
  {"T", "Device temperature (K)",TR_DOUBLE,false},
  {"RD", "Drain resistance (ohm)", TR_DOUBLE, false},
  {"RS", "Source resistance (ohm)", TR_DOUBLE, false},
  {"CGS0", "Gate-source capacitance (F)", TR_DOUBLE, false},
  {"CGD0", "Gate-drain capacitance (F)", TR_DOUBLE, false},
  {"CSC", "Maximum space-charge capacitance (F/m^2)", TR_DOUBLE, false},
  {"FREQ", "signal frequency (Hz)",TR_DOUBLE, false},
  {"FEFF", "Frequenct effect compensation", TR_DOUBLE, false},
  {"TAU", "Relaxation time constant (sec)", TR_DOUBLE, false},
  {"NU" , "First Order temperature gradient",TR_DOUBLE, false},
  {"CHI" , "Temperature exponential part",TR_DOUBLE, false},
  {"K2", "Temperature exponential part ",TR_DOUBLE, false}
};

Mosntft::Mosntft(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(L = 11e-6);
  paramvalue[1] = &(LD = 0.0);
  paramvalue[2] = &(W = 41e-6);
  paramvalue[3] = &(WD = 0.0);
  paramvalue[4] = &(U0 = 0.450);
  paramvalue[5] = &(VT0 = 1.699);
  paramvalue[6] = &(PHI = 0.620);
  paramvalue[7] = &(NFS = 1.925E21 );
  paramvalue[8] = &(NSS = 0.0 );
  paramvalue[9] = &(T1 = 300e-9);
  paramvalue[10] = &(T2 = 0.0);
  paramvalue[11] = &(E1 = 3.9);
  paramvalue[12] = &(E2 = 0.0);
  paramvalue[13] = &(THETA = 17.8e-3);
  paramvalue[14] = &(ETA = 410.8e-6);
  paramvalue[15] = &(VMAX = 2783);
  paramvalue[16] = &(G0 = 9.728e-15);
  paramvalue[17] = &(DEFF = 1.968);
  paramvalue[18] = &(VX = 0.033);
  paramvalue[19] = &(TVST = 100e-3);
  paramvalue[20] = &(PSI = 0.2);
  paramvalue[21] = &(GAMMA = 0.008);
  paramvalue[22] = &(VTIME = 10E-3);
  paramvalue[23] = &(TREF = 1.5);
  paramvalue[24] = &(T = 300.15);
  paramvalue[25] = &(RD = 8030);
  paramvalue[26] = &(RS = 8030);
  paramvalue[27] = &(CGD0 = 42.03e-15);
  paramvalue[28] = &(CGS0 = 52.21e-15);
  paramvalue[29] = &(CSC = 158.8e-6);
  paramvalue[30] = &(FREQ = 100e3);
  paramvalue[31] = &(FEFF = 0.5);
  paramvalue[32] = &(TAU = 10.7e-9);
  paramvalue[33] = &(NU = 0.0);
  paramvalue[34] = &(CHI = 0.5);
  paramvalue[35] = &(K2 = 2.0);
  // Set the number of terminals
  setNumTerms(3);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(2);
}


void Mosntft::init() throw(string&)
{
	// constants
  leff = L - 2.0 * LD;
  weff = W - 2.0 * WD;

	// Total insulator film thickness
  tfm = T1 + T2;

	// film capacitance
  cfm = (E0 * E1 * E2)/((T2 * E1) + (T1 * E2));
  if (E1 == 0.0)
		cfm = E0 * E2 / T2;
  if (E2 == 0.0)
		cfm = E0 * E1 / T1;
	// vt
	vt = (kBoltzman * T) / q;

  DenseIntVector var(2);
  var[0] = 0;
  var[1] = 1;
  //var[2] = 2;
  DenseIntVector dvar(2);
  dvar[0] = 0;
  dvar[1] = 1;
  //dvar[2] = 2;
  DenseIntVector dvar2;
  DenseIntVector tvar;
  DenseDoubleVector nodelay;
  initializeAD(var,dvar,dvar2,tvar,nodelay);
}


void Mosntft::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: vds
  // x[1]: vgs
  // x[2]: vds - vgs
  // x[3]: dvds/dt
  // x[4]: dvgs/dt
  // x[5]: d(vgs - vds) / dt
  // effort[0]: ugs , flow[0]: id
  // effort[1]: uds , flow[1]: ig

  AD MODE;

  //other parameters
  AD vto, eg, vbi,phi,uo, RATIO;
  AD vfb, vdsat, beta, von;
  AD mueff, Ids, Idssub,Igs=0,Icgs,Icgd;

  //capacitance
  AD Cgdi, Cgsi, ids, Cgs, Cgd;
  AD vgsx, vdsc, cfmmis, vdsx, epsfm,cfmlw,fval;
  AD cd1, cdnorm, fdrain;
  AD vgs,vds,vbs;
  AD vth,kp,xn, vtheff;

  vto = VT0;
  vbi = vto;
  phi = PHI;
  Cgdi = 0.0;
  Cgsi = 0.0;

  vth = vto + (ETA * x[0]);
  uo = U0 * 1e-4;
  mueff = (uo) / (1.0 + (THETA * (x[1] - vth)));
  beta = (weff / leff) * (cfm) ;

  vdsc = leff * VMAX / mueff;
  von = vth;
  if(x[1]<=0 || x[0] <=0)
    vtheff=vth;
  else
    vtheff = vth + GAMMA * (pow((x[1]),PSI) - (pow((x[0]),PSI)/2.0)) * exp((-q * VX/(kBoltzman * T)) * log(TVST));
  if (NFS != 0.0)
	{
		xn = 1.0 + ((q * NFS * 1e4 * weff * leff) / cfm);
    von = vth + (ETA*x[0]) + xn * vt;
	}

  if (x[1] > vth)
    vgsx = x[1];
  else
    vgsx = von;
  vdsat = (vgsx - vth) + vdsc - sqrt((((vgsx - vtheff) * (vgsx - vtheff)) + (vdsc * vdsc)));

  if (x[0] < vdsat)
    vdsx = x[0];
  else
    vdsx = vdsat;

  //general current equations
  Ids = (beta * mueff ) * (vdsx)*(vgsx - vth - (vdsx/2.0));
  Idssub = Ids * exp(((x[1] - von) / (xn * vt)) + 1.0);

  //off region
  if (x[1] < 0)
    Ids = (Idssub+ G0 * (x[1] + DEFF * x[0]));

  //subthreshold
  else if (x[1] < vth)
	{
		Ids = Idssub +G0 * (x[1] + DEFF * x[0]);
	}

  //saturation
  else if (x[1] <= (x[0] + vth))
	{
		//Ids = (beta * uo)*vdsc/(1+vdsx)*(((vgsx-vth)*vdsx)-((vdsx* vdsx)/2.0)) +G0*(x[0]+DEFF*x[1]);
		Ids= beta * mueff*vdsx*(vgsx - vth - (vdsx/2.0));
		Ids=Ids;
		//*vdsc/(1+vdsx)+G0*(x[1]+DEFF*x[0]);
		//      Ids = Ids +G0*(x[0]+DEFF*x[1]);
	}

  //capacitor calculation
  cfmmis = (cfm * CSC) / (cfm + CSC);
  cfmlw = cfmmis*L*W;
  epsfm = cfm *(T2+T1)/E0;
  if (x[1] >= vth)
	{
		if(x[0] < vdsat)
		{
			Cgsi = cfmlw * (1.0 - (pow(x[0],2)/(12.0 * pow((2.0 * x[1] - 2.0 * vtheff - x[0]),2))));
			Cgdi = cfmlw * (0.5 - ((2.0 * x[0] * (2.0 * x[1] - 2.0 * vtheff - x[0]) +
			pow(x[0],2))/(6.0 * pow((2.0 * x[1] - 2.0 * vtheff - x[0]),2))));
			Cgdi= Cgdi;
		}
		if(x[0] >= vdsat)
		{
			Cgsi = 0.5*cfmlw;
			Cgdi = cfmlw;
		}
	}
  else
	{
		Cgsi =0.0;
		Cgdi=0.0;
	}
  Cgd =Cgdi+CGD0;
  Cgs =Cgsi+CGS0;

  Icgd = Cgd * (x[3]-x[2]);
  Icgs = Cgs * x[3] ;

  effort[0] = x[0];
  effort[1] = x[1];
  flow[0] =Ids-Icgd;
  flow[1] = Icgs+Icgd;
}
