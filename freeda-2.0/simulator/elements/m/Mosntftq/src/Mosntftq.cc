// Based on the Lueder model
#include "Mosntftq.h"

// Static members
const unsigned Mosntftq::n_par = 16;
const double E0 = 8.854194e-12;
const double q = 1.6e-19;
int counter=0;
int flag=0;

// Element information
ItemInfo Mosntftq::einfo =
{
  "mosntftq",
  "a-Si TFT MODEl",
  "ECe718 Class Project",
  DEFAULT_ADDRESS"transistor>mosfet",
  "2003_05_15"
};

// Parameter information
ParmInfo Mosntftq::pinfo[] =
{
  {"l", "Channel length (m)",TR_DOUBLE, false},
  {"w", "Channel width (m)", TR_DOUBLE, false},
  {"u", "Surface mobility (cm^2/V-s)", TR_DOUBLE, false},
  {"vto", "Zero-bias threshold voltage (V)", TR_DOUBLE, false},
  {"tts", "Device temperature (K)",TR_DOUBLE,false},
  {"ctr", "Mobility modulation (V^-1)", TR_DOUBLE, false},
  {"rd", "Static feedback on threshold voltage (V^-1) ", TR_DOUBLE, false},
  {"vthm", "Saturated velocity (m/sec)", TR_DOUBLE, false},
  {"beta", "Conductance (ohm-1)", TR_DOUBLE, false},
  {"roff", "Coefficient of drain leakage", TR_DOUBLE, false},
  {"cgso", "Gate-source capacitance (F)", TR_DOUBLE, false},
  {"cgdo", "Gate-drain capacitance (F)", TR_DOUBLE, false},
  {"cl", "Maximum space-charge capacitance (F/m^2)", TR_DOUBLE, false},
  {"lambda", "signal frequency (Hz)",TR_DOUBLE, false},
  {"rs", "Frequenct effect compensation", TR_DOUBLE, false},
  {"tnom", "Frequenct effect compensation", TR_DOUBLE, false}
};

Mosntftq::Mosntftq(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(l = 5.5e-6);
  paramvalue[1] = &(w = 80e-6);
  paramvalue[2] = &(u = 0.12e-4);
  paramvalue[3] = &(vto = 3.2);
  paramvalue[4] = &(tts = 232);
  paramvalue[5] = &(ctr = 0.06);
  paramvalue[6] = &(rd = 83e3);
  paramvalue[7] = &(vthm = 4.2 );
  paramvalue[8] = &(beta = 3);
  paramvalue[9] = &(roff = 4e15);
  paramvalue[10] = &(cgso = 80e-15);
  paramvalue[11] = &(cgdo = 80e-15);
  paramvalue[12] = &(cl = 400e-6);
  paramvalue[13] = &(lambda = 450);
  paramvalue[14] = &(rs= 83e3);
  paramvalue[15] = &(tnom= 300.15);

  // Set the number of terminals
  setNumTerms(4);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
	setNumberOfStates(3);
}

void Mosntftq::init() throw(string&)
{
  // constants
  A = exp(-(((log(1/ctr))/(vthm-vto))*vto));
  D = ((log(1/ctr))/(vthm-vto));

  // vt
	vt = (kBoltzman * tnom) / q;

  DenseIntVector var(3);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  initializeAD(var, var);
}

void Mosntftq::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: vds
  // x[1]: vgs
  // x[2]: vds - vgs
  // x[3]: dvds/dt
  // x[4]: dvgs/dt
  // x[5]: d(vgs - vds) / dt
  // vp[0]: ugs , ip[0]: id
  // vp[1]: uds , ip[1]: ig
  // vp[2]: usb , ip[2]:is
  AD MODE;

  AD  Ids;
  AD vgs,vds,vbs,qdrn,qgate,qsrc;
  AD kp,vtheff,kt,xs;

  vtheff = vthm-vto;
  kp = (w/l)*u*cl;
  kt = -4.0/15.0*pow((tnom/tts),2.0)+tnom/tts+1.0/15.0;
  xs = (pow((1-x[0]/lambda),-1))/2.0;

	//Subtreshold region
  if(x[1]<vto)
	{

		if(x[0]<vtheff)
			Ids = kp *(pow(vtheff,2.0/kt) - pow((vtheff-x[0]),2.0/kt)) * exp(beta*(x[1]-vto)) * xs * ctr;

		if(x[0]>vtheff)
			Ids = kp * (pow(vtheff, 2.0/kt)) * exp(beta*(x[1]-vto)) * xs * ctr;
	}

  //Transition region
  if((x[1]>=vto)&&(x[1]<vthm))
	{
		if(x[0]<vtheff)
			Ids = kp *(pow(vtheff,2.0/kt) - pow((vtheff-x[0]),2.0/kt)) * exp(D*x[1]) * xs * ctr * A;
		if(x[0]>vtheff)
			Ids = kp *(pow(vtheff,2.0/kt)) * exp(D*x[1]) * xs * ctr * A;
	}

  //above threshold region
  if(x[1]>=vthm)
	{
		if(x[0]<(x[1]-vto))
			Ids = kp *(pow((x[1]-vto),2.0/kt) - pow((x[1]-vto-x[0]),2.0/kt)) * xs ;
		if(x[0]>vtheff)
			Ids = kp *(pow((x[1]-vto),2.0/kt)) * xs ;
	}

  AD ck, a;
  a= (x[1]-vto)/(x[0]-vto);
  ck = w*l*cl;

  qdrn =0;
  qgate=0;
  qsrc=0;
  AD temp=0;
  temp= 1.0+2.0/kt;
  if((a>=0) && (a!=1))
	{
		qgate = ck*(x[1]-vto)*(1.0-pow(a,temp));
		qgate = qgate/((1.0+kt/2.0)*(1-pow(a,2.0/kt)));
		temp = 1.0+4.0/kt;
		qdrn = -qgate/(1-pow(a,2.0/kt)) + ck*(x[1]-vto)*(1.0-pow(a,temp))/((2.0+kt/2.0)*(1-pow(a,2.0/kt))*(1-pow(a,2.0/kt)));
		qsrc = qgate*pow(a,2.0/kt)/(1-pow(a,2.0/kt)) -ck*(x[1]-vto)*(1.0-pow(a,temp))/((2.0+kt/2.0)*(1-pow(a,2.0/kt))*(1-pow(a,2.0/kt)));
	}

  // Calculate dynamic current contributions due to charge
  // iqd is dynamic contribution to drain current
  // iqg is dynamic contribution to gate current
  // iqs is dynamic contribution to source current
  AD iqd, iqg, iqs;
  iqd = qdrn.fastAccessDx(0)*x[3] + qdrn.fastAccessDx(1)*x[4] + qdrn.fastAccessDx(2)*x[5];
  iqg = qgate.fastAccessDx(0)*x[3] + qgate.fastAccessDx(1)*x[4] + qgate.fastAccessDx(2)*x[5];
  iqs = qsrc.fastAccessDx(0)*x[3] + qsrc.fastAccessDx(1)*x[4] + qsrc.fastAccessDx(2)*x[5];

  flow[0] =  Ids + iqd; //DC Drain current
  flow[1] = iqg; //DC Gate current
  flow[2] =  -Ids - iqs; //DC Source current

  effort[0] = x[0];  //Vdb
  effort[1] = x[1]; //Vgb
  effort[2] = x[2];//Vsb
}
