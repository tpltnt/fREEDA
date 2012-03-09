#include "Mospekv2.h"

const unsigned Mospekv2::n_par = 40;

// Element information
ItemInfo Mospekv2::einfo =
{
  "mospekv2",
  "EPFL-EKV MOST Model",
  "William Cornwell & Mohamed Guled & Jordan Novgrod & Mac Perry & Ryon Stewart",
  DEFAULT_ADDRESS"transistor>mosfet",
  "2003_05_10"
};

// Parameter information
ParmInfo Mospekv2::pinfo[] =
{
	// Device input variables
	{"l","channel length (m)", TR_DOUBLE, false},  // not sure about (M)
	{"w","channel width (m)", TR_DOUBLE, false},
	{"np","parallel multiple device number", TR_DOUBLE, false},
	{"ns","series multiple device number", TR_DOUBLE, false},
	// EKV intrinsic model process related parameters
	{"cox","gate oxide capacitance per unit area (F/m^2)", TR_DOUBLE, false},
	{"xj","junction depth (m)", TR_DOUBLE, false},
	{"dw","channel width correction (m)", TR_DOUBLE, false},
	{"dl","channel length correction (m)", TR_DOUBLE, false},
	// EKV basic intrinsic model parameters
	{"vto","long-channel threshold voltage (V)", TR_DOUBLE, false},
	{"gamma","body effect parameter (sqrtV)", TR_DOUBLE, false},
	{"phi","bulk Fermi potential (V)", TR_DOUBLE, false},
	{"kp","trasconductance parameter (A/V^2)", TR_DOUBLE, false},
	{"e0","mobility reduction coefficient (V/m)", TR_DOUBLE, false},
	{"ucrit","longitudinal critical field (V/m)", TR_DOUBLE, false},
	// Optional EKV intrinsic model parameters
	{"tox","oxide thickness (m)", TR_DOUBLE, false},
	{"nsub","channel doping (cm^-3)", TR_DOUBLE, false},
	{"vfb","flat-band voltage (V)", TR_DOUBLE, false},
	{"uo","low-field mobility (cm^2/Vs)", TR_DOUBLE, false},
	{"vmax","saturation velocity (m/s)", TR_DOUBLE, false},
	{"theta","mobility reduction coefficient (1/V)", TR_DOUBLE, false},
	// Channel length modulation and charge sharing parameters
	{"lambda","depletion length coefficient", TR_DOUBLE, false},
	{"weta","narrow-channel effect coefficient", TR_DOUBLE, false},
	{"leta","short-channel effect coefficient", TR_DOUBLE, false},
	// Reverse short-channel effect parameters
	{"q0","reverse short channel effect peak charge density (A*s/m^2)", TR_DOUBLE, false},
	{"lk","reverse short channel effect coefficient (m)", TR_DOUBLE, false},
	// Impact ionization related parameters
	{"iba","first impact ionization coefficient (1/m)", TR_DOUBLE, false},
	{"ibb","second impact ionization coefficient (V/m)", TR_DOUBLE, false},
	{"ibn","saturation voltage factor for impact ionization", TR_DOUBLE, false},
	// Intrinsic model temperature parameters
	{"tcv","threshold voltage temperature coefficient (V/K)", TR_DOUBLE, false},
	{"bex","mobility temperature exponent", TR_DOUBLE, false},
	{"ucex","longitudinal critical field temprature exponent", TR_DOUBLE, false},
	{"ibbt","temperature coefficient for ibb (1/K)", TR_DOUBLE, false},
	// Matching parameters
	{"avto","area related threshold voltage mismatch parameter (Vm)", TR_DOUBLE, false},
	{"akp","area related gain mismatch parameter (m)", TR_DOUBLE, false},
	{"agamma","area related body effect mismatch parameter (sqrtV*m)", TR_DOUBLE, false},
	// Flicker noise parameters
	{"kf","flicker noise coefficient", TR_DOUBLE, false},
	{"af","flicker noise exponent", TR_DOUBLE, false},
	// Setubp parameters
	{"nqs","Non-Quasi-Static operation switch", TR_DOUBLE, false},
	{"satlim","ratio defining the saturation limit if/ir", TR_DOUBLE, false},
	{"xqc","charge/capacitance model selector", TR_DOUBLE, false},
};

Mospekv2::Mospekv2(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values which come from ekv documentation
	paramvalue[0] = &(l=1e-6);
	paramvalue[1] = &(w=1e-6);
	paramvalue[2] = &(np=1.0);
	paramvalue[3] = &(ns=1.0);
	paramvalue[4] = &(cox=0.7e-3);
	paramvalue[5] = &(xj=0.1e-6);
	paramvalue[6] = &(dw=0);
	paramvalue[7] = &(dl=0);
	paramvalue[8] = &(vto=-0.5);
	paramvalue[9] = &(gamma=1.0);
	paramvalue[10] = &(phi=0.7);
	paramvalue[11] = &(kp=50.0e-6);
	paramvalue[12] = &(e0=1.0e12);
	paramvalue[13] = &(ucrit=2.0e6);
	paramvalue[14] = &(tox=0);
	paramvalue[15] = &(nsub=0);
	paramvalue[16] = &(vfb=0);
	paramvalue[17] = &(uo=0);
	paramvalue[18] = &(vmax=0);
	paramvalue[19] = &(theta=0);
	paramvalue[20] = &(lambda=0.5);
	paramvalue[21] = &(weta=0.25);
	paramvalue[22] = &(leta=0.1);
	paramvalue[23] = &(q0=0);
	paramvalue[24] = &(lk=0.29e-6);
	paramvalue[25] = &(iba=0);
	paramvalue[26] = &(ibb=3.0e8);
	paramvalue[27] = &(ibn=1);
	paramvalue[28] = &(tcv=-1.0e-3);
	paramvalue[29] = &(bex=-1.5);
	paramvalue[30] = &(ucex=0.8);
	paramvalue[31] = &(ibbt=9.0e-4);
	paramvalue[32] = &(avto=0);
	paramvalue[33] = &(akp=0);
	paramvalue[34] = &(agamma=0);
	paramvalue[35] = &(kf=0);
	paramvalue[36] = &(af=1);
	paramvalue[37] = &(nqs=0);
	paramvalue[38] = &(satlim=54.59815003314424);
	paramvalue[39] = &(xqc=0.4);

  // Set the number of terminals
  setNumTerms(4);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(3);
};

double vtp(double temper)
{
	return 1.3807e-23*temper/1.602e-19; 			// (7) thermal voltage
	//return k*temper/q;
}

double egp(double temper)
{
	return 1.16-(7.02e-4*temper*temper)/(temper + 1108.0); 		// (8) [eV] energy gap
}

double nip(double temper, double trefer) 	// Intrinsic carrier consentration
{
	return 1.45e16*(temper/trefer)*exp(0.5*egp(trefer)/vtp(trefer)
         -0.5*egp(temper)/(vtp(temper)));
         // [m^-3] Intrinsic carrier concentration (9)
}

void Mospekv2::init() throw(string&)
{
  // Set Constants
	q = 1.602e-19;  			// Magnitude of electron charge
	k = 1.3807e-23; 			// Boltzmann's constant
	tref = 300.15;  			// Reference Temperature

	scale=1; //scale factor set to default of 1
	epsilonsi=scale*104.5e-12;			// (4) [f_/m_] p_ermittivity of silicon
	epsilonox=scale*34.5e-12;			// (5) [f_/m_] p_ermittivity of silicon dioxide
 	tnom=300.15; 			// set nominal temp to tref to simplify equations.
	temp=tnom;   			// Could be further implemented. since we have used it throughout

	// Intrinsic prameters initializtion
	// all parameters are set to default unless otherwise spcified
	if(tox>0)
	{ 						// (10)
		cox=epsilonox/tox;
	}

	if(nsub>0)
	{
		gamma=sqrt(2*q*epsilonsi*(nsub*1e6))/cox; 	// (11)
	}

	if(nsub>0)
	{
		phi=2*vtp(tnom)*log((nsub*1e6)/nip(tnom, tref)); 	// (12)
	}

	if(vfb!=0)
	{
		vto=vfb+phi+gamma+sqrt(phi); 			// (13)
	}

	if(uo>0)
	{
		kp=(uo*1e-4)*cox; 				// (14)
	}

	if(vmax>0 && uo>0)
	{
		ucrit=vmax/(uo*1e-4); 				// (15)
	}

	if(theta!=0)
	{
		e0=0;						// (16)
	}

	//Intrinsic parameters temperature dependence
	// Yet to be implemented.
	/*
	FOR TEMPERATURE VARIAnCES

	To be implemented at a later date

	double vto(double vto, double tcv, double temp, double tnom) // (17)
	{
		return vto-tcv*(temp-tnom);
	}

	double kp(double kp, double temp, double tnom, double bex) // (18)
	{
		return kp*pow((temp/tnom),bex);
	}

	double ucrit(double ucrit, double temp, double tnom, double ucex) // (19)
	{
		return ucrit*((temp/tnom),ucex);
	}

	double phi(double phi, double temp, double tnom, double vtp) // (20)
	{
		return phi*(temp/tnom)-3*vtp*log(temp/tnom)-egp(tnom)*(temp/tnom)+egp(temp);
	}

	double ibb(double ibb, double temp, double tnom, double ibbt) // (21)
	{
		return ibb*(1+ibbt*(temp-tnom));
	}
	*/

	// Effective channel length and width
	weff=w+dw; 						// (25)
	leff=l+dl; 						// (26)

	// Short distance matching
	vtoa=vto+avto/sqrt(np*weff*ns*leff); 			// (27)
	kpa=kp*(1+akp/sqrt(np*weff*ns*leff)); 			// (28)
	gammaa=gamma+agamma/sqrt(np*weff*ns*leff); 		// (29)

	// Reverse short-channel effect (RSCE)
	ce=4*(22e-3)*(22e-3); 					// (30)
	ca=0.028;
	xi=ca*(10*leff/lk-1); 					// (31)
	deltavrsce=(2*q0)/cox*1/((1+0.5*(xi+sqrt(xi*xi+ce)))*(1+0.5*(xi+sqrt(xi*xi+ce)))); // (32)

  // For automatic differentiation
  DenseIntVector var(3);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  initializeAD(var, var);
}

void Mospekv2::eval(AD * x, AD * effort, AD * flow)
{
  // Input Matrix
  // x[0]: vds
  // x[1]: vgs
  // x[2]: vbs
  // x[3]: dvds/dt
  // x[4]: dvgs/dt
  // x[5]: dvbs/dt

  // Output Matrix
  // effort[0]: vdb , flow[0]: id
  // effort[1]: vgb , flow[1]: ig
  // effort[2]: vsb , flow[2]: is

  AD vg, vs, vd, vgprime, vp0, vprimes, vprimed,gammadegree, gammaprime, vp, n;
  AD isubf, vc, vdss, vdssprime, deltav, vds2, vip, lc, deltal;
  AD lprime, lmin, leq, virprime, irprime, ir, beta0, nu, qb0;
  AD beta0prime, vpprime, beta, isp, ids, idb, vib;
  AD verytemp, vpdt;
  AD nq, xf, xr, qd, qs, qi, qb;
  AD qid, qig, qis;

  ////////////////INPUT Variables///////////////////
  vg = -(x[1] - x[2]); 	// vgs-vbs   			   (22)
  vs = -(-x[2]);		// -vbs				   (23)
  vd = -(x[0] - x[2]);	// vds-vbs			   (24)

  //  Effective gate voltage including RSCE
  vgprime=vg-vtoa-deltavrsce+phi+gammaa*sqrt(phi); 	// (33)

  // Effective substrate factor including charge-sharing for short and narrow channels
  if (vgprime > 0.0)
    vp0 = vgprime-phi-gammaa*(sqrt(vgprime+(gammaa/2)*(gammaa/2))-gammaa/2);
  else
    vp0 = -phi;

  vprimes = 0.5*(vs+phi+sqrt((vs+phi)*(vs+phi)+(4*vtp(temp))*(4*vtp(temp)))); // (35)
  vprimed = 0.5*(vd+phi+sqrt((vd+phi)*(vd+phi)+(4*vtp(temp))*(4*vtp(temp)))); // (35b)

  gammadegree=gammaa-epsilonsi/cox*(leta/leff*(sqrt(vprimes)+sqrt(vprimed))-(3*weta)/weff*sqrt(vp0+phi)); //(36)

  gammaprime=0.5*(gammadegree+sqrt(gammadegree*gammadegree+0.1*vtp(temp))); 	// (37)

  // Pinch-off voltage including short- and narrow-channel effects
  if (vgprime > 0.0)
    vp = vgprime-phi-gammaprime*(sqrt(vgprime+(gammaprime/2)*(gammaprime/2))-gammaprime/2);
  else
    vp = -phi;

  // Slope Factor

  n=1+gammaa/(2*sqrt(vp+phi+4*vtp(temp))); 			// (39)

  // Forward normalized current

	verytemp=(vp-vs)/vtp(temp);

  //isubf=funct(verytemp);				// (44)
  // Uses smoothing function above called funct

  isubf=log(1+exp(verytemp/2))*log(1+exp(verytemp/2));   //start with assigning value
	/*
	if ((verytemp >-40) && (verytemp < 40))
	{
		isubf=log(1+exp(verytemp/2))*log(1+exp(verytemp/2));
	}
	*/

	if(verytemp>40)
		isubf=(verytemp/2)*(verytemp/2);

  if (verytemp<=-40)
		isubf=exp(verytemp);


	// Velocity saturation voltage
  vc=ucrit*ns*leff; 					// (45)

  vdss = vc*(sqrt(.25 + (vtp(temp)/vc)*sqrt(isubf)) - .5); 	// (46)

  // Drain-to-source saturation voltage for reverse normalized current
  vdssprime = vc*(sqrt(.25 + (vtp(temp)/vc)*(sqrt(isubf) - .75*log(isubf))) - .5) + vtp(temp)*(log(vc/(2*vtp(temp))) - .6); //(47)

  // Channel-length Modulation
  deltav = 4*vtp(temp)*sqrt(lambda*(sqrt(isubf) - vdss/vtp(temp)) + .015625); // (48)

  vds2 = (vd - vs)/2; 					// (49)

  vip = sqrt(vdss*vdss + deltav*deltav) - sqrt((vds2 - vdss)*(vds2 - vdss) + deltav*deltav); //(50)

  lc = sqrt((epsilonsi/cox)*xj); 			// (51)

  deltal = lambda*lc*log(1 + (vds2 - vip)/(lc*ucrit)); 	// (52)

  //Equivalent channel length including channel-length modulation and velocity saturation
  lprime = ns*leff - deltal + (vds2 + vip)/ucrit; 	// (53)
  lmin = ns*leff/10; 					// (54)
  leq = .5*(lprime + sqrt(lprime*lprime + lmin*lmin)); 	// (55)

  // Reverse normalized current
  // Uses smoothing function above called funct
  virprime = ((vp - vds2 - vs - sqrt(vdssprime*vdssprime + deltav*deltav)
  + sqrt((vds2 -vdssprime)*(vds2 - vdssprime) + deltav*deltav))/(vtp(temp))); //(56a)

  irprime=log(1+exp(virprime/2))*log(1+exp(virprime/2));  //(56b)
	/*

	if ((virprime >-40) && (virprime < 40))
	{
		irprime = log(1+exp(virprime/2))*log(1+exp(virprime/2));
	}

	*/

	if (virprime>=40)
		irprime = (virprime/2)*(virprime/2);

	if (virprime<=-40)
		irprime = exp(virprime);

  vpdt=(vp-vd)/vtp(temp);
  // Uses smoothing function above called funct

  ir=log(1+exp(vpdt/2))*log(1+exp(vpdt/2));
	if (vpdt>=40)
		ir=(vpdt/2)*(vpdt/2);
	if (vpdt<=-40)
		ir=exp(vpdt);
	if ((vpdt >-40) && (vpdt < 40))
		ir=log(1+exp(vpdt/2))*log(1+exp(vpdt/2));

	// Transconductance factor and mobility reduction due to vertical field
	beta0 = kpa*((np*weff)/leq); 				//(58)

	nu = 1/3; //FOR PMOS

  //Quasi-static model equations... Yet to be implemented
  /*
	nq = 1 + gammaa/(2*sqrt(vp+phi+0.000001));  	// (69)

	xf = sqrt(1/4 + isubf);				// (70)

	xr = sqrt(1/4 + ir);				// (71)

	qd = -nq*((4/15)*(3*xr*xr*xr + 6*xr*xr*xf + 4*xr*xf*xf + 2*xf*xf*xf)/((xf + xr)*(xr + xf)) - .5); //(72)

	qs = -nq*((4/15)*(3*xf*xf*xf + 6*xf*xf*xr + 4*xf*xr*xr + 2*xr*xr*xr)/((xf + xr)*(xr + xf)) - .5); //(73)

	qi = qs + qd; 					// (74)

	if (vgprime > 0)
	{
		qb = (-gammaa*sqrt(vp + phi + 0.000001)*(1/vtp(temp)) - ((nq - 1)/nq)*qi;
	}
	else
	{
		qb = -vgprime*(1/vtp(temp));
	} 						// (75)

  */
  // Transconductance factor and mobility reduction due to vertical field cont'd
  // more complicated way of calculating beta and we decided to go with simplified model provided below
  /*
	qb0 = gammaa*sqrt(phi); 				// (60)

	beta0prime = beta0*(1 + cox/(e0*epsilonsi)*qb0); 	// (61)

	beta = beta0prime/(1 + cox/(e0*epsilonsi)*vtp(temp)*sqrt((qb + nu*qi)*(qb + nu*qi))); // (62)
  */

  // Transconductance factor and mobility reduction due to vertical field cont'd
  //Easy way of calculating beta

  vpprime = .5*(vp + sqrt(vp*vp + 2*vtp(temp)*vtp(temp))); 		// (63)

  beta = beta0/(1 + theta*vpprime); 			// (64)

  // Specific current

  isp = 2*n*beta*vtp(temp)*vtp(temp); 				// (65)

  // Drain-to-source current

  ids = isp*(isubf - irprime); 				// (66)

  // Impact ionization current

  vib= vd-vs-ibn*2*vdss;  				//(67)

	if (vib > 0)
		idb = ids*(iba/ibb)*vib*exp((-ibb*lc)/vib);
	else
		idb = 0.0; // (68)

  //////////////OUTPUT Variables/////////////////
  effort[0] = vd;
  effort[1] = vg;
  effort[2] = vs;

  flow[0] = -(ids + idb);	// Total drain current
  flow[1] = 0.0;		// gate current
  flow[2] = ids;		// total souce current
}
