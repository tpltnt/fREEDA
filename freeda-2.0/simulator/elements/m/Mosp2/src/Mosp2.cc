#include "Mosp2.h"

// Static members
const unsigned Mosp2::n_par = 31;

// Element information
ItemInfo Mosp2::einfo =
{
  "mosp2",
  "Level 2 Grove-Frohman Model",
  "Aaron Walker",
  DEFAULT_ADDRESS"transistor>mosfet",
  "2002_10_10"
};

// Parameter information
ParmInfo Mosp2::pinfo[] =
{
  {"vt0",    "Zero-bias threshold voltage (V)", TR_DOUBLE, false},
  {"kp",     "Transconductance parameter (A/V^2)", TR_DOUBLE, false},
  {"gamma",  "Bulk threshold parameter (V^1/2)",TR_DOUBLE, false},
  {"phi",    "Surface inversion potential (V)", TR_DOUBLE, false},
  {"lambda", "Channel-length modulation (1/V)", TR_DOUBLE, false},
  {"rd",     "Drain ohmic resistance", TR_DOUBLE, false},
  {"rs",     "Source ohmic resistance", TR_DOUBLE, false},
  {"is",     "Bulk junction saturation current", TR_DOUBLE, false},
  {"pb",     "Bulk junction potential", TR_DOUBLE, false},
  {"js",     "Bulk junction saturation current density", TR_DOUBLE, false},
  {"tox",    "Oxide thickness (m)", TR_DOUBLE, false},
  {"ld",     "Lateral diffusion (m)", TR_DOUBLE, false},
  {"u0",     "Surface mobility (cm^2/V-s)", TR_DOUBLE, false},
  {"fc",     "Forward bias junction fit parameter", TR_DOUBLE, false},
  {"nsub",   "Substrate doping (1/cm^3)", TR_DOUBLE, false},
  {"tpg",    "Type of gate material",TR_DOUBLE, false},
  {"nss",    "Surface state density.", TR_DOUBLE, false},
  {"delta",  "Width effect on threshold", TR_DOUBLE, false},
  {"uexp",   "Crit. field exp for mob. deg.",TR_DOUBLE, false},
  {"ucrit",  "Crit. field for mob. degradation",TR_DOUBLE, false},
  {"vmax",   "Maximum carrier drift velocity",TR_DOUBLE, false},
  {"xj",     "Junction depth",TR_DOUBLE, false},
  {"neff",   "Total channel charge coeff.",TR_DOUBLE, false},
  {"nfs",    "Fast surface state density",TR_DOUBLE, false},
  {"tnom",   "Nominal device temperature (C)", TR_DOUBLE, false},
  {"kf",     "Flicker noise coefficient", TR_DOUBLE, false},
  {"af",     "Flicker noise exponent", TR_DOUBLE, false},
  {"t",      "Device temperature (C)", TR_DOUBLE, false},
  {"l",      "Gate length (m)", TR_DOUBLE, false},
  {"w",      "Gate width (m)", TR_DOUBLE, false},
  {"alpha",  "Impact ionization current coefficient (1/V)", TR_DOUBLE, false},
};


Mosp2::Mosp2(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0]  = &(vt0 = 0.0);
  paramvalue[1]  = &(kp = 2.0e-5);
  paramvalue[2]  = &(gamma = 0.0);
  paramvalue[3]  = &(phi = 0.6);
  paramvalue[4]  = &(lambda = 0.0);
  paramvalue[5]  = &(rd = 0.0);
  paramvalue[6]  = &(rs = 0.0);
  paramvalue[7]  = &(is = 1e-14);
  paramvalue[8]  = &(pb = 0.8);
  paramvalue[9]  = &(js = 0.0);
  paramvalue[10] = &(tox = 1.0e-7);
  paramvalue[11] = &(ld = 0.0);
  paramvalue[12] = &(u0 = 600.0);
  paramvalue[13] = &(fc = 0.5);
  paramvalue[14] = &(nsub = 1e15);
  paramvalue[15] = &(tpg = 1.0);
  paramvalue[16] = &(nss = 0.0);
  paramvalue[17] = &(delta = 0);	// MOS2narrowFactor, default=0
  paramvalue[18] = &(uexp = 0); 	// MOS2critFieldExp, default=0
  paramvalue[19] = &(ucrit = 1e4);	// MOS2critField, default=1e4
  paramvalue[20] = &(vmax = 0); 	// MOS2maxDriftVel, default=0
  paramvalue[21] = &(xj = 0);		// MOS2junctionDepth, default=0
  paramvalue[22] = &(neff = 1); 	// MOS2channelCharge, default=1
  paramvalue[23] = &(nfs = 0);  	// MOS2fastSurfaceStateDensity, default=0
  paramvalue[24] = &(tnom = 27);
  paramvalue[25] = &(kf = 0.0);
  paramvalue[26] = &(af = 1);
  paramvalue[27] = &(t = 27);
  paramvalue[28] = &(l = 2.0e-6);
  paramvalue[29] = &(w = 50e-6);
  paramvalue[30] = &(alpha = 0.0);

  // Set the number of terminals
  setNumTerms(4);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(3);
}

void Mosp2::init() throw(string&)
{
  // Set Constants
  pi = 3.14159265358979323846;
  epsilon0 = 8.854194e-12;
  epsilons = 11.7 * epsilon0;
  Kox = 3.9;
  tnom = tnom + 273.15;
  tref = tnom;			// For now
  tdev = t + 273.15;
  type = -1.0;			// Hardcoded to PMOS
  nssm = nss * 1e4;		// Converted to 1/m^2

  // Model constants
  vtnom = (kBoltzman * tnom)/eCharge;
  leff = l - 2.0 * ld;

  // Calculate any missing parameters from user-defined settings
  // COX (F/m^2)
  if(cox == 0.0)
	{
    if (tox != 0.0)
      cox = (Kox * epsilon0) / tox;
    else
      cox = 0.0;
	}

  // KP
  if(kp == 0.0)
    if((u0 != 0.0) && (tox != 0.0))
      kp = u0 * cox * 1e-4 /*(m^2/cm^2)*/;
	else
		kp = 2.0718e-5;

  // PHI
  // First calculate the energy gap
  egfet1 = 1.16 - 7.02e-4 * tnom * tnom / (tnom + 1108);
  // The intrinsic carrier concentration at 300 K is 1.45e16 m^-3
  ni = 1.45e16;
  if(phi == 0.0)
    phi = 2.0 * vtnom * log(nsub * 1e6 / ni);

  // Calculate common square root of phi term
  sqrphi = sqrt(phi);

  // GAMMA
  if(gamma == 0.0)
    gamma_cmp = sqrt(2.0 * eCharge * epsilons * nsub * 1e6) / cox;
  else
    gamma_cmp = gamma;

  // VTO, assume poly gate.
  // Calculate channel Fermi voltage
  fermis = type * 0.5 * phi;
  // The work function voltage for the gate is given by
  wkfng = 3.2;
  // Test the gate type
  if(tpg != 0 || tpg != 1 || tpg || -1)
    tpg = 1;
  if(tpg != 0)
	{
    fermig = type * tpg * 0.5 * egfet1;
    wkfng = 3.25 + .5 * egfet1 - fermig;
  }
  wkfngs =  wkfng - (3.25 + .5 * egfet1 + fermis);

  // Now calculate vfb
  vfb = wkfngs - eCharge * nssm / cox;

  // Zero-bias threshold voltage
  if (vt0 == 0.0)
    vt0 = vfb + type * (gamma_cmp * sqrphi + phi);

  // Calculate temperature dependent values that are not
  // functions of state variables, assume temperature is constant
  // for entire simulation.  May change if thermal ports added
  // to model later.  Temperature model used is equivalent to
  // TLEV=0 in HSPICE.

  // Factors needed in temperature equations
  ratio = tdev / tnom;
  fact1 = tnom / tref;
  fact2 = tdev / tref;		// Should be tdev / tref, need to check
  ratio4 = ratio * sqrt(ratio);
  kt = tdev * kBoltzman;
  vt = kt / eCharge;
  ktnom = tnom * kBoltzman;
  egfet = 1.16 - (7.02e-4 * tdev * tdev) / (tdev + 1108);
  arg = -egfet / (kt + kt) + 1.1150877/(ktnom + ktnom);
  pbfact = -2.0*vt *(1.5*log(fact2) + eCharge * arg);
  arg1 = -egfet1/(ktnom+ktnom)+1.1150877/(kBoltzman*(tref+tref));
  pbfact1 = -2.0*vtnom * (1.5*log(fact1) + eCharge*arg1);
  pbo = (pb - pbfact1)/fact1;

  // Bulk Potential
  tpb = fact2 * pbo + pbfact;

  // Transconductance kp
  tkp = kp / ratio4;

  // Surface mobility
  tu0 = u0 / ratio4;

  // PHI
  phio = (phi - pbfact1) / fact1;
  tphi = fact2 * phio + pbfact;
  sqrtphi = sqrt(tphi);
  stphi3 = tphi * sqrtphi;

  // Built in voltage vbi
  tvbi = vt0-type*(gamma*sqrphi)+0.5*(egfet1-egfet)+type*0.5*(tphi-phi);

  // Zero-bias threshold voltage
  tvt0 = tvbi + type * gamma_cmp * sqrtphi;

  // Miscellaneous terms
  beta = tkp * w / leff;
  xd = sqrt((epsilons+epsilons)/(eCharge*nsub*1e6));
  sbiarg = sqrt(tpb);

  DenseIntVector var(3);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  initializeAD(var, var);
}

void Mosp2::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: vds
  // x[1]: vgs
  // x[2]: vbs
  // x[3]: dvds/dt
  // x[4]: dvgs/dt
  // x[5]: dvbs/dt
  // effort[0]: vdb , flow[0]: id
  // effort[1]: vgb , flow[1]: ig
  // effort[2]: vsb , flow[2]: is

  AD vth, veff, Ids, Ion, von;
  AD arg, sarg, sarg1, sarg2, sargtemp, vgst;
  AD phiMinVbs, Xs, X1, Xb, Xbs, X2term, X2, vbin;
  AD args, argd, Fdd, Fsd, dbxwd, dbxws, temp;
  AD dbargs, dbargd, gammasd, vdsat, cfs, xn;
  AD argg, Xs3, body, ufact, ueff, vgsx, gammad;
  AD gammad2, vdsat1, bodys, argv, sargv;
  AD xlfact, xlambda, xdv, xlv, xls;
  AD xwb, xld, clfact, xleff, deltal, beta1, vdson;
  AD bodysub, Isub1, Isub2, Ilin, Isat, Ion1, Ion2;
  AD alphax, term1, term2, qdrain, qsource, qbulk, qgate;

  // DC Current calculation, follows Berkely SPICE formulation
  // No test for inverse connection of source and drain nodes from node convention
  // m1: drain gate source bulk

  // Compute some needed factors
  phiMinVbs = tphi - type*(x[2]);
  xlambda = lambda;

  // Xs, Xb, X1 computaitions
  if (type*(x[2]+0.0) > 0.0)
    Xs = sqrtphi/(1.0+0.5*type*x[2]/tphi);
  else
    Xs = sqrt(phiMinVbs);

  if (type*(x[2]+0.0) > 0.0)
    X1 = -0.5*Xs*Xs/stphi3;
  else
    X1 = -0.5/Xs;

  if (type*(x[2]- x[0]) > 0.0)
    Xb = sqrtphi/(1.0+0.5*type*(x[2]-x[0])/tphi);
  else
    Xb = sqrt(phiMinVbs + type*x[0]);

  if (type*(x[2]-x[0]) > 0.0)
    X2term = -0.5*Xb*Xb/stphi3;
  else
    X2term = -0.5/Xb;

  // Calculate threshold voltage narrow channel effect
  factor = 0.25 * delta * pi * epsilons/(cox * w);
  eta = 1.0 + factor;
  vbin = tvbi*type + factor*phiMinVbs;

	//  printf("\n\ntvbi = %0e\n", tvbi);
  // Calculate effective body effect gamma
  Fdd = 0.0;
  Fsd = 0.0;
  dbargs = 0.0;
  dbargd = 0.0;
  if(gamma > 0.0 || nsub > 0.0)
	{
    if(xj > 0)
		{
      args = sqrt(1+2*xd*Xs/xj);
      argd = sqrt(1+2*xd*Xb/xj);
      Fdd = 0.5*xj/leff*( argd - 1.0);
      Fsd = 0.5*xj/leff*( args - 1.0);
    }
    dbxwd = xd*X2term;
    dbxws = xd*X1;
    temp = 0.5/leff;
    dbargs = temp*dbxws/args;
    dbargd = temp*dbxwd/argd;
    X2 = -gamma*(dbargs + dbargd);
  }
  else
	{
		X2 = 0.0;
  }
  gammasd = gamma*(1-Fdd-Fsd);

  // Initial threshold voltage
  von = (vbin+gammasd*Xs);
  vth = von;
  vdsat = 0.0;

  // Threshold voltage modification for fast surface state density, nfs
  if(nfs != 0.0 && cox != 0.0)
	{
    cfs = eCharge*nfs*1e4;
    xn = 1.0 + factor - gammasd*X1 - X2*Xs + w*leff*cfs/cox;
    temp = xn * vt;
    von = von + temp;
    argg = 1.0/temp;
  }

  // Calculate effective voltage
  veff = type*(x[1] - von);

  // Compute constants
  Xs3 = Xs*Xs*Xs;
  body = Xb*Xb*Xb - Xs3;

  // Compute mobility modification factor
  // Default value
  ufact = 1.0;
  ueff = u0*1e-4;
  if(nfs == 0.0)
	if(cox != 0.0)
	{
		temp = ucrit * 100.0 * epsilons / cox;
		// If temp > veff ufact = 1.0 else ufact = exp(uexp*log(temp/veff))
    if (temp > veff)
      ufact = 1.0;
    else
      ufact = exp(uexp*log(temp/veff));

    if (temp > veff)
      ueff = u0*1e-4;
    else
      ueff = u0*1e-4*ufact;
	}

  // Compute saturation voltage using the Grove-Frohman equation
  vgsx = type*x[1];
  gammad = gammasd/eta;
  if( (nfs != 0.0) && (cox != 0.0))
	{
    if (type*(von-x[1]) > 0.0)
      vgsx = von;
    else
      vgsx = type*x[1]+0.0;
  }

  if(gammad > 0)
	{
    gammad2 = gammad*gammad;
    argv = (vgsx - vbin)/eta + phiMinVbs;
    arg = sqrt(1.0+4.0*argv/gammad2);
    if (argv > 0.0)
      vdsat1 = (vgsx-vbin)/eta+gammad2*(1.0-arg)/2.0;
    else
      vdsat1 = 0.0;
  }
  else
    vdsat1 = (vgsx-vbin)/eta;

  // Assure that computed saturation voltage is > 0
  if (vdsat1 > 0.0)
    vdsat = vdsat1;
  else
    vdsat = 0.0;

  // Evaluate effective channel length
  // Compute constants with the new vdsat
  if (type*x[2] > vdsat)
    Xbs = sqrtphi/(1.0+0.5*(type*x[2]-vdsat)/tphi);
  else
    Xbs = sqrt(phiMinVbs + vdsat);
  bodys = Xbs*Xbs*Xbs - Xs3;
  gammad = gammasd;

  // Only compute the lambda parameter if it is not supplied by the user.
  // Also if vds = 0.0, there is no channel length modulations so the
  // parameter is not calculated.
  if(nsub != 0.0)
	if (xlambda <= 0.0)
	{
		if(vmax <= 0.0)
		{
			argv = (type*x[0]-vdsat)/4.0;
			sargv = sqrt(1.0+argv*argv);
			arg = sqrt(argv+sargv);
      if (type*x[0]+0.0 > 0.0)
        xlfact = xd/(leff*type*x[0]);
      else
        xlfact = 0.0;
			xlambda = xlfact*arg;
		}
		else
		{
			argv = (vgsx-vbin)/eta-vdsat;
			xdv = xd/sqrt(neff);
			xlv = vmax*xdv/(2.0*ueff);
      if (type*x[0] > vdsat)
        argv = type * x[0] - vdsat;
      else
        argv = 0.0;
			xls = sqrt(xlv*xlv+argv);
      if (x[0] + 0.0 > 0.0)
        xlfact = xdv/(leff*type*x[0]);
      else
        xlfact = 0.0;
			xlambda = xlfact*(xls-xlv);
		}
  }

  // Test for vds == 0.0, if xlambda was not computed above, it will
  // equal the user-supplied lambda.
  if (x[0]*x[0] <= 0.0)
    xlambda = lambda;

  // limit channel shortening at punch-through
  xwb = xd*sbiarg;
  xld = leff-xwb;
  clfact = 1.0-xlambda*type*x[0];
  xleff = leff*clfact;
  deltal = xlambda*type*x[0]*leff;
  if (nsub == 0.0)
	{
    xwb = 0.25e-6;
  }
  if (xwb > xleff)
    xleff = xwb/(1.0+(deltal-xld)/xwb);
  clfact = xleff/leff;

  // All voltages are retained in the same form between NMOS and PMOS.  The
  // difference is the computed currents are negative for the PMOS case and
  // the definitions of saturation and linear regions are reversed.

  // Compute effective beta
  beta1 = beta * ufact/clfact;

	// Find vdson, minimum of vds and vdsat
  if (type*x[0] > vdsat)
    vdson = vdsat;
  else
    vdson = type*x[0] + 0.0;

  // Compute subthreshold current, if vdsat <= 0 then the subthreshold
  // current is zero.
  // body term depends on vds
  if (type*x[0] > vdsat)
    bodysub = bodys;
  else
    bodysub = body;
  Isub1 = beta1*((von-vbin-eta*vdson/2.0)*vdson-gammad*bodysub/1.5);
  Isub1 = Isub1*exp(argg*(type*x[1]-von));
  if (vdsat > 0.0)
    Isub2 = Isub1;
  else
    Isub2 = 0.0;

  // Compute linear current
  Ilin = beta1*((type*x[1]-vbin-eta*type*x[0]/2.0)*type*x[0]-gammad*body/1.5);

  // Compute saturation current
  Isat = beta1*((type*x[1]-vbin-eta*vdsat/2.0)*vdsat-gammad*bodys/1.5);

  // Decide between linear and saturation regions for intermediate
  // current assignment.  If vdsat > x[0] => linear, else => satur.
  if (vdsat > type*x[0])
    Ion1 = Ilin;
  else
    Ion1 = Isat;

  // Decide between on and subthreshold for intermediate
  // current assignment. If x[1] > von, then FET is on, otw off/subth.
  if (type*x[1] > von)
    Ion2 = Ion1;
  else
    Ion2 = Isub2;

  // If vds > 0 then the FET is on, or in subthreshold, otherwise
  // it is cutoff.
  if (type*x[0] < 1e-10)
    Ids = 0.0;
  else
    Ids = Ion2;

  // The capacitance equations have been removed and a charge
  // conservation model has been implemented.  This model ensures
  // that charge is conserved in transient analysis which is not
  // the case with the Meyer capacitance model.  The charge model
  // is the Yang-Chatterjee model.

  // Compute some common coefficients
  // Short channel effects coefficient
  vth = type*vth;
  alphax = alpha + gamma_cmp * type * (x[1] - vth);
  coxwl = cox * w * leff;
  gammasqr = gamma_cmp * gamma_cmp;

  term1 = alphax * x[0] * x[0] / (4.0 * type * (x[1] - vth - alphax * x[0] / 2.0));
  term2 = type * (x[1] - vth) / 2.0;

  // Accumulation region
  if((type*x[1]) < (type*(vfb + x[2])) )
	{
    qdrain = 0;
    qsource = 0;
    qbulk = -coxwl*type*(x[1] - vfb - x[2]);
    qgate = -(qdrain + qsource + qbulk);
  }

  // Subthreshold region
  else if(( (type*x[1]) >= (type*(vfb + x[2])) ) && ( (type*x[1]) < (type*vth) ))
	{
    qdrain = 0;
    qsource = 0;
    qbulk = -coxwl*gammasqr/2.0*(-1.0 + sqrt((1.0 + 4.0*type*(x[1] - vfb - x[2]))/gammasqr));
    qgate = -(qdrain + qsource + qbulk);
  }

  // Linear region
  else if( (type*(x[1] - vth)) > (type*x[0]) )
	{
    qdrain = -coxwl*(term2 - 3.0 * alphax * type * x[0]/4.0 + alphax*term1/2.0);
    qsource = -coxwl*(term2 + alphax * type * x[0]/4.0 - alphax*term1/6.0);
    qbulk = coxwl*( (type * (vfb  - vth + (1.0-alphax)*x[0]/2) + tphi) - (1-alphax)*term1/3.0);
    qgate = -(qdrain + qsource + qbulk);
  }

  // Saturation region
  else
	{
    qdrain = 0;
    qsource = -coxwl *(2.0/3.0 * type * (x[1] - vth));
    qbulk = coxwl * (type * (vfb - vth + (1.0-alphax)*(x[1] - vth)/(3.0*alphax)) + tphi) ;
    qgate = -(qdrain + qsource + qbulk);
	}

  // Calculate dynamic current contributions due to charge
  // iqd is dynamic contribution to drain current
  // iqg is dynamic contribution to gate current
  // iqs is dynamic contribution to source current
  AD iqd, iqg, iqs;
  iqd = type * (qdrain.fastAccessDx(0)*x[3] + qdrain.fastAccessDx(1)*x[4] + qdrain.fastAccessDx(2)*x[5]);
  iqg = type * (qgate.fastAccessDx(0)*x[3] + qgate.fastAccessDx(1)*x[4] + qgate.fastAccessDx(2)*x[5]);
  iqs = type * (qsource.fastAccessDx(0)*x[3] + qsource.fastAccessDx(1)*x[4] + qsource.fastAccessDx(2)*x[5]);

  flow[0] = type * Ids + iqd; // DC + dynamic drain current
  flow[1] = iqg; // DC + dynamic gate current
  flow[2] = -type *Ids - iqs; 	// DC + dynamic source current

  effort[0] = x[0] - x[2];  // vdb = vds - vbs
  effort[1] = x[1] - x[2];	// vgb = vgs - vbs
  effort[2] = -x[2];	// vsb = -vbs
}
