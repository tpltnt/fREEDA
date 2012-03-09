#include "Mosp1.h"

// Static members
const unsigned Mosp1::n_par = 34;

// Element information
ItemInfo Mosp1::einfo =
{
  "mosp1",
  "Schichman-Hodges Mosfet model",
  "Aaron Walker",
  DEFAULT_ADDRESS"transistor>mosfet",
  "2002_10_10"
};

// Parameter information
ParmInfo Mosp1::pinfo[] =
{
  {"vt0", "Zero-bias threshold voltage (V)", TR_DOUBLE, false},
  {"kp", "Transconductance parameter (A/V^2)", TR_DOUBLE, false},
  {"gamma", "Bulk threshold parameter (V^1/2)",TR_DOUBLE, false},
  {"phi", "Surface inversion potential (V)", TR_DOUBLE, false},
  {"lambda", "Channel-length modulation (1/V)", TR_DOUBLE, false},
  {"rd", "Drain ohmic resistance", TR_DOUBLE, false},
  {"rs", "Source ohmic resistance", TR_DOUBLE, false},
  {"cbd", "B-D junction capacitance", TR_DOUBLE, false},
  {"cbs", "B-S junction capacitance", TR_DOUBLE, false},
  {"is", "Bulk junction saturation current", TR_DOUBLE, false},
  {"pb", "Bulk junction potential", TR_DOUBLE, false},
  {"cgso", "Gate-source overlap capacitance", TR_DOUBLE, false},
  {"cgdo", "Gate-drain overlap capacitance", TR_DOUBLE, false},
  {"cgbo", "Gate-bulk overlap capacitance", TR_DOUBLE, false},
  {"rsh", "Sheet resistance", TR_DOUBLE, false},
  {"cj", "Bottom junction capacitance per area", TR_DOUBLE, false},
  {"mj", "Bottom grading coefficient", TR_DOUBLE, false},
  {"cjsw", "Side junction capacitance per area", TR_DOUBLE, false},
  {"mjsw", "Side grading coefficient", TR_DOUBLE, false},
  {"js", "Bulk junction saturation current density", TR_DOUBLE, false},
  {"tox", "Oxide thickness (m)", TR_DOUBLE, false},
  {"ld", "Lateral diffusion (m)", TR_DOUBLE, false},
  {"u0", "Surface mobility (cm^2/V-s)", TR_DOUBLE, false},
  {"fc", "Forward bias junction fit parameter", TR_DOUBLE, false},
  {"nsub", "Substrate doping (1/m^3)", TR_DOUBLE, false},
  {"tpg", "Type of gate material",TR_DOUBLE, false},
  {"nss", "Surface state density.", TR_DOUBLE, false},
  {"tnom", "Nominal device temperature (C)", TR_DOUBLE, false},
  {"kf", "Flicker noise coefficient", TR_DOUBLE, false},
  {"af", "Flicker noise exponent", TR_DOUBLE, false},
  {"t", "Device temperature (C)", TR_DOUBLE, false},
  {"l", "Gate length (m)", TR_DOUBLE, false},
  {"w", "Gate width (m)", TR_DOUBLE, false},
  {"alpha", "Impact ionization current coefficient (1/V)", TR_DOUBLE, false}
};

Mosp1::Mosp1(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0]  = &(vt0 = 0.0);
  paramvalue[1]  = &(kp = 2.0e-5);
  paramvalue[2]  = &(gamma = 0.0);
  paramvalue[3]  = &(phi = 0.6);
  paramvalue[4]  = &(lambda = 0.0);
  paramvalue[5]  = &(rd = 0.0);
  paramvalue[6]  = &(rs = 0.0);
  paramvalue[7]  = &(cbd = 0.0);
  paramvalue[8]  = &(cbs = 0.0);
  paramvalue[9]  = &(is = 1e-14);
  paramvalue[10] = &(pb = 0.8);
  paramvalue[11] = &(cgso = 0.0);
  paramvalue[12] = &(cgdo = 0.0);
  paramvalue[13] = &(cgbo = 0.0);
  paramvalue[14] = &(rsh = 0.0);
  paramvalue[15] = &(cj = 0.0);
  paramvalue[16] = &(mj = 0.5);
  paramvalue[17] = &(cjsw = 0.0);
  paramvalue[18] = &(mjsw = 0.5);
  paramvalue[19] = &(js = 0.0);
  paramvalue[20] = &(tox = 1.0e-7);
  paramvalue[21] = &(ld = 0.0);
  paramvalue[22] = &(u0 = 600.0);
  paramvalue[23] = &(fc = 0.5);
  paramvalue[24] = &(nsub = 1e15);
  paramvalue[25] = &(tpg = 1.0);
  paramvalue[26] = &(nss = 0.0);
  paramvalue[27] = &(tnom = 27);
  paramvalue[28] = &(kf = 0.0);
  paramvalue[29] = &(af = 1);
  paramvalue[30] = &(t = 27);
  paramvalue[31] = &(l = 2.0e-6);
  paramvalue[32] = &(w = 50e-6);
  paramvalue[33] = &(alpha = 0.0);

  // Set the number of terminals
  setNumTerms(4);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(3);
}

void Mosp1::init() throw(string&)
{
  // Set Constants
  pi = 3.14159265358979323846;
  epsilon0 = 8.854194e-12;
  epsilons = 11.7 * epsilon0;
  Kox = 3.9;
  tnom = tnom + 273.15;
  tref = tnom;			// For now
  tdev = t + 273.15;
  type = -1.0;  		// Hardcoded to PMOS
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
    gamma = sqrt(2.0 * eCharge * epsilons * nsub * 1e6) / cox;

  // VTO, assume poly gate.

  // Calculate channel Fermi voltage
  fermis = type * 0.5 * phi;
  // The work function voltage for the gate is given by
  wkfng = 3.2;
  // Test the gate type
  if(tpg != 0 || tpg != 1 || tpg || -1)
    tpg = 1;
  if(tpg != 0){
    fermig = type * tpg * 0.5 * egfet1;
    wkfng = 3.25 + .5 * egfet1 - fermig;
  }
  wkfngs =  wkfng - (3.25 + .5 * egfet1 + fermis);

  // Now calculate vfb
  vfb = wkfngs - eCharge * nssm / cox;

  // Zero-bias threshold voltage
  if (vt0 == 0.0)
    vt0 = vfb + type * (gamma * sqrphi + phi);

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

  // Transconductance kp
  tkp = kp / ratio4;

  // Surface mobility
  tu0 = u0 / ratio4;

  // PHI
  phio = (phi - pbfact1) / fact1;
  tphi = fact2 * phio + pbfact;
  sqrtphi = sqrt(tphi);

  // Built in voltage vbi
  tvbi = vt0-type*(gamma*sqrphi)+0.5*(egfet1-egfet)+type*0.5*(tphi-phi);

  // Zero-bias threshold voltage
  tvt0 = tvbi + type * gamma * sqrtphi;

  // Miscellaneous terms
  beta = tkp * w / leff;

  DenseIntVector var(3);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  initializeAD(var, var);
}


void Mosp1::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: vds
  // x[1]: vgs
  // x[2]: vbs
  // x[3]: dvds/dt
  // x[4]: dvgs/dt
  // x[5]: dvbs/dt
  // effort[3]: vdb , flow[3]: id
  // effort[4]: vgb , flow[4]: ig
  // effort[5]: vsb , flow[5]: is

  //All the active variables for ADOL-C "adoubles" must be initialiazed here
  AD vth, veff, Ids, Ion ;
  AD arg, betap, sarg, sarg1, sarg2, sargtemp, vgst;
  AD alphax, term1, term2, qdrain, qsource, qbulk, qgate;

  // DC Current calculation, follows Berkely SPICE formulation
  // No test for inverse connection of source and drain nodes from node convention
  // m1 drain gate source bulk
  sarg1 = sqrt(tphi - x[2]);
  sargtemp = sqrt(tphi);
  sargtemp = sargtemp - x[2] / (sargtemp + sargtemp);
  if (sargtemp > 0.0)
    sarg2 = sargtemp;
  else
    sarg2 = 0.0;
  // If vbs > 0 then sarg = sarg2 else sarg = sarg1
  if ( x[2] > 0.0)
    sarg = sarg2;
  else
    sarg = sarg1;

  // All voltages are retained in the same form between NMOS and PMOS.  The
  // difference is the computed currents are negative for the PMOS case and
  // the definitions of saturation and linear regions are reversed.

  // Compute threshold voltage, ensure that it is negative for PMOS
  vth = type * (tvbi * type + gamma *sarg);

  // Compute effective voltage
  veff = x[1] - vth;

  // Compute on current
  betap = beta * (1.0 + lambda * x[0] * type);
  if (type*(veff - x[0]) > 0.0)
    Ion = betap * (veff - x[0] * 0.5) * x[0];
  else
    Ion = 0.5 * betap * veff * veff;

  // Depending on the sign of veff, assign the channel current Ids
  if (type*veff > 0.0)
    Ids = Ion;
  else
    Ids = 0.0;

  // The capacitance equations have been removed and a charge
  // conservation model has been implemented.  This model ensures
  // that charge is conserved in transient analysis which is not
  // the case with the Meyer capacitance model.  The charge model
  // is the Yang-Chatterjee model.

  // Compute some common coefficients
  // Short channel effects coefficient
  alphax = alpha + gamma * type * (x[1] - vth);
  coxwl = cox * w * leff;
  gammasqr = gamma * gamma;

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
