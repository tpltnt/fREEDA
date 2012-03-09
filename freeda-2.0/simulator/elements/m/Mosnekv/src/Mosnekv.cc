#include "Mosnekv.h"

// Static members
const unsigned Mosnekv::n_par = 44;

// Element information
ItemInfo Mosnekv::einfo =
{
  "mosnekv",
  "EPFL EKV MOSFET model",
  "Wonhoon Jang",
  DEFAULT_ADDRESS"transistor>mosfet",
  "2003_05_15"
};

// Parameter information
ParmInfo Mosnekv::pinfo[] =
{
  {"type", "N-channel or P-channel MOS", TR_DOUBLE, false},
  // Device input variables
  {"l", "Gate length (m)", TR_DOUBLE, false},
  {"w", "Gate width (m)", TR_DOUBLE, false},
  {"np", "Parallel multiple device number", TR_DOUBLE, false},
  {"ns", "Serial multiple device number", TR_DOUBLE, false},
	// Process related parameters
  {"cox", "Gate oxide capacitance per area (F/m^2)", TR_DOUBLE, false},
  {"xj", "Junction depth (m)", TR_DOUBLE, false},
  {"dw", "Channel width correction (m)", TR_DOUBLE, false},
  {"dl", "Channel length correction (m)", TR_DOUBLE, false},
	// Basic intrinsic model parameters
  {"vto", "Long_channel threshold voltage (V)", TR_DOUBLE, false},
  {"gamma", "Body effect parameter (V^1/2)",TR_DOUBLE, false},
  {"phi", "Bulk Fermi potential (V)", TR_DOUBLE, false},
  {"kp", "Transconductance parameter (A/V^2)", TR_DOUBLE, false},
  {"eo", "Mobility reduction coefficient (V/m)",TR_DOUBLE, false},
  {"ucrit", "Longitudinal critical field (V/m)", TR_DOUBLE, false},
	// Operational parameters
  {"tox", "Oxide thickness (m)", TR_DOUBLE, false},
  {"nsub", "Channel doping (1/cm^3)", TR_DOUBLE, false},
  {"vfb", "Flat-band voltage (V)", TR_DOUBLE, false},
  {"uo", "Low-field mobility(cm^2/(V.s))", TR_DOUBLE, false},
  {"vmax", "Saturation velocity (m/s)", TR_DOUBLE, false},
  {"theta", "Mobility recuction coefficient (1/V)", TR_DOUBLE, false},
	// Channel length modulation and charge sharing parameters
  {"lambda", "Channel-length modulation", TR_DOUBLE, false},
  {"weta", "Narrow-channel effect coefficient", TR_DOUBLE, false},
  {"leta", "Short-channel effect coefficient", TR_DOUBLE, false},
	// Reverse short-channel effect parameters
  {"qo", "Reverse short channel effect peak charge density (A.s/m^2)", TR_DOUBLE, false},
  {"lk", "Reverse short channel effect characteristic length (m)", TR_DOUBLE, false},
	// Impact ionization related parameters
  {"iba", "First impact ionization coefficient (1/m)", TR_DOUBLE, false},
  {"ibb", "Second impact ionization coefficient (V/m)", TR_DOUBLE, false},
  {"ibn", "Saturation voltage factor for impact ionization", TR_DOUBLE, false},
	// Intrinsic model temperature parameters
  {"tcv", "Threshold voltage temperature coefficient (V/K)", TR_DOUBLE, false},
  {"bex", "Mobility temperature exponent", TR_DOUBLE, false},
  {"ucex", "Longitudinal critical field temperature exponent", TR_DOUBLE, false},
  {"ibbt", "Temperature coefficient for IBB (1/K)", TR_DOUBLE, false},
	// Matching parameters
  {"avto", "Area related threshold voltage mismatch parameter (Vm)", TR_DOUBLE, false},
  {"akp", "Area related gain mismatch parameter (m)", TR_DOUBLE, false},
  {"agamma", "Area related body effect mismatch parameter (V^(1/2)m)", TR_DOUBLE, false},
	// Flicker noise parameters
  {"kf", "Flicker noise coefficient", TR_DOUBLE, false},
  {"af", "Flicker noise exponent", TR_DOUBLE, false},
	// Setup parameters
  {"nqs", "Non-Quasi-Static operation switch", TR_DOUBLE, false},
  {"satlim", "Ratio defining the saturation limit if/ir", TR_DOUBLE, false},
  {"xqc", "Charge/capacitance model selector", TR_DOUBLE, false},
  {"scale", "Scale parameter", TR_DOUBLE, false},
  {"tnom", "Nominal temperature of model parameters (K)", TR_DOUBLE, false},
  {"tmp", "Model simulation temperature (K)", TR_DOUBLE, false},
};

Mosnekv::Mosnekv(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0]  = &(type = NMOS);
  paramvalue[1]  = &(l = 1.0e-6);
  paramvalue[2]  = &(w = 1.0e-6);
  paramvalue[3]  = &(np = 1.0);
  paramvalue[4]  = &(ns = 1.0);
  paramvalue[5]  = &(cox = 0.0);
  paramvalue[6]  = &(xj = 1.0e-7);
  paramvalue[7]  = &(dw = 0.0);
  paramvalue[8]  = &(dl = 0.0);
  paramvalue[9]  = &(vto = 0.0);
  paramvalue[10]  = &(gamma = 0.0);
  paramvalue[11]  = &(phi = 0.0);
  paramvalue[12] = &(kp = 0.0);
  paramvalue[13] = &(eo = 0.0);
  paramvalue[14] = &(ucrit = 0.0);
  paramvalue[15] = &(tox = 0.0);
  paramvalue[16] = &(nsub = 0.0);
  paramvalue[17] = &(vfb = -2003.0);
  paramvalue[18] = &(uo = 0.0);
  paramvalue[19] = &(vmax = 0.0);
  paramvalue[20] = &(theta = -1.0);
  paramvalue[21] = &(lambda = 0.5);
  paramvalue[22] = &(weta = 0.25);
  paramvalue[23] = &(leta = 0.1);
  paramvalue[24] = &(qo = 0.0);
  paramvalue[25] = &(lk = 2.9e-7);
  paramvalue[26] = &(iba = 0.0);
  paramvalue[27] = &(ibb = 3.0e8);
  paramvalue[28] = &(ibn = 1.0);
  paramvalue[29] = &(tcv = 1.0e-3);
  paramvalue[30] = &(bex = -1.5);
  paramvalue[31] = &(ucex = 0.8);
  paramvalue[32] = &(ibbt = 9.0e-4);
  paramvalue[33] = &(avto = 0.0);
  paramvalue[34] = &(akp = 0.0);
  paramvalue[35] = &(agamma = 0.0);
  paramvalue[36] = &(kf = 0.0);
  paramvalue[37] = &(af = 1.0);
  paramvalue[38] = &(nqs = 0.0);
  paramvalue[39] = &(satlim = 54.59815003314424);
  paramvalue[40] = &(xqc = 0.4);
  paramvalue[41] = &(scale = 1.0);
  paramvalue[42] = &(tnom = 300.15);
  paramvalue[43] = &(tmp = 300.15);

  // Set the number of terminals
  setNumTerms(4);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(3);
}

void Mosnekv::init() throw(string&)
{
  // Set Constants
  epsilonox = scale*34.5e-12;
  epsilonsi = scale*104.5e-12;
  q=1.602e-19;
  k=1.3807e-23;
  Tref=300.15;
  vtT = (k*tmp)/q;	// Thermal voltage
  vtTnom = (k*tnom)/q;
  vtTref = (k*Tref)/q;
  egT = 1.16-0.000702*tmp*tmp/(tmp+1108.0);		// Energy gap
  egTnom = 1.16-0.000702*tnom*tnom/(tnom+1108.0);
  egTref = 1.16-0.000702*Tref*Tref/(Tref+1108.0);
  niT = 1.45e16*(tmp/Tref)*exp(egTref/(2.0*vtTref)-egT/(2.0*vtT));
  niTnom = 1.45e16*(tnom/Tref)*exp(egTref/(2.0*vtTref)-egTnom/(2.0*vtTnom));

  // For P-channel devices
  if(type == -1)
	{
    tcv = type * tcv;
    vto = type * vto;
  }

  // Calculate any missing parameters from user-defined settings
  // COX
  if(cox == 0.0)
	{
    if (tox > 0.0)
      cox = epsilonox / tox;
    else
      cox = 0.7e-3;
  }
  // GAMMA
  if(gamma == 0.0)
	{
    if (nsub > 0.0)
      gamma = sqrt(2.0 * epsilonsi * nsub * 1.0e6) / cox;
    else
      gamma = 1.0;
  }
  // PHI
  if(phi == 0.0)
	{
    if (nsub > 0.0)
      phi = 2.0 * vtTnom * log(nsub * 1.0e6 / niTnom);
    else
      phi = 0.7;
  }
  // VTO
  if(vto == 0.0)
	{
    if(vfb != -2003.0)
      vto = vfb + phi + gamma * sqrt(phi);
    else
      vto = 0.5;
  }
  // KP
  if(kp == 0.0)
	{
    if(uo > 0.0)
      kp = uo * cox * 1.0e-4 /*(m^2/cm^2)*/;
    else
      kp = 5.0e-5;
  }
  // UCRIT
  if(ucrit == 0.0)
	{
    if((vmax > 0) && (uo > 0))
      ucrit = vmax / (uo * 1.0e-4);
    else
      ucrit = 2.0e6;
  }
  // EO
  if(eo == 0.0)
	{
    if(theta >=0)
      eo = 0.0;
    else
      eo = 1.0e12;
  }

  // Intrinsic parameters temperature dependence
  vtoT=vto-tcv*(tmp-tnom);
  kpT=kp*pow(tmp/tnom,bex);
  ucritT=ucrit*pow(tmp/tnom,ucex);
  phiT=phi*tmp/tnom-3.0*vtT*log(tmp/tnom)-egTnom*tmp/tnom+egT;
  ibbT=ibb*(1.0+ibbt*(tmp-tnom));

  // create tape
  DenseIntVector var(3);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  initializeAD(var, var);
}

void Mosnekv::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: vdb
  // x[1]: vgb
  // x[2]: vsb
  // x[3]: dvdb/dt
  // x[4]: dvgb/dt
  // x[5]: dvsb/dt

  double weff = w + dw;
  double leff = l + dl;
  double vtoa = vtoT + avto / sqrt(np * weff * ns * leff);
  double kpa = kpT * (1 + akp / sqrt(np * weff * ns * leff));
  double gammaa = gamma + agamma / sqrt(np * weff * ns * leff);
  double cEps = 4.0 * pow(22.0e-3,2);
  double cA = 0.028;
  double xi  = cA * (10.0 * leff / lk -1.0);
  double deltavRSCE = 2.0 * qo / (cox * pow(1.0 + 0.5 * (xi + sqrt(xi*xi + cEps)), 2));
  double vc = ucritT * ns * leff;
  double lc = sqrt(epsilonsi * xj / cox);
  double lmin = ns * leff / 10.0;
  if(type == 1)
    eta = 0.5;
  else
    eta = 0.3333333333333;
  double qbo = gammaa * sqrt(phiT);
  double qox = 0.0;
  double C_ox = cox * np * weff * ns * leff;

  AD vgprime, vpo1, vpo, vsprime, vdprime, gammao, gammaprime, vp1, vp;
  AD n, i_f, vdss, vdssprime, deltav, vds, vip, deltal, lprime, leq;
  AD irprime, ir, betao, betaoprime, beta, is, ids, vib; //vpprime;
  AD idb1, idb, id, nq, xf, xr, qd, qs, qi, qb1, qb, qg, QI, QB, QD;
  AD QS, QG;
  AD qid, qig, qis;

  // Effective gate voltage in cluding reverse short channel effect
  vgprime =  type*x[1] - vtoa - deltavRSCE + phiT + gammaa * sqrt(phiT);

  // Effective substrate factor including charge-sharing for short and narrow
  // channels
  // Pinch-off voltage for narrow-channel effect
  vpo1 = vgprime - phiT - gammaa * (sqrt(vgprime + gammaa*gammaa / 4.0) - gammaa
  / 2.0);
  if (vgprime > 0.0)
    vpo = vpo1;
  else
    vpo = -phiT;

  // Effective substrate factor accounting for charge-sharing
  vsprime = 0.5 * (type*x[2] + phiT + sqrt(pow(type*x[2] + phiT,2) + 16.0 * vtT*vtT));
  vdprime = 0.5 * ( type*x[0] + phiT + sqrt(pow( type*x[0] + phiT,2) + 16.0 * vtT*vtT));

  // Pinch-off voltage including short- and narrow-channel effect
  gammao = gammaa - epsilonsi * (leta * (sqrt(vsprime) + sqrt(vdprime))
  / leff - 3.0 * weta * sqrt(vpo + phiT) / weff) / cox;
  gammaprime = 0.5 * (gammao + sqrt(gammao*gammao + 0.1 * vtT));
  vp1 = vgprime - phiT - gammaprime * (sqrt(vgprime+pow(gammaprime / 2.0,2)) -
  gammaprime / 2.0);
  if (vgprime > 0.0)
    vp = vp1;
  else
    vp = -phiT;

  // Slope factor
  n = 1.0 + gammaa / (2.0 * sqrt(vp + phiT + 4.0 * vtT));

  // Forward normalized current
  i_f = log(1.0 + exp((vp - type*x[2]*vtT) / (2.0*vtT)))*log(1.0 + exp((vp - type*x[2]) / (2.0*vtT)));

  // Velocity saturation voltage
  vdss = vc * (sqrt(0.25 + vtT * sqrt(i_f) / vc) - 0.5);

  // Drain-to-source saturation voltage for reverse normalized current
  vdssprime = vc * (sqrt(0.25 + vtT * (sqrt(i_f) - 0.75 * log(i_f))/vc) - 0.5) +
  vtT * (log(vc / (2.0 * vtT)) - 0.6);

  // Channel-length modulation
  deltav = 4.0 * vtT * sqrt(lambda * (sqrt(i_f) - vdss / vtT) + 1.0/64.0);
  vds = (type*x[0] - type*x[2]) / 2.0;
  vip = sqrt(vdss*vdss + deltav*deltav) - sqrt(pow(vds - vdss,2) + deltav*deltav);
  deltal = lambda * lc * log(1.0 + (vds - vip) / (lc * ucritT));

  // Equivalent channel length including channel-length moculation and velocity
  // saturation
  lprime = ns * leff - deltal + (vds + vip) / ucritT;
  leq = 0.5 * (lprime + sqrt(lprime*lprime + lmin*lmin));

  // Reverse normalized current
  irprime = log(1.0 + exp(((vp - vds - type*x[2] - sqrt(vdssprime*vdssprime +
  deltav*deltav) + sqrt((vds-vdssprime)*(vds-vdssprime) + deltav*deltav)) /
  vtT)/ 2.0))*log(1.0 + exp(((vp - vds - type*x[2] - sqrt(vdssprime*vdssprime +
  deltav*deltav) + sqrt((vds-vdssprime)*(vds-vdssprime) + deltav*deltav)) / vtT)/ 2.0));

  // Reverse normalized currect for mobility model, intrinsic
  //charges/capacitances, thermal noise model and NQS time-constant
  ir = log(1.0 + exp((vp -  type*x[0]) / (2.0*vtT)))*log(1.0 + exp((vp -  type*x[0]) / (2.0*vtT)));

  // Transconductance factor and mobility reduction due to vertical field
  betao = kpa * np * weff / leq;
  betaoprime = betao * (1.0 + cox * qbo / (eo * epsilonsi));

  // Quasi-static model equations
  // Dynamic model for the intrinsic node charges
  nq = 1.0 + gammaa / (2.0 * sqrt(vp + phiT + 1.0e-6));

  // Normalized intrinsic node charges
  xf = sqrt(0.25 + i_f);
  xr = sqrt(0.25 + ir);
  qd = -nq * (4.0 * (3.0 * xr*xr*xr + 6.0 * xr*xr * xf + 4.0 * xr * xf*xf + 2.0 * xf*xf*xf) / (15.0 *
  pow(xf + xr,2)) - 0.5);
  qs = -nq * (4.0 * (3.0 * xf*xf*xf + 6.0 * xf*xf * xr + 4.0 * xf * xr*xr + 2.0 * xr*xr*xr) / (15.0 *
  pow(xf + xr,2)) - 0.5);
  qi = qs + qd;
  qb1 = -gammaa * sqrt(vp + phiT + 1.0e-6) / vtT - (nq - 1.0) * qi / nq;
  if (vgprime > 0.0)
    qb = qb1;
  else
    qb = -vgprime / vtT;
  qg = -qi - qox - qb;

  // Rigorous mobility reduction model
  AD qtemp = qb + eta*qi;
  if (qtemp < 0.0)
    qtemp = -qtemp;
  beta = betaoprime / (1.0 + cox * vtT * qtemp / (eo * epsilonsi));

  // Specific current
  is = 2.0 * n * beta * vtT*vtT;

  // Drain-to-source current
  ids = type*is * (i_f - irprime);

  // Impact ionization current
  vib =  type*x[0] - type*x[2] - 2.0 * ibn * vdss;
  idb1 = ids * iba * vib * exp(-ibbT * lc / vib) / ibbT;
  if (vib > 0.0)
    idb = idb1;
  else
    idb = 0.0;
  id = ids + idb;

  // Total charges
  QI = C_ox * vtT * qi;
  QB = C_ox * vtT * qb;
  QD = C_ox * vtT * qd;
  QS = C_ox * vtT * qs;
  QG = C_ox * vtT * qg;

  // The capacitors considered in this model are as shown in
  // the equivalent circuit in the EKV model documentation.
  // In general, Cxy = +/-dQx/dVy * dVy/dt, where the "+" sign is
  // used if x == y, else "-" sign is used

  // Cds contribution
  qid = - QD.fastAccessDx(2) * x[5];

  // Cgd + Cgs contribution
  qig = -QG.fastAccessDx(0) * x[3] - QG.fastAccessDx(2) * x[5];

  // Csd + Csg contribution
  qis = -QS.fastAccessDx(0) * x[3] - QS.fastAccessDx(1) * x[4];

  // Assign DC currents
  flow[0] = id + qid;    //DC Drain current
  flow[1] = qig;   //DC Gate current
  flow[2] = -id + qis;   //DC Source current

  // Assign known output voltages
  effort[0] = x[0];  // vdb
  effort[1] = x[1];	 // vgb
  effort[2] = x[2];	 // vsb
}
