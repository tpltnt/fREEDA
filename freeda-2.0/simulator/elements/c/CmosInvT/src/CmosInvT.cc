#include "../../../../analysis/TimeDomainSV.h"
#include "CmosInvT.h"

// Define the number of parameters
const unsigned CmosInvT :: n_par = 72;

// Element Information
ItemInfo CmosInvT::einfo =
{
  "cmosinvt",
  "CMOS Inverter-Thermal",
  "Shivam Priyadarshi",
  DEFAULT_ADDRESS"category:thermal-cmos,digital"
};

// Parameter Information
ParmInfo CmosInvT::pinfo[] =
{       // MOSFET Device Input Variables
  {"ln","Channel Length of NMOS (m)", TR_DOUBLE, false},
  {"wn","Channel Width of NMOS (m)", TR_DOUBLE, false},
  {"lp","Channel Length of PMOS (m)", TR_DOUBLE, false},
  {"wp","Channel Width of PMOS (m)", TR_DOUBLE, false},
  {"np_n", "Parallel multiple device number of NMOS ", TR_DOUBLE, false},
  {"np_p", "Parallel multiple device number of PMOS ", TR_DOUBLE, false},
  {"ns_n", "Serial multiple device number of NMOS", TR_DOUBLE, false},
  {"ns_p", "Serial multiple device number of PMOS", TR_DOUBLE, false},
  // Process related parameters
  {"cox_n", "Gate oxide capacitance per area of NMOS (F/m^2)", TR_DOUBLE, false},
  {"cox_p", "Gate oxide capacitance per area of PMOS (F/m^2)", TR_DOUBLE, false},
  {"xj_n", "Junction depth of NMOS (m)", TR_DOUBLE, false},
  {"xj_p", "Junction depth of PMOS (m)", TR_DOUBLE, false},
  {"dw_n", "Channel width correction of NMOS (m)", TR_DOUBLE, false},
  {"dw_p", "Channel width correction of PMOS (m)", TR_DOUBLE, false},
  {"dl_n", "Channel length correction of NMOS (m)", TR_DOUBLE, false},
  {"dl_p", "Channel length correction of PMOS (m)", TR_DOUBLE, false},
  // Basic intrinsic model parameters
  {"vto_n", "Long_channel threshold voltage of NMOS (V)", TR_DOUBLE, false},
  {"vto_p", "Long_channel threshold voltage of PMOS (V)", TR_DOUBLE, false},
  {"gamma_n", "Body effect parameter of NMOS (V^1/2)",TR_DOUBLE, false},
  {"gamma_p", "Body effect parameter of PMOS (V^1/2)",TR_DOUBLE, false},
  {"phi_n", "Bulk Fermi potential of NMOS (V)", TR_DOUBLE, false},
  {"phi_p", "Bulk Fermi potential of PMOS (V)", TR_DOUBLE, false},
  {"kp_n", "Transconductance parameter of NMOS (A/V^2)", TR_DOUBLE, false},
  {"kp_p", "Transconductance parameter of PMOS (A/V^2)", TR_DOUBLE, false},
  {"eo_n", "Mobility reduction coefficient of NMOS (V/m)",TR_DOUBLE, false},
  {"eo_p", "Mobility reduction coefficient of PMOS (V/m)",TR_DOUBLE, false},
  {"ucrit_n", "Longitudinal critical field of NMOS (V/m)", TR_DOUBLE, false},
  {"ucrit_p", "Longitudinal critical field of PMOS (V/m)", TR_DOUBLE, false},
  // Operational parameters
  {"tox_n", "Oxide thickness of NMOS (m)", TR_DOUBLE, false},
  {"tox_p", "Oxide thickness of PMOS (m)", TR_DOUBLE, false},
  {"nsub", "Channel doping of NMOS (1/cm^3)", TR_DOUBLE, false},
  {"psub", "Channel doping of PMOS (1/cm^3)", TR_DOUBLE, false},
  {"vfb_n", "Flat-band voltage of NMOS (V)", TR_DOUBLE, false},
  {"vfb_p", "Flat-band voltage of PMOS (V)", TR_DOUBLE, false},
  {"uo_n", "Low-field mobility of NMOS(cm^2/(V.s))", TR_DOUBLE, false},
  {"uo_p", "Low-field mobility of PMOS(cm^2/(V.s))", TR_DOUBLE, false},
  {"vmax_n", "Saturation velocity of NMOS (m/s)", TR_DOUBLE, false},
  {"vmax_p", "Saturation velocity of PMOS (m/s)", TR_DOUBLE, false},
  // Channel length modulation and charge sharing parameters
  {"lambda_n", "Channel-length modulation of NMOS", TR_DOUBLE, false},
  {"lambda_p", "Channel-length modulation of PMOS", TR_DOUBLE, false},
  {"weta_n", "Narrow-channel effect coefficient of NMOS", TR_DOUBLE, false},
  {"weta_p", "Narrow-channel effect coefficient of PMOS", TR_DOUBLE, false},
  {"leta_n", "Short-channel effect coefficient of NMOS", TR_DOUBLE, false},
  {"leta_p", "Short-channel effect coefficient of PMOS", TR_DOUBLE, false},
  // Reverse short-channel effect parameters
  {"qo_n", "Reverse short channel effect peak charge density of NMOS (A.s/m^2)", TR_DOUBLE, false},
  {"qo_p", "Reverse short channel effect peak charge density of PMOS (A.s/m^2)", TR_DOUBLE, false},
  {"lk_n", "Reverse short channel effect characteristic length of NMOS (m)", TR_DOUBLE, false},
  {"lk_p", "Reverse short channel effect characteristic length of PMOS (m)", TR_DOUBLE, false},
  // Impact ionization related parameters
  {"iba_n", "First impact ionization coefficient of NMOS (1/m)", TR_DOUBLE, false},
  {"iba_p", "First impact ionization coefficient of PMOS (1/m)", TR_DOUBLE, false},
  {"ibb_n", "Second impact ionization coefficient of NMOS (V/m)", TR_DOUBLE, false},
  {"ibb_p", "Second impact ionization coefficient of PMOS (V/m)", TR_DOUBLE, false},
  {"ibn_n", "Saturation voltage factor for impact ionization of NMOS", TR_DOUBLE, false},
  {"ibn_p", "Saturation voltage factor for impact ionization of PMOS", TR_DOUBLE, false},
  // Intrinsic model temperature parameters
  {"tcv_n", "Threshold voltage temperature coefficient of NMOS (V/K)", TR_DOUBLE, false},
  {"tcv_p", "Threshold voltage temperature coefficient of PMOS(V/K)", TR_DOUBLE, false},
  {"bex_n", "Mobility temperature exponent of NMOS", TR_DOUBLE, false},
  {"bex_p", "Mobility temperature exponent of PMOS", TR_DOUBLE, false},
  {"ucex_n", "Longitudinal critical field temperature exponent of NMOS", TR_DOUBLE, false},
  {"ucex_p", "Longitudinal critical field temperature exponent of PMOS", TR_DOUBLE, false},
  {"ibbt_n", "Temperature coefficient for IBB of NMOS (1/K)", TR_DOUBLE, false},
  {"ibbt_p", "Temperature coefficient for IBB of PMOS (1/K)", TR_DOUBLE, false},
  // Matching parameters
  {"avto_n", "Area related threshold voltage mismatch parameter of NMOS (Vm)", TR_DOUBLE, false},
  {"avto_p", "Area related threshold voltage mismatch parameter of PMOS (Vm)", TR_DOUBLE, false},
  {"akp_n", "Area related gain mismatch parameter of NMOS (m)", TR_DOUBLE, false},
  {"akp_p", "Area related gain mismatch parameter of PMOS (m)", TR_DOUBLE, false},
  {"agamma_n", "Area related body effect mismatch parameter of NMOS (V^(1/2)m)", TR_DOUBLE, false},
  {"agamma_p", "Area related body effect mismatch parameter of PMOS (V^(1/2)m)", TR_DOUBLE, false},
  // Setup Parameters
  {"scale", "Scale parameter", TR_DOUBLE, false},
  {"tnom", "Nominal temperature of model parameters (K)", TR_DOUBLE, false},
  {"tmp", "Model simulation temperature (K)", TR_DOUBLE, false},
  {"rth", "Thermal Resistance(ohm)", TR_DOUBLE, false},
};

// Constructor
  CmosInvT::CmosInvT(const string& iname)
:ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(ln = 1.0e-6);
  paramvalue[1] = &(wn = 1.0e-6);
  paramvalue[2] = &(lp = 1.0e-6);
  paramvalue[3] = &(wp = 1.0e-6);
  paramvalue[4] = &(np_n = 1.0);
  paramvalue[5] = &(np_p = 1.0);
  paramvalue[6] = &(ns_n = 1.0);
  paramvalue[7] = &(ns_p = 1.0);
  paramvalue[8] = &(cox_n = 0.0);
  paramvalue[9] = &(cox_p = 0.0);
  paramvalue[10] = &(xj_n = 1.0e-7);
  paramvalue[11] = &(xj_p = 1.0e-7);
  paramvalue[12] = &(dw_n = 0.0);
  paramvalue[13]= &(dw_p = 0.0);
  paramvalue[14]= &(dl_n = 0.0);
  paramvalue[15]= &(dl_p = 0.0);
  paramvalue[16]= &(vto_n = 0.0);
  paramvalue[17]= &(vto_p = 0.0);
  paramvalue[18]= &(gamma_n = 0.0);
  paramvalue[19]= &(gamma_p = 0.0);
  paramvalue[20]= &(phi_n = 0.0);
  paramvalue[21]= &(phi_p = 0.0);
  paramvalue[22]= &(kp_n = 0.0);
  paramvalue[23]= &(kp_p = 0.0);
  paramvalue[24] = &(eo_n = 0.0);
  paramvalue[25] = &(eo_p = 0.0);
  paramvalue[26] = &(ucrit_n = 0.0);
  paramvalue[27] = &(ucrit_p = 0.0);
  paramvalue[28] = &(tox_n = 0.0);
  paramvalue[29] = &(tox_p = 0.0);
  paramvalue[30] = &(nsub = 0.0);
  paramvalue[31] = &(psub = 0.0);
  paramvalue[32] = &(vfb_n = -2003.0);
  paramvalue[33] = &(vfb_p = 2003.0);
  paramvalue[34] = &(uo_n = 500);
  paramvalue[35] = &(uo_p = 200);
  paramvalue[36] = &(vmax_n = 0.0);
  paramvalue[37] = &(vmax_p = 0.0);
  paramvalue[38] = &(lambda_n = 0.5);
  paramvalue[39] = &(lambda_p = 0.5);
  paramvalue[40] = &(weta_n = 0.25);
  paramvalue[41] = &(weta_p = 0.25);
  paramvalue[42] = &(leta_n = 0.1);
  paramvalue[43] = &(leta_p = 0.1);
  paramvalue[44] = &(qo_n = 0.0);
  paramvalue[45] = &(qo_p = 0.0);
  paramvalue[46] = &(lk_n = 2.9e-7);
  paramvalue[47] = &(lk_p = 2.9e-7);
  paramvalue[48] = &(iba_n = 0.0);
  paramvalue[49] = &(iba_p = 0.0);
  paramvalue[50] = &(ibb_n = 3.0e8);
  paramvalue[51] = &(ibb_p = 3.0e8);
  paramvalue[52] = &(ibn_n = 1.0);
  paramvalue[53] = &(ibn_p = 1.0);
  paramvalue[54] = &(tcv_n = 1.0e-3);
  paramvalue[55] = &(tcv_p = 1.0e-3);
  paramvalue[56] = &(bex_n = -1.5);
  paramvalue[57] = &(bex_p = -1.5);
  paramvalue[58] = &(ucex_n = 0.8);
  paramvalue[59] = &(ucex_p = 0.8);
  paramvalue[60] = &(ibbt_n = 9.0e-4);
  paramvalue[61] = &(ibbt_p = 9.0e-4);
  paramvalue[62] = &(avto_n = 0.0);
  paramvalue[63] = &(avto_p = 0.0);
  paramvalue[64] = &(akp_n = 0.0);
  paramvalue[65] = &(akp_p = 0.0);
  paramvalue[66] = &(agamma_n = 0.0);
  paramvalue[67] = &(agamma_p = 0.0);
  paramvalue[68] = &(scale = 1.0);
  paramvalue[69] = &(tnom = 300.15);
  paramvalue[70] = &(tmp = 300.15);
  paramvalue[71] = &(rth = 1);

  // Set Number of Terminals
  setNumTerms(6);   // two extra terminals from thermal resistor

  // Set Flags
  setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);

  // Set Number of States
  setNumberOfStates(4);  // Temperature is also a state
}

// Initialization Function
void CmosInvT::init() throw(string&)
{
  // Set Constants
  scale = 1.0;
  epsilonox = scale*34.5e-12;
  epsilonsi = scale*104.5e-12;
  q=1.602e-19;
  k=1.3807e-23;
  Tref=300.15;

  // Thermal Parameter Calculation
  vtTnom = (k*tnom)/q;
  vtTref = (k*Tref)/q;
  egTnom = 1.16-0.000702*tnom*tnom/(tnom+1108.0);
  egTref = 1.16-0.000702*Tref*Tref/(Tref+1108.0);
  niTnom = 1.45e16*(tnom/Tref)*exp(egTref/(2.0*vtTref)-egTnom/(2.0*vtTnom));

  // Distinguishing between pmos and nmos
  type_n = 1;
  type_p = -1;

  // Calculate any missing parameters from user-defined settings
  // COX_N
  if(cox_n == 0.0)
  {
    if (tox_n > 0.0)
      cox_n = epsilonox / tox_n;
    else
      cox_n = 0.7e-3;
  }

  // COX_P
  if(cox_p == 0.0)
  {
    if (tox_p > 0.0)
      cox_p = epsilonox / tox_p;
    else
      cox_p = 0.7e-3;
  }

  //  cout << "\nCOX ="  << cox  << "\n";

  // GAMMA_N
  if(gamma_n == 0.0)
  {
    if (nsub > 0.0)
      gamma_n = sqrt(2.0 * q * epsilonsi * nsub * 1.0e6) / cox_n;
    else
      gamma_n = 1.0;
  }

  // cout << "\n GAMMAN=" << gamma_n << "\n";

  // GAMMA_P
  if(gamma_p == 0.0)
  {
    if (psub > 0.0)
      gamma_p = sqrt(2.0 * q * epsilonsi * psub * 1.0e6) / cox_p;
    else
      gamma_p = 1.0;
  }
  // cout << "\n GAMMAP=" << gamma_p << "\n";

  // PHI_N
  if(phi_n == 0.0)
  {
    if (nsub > 0.0)
      phi_n = 2.0 * vtTnom * log(nsub * 1.0e6 / niTnom);
    else
      phi_n = 0.7;
  }

  // PHI_P
  if(phi_p == 0.0)
  {
    if (psub > 0.0)
      phi_p = 2.0 * vtTnom * log(psub * 1.0e6 / niTnom);
    else
      phi_p = 0.7;
  }

  //VTO_N
  if(vto_n == 0.0)
  {
    if(vfb_n != -2003.0)
      vto_n = -1*vfb_n + phi_n + gamma_n * sqrt(phi_n);  // VFB negative for NMOS
    else
      vto_n = 0.5*type_n;
  }
  else
    vto_n = type_n*vto_n;


  //VTO_P
  if(vto_p == 0.0)
  {
    if(vfb_p != 2003.0)
      vto_p =  1*vfb_p + phi_p + gamma_p * sqrt(phi_p);  // VFB positive for PMOS
    else
      vto_p = -1*0.5*type_p;
  }
  else
    vto_p = type_p*vto_p;

  // KP_n
  if(kp_n == 0.0)
  {
    if(uo_n > 0.0)
      kp_n = uo_n * cox_n * 1.0e-4 /*(m^2/cm^2)*/;
    else
      kp_n = 5.0e-5;
  }

  // KP_p
  if(kp_p == 0.0)
  {
    if(uo_p > 0.0)
      kp_p = uo_p * cox_p * 1.0e-4 /*(m^2/cm^2)*/;
    else
      kp_p = 5.0e-5;
  }

  // UCRIT_N
  if(ucrit_n == 0.0)
  {
    if((vmax_n > 0) && (uo_n > 0))
      ucrit_n = vmax_n / (uo_n * 1.0e-4);
    else
      ucrit_n = 2.0e6;
  }
  // UCRIT_P
  if(ucrit_p == 0.0)
  {
    if((vmax_p > 0) && (uo_p > 0))
      ucrit_p = vmax_p / (uo_p * 1.0e-4);
    else
      ucrit_p = 2.0e6;
  }

  tcv_n = type_n*tcv_n;
  tcv_p = type_p*tcv_p;
  if(eo_n==0) {eo_n = 1.0e12;}
  if(eo_p==0) {eo_p = 1.0e12;}

  //create tape
  DenseIntVector var(4);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  var[3] = 3;   //Temperature
  initializeAD(var,var);
}

void CmosInvT::getLocalRefIdx(UnsignedVector& local_ref_vec,TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));  //Vdd
  term_list.push_back(getTerminal(1));  //Vin
  term_list.push_back(getTerminal(2));  //Vout
  term_list.push_back(getTerminal(3));  // Local reference terminal (Gnd)
  term_list.push_back(getTerminal(4));
  term_list.push_back(getTerminal(5)); // Local reference terminal
  local_ref_vec.push_back(3); // Local reference index
  local_ref_vec.push_back(5); // Local reference index
}

// Evaluate Function
void CmosInvT::eval(AD * x, AD * effort, AD * flow)
{

  //EKV Pmos and nmos models are used for building inverter which has three state variable :
  //Vdb-x[0],Vgb-x[1],Vsb-x[2]
  //For Simplicity it is assumed : Bulk & Source of Nmos = 0, Bulk & Source of Pmos = Vdd
  //Pmos and nmos state varibles are found in terms of inverter state variable:
  //x[0]  : Vdd    x[4] : d(Vdd)/dt   Vdb of Nmos = x[2]  Vdb of Pmos = x[2]-x[0]
  //x[1]  : Vin    x[5] : d(Vin)/dt   Vgb of Nmos = x[1]  Vgb of Pmos = x[1]-x[0]
  //x[2]  : Vout   x[6] : d(Vout)/dt  Vsb of Nmos = 0     Vsb  of Pmos = 0
  //x[3]  : Temperature_Rise = Temperature - tmp  : With respect to SINK temperature
  //x[7]  : d(Temperature_Rise)/dt
  //z1[0] : Current in Vdd
  //z1[1] : Current in Vin
  //z1[2] : Current in Vout
  //z1[3] : Heat Current
  //z1[4] : x[0] : Vdd
  //z1[5] : x[1] : Vin
  //z1[6] : x[2] : Vout
  //z1[7] : x[3] + tmp : Temperature = Temperature_Rise + tmp

  effort[3] = x[3] + tmp;  // Temperature

  AD vtoT_n,vtoT_p, kpT_n,kpT_p,ucritT_n,ucritT_p, phiT_n,phiT_p,ibbT_n,ibbT_p;
  AD vtT,egT,niT;


  // Thermal Parameter Calculation  : Calculated as function of Temp
  vtT = (k*effort[3])/q;    // Thermal voltage
  egT = 1.16-0.000702*effort[3]*effort[3]/(effort[3]+1108.0);               // Energy gap
  niT = 1.45e16*(effort[3]/Tref)*exp(egTref/(2.0*vtTref)-egT/(2.0*vtT)); // Temperature Dependency of Intrinsic Carrier Concentration

  // Intrinsic parameters temperature dependence
  vtoT_n=vto_n-tcv_n*(effort[3]-tnom);
  vtoT_p=vto_p-tcv_p*(effort[3]-tnom);

  kpT_n=kp_n*pow(effort[3]/tnom,bex_n);
  kpT_p=kp_p*pow(effort[3]/tnom,bex_p);
  ucritT_n=ucrit_n*pow(effort[3]/tnom,ucex_n);
  ucritT_p=ucrit_p*pow(effort[3]/tnom,ucex_p);
  phiT_n=phi_n*effort[3]/tnom-3.0*vtT*log(effort[3]/tnom)-egTnom*effort[3]/tnom+egT;
  phiT_p=phi_p*effort[3]/tnom-3.0*vtT*log(effort[3]/tnom)-egTnom*effort[3]/tnom+egT;
  ibbT_n=ibb_n*(1+ibbt_n*(effort[3]-tnom));
  ibbT_p=ibb_p*(1+ibbt_p*(effort[3]-tnom));



  double weff_n = wn + dw_n;
  double leff_n = ln + dl_n;
  double weff_p = wp + dw_p;
  double leff_p = lp + dl_p;


  AD vtoa_n = vtoT_n + avto_n / sqrt(np_n * weff_n * ns_n * leff_n);
  AD vtoa_p = vtoT_p + avto_p / sqrt(np_p * weff_p * ns_p * leff_p);
  AD kpa_n = kpT_n * (1 + akp_n / sqrt(np_n * weff_n * ns_n * leff_n));
  AD kpa_p = kpT_p * (1 + akp_p / sqrt(np_p * weff_p * ns_p * leff_p));
  double gammaa_n = gamma_n + agamma_n / sqrt(np_n * weff_n * ns_n * leff_n);
  double gammaa_p = gamma_p + agamma_p / sqrt(np_p * weff_p * ns_p * leff_p);
  double cEps = 4.0 * pow(22.0e-3,2);
  double cA = 0.028;
  double xi_n  = cA * (10.0 * leff_n / lk_n -1.0);
  double xi_p  = cA * (10.0 * leff_p / lk_p -1.0);
  double deltavRSCE_n = 2 * qo_n / (cox_n * pow(1 + 0.5 * (xi_n + sqrt(xi_n*xi_n + cEps)),2));
  double deltavRSCE_p = 2 * qo_p / (cox_p * pow(1 + 0.5 * (xi_p + sqrt(xi_p*xi_p + cEps)),2));
  AD vc_n = ucritT_n * ns_n * leff_n;
  AD vc_p = ucritT_p * ns_p * leff_p;
  double lc_n = sqrt(epsilonsi * xj_n / cox_n);
  double lc_p = sqrt(epsilonsi * xj_p / cox_p);
  double lmin_n = ns_n * leff_n / 10;
  double lmin_p = ns_p * leff_p / 10;
  AD qbo_n = gammaa_n * sqrt(phiT_n);
  AD qbo_p = gammaa_p * sqrt(phiT_p);
  double qox = 0.0;
  double C_ox_n = cox_n * np_n * weff_n * ns_n * leff_n;
  double C_ox_p = cox_p * np_p * weff_p * ns_p * leff_p;

  eta_n= 0.5;
  eta_p = 0.3333333333333;

  //All the active variables for ADOL-C "ADs" must be initialiazed here
  AD vgprime_n,vgprime_p, vpo1_n, vpo_n,vpo1_p,vpo_p, vsprime_n, vdprime_n,vsprime_p, vdprime_p, gammao_n, gammaprime_n, vp1_n, vp_n,gammao_p, gammaprime_p, vp1_p, vp_p;
  AD n_n,n_p,i_f_n,i_f_p,vdss_n,vdssprime_n,vdss_p,vdssprime_p, deltav_n, vds_n, vip_n, deltal_n,deltav_p, vds_p, vip_p, deltal_p,lprime_n,leq_n,lprime_p,leq_p;
  AD irprime_n,irprime_p, ir_n,ir_p, betao_n, betaoprime_n, beta_n,betao_p, betaoprime_p, beta_p, is_n, ids_n, vib_n,is_p, ids_p, vib_p; //vpprime;
  AD idb1_n, idb_n, id_n,idb1_p, idb_p, id_p, nq_n,nq_p,xf_n, xr_n, qd_n, qs_n, qi_n, qb1_n, qb_n, qg_n, QI_n, QB_n, QD_n,xf_p,xr_p, qd_p, qs_p, qi_p, qb1_p, qb_p, qg_p, QI_p, QB_p, QD_p;
  AD QS_n, QG_n,QS_p,QG_p;

  // Effective gate voltage including reverse short channel effect

  vgprime_n =  type_n*x[1] - vtoa_n - deltavRSCE_n + phiT_n + gammaa_n * sqrt(phiT_n); // Vgb of nmos : x[1]
  vgprime_p =  type_p*(x[1]-x[0]) - vtoa_p - deltavRSCE_p + phiT_p + gammaa_p * sqrt(phiT_p); // Vgb of pmos : x[1]-x[0]

  // Effective substrate factor including charge-sharing for short and narrow
  // channels
  // Pinch-off voltage for narrow-channel effect
  vpo1_n = vgprime_n - phiT_n - gammaa_n * (sqrt(vgprime_n + gammaa_n*gammaa_n / 4.0) - gammaa_n/ 2.0);
  //condassign(vpo_n, vgprime_n, vpo1_n, -phiT_n);
  if (vgprime_n > 0)
    vpo_n = vpo1_n;
  else
    vpo_n = -phiT_n;

  vpo1_p = vgprime_p - phiT_p - gammaa_p * (sqrt(vgprime_p + gammaa_p*gammaa_p / 4) - gammaa_p/ 2);
  //  condassign(vpo_p, vgprime_p, vpo1_p, -phiT_p);
  if (vgprime_p > zero)
    vpo_p = vpo1_p;
  else
    vpo_p = -phiT_p;



  // Effective substrate factor accounting for charge-sharing
  vsprime_n = 0.5 * (type_n*0 + phiT_n + sqrt(pow(type_n*0 + phiT_n,2) + 16.0 * vtT*vtT)); // vsb for nmos =0
  vdprime_n = 0.5 * ( type_n*x[2] + phiT_n + sqrt(pow( type_n*x[2] + phiT_n,2) + 16.0 * vtT*vtT)); // Vdb for nmos : x[2]

  vsprime_p = 0.5 * (type_p*0 + phiT_p + sqrt(pow(type_p*0 + phiT_p,2) + 16.0 * vtT*vtT)); // vsb for pmos =0
  vdprime_p = 0.5 * ( type_p*(x[2]-x[0]) + phiT_p + sqrt(pow( type_p*(x[2]-x[0]) + phiT_p,2) + 16.0 * vtT*vtT)); // Vdb for pmos : x[2]-x[0]

  // Pinch-off voltage including short- and narrow-channel effect
  gammao_n = gammaa_n - epsilonsi * (leta_n * (sqrt(vsprime_n) + sqrt(vdprime_n))/ leff_n - 3.0 * weta_n * sqrt(vpo_n + phiT_n) / weff_n) / cox_n;
  gammaprime_n = 0.5 * (gammao_n + sqrt(gammao_n*gammao_n + 0.1 * vtT));
  vp1_n = vgprime_n - phiT_n - gammaprime_n * (sqrt(vgprime_n+pow(gammaprime_n / 2.0,2)) -gammaprime_n / 2.0);
  //  condassign(vp_n, vgprime_n, vp1_n, -phiT_n);
  if (vgprime_n > zero)
    vp_n = vp1_n;
  else
    vp_n = -phiT_n;

  gammao_p = gammaa_p - epsilonsi * (leta_p * (sqrt(vsprime_p) + sqrt(vdprime_p))/ leff_p - 3.0 * weta_p * sqrt(vpo_p + phiT_p) / weff_p) / cox_p;
  gammaprime_p = 0.5 * (gammao_p + sqrt(gammao_p*gammao_p + 0.1 * vtT));
  vp1_p = vgprime_p - phiT_p - gammaprime_p * (sqrt(vgprime_p+pow(gammaprime_p / 2.0,2)) -gammaprime_p / 2.0);
  // condassign(vp_p, vgprime_p, vp1_p, -phiT_p);
  if (vgprime_p > zero)
    vp_p = vp1_p;
  else
    vp_p = -phiT_p;


  // Slope factor
  n_n = 1 + gammaa_n / (2.0 * sqrt(vp_n + phiT_n + 4.0 * vtT));
  n_p = 1 + gammaa_p / (2.0 * sqrt(vp_p + phiT_p + 4.0 * vtT));

  // Forward normalized current
  i_f_n = log(1.0 + exp((vp_n - type_n*0) / (2.0*vtT)))*log(1.0 + exp((vp_n - type_n*0) / (2.0*vtT))); //Vsb of nmos =0

  i_f_p = log(1 + exp((vp_p - type_p*0) / (2 * vtT)))*log(1 + exp((vp_p - type_p*0) / (2 * vtT))); // Vsb for pmos =0


  // Velocity saturation voltage
  vdss_n = vc_n * (sqrt(0.25 + vtT * sqrt(i_f_n) / vc_n) - 0.5);

  vdss_p = vc_p * (sqrt(0.25 + vtT * sqrt(i_f_p) / vc_p) - 0.5);

  // Drain-to-source saturation voltage for reverse normalized current
  vdssprime_n = vc_n * (sqrt(0.25 + vtT * (sqrt(i_f_n) - 0.75 * log(i_f_n))/vc_n) - 0.5) +vtT * (log(vc_n / (2.0 * vtT)) - 0.6);

  vdssprime_p = vc_p * (sqrt(0.25 + vtT * (sqrt(i_f_p) - 0.75 * log(i_f_p))/vc_p) - 0.5) +vtT * (log(vc_p / (2 * vtT)) - 0.6);

  // Channel-length modulation
  deltav_n = 4.0 * vtT * sqrt(lambda_n * (sqrt(i_f_n) - vdss_n / vtT) + 1.0/64.0);
  vds_n = (type_n*x[2] - type_n*0) / 2.0; // Vdb = x[2] and Vsb= 0 for nmos
  vip_n = sqrt(vdss_n*vdss_n + deltav_n*deltav_n) - sqrt(pow(vds_n - vdss_n,2) + deltav_n*deltav_n);
  deltal_n = lambda_n * lc_n * log(1.0 + (vds_n - vip_n) / (lc_n * ucritT_n));

  deltav_p = 4 * vtT * sqrt(lambda_p * (sqrt(i_f_p) - vdss_p / vtT) + 1 /64);
  vds_p = ( type_p*(x[2]-x[0]) - type_p*0) / 2; //Vdb =x[2]-x[0] & Vsb=0 pmos
  vip_p = sqrt(vdss_p*vdss_p + deltav_p*deltav_p) - sqrt(pow(vds_p - vdss_p,2) + deltav_p*deltav_p);
  deltal_p = lambda_p * lc_p * log(1 + (vds_p - vip_p) / (lc_p * ucritT_p));

  //cout << "VDS_n:" << vds_n << endl;
  //cout << "VDS_p:" << vds_p << endl;


  // Equivalent channel length including channel-length moculation and velocity
  // saturation
  lprime_n = ns_n * leff_n - deltal_n + (vds_n + vip_n) / ucritT_n;
  leq_n = 0.5 * (lprime_n + sqrt(lprime_n*lprime_n + lmin_n*lmin_n));

  lprime_p = ns_p * leff_p - deltal_p + (vds_p + vip_p) / ucritT_p;
  leq_p = 0.5 * (lprime_p + sqrt(lprime_p*lprime_p + lmin_p*lmin_p));

  // Reverse normalized current
  irprime_n = log(1.0 + exp(((vp_n - vds_n - type_n*0 - sqrt(vdssprime_n*vdssprime_n +deltav_n*deltav_n) + sqrt((vds_n-vdssprime_n)*(vds_n-vdssprime_n) + deltav_n*deltav_n)) /vtT)/ 2.0))*log(1.0 + exp(((vp_n - vds_n - type_n*0 - sqrt(vdssprime_n*vdssprime_n +deltav_n*deltav_n) + sqrt((vds_n-vdssprime_n)*(vds_n-vdssprime_n) + deltav_n*deltav_n)) / vtT)/ 2.0));

  irprime_p = log(1 + exp(((vp_p - vds_p - type_p*0 - sqrt(vdssprime_p*vdssprime_p +deltav_p*deltav_p) + sqrt((vds_p-vdssprime_p)*(vds_p-vdssprime_p) + deltav_p*deltav_p)) /vtT)/ 2))*log(1 + exp(((vp_p - vds_p - type_p*0 - sqrt(vdssprime_p*vdssprime_p +deltav_p*deltav_p) + sqrt((vds_p-vdssprime_p)*(vds_p-vdssprime_p) + deltav_p*deltav_p)) / vtT)/ 2));

  // Reverse normalized currect for mobility model, intrinsic
  //charges/capacitances, thermal noise model and NQS time-constant
  ir_n = log(1.0 + exp((vp_n -  type_n*x[2]) / (2.0*vtT)))*log(1.0 + exp((vp_n -  type_n*x[2]) / (2.0*vtT))); //Vdb of nmos x[2]

  ir_p = log(1 + exp((vp_p -  type_p*(x[2]-x[0])) / (2 * vtT)))*log(1 + exp((vp_p -  type_p*(x[2]-x[0])) / (2 * vtT))); //Vdb of pmos x[2]-x[0]

  // Transconductance factor and mobility reduction due to vertical field
  betao_n = kpa_n * np_n * weff_n / leq_n;
  betaoprime_n = betao_n * (1.0 + cox_n * qbo_n / (eo_n * epsilonsi));

  betao_p = kpa_p * np_p * weff_p / leq_p;
  betaoprime_p = betao_p * (1 + cox_p * qbo_p / (eo_p * epsilonsi));

  //cout << "\nBetaPrimeP ="  << betaoprime_p  << "\n";
  //cout << "\nBetaPrimeN ="  << betaoprime_n  << "\n";

  // Dynamic model for the intrinsic node charges
  nq_n = 1.0 + gammaa_n / (2.0 * sqrt(vp_n + phiT_n +1.0e-6));

  nq_p = 1 + gammaa_p / (2 * sqrt(vp_p + phiT_p +1e-6));

  // Normalized intrinsic node charges
  xf_n = sqrt(0.25 + i_f_n);
  xr_n = sqrt(0.25 + ir_n);
  qd_n = -nq_n * (4.0 * (3.0 * xr_n*xr_n*xr_n + 6.0 * xr_n*xr_n * xf_n + 4.0 * xr_n * xf_n*xf_n + 2.0 * xf_n*xf_n*xf_n) / (15.0 *pow(xf_n + xr_n,2)) - 0.5);
  qs_n = -nq_n * (4.0 * (3.0 * xf_n*xf_n*xf_n + 6.0 * xf_n*xf_n * xr_n + 4.0 * xf_n * xr_n*xr_n + 2.0 * xr_n*xr_n*xr_n) / (15.0 * pow(xf_n + xr_n,2)) - 0.5);
  qi_n = qs_n + qd_n;
  qb1_n = -gammaa_n * sqrt(vp_n + phiT_n + 1.0e-6) / vtT - (nq_n - 1.0) * qi_n / nq_n;
  //  condassign(qb_n, vgprime_n, qb1_n, -vgprime_n / vtT);
  if(vgprime_n > zero)
    qb_n = qb1_n;
  else
    qb_n = -vgprime_n / vtT;

  qg_n = -qi_n - qox - qb_n;

  xf_p = sqrt(0.25 + i_f_p);
  xr_p = sqrt(0.25 + ir_p);
  qd_p = -nq_p * (4 * (3 * xr_p*xr_p*xr_p + 6 * xr_p*xr_p * xf_p + 4 * xr_p * xf_p*xf_p + 2 * xf_p*xf_p*xf_p) / (15 *pow(xf_p + xr_p,2)) - 0.5);
  qs_p = -nq_p * (4 * (3 * xf_p*xf_p*xf_p + 6 * xf_p*xf_p * xr_p + 4 * xf_p * xr_p*xr_p + 2 * xr_p*xr_p*xr_p) / (15 *pow(xf_p + xr_p,2)) - 0.5);
  qi_p = qs_p + qd_p;
  qb1_p = -gammaa_p * sqrt(vp_p + phiT_p + 1e-6) / vtT - (nq_p - 1) * qi_p / nq_p;
  //  condassign(qb_p, vgprime_p, qb1_p, -vgprime_p / vtT);
  if(vgprime_p > zero)
    qb_p = qb1_p;
  else
    qb_p = -vgprime_p / vtT;

  qg_p = -qi_p - qox - qb_p;

  // Rigorous mobility reduction model
  beta_n = betaoprime_n / (1.0 + cox_n * vtT * fabs(qb_n + eta_n*qi_n) / (eo_n * epsilonsi));
  beta_p = betaoprime_p / (1 + cox_p * vtT * fabs(qb_p + eta_p*qi_p) / (eo_p * epsilonsi));

  // Specific current
  is_n = 2.0 * n_n * beta_n * vtT*vtT;
  is_p = 2 * n_p * beta_p * vtT*vtT;

  // Drain-to-source current
  ids_n = type_n*is_n * (i_f_n - irprime_n);
  ids_p = type_p*is_p * (i_f_p - irprime_p);

  // Impact ionization current
  vib_n =  type_n*x[2] - type_n*0 - 2.0 * ibn_n * vdss_n;
  idb1_n = ids_n * iba_n * vib_n * exp(-ibbT_n * lc_n / vib_n) / ibbT_n;
  // condassign(idb_n, vib_n, idb1_n, 0.0);
  if (vib_n > zero)
    idb_n = idb1_n;
  else
    idb_n = zero;

  id_n = ids_n + idb_n;

  vib_p =  type_p*(x[2]-x[0]) - type_p*0 - 2 * ibn_p * vdss_p;
  idb1_p = ids_p * iba_p * vib_p * exp(-ibbT_p * lc_p / vib_p) / ibbT_p;
  //  condassign(idb_p, vib_p, idb1_p, 0);
  if (vib_p > zero)
    idb_p = idb1_p;
  else
    idb_p = zero;

  id_p = ids_p + idb_p;

  //Total Charges

  QI_n = C_ox_n * vtT * qi_n;
  QB_n = C_ox_n * vtT * qb_n;
  QD_n = C_ox_n * vtT * qd_n;
  QS_n = C_ox_n * vtT * qs_n;
  QG_n = C_ox_n * vtT * qg_n;

  QI_p = C_ox_p * vtT * qi_p;
  QB_p = C_ox_p * vtT * qb_p;
  QD_p = C_ox_p * vtT * qd_p;
  QS_p = C_ox_p * vtT * qs_p;
  QG_p = C_ox_p * vtT * qg_p;


  //                 New Code According to ADInterface Package for Differentiation

  //...... Current Entering in VDD = AC(Source Charge Derivative) + DC current = d(QS_p)/dt - id_p
  // d(QS_p) / dt = f(Vsb_p,Vdb_p,Vgb_p) = f(0,x[2]-x[0],x[1]-x[0])
  // Csd + Csg considered
  AD iqs_p = type_p*(-QS_p.fastAccessDx(0) * x[4] - QS_p.fastAccessDx(1) * x[5] - QS_p.fastAccessDx(2) * x[6]);


  // .....Current Entering in Gate = AC (Pmos and Nmos Gate Charge derivative) + DC(0) = d(QG_p)/dt + d(QG_n)/dt
  // d(QG_p) / dt = f(Vsb_p,Vdb_p,Vgb_p) = f(0,x[2]-x[0],x[1]-x[0])
  // d(QG_n) / dt = f(Vsb_p,Vdb_p,Vgb_p) = f(0,x[2],x[1])

  // Cgd + Cgs considered
  //AD iqg_p =  type_p*(QG_p.fastAccessDx(0) * x[3] + QG_p.fastAccessDx(1) * x[4] + QG_p.fastAccessDx(2) * x[5]);
  AD iqg_p =  type_p*(-QG_p.fastAccessDx(0) * x[4] - QG_p.fastAccessDx(2) * x[6]);

  //AD iqg_n = type_n*(QG_n.fastAccessDx(1) * x[4] + QG_n.fastAccessDx(2) * x[5]);
  AD iqg_n = type_n*(-QG_n.fastAccessDx(2) * x[6]);

  // .....Current Entering in Output AC (Drain Charge Derivative) + DC = dQD_p/dt + dQD_n/dt + id_p + id_n
  // d(QD_p) / dt = f(Vsb_p,Vdb_p,Vgb_p) = f(0,x[2]-x[0],x[1]-x[0])
  // d(QD_n) / dt = f(Vsb_p,Vdb_p,Vgb_p) = f(0,x[2],x[1])

  // Cds + Cdg considered
  //AD iqd_p =  type_p*(QD_p.fastAccessDx(0) * x[3] + QD_p.fastAccessDx(1) * x[4] + QD_p.fastAccessDx(2) * x[5]);
  AD iqd_p =  type_p*(-QD_p.fastAccessDx(0) * x[4] - QD_p.fastAccessDx(1) * x[5]);

  //AD iqd_n = type_n*(QD_n.fastAccessDx(1) * x[4] + QD_n.fastAccessDx(2) * x[5]); ;
  AD iqd_n = type_n*(-QD_n.fastAccessDx(1) * x[5]);

  // Assign DC currents
  flow[0] = -(id_p+iqs_p); // DC current entering in VDD
  flow[1] = iqg_p + iqg_n ;   // DC gate current
  flow[2] = id_n + id_p + iqd_n + iqd_p; // DC output current
  //flow[3] = -1*rth*(type_n*(id_n+iqd_n)*vds_n + type_p*(id_p+iqd_p)*vds_p); // IHeat = -(Ids + Idb)*Vds
  flow[3] = -1*rth*(type_n*id_n*vds_n + type_p*id_p*vds_p); // IHeat = -(Ids + Idb)*Vds

  // Assign known output voltages
  effort[0] = x[0];  // Vdd
  effort[1] = x[1]; // Vin
  effort[2] = x[2]; // Vout
}

