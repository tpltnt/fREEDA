#include "Cmos2NandT.h"

// Define the number of parameters
const unsigned Cmos2NandT :: n_par = 72;

// Element Information
ItemInfo Cmos2NandT::einfo =
{
  "cmos2nandt",
  "CMOS 2 Input Thermal Nand Gate",
  "Shivam Priyadarshi",
  DEFAULT_ADDRESS"category:cmos,digital"
};

// Parameter Information
ParmInfo Cmos2NandT::pinfo[] =
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
Cmos2NandT::Cmos2NandT(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
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
  setNumTerms(7);

  // Set Flags
  setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);

  // Set Number of States
  setNumberOfStates(5);
}

// Initialization Function
void Cmos2NandT::init() throw(string&)
{
  // Set Constants
  scale=1.0;
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
  // cout << "\n PHIN=" << phi_n << "\n";

  // PHI_P
  if(phi_p == 0.0)
  {
    if (psub > 0.0)
      phi_p = 2.0 * vtTnom * log(psub * 1.0e6 / niTnom);
    else
      phi_p = 0.7;
  }
  // cout << "\n PHIP=" << phi_p << "\n";
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

  cout << "\n VTON=" << vto_n << "\n";

  //VTO_P
  if(vto_p == 0.0)
  {
    if(vfb_p != 2003.0)
      vto_p =  1*vfb_p + phi_p + gamma_p * sqrt(phi_p);  // VFB positive for PMOS
    else
      vto_p = -1*0.5*type_p;  // Pmos threshold is -ve
  }
  else
    vto_p = type_p*vto_p;
  cout << "\n VTOP=" << vto_p << "\n";

  // KP_n
  if(kp_n == 0.0)
  {
    if(uo_n > 0.0)
      kp_n = uo_n * cox_n * 1.0e-4 /*(m^2/cm^2)*/;
    else
      kp_n = 5.0e-5;
  }
  // cout << "\n KPN=" << kp_n << "\n";

  // KP_p
  if(kp_p == 0.0)
  {
    if(uo_p > 0.0)
      kp_p = uo_p * cox_p * 1.0e-4 /*(m^2/cm^2)*/;
    else
      kp_p = 5.0e-5;
  }
  // cout << "\n KPP=" << kp_p << "\n";

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

  /*cout << "\n UCRITN=" << ucrit_n << "\n";
    cout << "\n UCRITP=" << ucrit_p << "\n"; */

  tcv_n = type_n*tcv_n;
  tcv_p = type_p*tcv_p;
  if(eo_n==0) {eo_n = 1.0e12;}
  if(eo_p==0) {eo_p = 1.0e12;}

  //create tape
  DenseIntVector var(5);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  var[3] = 3;
  var[4] = 4;
  initializeAD(var,var);
}

void Cmos2NandT::getLocalRefIdx(UnsignedVector& local_ref_vec,TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));  //Vdd
  term_list.push_back(getTerminal(1));  //Vin1
  term_list.push_back(getTerminal(2));  //Vin2
  term_list.push_back(getTerminal(3));  //Vout
  term_list.push_back(getTerminal(4));  // Local reference terminal (Gnd)
  term_list.push_back(getTerminal(5));
  term_list.push_back(getTerminal(6)); // Local reference terminal (Tref)
  local_ref_vec.push_back(4); // Local reference index
  local_ref_vec.push_back(6); // Local reference index
}

// Evaluate Function
void Cmos2NandT::eval(AD * x, AD * effort, AD * flow)
{

  //EKV Pmos and nmos models are used for building nand gate  which has four state variable :
  //x[0]- Vdd, x[1] - Vin1, x[2] - Vin2, x[3] - Vout
  //For Simplicity it is assumed : Bulk & Source of Nmos = 0, Bulk & Source of Pmos = Vdd
  //Pmos and nmos state varibles are found in terms of nand state variable:
  //MP1 : Vdb- x[3]-x[0]  Vgb- x[1]-x[0] Vsb=0
  //MP2 : Vdb- x[3]-x[0]  Vgb- x[2]-x[0] Vsb=0
  //MN1 : Vdb- x[3] Vgb- x[1] Vsb=0;
  //MN2 : Vdb- x[3] Vgb- x[2] Vsb=0;
  //x[0]:VDD
  //x[1]:Vin1
  //x[2]:Vin2
  //x[3]:Vout
  //x[4]:Temperature-tmp=Temp_rise->With respect to tmp(Model Simulation Temperature i.e Ambient/HeatSink Temperature)
  //x[5]:d(VDD)/dt
  //x[6]:d(Vin1)/dt
  //x[7]:d(Vin2)/dt
  //x[8]:d(Vout)/dt
  //x[9]:d(Temp_rise)/dt
  //z1[0] : Current in Vdd
  //z1[1] : Current in Vin1
  //z1[2] : Current in Vin2
  //z1[3] : Current in Vout
  //z1[4] : Heat Current
  //z1[5] : x[0] : VDD
  //z1[6] : x[1] : Vin1
  //z1[7] : x[2] : Vin2
  //z1[8] : x[3] : Vout
  //z1[9] : x[4] + tmp : Temp_rise + Temp_Ambient/Sink = Total Temperature

  effort[4] = x[4] + tmp;  // Total Temperature

  AD vtoT_n,vtoT_p, kpT_n,kpT_p,ucritT_n,ucritT_p, phiT_n,phiT_p,ibbT_n,ibbT_p;
  AD vtT,egT,niT;


  // Thermal Parameter Calculation  : Calculated as function of Temp
  vtT = (k*effort[4])/q;    // Thermal voltage
  egT = 1.16-0.000702*effort[4]*effort[4]/(effort[4]+1108.0);               // Energy gap
  niT = 1.45e16*(effort[4]/Tref)*exp(egTref/(2.0*vtTref)-egT/(2.0*vtT)); // Temperature Dependency of Intrinsic Carrier Concentration

  // Intrinsic parameters temperature dependence
  vtoT_n=vto_n-tcv_n*(effort[4]-tnom);
  vtoT_p=vto_p-tcv_p*(effort[4]-tnom);

  kpT_n=kp_n*pow(effort[4]/tnom,bex_n);
  kpT_p=kp_p*pow(effort[4]/tnom,bex_p);
  ucritT_n=ucrit_n*pow(effort[4]/tnom,ucex_n);
  ucritT_p=ucrit_p*pow(effort[4]/tnom,ucex_p);
  phiT_n=phi_n*effort[4]/tnom-3.0*vtT*log(effort[4]/tnom)-egTnom*effort[4]/tnom+egT;
  phiT_p=phi_p*effort[4]/tnom-3.0*vtT*log(effort[4]/tnom)-egTnom*effort[4]/tnom+egT;
  ibbT_n=ibb_n*(1+ibbt_n*(effort[4]-tnom));
  ibbT_p=ibb_p*(1+ibbt_p*(effort[4]-tnom));


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
  AD vgprime_n1,vgprime_p1,vgprime_n2,vgprime_p2, vpo1_n1, vpo_n1,vpo1_p1,vpo_p1,vpo1_n2, vpo_n2,vpo1_p2,vpo_p2, vsprime_n1, vdprime_n1,vsprime_p1, vdprime_p1,vsprime_n2, vdprime_n2,vsprime_p2,vdprime_p2;

  AD gammao_n1, gammaprime_n1, vp1_n1, vp_n1,gammao_p1, gammaprime_p1, vp1_p1, vp_p1, gammao_n2, gammaprime_n2, vp1_n2, vp_n2,gammao_p2, gammaprime_p2, vp1_p2, vp_p2 ;

  AD n_n1,n_p1,n_n2,n_p2,i_f_n1,i_f_p1, i_f_n2,i_f_p2, vdss_n1,vdssprime_n1,vdss_p1,vdssprime_p1,vdss_n2,vdssprime_n2,vdss_p2,vdssprime_p2;

  AD deltav_n1, vds_n1, vip_n1, deltal_n1,deltav_p1, vds_p1, vip_p1, deltal_p1 ,deltav_n2, vds_n2, vip_n2, deltal_n2,deltav_p2, vds_p2, vip_p2, deltal_p2;

  AD lprime_n1,leq_n1,lprime_p1,leq_p1,lprime_n2,leq_n2,lprime_p2,leq_p2,irprime_n1,irprime_p1,irprime_n2,irprime_p2,ir_n1,ir_p1,ir_n2,ir_p2;

  AD  betao_n1, betaoprime_n1, beta_n1,betao_p1, betaoprime_p1, beta_p1,betao_n2,betaoprime_n2, beta_n2,betao_p2, betaoprime_p2, beta_p2, is_n1, ids_n1, vib_n1,is_p1, ids_p1, vib_p1,is_n2,ids_n2, vib_n2,is_p2, ids_p2, vib_p2;

  AD idb1_n1, idb_n1, id_n1,idb1_p1, idb_p1, id_p1, nq_n1,nq_p1,nq_n2,nq_p2, idb1_n2, idb_n2, id_n2,idb1_p2, idb_p2, id_p2;

  AD xf_n1, xr_n1, qd_n1, qs_n1, qi_n1, qb1_n1, qb_n1, qg_n1, QI_n1, QB_n1, QD_n1,xf_p1,xr_p1, qd_p1, qs_p1, qi_p1, qb1_p1, qb_p1, qg_p1, QI_p1, QB_p1, QD_p1;

  AD xf_n2, xr_n2, qd_n2, qs_n2, qi_n2, qb1_n2, qb_n2, qg_n2, QI_n2, QB_n2, QD_n2,xf_p2,xr_p2, qd_p2, qs_p2, qi_p2, qb1_p2, qb_p2, qg_p2, QI_p2, QB_p2, QD_p2;

  AD QS_n1, QG_n1,QS_p1,QG_p1,QS_n2,QG_n2,QS_p2,QG_p2,QD_O,in_O,QS_O,vdsn_O;

  // Effective gate voltage including reverse short channel effect

  vgprime_n1 =  type_n*x[1] - vtoa_n - deltavRSCE_n + phiT_n + gammaa_n * sqrt(phiT_n); // Vgb of MN1 : x[1]
  vgprime_n2 =  type_n*x[2] - vtoa_n - deltavRSCE_n + phiT_n + gammaa_n * sqrt(phiT_n); // Vgb of MN2 : x[2]

  vgprime_p1 =  type_p*(x[1]-x[0]) - vtoa_p - deltavRSCE_p + phiT_p + gammaa_p * sqrt(phiT_p); // Vgb of MP1 : x[1]-x[0]
  vgprime_p2 =  type_p*(x[2]-x[0]) - vtoa_p - deltavRSCE_p + phiT_p + gammaa_p * sqrt(phiT_p); // Vgb of MP2 : x[2]-x[0]

  // Effective substrate factor including charge-sharing for short and narrow
  // channels
  // Pinch-off voltage for narrow-channel effect
  vpo1_n1 = vgprime_n1 - phiT_n - gammaa_n * (sqrt(vgprime_n1 + gammaa_n*gammaa_n / 4.0) - gammaa_n/ 2.0);
  //condassign(vpo_n1, vgprime_n1, vpo1_n1, -phiT_n);
  if(vgprime_n1 > zero)
    vpo_n1 = vpo1_n1;
  else
    vpo_n1 = -phiT_n;

  vpo1_n2 = vgprime_n2 - phiT_n - gammaa_n * (sqrt(vgprime_n2 + gammaa_n*gammaa_n / 4.0) - gammaa_n/ 2.0);
  // condassign(vpo_n2, vgprime_n2, vpo1_n2, -phiT_n);
  if(vgprime_n2 > zero)
    vpo_n2 = vpo1_n2;
  else
    vpo_n2 = -phiT_n;

  vpo1_p1 = vgprime_p1 - phiT_p - gammaa_p * (sqrt(vgprime_p1 + gammaa_p*gammaa_p / 4) - gammaa_p/ 2);
  //  condassign(vpo_p1, vgprime_p1, vpo1_p1, -phiT_p);
  if(vgprime_p1 > zero)
    vpo_p1 = vpo1_p1;
  else
    vpo_p1 = -phiT_p;

  vpo1_p2 = vgprime_p2 - phiT_p - gammaa_p * (sqrt(vgprime_p2 + gammaa_p*gammaa_p / 4) - gammaa_p/ 2);
  //  condassign(vpo_p2, vgprime_p2, vpo1_p2, -phiT_p);
  if(vgprime_p2 > zero)
    vpo_p2 = vpo1_p2;
  else
    vpo_p2 = -phiT_p;




  // Effective substrate factor accounting for charge-sharing
  vsprime_n1 = 0.5 * (type_n*0 + phiT_n + sqrt(pow(type_n*0 + phiT_n,2) + 16.0 * vtT*vtT)); // vsb for MN1 =0
  vdprime_n1 = 0.5 * ( type_n*x[3] + phiT_n + sqrt(pow( type_n*x[3] + phiT_n,2) + 16.0 * vtT*vtT)); // Vdb for MN1 : x[3]
  vsprime_n2 = 0.5 * (type_n*0 + phiT_n + sqrt(pow(type_n*0 + phiT_n,2) + 16.0 * vtT*vtT)); // vsb for MN2 =0
  vdprime_n2 = 0.5 * ( type_n*x[3] + phiT_n + sqrt(pow( type_n*x[3] + phiT_n,2) + 16.0 * vtT*vtT)); // Vdb for MN2 : x[3]


  vsprime_p1 = 0.5 * (type_p*0 + phiT_p + sqrt(pow(type_p*0 + phiT_p,2) + 16.0 * vtT*vtT)); // vsb for MP1 =0
  vdprime_p1 = 0.5 * ( type_p*(x[3]-x[0]) + phiT_p + sqrt(pow( type_p*(x[3]-x[0]) + phiT_p,2) + 16.0 * vtT*vtT)); // Vdb for MP1 : x[3]-x[0]
  vsprime_p2 = 0.5 * (type_p*0 + phiT_p + sqrt(pow(type_p*0 + phiT_p,2) + 16.0 * vtT*vtT)); // vsb for MP2 =0
  vdprime_p2 = 0.5 * ( type_p*(x[3]-x[0]) + phiT_p + sqrt(pow( type_p*(x[3]-x[0]) + phiT_p,2) + 16.0 * vtT*vtT)); // Vdb for MP2 : x[3]-x[0]

  // Pinch-off voltage including short- and narrow-channel effect
  gammao_n1 = gammaa_n - epsilonsi * (leta_n * (sqrt(vsprime_n1) + sqrt(vdprime_n1))/ leff_n - 3.0 * weta_n * sqrt(vpo_n1 + phiT_n) / weff_n) / cox_n;
  gammaprime_n1 = 0.5 * (gammao_n1 + sqrt(gammao_n1*gammao_n1 + 0.1 * vtT));
  gammao_n2 = gammaa_n - epsilonsi * (leta_n * (sqrt(vsprime_n2) + sqrt(vdprime_n2))/ leff_n - 3.0 * weta_n * sqrt(vpo_n2 + phiT_n) / weff_n) / cox_n;
  gammaprime_n2 = 0.5 * (gammao_n2 + sqrt(gammao_n2*gammao_n2 + 0.1 * vtT));

  vp1_n1 = vgprime_n1 - phiT_n - gammaprime_n1 * (sqrt(vgprime_n1+pow(gammaprime_n1 / 2.0,2)) -gammaprime_n1 / 2.0);
  //  condassign(vp_n1, vgprime_n1, vp1_n1, -phiT_n);
  if(vgprime_n1 > zero)
    vp_n1 = vp1_n1;
  else
    vp_n1 = -phiT_n;

  vp1_n2 = vgprime_n2 - phiT_n - gammaprime_n2 * (sqrt(vgprime_n2+pow(gammaprime_n2 / 2.0,2)) -gammaprime_n2 / 2.0);
  //  condassign(vp_n2, vgprime_n2, vp1_n2, -phiT_n);
  if(vgprime_n2 > zero)
    vp_n2 = vp1_n2;
  else
    vp_n2 = -phiT_n;

  gammao_p1 = gammaa_p - epsilonsi * (leta_p * (sqrt(vsprime_p1) + sqrt(vdprime_p1))/ leff_p - 3.0 * weta_p * sqrt(vpo_p1 + phiT_p) / weff_p) / cox_p;
  gammaprime_p1 = 0.5 * (gammao_p1 + sqrt(gammao_p1*gammao_p1 + 0.1 * vtT));
  gammao_p2 = gammaa_p - epsilonsi * (leta_p * (sqrt(vsprime_p2) + sqrt(vdprime_p2))/ leff_p - 3.0 * weta_p * sqrt(vpo_p2 + phiT_p) / weff_p) / cox_p;
  gammaprime_p2 = 0.5 * (gammao_p2 + sqrt(gammao_p2*gammao_p2 + 0.1 * vtT));

  vp1_p1 = vgprime_p1 - phiT_p - gammaprime_p1 * (sqrt(vgprime_p1+pow(gammaprime_p1 / 2.0,2)) -gammaprime_p1 / 2.0);
  //  condassign(vp_p1, vgprime_p1, vp1_p1, -phiT_p);
  if(vgprime_p1 > zero)
    vp_p1 = vp1_p1;
  else
    vp_p1 = -phiT_p;

  vp1_p2 = vgprime_p2 - phiT_p - gammaprime_p2 * (sqrt(vgprime_p2+pow(gammaprime_p2 / 2.0,2)) -gammaprime_p2 / 2.0);
  //  condassign(vp_p2, vgprime_p2, vp1_p2, -phiT_p);
  if(vgprime_p2 > zero)
    vp_p2 = vp1_p2;
  else
    vp_p2 = -phiT_p;

  // Slope factor
  n_n1 = 1 + gammaa_n / (2.0 * sqrt(vp_n1 + phiT_n + 4.0 * vtT));
  n_n2 = 1 + gammaa_n / (2.0 * sqrt(vp_n2 + phiT_n + 4.0 * vtT));

  n_p1 = 1 + gammaa_p / (2.0 * sqrt(vp_p1 + phiT_p + 4.0 * vtT));
  n_p2 = 1 + gammaa_p / (2.0 * sqrt(vp_p2 + phiT_p + 4.0 * vtT));

  // Forward normalized current
  i_f_n1 = log(1.0 + exp((vp_n1 - type_n*0) / (2.0*vtT)))*log(1.0 + exp((vp_n1 - type_n*0) / (2.0*vtT))); //Vsb of nmos =0
  i_f_n2 = log(1.0 + exp((vp_n2 - type_n*0) / (2.0*vtT)))*log(1.0 + exp((vp_n2 - type_n*0) / (2.0*vtT))); //Vsb of nmos =0

  i_f_p1 = log(1 + exp((vp_p1 - type_p*0) / (2 * vtT)))*log(1 + exp((vp_p1 - type_p*0) / (2 * vtT))); // Vsb for pmos =0
  i_f_p2 = log(1 + exp((vp_p2 - type_p*0) / (2 * vtT)))*log(1 + exp((vp_p2 - type_p*0) / (2 * vtT))); // Vsb for pmos =0



  // Velocity saturation voltage
  vdss_n1 = vc_n * (sqrt(0.25 + vtT * sqrt(i_f_n1) / vc_n) - 0.5);
  vdss_n2 = vc_n * (sqrt(0.25 + vtT * sqrt(i_f_n2) / vc_n) - 0.5);

  vdss_p1 = vc_p * (sqrt(0.25 + vtT * sqrt(i_f_p1) / vc_p) - 0.5);
  vdss_p2 = vc_p * (sqrt(0.25 + vtT * sqrt(i_f_p2) / vc_p) - 0.5);


  // Drain-to-source saturation voltage for reverse normalized current
  vdssprime_n1 = vc_n * (sqrt(0.25 + vtT * (sqrt(i_f_n1) - 0.75 * log(i_f_n1))/vc_n) - 0.5) +vtT * (log(vc_n / (2.0 * vtT)) - 0.6);
  vdssprime_n2 = vc_n * (sqrt(0.25 + vtT * (sqrt(i_f_n2) - 0.75 * log(i_f_n2))/vc_n) - 0.5) +vtT * (log(vc_n / (2.0 * vtT)) - 0.6);

  vdssprime_p1 = vc_p * (sqrt(0.25 + vtT * (sqrt(i_f_p1) - 0.75 * log(i_f_p1))/vc_p) - 0.5) +vtT * (log(vc_p / (2 * vtT)) - 0.6);
  vdssprime_p2 = vc_p * (sqrt(0.25 + vtT * (sqrt(i_f_p2) - 0.75 * log(i_f_p2))/vc_p) - 0.5) +vtT * (log(vc_p / (2 * vtT)) - 0.6);


  // Channel-length modulation
  deltav_n1 = 4.0 * vtT * sqrt(lambda_n * (sqrt(i_f_n1) - vdss_n1 / vtT) + 1.0/64.0);
  vds_n1 = (type_n*x[3] - type_n*0) / 2.0; // Vdb = x[3] and Vsb= 0 for both nmos
  vip_n1 = sqrt(vdss_n1*vdss_n1 + deltav_n1*deltav_n1) - sqrt(pow(vds_n1 - vdss_n1,2) + deltav_n1*deltav_n1);
  deltal_n1 = lambda_n * lc_n * log(1.0 + (vds_n1 - vip_n1) / (lc_n * ucritT_n));

  deltav_n2 = 4.0 * vtT * sqrt(lambda_n * (sqrt(i_f_n2) - vdss_n2 / vtT) + 1.0/64.0);
  vds_n2 = (type_n*x[3] - type_n*0) / 2.0; // Vdb = x[3] and Vsb= 0 for both nmos
  vip_n2 = sqrt(vdss_n2*vdss_n2 + deltav_n2*deltav_n2) - sqrt(pow(vds_n2 - vdss_n2,2) + deltav_n2*deltav_n2);
  deltal_n2 = lambda_n * lc_n * log(1.0 + (vds_n2 - vip_n2) / (lc_n * ucritT_n));


  deltav_p1 = 4 * vtT * sqrt(lambda_p * (sqrt(i_f_p1) - vdss_p1 / vtT) + 1 /64);
  vds_p1 = ( type_p*(x[3]-x[0]) - type_p*0) / 2; //Vdb =x[3]-x[0] & Vsb=0 for both  pmos
  vip_p1 = sqrt(vdss_p1*vdss_p1 + deltav_p1*deltav_p1) - sqrt(pow(vds_p1 - vdss_p1,2) + deltav_p1*deltav_p1);
  deltal_p1 = lambda_p * lc_p * log(1 + (vds_p1 - vip_p1) / (lc_p * ucritT_p));

  deltav_p2 = 4 * vtT * sqrt(lambda_p * (sqrt(i_f_p2) - vdss_p2 / vtT) + 1 /64);
  vds_p2 = ( type_p*(x[3]-x[0]) - type_p*0) / 2; //Vdb =x[3]-x[0] & Vsb=0 for both  pmos
  vip_p2 = sqrt(vdss_p2*vdss_p2 + deltav_p2*deltav_p2) - sqrt(pow(vds_p2 - vdss_p2,2) + deltav_p2*deltav_p2);
  deltal_p2 = lambda_p * lc_p * log(1 + (vds_p2 - vip_p2) / (lc_p * ucritT_p));


  // Equivalent channel length including channel-length moculation and velocity
  // saturation
  lprime_n1 = ns_n * leff_n - deltal_n1 + (vds_n1 + vip_n1) / ucritT_n;
  leq_n1 = 0.5 * (lprime_n1 + sqrt(lprime_n1*lprime_n1 + lmin_n*lmin_n));
  lprime_n2 = ns_n * leff_n - deltal_n2 + (vds_n2 + vip_n2) / ucritT_n;
  leq_n2 = 0.5 * (lprime_n2 + sqrt(lprime_n2*lprime_n2 + lmin_n*lmin_n));

  lprime_p1 = ns_p * leff_p - deltal_p1 + (vds_p1 + vip_p1) / ucritT_p;
  leq_p1 = 0.5 * (lprime_p1 + sqrt(lprime_p1*lprime_p1 + lmin_p*lmin_p));
  lprime_p2 = ns_p * leff_p - deltal_p2 + (vds_p2 + vip_p2) / ucritT_p;
  leq_p2 = 0.5 * (lprime_p2 + sqrt(lprime_p2*lprime_p2 + lmin_p*lmin_p));


  // Reverse normalized current
  irprime_n1 = log(1.0 + exp(((vp_n1 - vds_n1 - type_n*0 - sqrt(vdssprime_n1*vdssprime_n1 +deltav_n1*deltav_n1) + sqrt((vds_n1-vdssprime_n1)*(vds_n1-vdssprime_n1) + deltav_n1*deltav_n1)) /vtT)/ 2.0))*log(1.0 + exp(((vp_n1 - vds_n1 - type_n*0 - sqrt(vdssprime_n1*vdssprime_n1 +deltav_n1*deltav_n1) + sqrt((vds_n1-vdssprime_n1)*(vds_n1-vdssprime_n1) + deltav_n1*deltav_n1)) / vtT)/ 2.0));

  irprime_n2 = log(1.0 + exp(((vp_n2 - vds_n2 - type_n*0 - sqrt(vdssprime_n2*vdssprime_n2 +deltav_n2*deltav_n2) + sqrt((vds_n2-vdssprime_n2)*(vds_n2-vdssprime_n2) + deltav_n2*deltav_n2)) /vtT)/ 2.0))*log(1.0 + exp(((vp_n2 - vds_n2 - type_n*0 - sqrt(vdssprime_n2*vdssprime_n2 +deltav_n2*deltav_n2) + sqrt((vds_n2-vdssprime_n2)*(vds_n2-vdssprime_n2) + deltav_n2*deltav_n2)) / vtT)/ 2.0));

  irprime_p1 = log(1 + exp(((vp_p1 - vds_p1 - type_p*0 - sqrt(vdssprime_p1*vdssprime_p1 +deltav_p1*deltav_p1) + sqrt((vds_p1-vdssprime_p1)*(vds_p1-vdssprime_p1) + deltav_p1*deltav_p1)) /vtT)/ 2))*log(1 + exp(((vp_p1 - vds_p1 - type_p*0 - sqrt(vdssprime_p1*vdssprime_p1 +deltav_p1*deltav_p1) + sqrt((vds_p1-vdssprime_p1)*(vds_p1-vdssprime_p1) + deltav_p1*deltav_p1)) / vtT)/ 2));

  irprime_p2 = log(1 + exp(((vp_p2 - vds_p2 - type_p*0 - sqrt(vdssprime_p2*vdssprime_p2 +deltav_p2*deltav_p2) + sqrt((vds_p2-vdssprime_p2)*(vds_p2-vdssprime_p2) + deltav_p2*deltav_p2)) /vtT)/ 2))*log(1 + exp(((vp_p2 - vds_p2 - type_p*0 - sqrt(vdssprime_p2*vdssprime_p2 +deltav_p2*deltav_p2) + sqrt((vds_p2-vdssprime_p2)*(vds_p2-vdssprime_p2) + deltav_p2*deltav_p2)) / vtT)/ 2));


  // Reverse normalized currect for mobility model, intrinsic
  //charges/capacitances, thermal noise model and NQS time-constant
  ir_n1 = log(1.0 + exp((vp_n1 -  type_n*x[3]) / (2.0*vtT)))*log(1.0 + exp((vp_n1 -  type_n*x[3]) / (2.0*vtT))); //Vdb of both nmos x[3]
  ir_n2 = log(1.0 + exp((vp_n2 -  type_n*x[3]) / (2.0*vtT)))*log(1.0 + exp((vp_n2 -  type_n*x[3]) / (2.0*vtT))); //Vdb of both nmos x[3]

  ir_p1 = log(1 + exp((vp_p1 -  type_p*(x[3]-x[0])) / (2 * vtT)))*log(1 + exp((vp_p1 -  type_p*(x[3]-x[0])) / (2 * vtT))); //Vdb of pmos x[3]-x[0]
  ir_p2 = log(1 + exp((vp_p2 -  type_p*(x[3]-x[0])) / (2 * vtT)))*log(1 + exp((vp_p2 -  type_p*(x[3]-x[0])) / (2 * vtT))); //Vdb of pmos x[3]-x[0]

  // Transconductance factor and mobility reduction due to vertical field
  betao_n1 = kpa_n * np_n * weff_n / leq_n1;
  betaoprime_n1 = betao_n1 * (1.0 + cox_n * qbo_n / (eo_n * epsilonsi));
  betao_n2 = kpa_n * np_n * weff_n / leq_n2;
  betaoprime_n2 = betao_n2 * (1.0 + cox_n * qbo_n / (eo_n * epsilonsi));

  betao_p1 = kpa_p * np_p * weff_p / leq_p1;
  betaoprime_p1 = betao_p1 * (1 + cox_p * qbo_p / (eo_p * epsilonsi));
  betao_p2 = kpa_p * np_p * weff_p / leq_p2;
  betaoprime_p2 = betao_p2 * (1 + cox_p * qbo_p / (eo_p * epsilonsi));


  // Dynamic model for the intrinsic node charges
  nq_n1 = 1.0 + gammaa_n / (2.0 * sqrt(vp_n1 + phiT_n +1.0e-6));
  nq_n2 = 1.0 + gammaa_n / (2.0 * sqrt(vp_n2 + phiT_n +1.0e-6));

  nq_p1 = 1 + gammaa_p / (2 * sqrt(vp_p1 + phiT_p +1e-6));
  nq_p2 = 1 + gammaa_p / (2 * sqrt(vp_p2 + phiT_p +1e-6));

  // Normalized intrinsic node charges
  xf_n1 = sqrt(0.25 + i_f_n1);
  xr_n1 = sqrt(0.25 + ir_n1);
  qd_n1 = -nq_n1 * (4.0 * (3.0 * xr_n1*xr_n1*xr_n1 + 6.0 * xr_n1*xr_n1 * xf_n1 + 4.0 * xr_n1 * xf_n1*xf_n1 + 2.0 * xf_n1*xf_n1*xf_n1) / (15.0 *pow(xf_n1 + xr_n1,2)) - 0.5);
  qs_n1 = -nq_n1 * (4.0 * (3.0 * xf_n1*xf_n1*xf_n1 + 6.0 * xf_n1*xf_n1 * xr_n1 + 4.0 * xf_n1 * xr_n1*xr_n1 + 2.0 * xr_n1*xr_n1*xr_n1) / (15.0 * pow(xf_n1 + xr_n1,2)) - 0.5);
  qi_n1 = qs_n1 + qd_n1;
  qb1_n1 = -gammaa_n * sqrt(vp_n1 + phiT_n + 1.0e-6) / vtT - (nq_n1 - 1.0) * qi_n1 / nq_n1;
  // condassign(qb_n1, vgprime_n1, qb1_n1, -vgprime_n1 / vtT);
  if(vgprime_n1 > zero)
    qb_n1 = qb1_n1;
  else
    qb_n1 = -vgprime_n1 / vtT;

  qg_n1 = -qi_n1 - qox - qb_n1;

  xf_n2 = sqrt(0.25 + i_f_n2);
  xr_n2 = sqrt(0.25 + ir_n2);
  qd_n2 = -nq_n2 * (4.0 * (3.0 * xr_n2*xr_n2*xr_n2 + 6.0 * xr_n2*xr_n2 * xf_n2 + 4.0 * xr_n2 * xf_n2*xf_n2 + 2.0 * xf_n2*xf_n2*xf_n2) / (15.0 *pow(xf_n2 + xr_n2,2)) - 0.5);
  qs_n2 = -nq_n2 * (4.0 * (3.0 * xf_n2*xf_n2*xf_n2 + 6.0 * xf_n2*xf_n2 * xr_n2 + 4.0 * xf_n2 * xr_n2*xr_n2 + 2.0 * xr_n2*xr_n2*xr_n2) / (15.0 * pow(xf_n2 + xr_n2,2)) - 0.5);
  qi_n2 = qs_n2 + qd_n2;
  qb1_n2 = -gammaa_n * sqrt(vp_n2 + phiT_n + 1.0e-6) / vtT - (nq_n2 - 1.0) * qi_n2 / nq_n2;
  //  condassign(qb_n2, vgprime_n2, qb1_n2, -vgprime_n2 / vtT);
  if(vgprime_n2 > zero)
    qb_n2 = qb1_n2;
  else
    qb_n2 = -vgprime_n2 / vtT;
  qg_n2 = -qi_n2 - qox - qb_n2;


  xf_p1 = sqrt(0.25 + i_f_p1);
  xr_p1 = sqrt(0.25 + ir_p1);
  qd_p1 = -nq_p1 * (4 * (3 * xr_p1*xr_p1*xr_p1 + 6 * xr_p1*xr_p1 * xf_p1 + 4 * xr_p1 * xf_p1*xf_p1 + 2 * xf_p1*xf_p1*xf_p1) / (15 *pow(xf_p1 + xr_p1,2)) - 0.5);
  qs_p1 = -nq_p1 * (4 * (3 * xf_p1*xf_p1*xf_p1 + 6 * xf_p1*xf_p1 * xr_p1 + 4 * xf_p1 * xr_p1*xr_p1 + 2 * xr_p1*xr_p1*xr_p1) / (15 *pow(xf_p1 + xr_p1,2)) - 0.5);
  qi_p1 = qs_p1 + qd_p1;
  qb1_p1 = -gammaa_p * sqrt(vp_p1 + phiT_p + 1e-6) / vtT - (nq_p1 - 1) * qi_p1 / nq_p1;
  //  condassign(qb_p1, vgprime_p1, qb1_p1, -vgprime_p1 / vtT);
  if(vgprime_p1 > zero)
    qb_p1 = qb1_p1;
  else
    qb_p1 = -vgprime_p1 / vtT;

  qg_p1 = -qi_p1 - qox - qb_p1;

  xf_p2 = sqrt(0.25 + i_f_p2);
  xr_p2 = sqrt(0.25 + ir_p2);
  qd_p2 = -nq_p2 * (4 * (3 * xr_p2*xr_p2*xr_p2 + 6 * xr_p2*xr_p2 * xf_p2 + 4 * xr_p2 * xf_p2*xf_p2 + 2 * xf_p2*xf_p2*xf_p2) / (15 *pow(xf_p2 + xr_p2,2)) - 0.5);
  qs_p2 = -nq_p2 * (4 * (3 * xf_p2*xf_p2*xf_p2 + 6 * xf_p2*xf_p2 * xr_p2 + 4 * xf_p2 * xr_p2*xr_p2 + 2 * xr_p2*xr_p2*xr_p2) / (15 *pow(xf_p2 + xr_p2,2)) - 0.5);
  qi_p2 = qs_p2 + qd_p2;
  qb1_p2 = -gammaa_p * sqrt(vp_p2 + phiT_p + 1e-6) / vtT - (nq_p2 - 1) * qi_p2 / nq_p2;
  //  condassign(qb_p2, vgprime_p2, qb1_p2, -vgprime_p2 / vtT);
  if(vgprime_p2 > zero)
    qb_p2 = qb1_p2;
  else
    qb_p2 = -vgprime_p2 / vtT;
  qg_p2 = -qi_p2 - qox - qb_p2;


  // Rigorous mobility reduction model
  beta_n1 = betaoprime_n1 / (1.0 + cox_n * vtT * fabs(qb_n1 + eta_n*qi_n1) / (eo_n * epsilonsi));
  beta_n2 = betaoprime_n2 / (1.0 + cox_n * vtT * fabs(qb_n2 + eta_n*qi_n2) / (eo_n * epsilonsi));

  beta_p1 = betaoprime_p1 / (1 + cox_p * vtT * fabs(qb_p1 + eta_p*qi_p1) / (eo_p * epsilonsi));
  beta_p2 = betaoprime_p2 / (1 + cox_p * vtT * fabs(qb_p2 + eta_p*qi_p2) / (eo_p * epsilonsi));

  // Specific current
  is_n1 = 2.0 * n_n1 * beta_n1 * vtT*vtT;
  is_n2 = 2.0 * n_n2 * beta_n2 * vtT*vtT;

  is_p1 = 2 * n_p1 * beta_p1 * vtT*vtT;
  is_p2 = 2 * n_p2 * beta_p2 * vtT*vtT;

  // Drain-to-source current
  ids_n1 = type_n*is_n1 * (i_f_n1 - irprime_n1);
  ids_n2 = type_n*is_n2 * (i_f_n2 - irprime_n2);
  ids_p1 = type_p*is_p1 * (i_f_p1 - irprime_p1);
  ids_p2 = type_p*is_p2 * (i_f_p2 - irprime_p2);

  // Impact ionization current
  vib_n1 =  type_n*x[3] - type_n*0 - 2.0 * ibn_n * vdss_n1;
  idb1_n1 = ids_n1 * iba_n * vib_n1 * exp(-ibbT_n * lc_n / vib_n1) / ibbT_n;
  //  condassign(idb_n1, vib_n1, idb1_n1, 0.0);
  if(vib_n1 > zero)
    idb_n1 = idb1_n1;
  else
    idb_n1 = zero;

  id_n1 = ids_n1 + idb_n1;

  vib_n2 =  type_n*x[3] - type_n*0 - 2.0 * ibn_n * vdss_n2;
  idb1_n2 = ids_n2 * iba_n * vib_n2 * exp(-ibbT_n * lc_n / vib_n2) / ibbT_n;
  //  condassign(idb_n2, vib_n2, idb1_n2, 0.0);
  if(vib_n2 > zero)
    idb_n2 = idb1_n2;
  else
    idb_n2 = zero;

  id_n2 = ids_n2 + idb_n2;

  vib_p1 =  type_p*(x[3]-x[0]) - type_p*0 - 2 * ibn_p * vdss_p1;
  idb1_p1 = ids_p1 * iba_p * vib_p1 * exp(-ibbT_p * lc_p / vib_p1) / ibbT_p;
  //  condassign(idb_p1, vib_p1, idb1_p1, 0);
  if(vib_p1 > zero)
    idb_p1 = idb1_p1;
  else
    idb_p1 = zero;

  id_p1 = ids_p1 + idb_p1;

  vib_p2 =  type_p*(x[3]-x[0]) - type_p*0 - 2 * ibn_p * vdss_p2;
  idb1_p2 = ids_p2 * iba_p * vib_p2 * exp(-ibbT_p * lc_p / vib_p2) / ibbT_p;
  //  condassign(idb_p2, vib_p2, idb1_p2, 0);
  if(vib_p2 > zero)
    idb_p2 = idb1_p2;
  else
    idb_p2 = zero;

  id_p2 = ids_p2 + idb_p2;


  //Total Charges

  QI_n1 = C_ox_n * vtT * qi_n1;
  QB_n1 = C_ox_n * vtT * qb_n1;
  QD_n1 = C_ox_n * vtT * qd_n1;
  QS_n1 = C_ox_n * vtT * qs_n1;
  QG_n1 = C_ox_n * vtT * qg_n1;

  QI_n2 = C_ox_n * vtT * qi_n2;
  QB_n2 = C_ox_n * vtT * qb_n2;
  QD_n2 = C_ox_n * vtT * qd_n2;
  QS_n2 = C_ox_n * vtT * qs_n2;
  QG_n2 = C_ox_n * vtT * qg_n2;

  QI_p1 = C_ox_p * vtT * qi_p1;
  QB_p1 = C_ox_p * vtT * qb_p1;
  QD_p1 = C_ox_p * vtT * qd_p1;
  QS_p1 = C_ox_p * vtT * qs_p1;
  QG_p1 = C_ox_p * vtT * qg_p1;

  QI_p2 = C_ox_p * vtT * qi_p2;
  QB_p2 = C_ox_p * vtT * qb_p2;
  QD_p2 = C_ox_p * vtT * qd_p2;
  QS_p2 = C_ox_p * vtT * qs_p2;
  QG_p2 = C_ox_p * vtT * qg_p2;

  //  condassign(QD_O,(x[2]-x[1]),QD_n1,QD_n2);
  //  condassign(in_O,(x[2]-x[1]),id_n1,id_n2);
  //  condassign(vdsn_O,(x[2]-x[1]),vds_n2,vds_n1);

  if((x[2] - x[1]) > zero)
    QD_O = QD_n1;
  else
    QD_O = QD_n2;

  if((x[2] - x[1]) > zero)
    in_O = id_n1;
  else
    in_O = id_n2;

  if((x[2] - x[1]) > zero)
    vdsn_O = vds_n2;
  else
    vdsn_O = vds_n1;

  //------- Current Entering in VDD
  // The capacitors considered in this model are as shown in
  // the equivalent circuit in the EKV model documentation.
  // In general, Cxy = +/-dQx/dVy * dVy/dt, where the "+" sign is
  // used if x == y, else "-" sign is used

  // QS_p1 = f(Vsb_p1,Vgb_p1,Vdb_p1) = f(0,x[1]-x[0],x[3]-x[0]) Csd+Csg considered
  AD iqs_p1 = type_p*(-QS_p1.fastAccessDx(0) * x[5] - QS_p1.fastAccessDx(1) * x[6] - QS_p1.fastAccessDx(3) *x[8]);

  // QS_p2 = f(Vsb_p2,Vgb_p2,Vdb_p2) = f(0,x[2]-x[0],x[3]-x[0]) Csd+Csg considered
  AD iqs_p2 = type_p*(-QS_p2.fastAccessDx(0) * x[5] - QS_p2.fastAccessDx(2) * x[7] - QS_p2.fastAccessDx(3) * x[8]);

  // .....Current Entering in Gate Vp1 = AC (Pmos and Nmos Gate Charge derivative) + DC(0) = d(QG_p1)/dt + d(QG_n1)/dt
  // QG_p1 = f(Vsb_p1,Vgb_p1,Vdb_p1) = f(0,x[1]-x[0],x[3]-x[0])
  //AD iqg_p1 = type_p*(QG_p1.fastAccessDx(0) * x[4] + QG_p1.fastAccessDx(1) * x[5] + QG_p1.fastAccessDx(3) *x[7]);

  // Cgs + Cgd Considered
  AD iqg_p1 = type_p*(-QG_p1.fastAccessDx(0) * x[5] -  QG_p1.fastAccessDx(3) *x[8]);


  // QG_n1 = f(Vsb_n1,Vgb_n1,Vdb_n1) = f(0,x[1],x[3]
  //AD iqg_n1 = type_n*(QG_n1.fastAccessDx(1) * x[5] + QG_n1.fastAccessDx(3) * x[7]);
  // Cgs + Cgd Considered
  AD iqg_n1 = type_n*(-QG_n1.fastAccessDx(3) * x[8]);

  // .....Current Entering in Gate Vp2 = AC (Pmos and Nmos Gate Charge derivative) + DC(0) = d(QG_p2)/dt + d(QG_n2)/dt
  // QG_p2 = f(Vsb_p2,Vgb_p2,Vdb_p2) = f(0,x[2]-x[0],x[3]-x[0])
  //AD iqg_p2 = type_p*(QG_p2.fastAccessDx(0) * x[4] + QG_p2.fastAccessDx(2) * x[6] + QG_p2.fastAccessDx(3) *x[7]);
  //Cgs + Cgd Considered
  AD iqg_p2 = type_p*(-QG_p2.fastAccessDx(0) * x[5] - QG_p2.fastAccessDx(3) *x[8]);

  // QG_n2 = f(Vsb_n2,Vgb_n2,Vdb_n2) = f(0,x[2],x[3]
  //AD iqg_n2 = type_n*(QG_n2.fastAccessDx(2) * x[6] + QG_n2.fastAccessDx(3) * x[7]);
  // Cgs + Cgd Considered
  AD iqg_n2 = type_n*(-QG_n2.fastAccessDx(3) * x[8]);

  // .....Current Entering in Output AC (Drain Charge Derivative) + DC = dQD_p1/dt + dQD_p2/dt + dQD_O/dt +id_p1 + id_p2 + in_O
  // QD_p1 = f(Vsb_p1,Vgb_p1,Vdb_p1) = f(0,x[1]-x[0],x[3]-x[0])
  //AD iqd_p1 = type_p*(QD_p1.fastAccessDx(0) * x[4] + QD_p1.fastAccessDx(1) * x[5] + QD_p1.fastAccessDx(3) *x[7]);
  // Cdg + Cds considered
  AD iqd_p1 = type_p*(-QD_p1.fastAccessDx(0) * x[5] - QD_p1.fastAccessDx(1) * x[6]);

  // QD_p2 = f(Vsb_p2,Vgb_p2,Vdb_p2) = f(0,x[2]-x[0],x[3]-x[0])
  //AD iqd_p2 = type_p*(QD_p2.fastAccessDx(0) * x[4] + QD_p2.fastAccessDx(2) * x[6] + QD_p2.fastAccessDx(3) * x[7]);
  // Cdg + Cds considered
  AD iqd_p2 = type_p*(-QD_p2.fastAccessDx(0) * x[5] - QD_p2.fastAccessDx(2) * x[7]);

  AD iqd_O;
  if ((x[2] - x[1]) > 0)
  {
    //iqd_O = type_n*(QD_n1.fastAccessDx(1) * x[5] + QD_n1.fastAccessDx(3) * x[7]);
    iqd_O = type_n*(-QD_n1.fastAccessDx(1) * x[6]);  // Cdg +  Cds considered

  }
  else
  {
    //iqd_O =  type_n*(QD_n2.fastAccessDx(2) * x[6] + QD_n2.fastAccessDx(3) * x[7]);
    iqd_O =  type_n*(-QD_n2.fastAccessDx(2) * x[7]); // Cdg + Cds Considered

  }


  flow[0] = -(id_p1 + id_p2+ iqs_p1 + iqs_p2); // AC + DC current entering in VDD
  flow[1] = 0  + iqg_n1 + iqg_p1;   // AC + DC gate current Vin1
  flow[2] = 0  + iqg_n2 + iqg_p2;   // AC + DC gate current Vin2
  flow[3] = in_O + id_p1 + id_p2 + iqd_p1 + iqd_p2  + iqd_O ; // AC + DC output current
  //  flow[4] =  -1*rth*(type_n*in_O*vdsn_O + type_p*id_p1*vds_p1 + type_p*id_p2*vds_p2);

  flow[4] =  -1*rth*(type_n*(in_O+iqd_O)*vdsn_O + type_p*(id_p1+iqs_p1)*vds_p1 + type_p*(id_p2+iqs_p2)*vds_p2);

  //flow[4]  = -1*(type_n*in_O*vdsn_O + type_p*id_p1*vds_p1 + type_p*id_p2*vds_p2 + (type_n*in_O + type_p*id_p1 + type_p*id_p2)*(type_n*in_O + type_p*id_p1 + type_p*id_p2)*rth);

  effort[0] = x[0];  // Vdd
  effort[1] = x[1]; // Vin1
  effort[2] = x[2];  //Vin2
  effort[3] = x[3]; // Vout
}
