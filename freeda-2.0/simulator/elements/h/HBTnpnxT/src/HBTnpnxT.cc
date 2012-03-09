#include "HBTnpnxT.h"

const unsigned HBTnpnxT :: n_par = 104;

ItemInfo HBTnpnxT::einfo =
{
  "hbtnpnxt",
  "HBT Model, NPN, UCSD, Electro-thermal",
  "Jian Ding, Sonali Luniya",
  DEFAULT_ADDRESS"category:transistor>HBT",
  "2005_12_20"
};

//parameter Information
ParmInfo HBTnpnxT :: pinfo[] =
{
  {"type","Type, only NPN currently", TR_DOUBLE, false},
  {"afn","Flicker noise exponent for current", TR_DOUBLE, false},
  {"bf","Forward ideal current gain", TR_DOUBLE, false},
  {"bfn","BE flicker noise exponent for frequency", TR_DOUBLE, false},
  {"bkdn","Flag denoting that BC breakdown should be included [logic]", TR_BOOLEAN, false},
  {"br","Reverse ideal current gain", TR_DOUBLE, false},
  {"bvc","Collector-base breakdown voltage BVcbo [V]", TR_DOUBLE, false},
  {"ccmin","Minimum value of intrinsic BC Cj [F]", TR_DOUBLE, false},
  {"cemin","Minimum BE capacitance [F]", TR_DOUBLE, false},
  {"cjc","Intrinsic BC depletion capacitance at zero bias [F]", TR_DOUBLE, false},
  {"cjcx","Extrinsic BC depletion capacitance at zero bias [F]", TR_DOUBLE, false},
  {"cje","BE depletion capacitance at zero bias [F]", TR_DOUBLE, false},
  {"cjs","Collector-substrate depletion capacitance (0 bias) [F]", TR_DOUBLE, false},
  {"cth","Thermal capacitance of device [C/joule]", TR_DOUBLE, false},
  {"cxmin","Minimum extrinsic Cbc [F]", TR_DOUBLE, false},
  {"dtmax","Maximum expected temperature rise above heatsink [C]", TR_DOUBLE, false},
  {"eaa","Added activation energy for ISE temp dependence [V] ", TR_DOUBLE, false},
  {"eab","Added activation energy for ISC temp dependence [V] ", TR_DOUBLE, false},
  {"eac","Activation energy for ISB temperature dependence [V]", TR_DOUBLE, false},
  {"eae","Activation energy for ISA temperature dependence [V]", TR_DOUBLE, false},
  {"eax","Added activation energy for ISEX temp dependence [V]", TR_DOUBLE, false},
  {"eg","Activation energy for IS temperature dependence [V]", TR_DOUBLE, false},
  {"fa","Factor for specification of avalanche voltage", TR_DOUBLE, false},
  {"fc","Factor for start of high bias BC Cj approximation", TR_DOUBLE, false},
  {"fce","Factor for start of high bias BE Cj approximation", TR_DOUBLE, false},
  {"fex","Factor to determine excess phase", TR_DOUBLE, false},
  {"icrit0","Critical current for intrinsic Cj variation [A]", TR_DOUBLE, false},
  {"ics","Saturation value for collector-substrate current [A]", TR_DOUBLE, false},
  {"ik","Knee current for dc high injection effect [A]", TR_DOUBLE, false},
  {"ikrk","Characteristic current for Kirk effect [A]", TR_DOUBLE, false},
  {"is","Saturation value for forward collector current [A]", TR_DOUBLE, false},
  {"isa","Collector current EB barrier limiting current [A]", TR_DOUBLE, false},
  {"isb","Collector current BC barrier limiting current [A]", TR_DOUBLE, false},
  {"isc","Saturation value for intrinsic bc junction current [A]", TR_DOUBLE, false},
  {"iscx","Saturation value for extrinsic bc junction current [A]", TR_DOUBLE, false},
  {"ise","Saturation value for nonideal base current [A]", TR_DOUBLE, false},
  {"isex","Saturation value for emitter leakage diode [A]", TR_DOUBLE, false},
  {"itc","Characteristic current for TFC [A]", TR_DOUBLE, false},
  {"itc2","Characteristic current for TFC [A]", TR_DOUBLE, false},
  {"kfn","BE flicker noise constant", TR_DOUBLE, false},
  {"mjc","Exponent for voltage variation of Intrinsic BC Cj", TR_DOUBLE, false},
  {"mjcx","Exponent for voltage variation of Extrinsic BC Cj", TR_DOUBLE, false},
  {"mje","Exponent for voltage variation of BE Cj", TR_DOUBLE, false},
  {"mjs","Exponent for voltage variation of CS Cj", TR_DOUBLE, false},
  {"na","Collector current EB barrier ideality factor", TR_DOUBLE, false},
  {"nb","Collector current BC barrier ideality factor", TR_DOUBLE, false},
  {"nbc","Exponent for BC multiplication factor vs voltage", TR_DOUBLE, false},
  {"nc","Ideality factor for intrinsic bc junction current", TR_DOUBLE, false},
  {"ncs","Ideality factor for collector-substrate current", TR_DOUBLE, false},
  {"ncx","Ideality factor for extrinsic bc junction current", TR_DOUBLE, false},
  {"ne","Ideality factor for nonideal forward base current", TR_DOUBLE, false},
  {"nex","Ideality factor for emitter leakage diode", TR_DOUBLE, false},
  {"nf","Forward collector current ideality factor", TR_DOUBLE, false},
  {"nr","Reverse current ideality factor", TR_DOUBLE, false},
  {"rbi","Intrinsic base resistance [ohm]", TR_DOUBLE, false},
  {"rbx","Extrinsic base resistance [ohm]", TR_DOUBLE, false},
  {"rci","Intrinsic collector resistance [ohm]", TR_DOUBLE, false},
  {"rcx","Extrinsic collector resistance [ohm]", TR_DOUBLE, false},
  {"re","Emitter resistance [ohm]", TR_DOUBLE, false},
  {"rex","Extrinsic emitter leakage diode series resistance [ohm]", TR_DOUBLE, false},
  {"rth","Thermal resistance from device to thermal ground [C/W]", TR_DOUBLE, false},
  {"tbcxs","Excess BC heterojunction transit time [s]", TR_DOUBLE, false},
  {"tbexs","Excess BE heterojunction transit time [s]", TR_DOUBLE, false},
  {"tfb","Base transit time [s]", TR_DOUBLE, false},
  {"tfc0","Collector forward transit time [s]", TR_DOUBLE, false},
  {"tkrk","Forward transit time for Kirk effect [s]", TR_DOUBLE, false},
  {"tnc","Coefficient for NC temperature dependence", TR_DOUBLE, false},
  {"tne","Coefficient for NE temperature dependence", TR_DOUBLE, false},
  {"tnex","Coefficient for NEX temperature dependence", TR_DOUBLE, false},
  {"tnom","Temperature at which model parameters are given [K]", TR_DOUBLE, false},
  {"tr","Reverse charge storage time for intrinsic BC diode [s]", TR_DOUBLE, false},
  {"trx","Charge storage time for extrinsic BC diode [s]", TR_DOUBLE, false},
  {"tre","Charge storage time[s]", TR_DOUBLE, false},
  {"tvjc","Coefficient for VJC temperature dependence [V/C]", TR_DOUBLE, false},
  {"tvjci","Coefficient for VJCI temperature dependence [V/C]", TR_DOUBLE, false},
  {"tvjcx","Coefficient for VJCX temperature dependence [V/C]", TR_DOUBLE, false},
  {"tvje","Coefficient for VJE temperature dependence [V/C]", TR_DOUBLE, false},
  {"tvjs","Coefficient for VJS temperature dependence [V/C]", TR_DOUBLE, false},
  {"vaf","Forward Early voltage [V]", TR_DOUBLE, false},
  {"var","Reverse Early voltage [V]", TR_DOUBLE, false},
  {"vjc","Intrinsic BC diode builtin potential for Cj estimation [V]", TR_DOUBLE, false},
  {"vjci","Vjci [V]", TR_DOUBLE, false},
  {"vjcx","Extrinsic BC diode builtin potential for Cj estimation [V]", TR_DOUBLE, false},
  {"vje","BE diode builtin potential for Cj estimation [V]", TR_DOUBLE, false},
  {"vjs","CS diode builtin potential for Cj estimation [V]", TR_DOUBLE, false},
  {"vkrk","Characteristic Voltage for Kirk effect [V]", TR_DOUBLE, false},
  {"vtc","Characteristic voltage for TFC [V]", TR_DOUBLE, false},
  {"xcjc","Factor for partitioning extrinsic BC Cj", TR_DOUBLE, false},
  {"xrb","Exponent for RB temperature dependence", TR_DOUBLE, false},
  {"xrc","Exponent for RC temperature dependence", TR_DOUBLE, false},
  {"xre","Exponent for RE temperature dependence", TR_DOUBLE, false},
  {"xrex","Exponent for REX temperature dependence", TR_DOUBLE, false},
  {"xrt","Exponent for RTH temperature dependence", TR_DOUBLE, false},
  {"xtb","Exponent for beta temperature dependence", TR_DOUBLE, false},
  {"xti","Exponent for IS temperature dependence", TR_DOUBLE, false},
  {"xtikrk","Exponent for IKRK temperature dependence", TR_DOUBLE, false},
  {"xtitc","Exponent for ITC temperature dependence", TR_DOUBLE, false},
  {"xtitc2","Exponent for ITC2 temperature dependence", TR_DOUBLE, false},
  {"xttf","Exponent for TF temperature dependence", TR_DOUBLE, false},
  {"xttkrk","Exponent for TKRK temperature dependence", TR_DOUBLE, false},
  {"xtvkrk","Exponent for VKRK temperature dependence", TR_DOUBLE, false},
  {"kirchhoff","Use Kirchhoff transformation flag", TR_BOOLEAN, false},
  {"b","Thermal conductivity temperature exponent", TR_DOUBLE, false},
  {"ts","Kirchhoff transformation temperature (K)", TR_DOUBLE, false}
};

HBTnpnxT::HBTnpnxT(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  //set default parameter values
  paramvalue[0] = &(type = 1);
  paramvalue[1] = &(afn = 1.5);
  paramvalue[2] = &(bf = 200000);
  paramvalue[3] = &(bfn = 1);
  paramvalue[4] = &(bkdn = false);
  paramvalue[5] = &(br = 1000);
  paramvalue[6] = &(bvc = 28);
  paramvalue[7] = &(ccmin = 9.89E-15);
  paramvalue[8] = &(cemin = 1.09E-14);
  paramvalue[9] = &(cjc = 1.4E-14);
  paramvalue[10] = &(cjcx = 1.6E-14);
  paramvalue[11] = &(cje = 1.88E-14);
  paramvalue[12] = &(cjs = 0);
  paramvalue[13] = &(cth = 3E-10);
  paramvalue[14] = &(cxmin = 3.6E-14);
  paramvalue[15] = &(dtmax = 1000);
  paramvalue[16] = &(eaa = 0.105);
  paramvalue[17] = &(eab = 0);
  paramvalue[18] = &(eac = -0.1);
  paramvalue[19] = &(eae = -0.495);
  paramvalue[20] = &(eax = 0);
  paramvalue[21] = &(eg = 0.271);
  paramvalue[22] = &(fa = 0.995);
  paramvalue[23] = &(fc = 0.8);
  paramvalue[24] = &(fce = 0.9);
  paramvalue[25] = &(fex = 1);
  paramvalue[26] = &(icrit0 = 0.035);
  paramvalue[27] = &(ics = 0);
  paramvalue[28] = &(ik = 2);
  paramvalue[29] = &(ikrk = 0.023);
  paramvalue[30] = &(is = 5E-16);
  paramvalue[31] = &(isa = 1E10);
  paramvalue[32] = &(isb = 2E-13);
  paramvalue[33] = &(isc = 2.16E-15);
  paramvalue[34] = &(iscx = 3E-9);
  paramvalue[35] = &(ise = 1.05E-17);
  paramvalue[36] = &(isex = 3E-9);
  paramvalue[37] = &(itc = 0.06);
  paramvalue[38] = &(itc2 = 0.04);
  paramvalue[39] = &(kfn = 0);
  paramvalue[40] = &(mjc = 0.5);
  paramvalue[41] = &(mjcx = 0.3);
  paramvalue[42] = &(mje = 0.507);
  paramvalue[43] = &(mjs = 0.01);
  paramvalue[44] = &(na = 10);
  paramvalue[45] = &(nb = 3);
  paramvalue[46] = &(nbc = 6);
  paramvalue[47] = &(nc = 1);
  paramvalue[48] = &(ncs = 2);
  paramvalue[49] = &(ncx = 22);
  paramvalue[50] = &(ne = 1.15);
  paramvalue[51] = &(nex = 22);
  paramvalue[52] = &(nf = 1.15);
  paramvalue[53] = &(nr = 1.03);
  paramvalue[54] = &(rbi = 10.5);
  paramvalue[55] = &(rbx = 0);
  paramvalue[56] = &(rci = 0);
  paramvalue[57] = &(rcx = 1);
  paramvalue[58] = &(re = 2.5);
  paramvalue[59] = &(rex = 0);
  paramvalue[60] = &(rth = 2200);
  paramvalue[61] = &(tbcxs = 0);
  paramvalue[62] = &(tbexs = 0);
  paramvalue[63] = &(tfb = 2.5E-13);
  paramvalue[64] = &(tfc0 = 0);
  paramvalue[65] = &(tkrk = 5.5E-14);
  paramvalue[66] = &(tnc = 0);
  paramvalue[67] = &(tne = 0);
  paramvalue[68] = &(tnex = 0);
  paramvalue[69] = &(tnom = 300);
  paramvalue[70] = &(tr = 3.5E-10);
  paramvalue[71] = &(trx = 3.5E-10);
  paramvalue[72] = &(tre = 3.5E-10);
  paramvalue[73] = &(tvjc = -0.0015);
  paramvalue[74] = &(tvjci = -0.0015);
  paramvalue[75] = &(tvjcx = -0.0015);
  paramvalue[76] = &(tvje = -0.0015);
  paramvalue[77] = &(tvjs = -0.0015);
  paramvalue[78] = &(vaf = 300);
  paramvalue[79] = &(var = 100);
  paramvalue[80] = &(vjc = 0.242);
  paramvalue[81] = &(vjci = 0.242);
  paramvalue[82] = &(vjcx = 0.35);
  paramvalue[83] = &(vje = 1);
  paramvalue[84] = &(vjs = 1.4);
  paramvalue[85] = &(vkrk = 0.25);
  paramvalue[86] = &(vtc = 0.5);
  paramvalue[87] = &(xcjc = 1);
  paramvalue[88] = &(xrb = 0.5);
  paramvalue[89] = &(xrc = 0.5);
  paramvalue[90] = &(xre = 0.5);
  paramvalue[91] = &(xrex = 0.5);
  paramvalue[92] = &(xrt = 1.2);
  paramvalue[93] = &(xtb = -2.8);
  paramvalue[94] = &(xti = 2);
  paramvalue[95] = &(xtikrk = 0.6);
  paramvalue[96] = &(xtitc = 0);
  paramvalue[97] = &(xtitc2 = 0);
  paramvalue[98] = &(xttf = 0.75);
  paramvalue[99] = &(xttkrk = 0.6);
  paramvalue[100] = &(xtvkrk = 0.6);
  paramvalue[101] = &(kirchhoff = false);
  paramvalue[102] = &(b = 1.22);
  paramvalue[103] = &(ts = 300);

  //set the number of terminals
  setNumTerms(6);

  //set Flags
  setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);

  //set number of states
  setNumberOfStates(4);
}

void HBTnpnxT::init() throw(string&)
{
  //create tape
  DenseIntVector var(4);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  var[3] = 3;
  initializeAD(var, var);
}

void HBTnpnxT::getLocalRefIdx(UnsignedVector& local_ref_vec,
    TerminalVector& term_list)
{
  bm1 = b - one;

  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1));
  term_list.push_back(getTerminal(2));
  term_list.push_back(getTerminal(3)); // Local reference terminal
  term_list.push_back(getTerminal(4));
  term_list.push_back(getTerminal(5)); // Local reference terminal

  local_ref_vec.push_back(3); // Local reference index
  local_ref_vec.push_back(5); // Local reference index
}

void HBTnpnxT :: eval(AD * x, AD * effort, AD * flow)
{
  //x[vbei, vbci, vcs, deltaT]
  //x[0] : vbei
  //x[1] : vbci
  //x[2] : vcs
  //x[3] : T - tnom

  AD KT, vTd, rTd;
  AD delta, q1, q2, qb;
  AD Icc, Ibei, Ibci, Icf, Icr, Ibcx, Ibex, Ibk, Ics;
  AD Mf, g1, Mf1, Tfbt;
  AD Vbbi, Vbe, Vbc;
  AD Qbej, Qbcj, Qfdiff, Qkrk, CjcH, Qbcf, Qbcm, Vmin;
  AD Qbe, Qbci, Qbcx, Qbcxx, Qcs;
  AD ftt, qcc, qcc1, qcc2, Icrit;
  AD temp;

//#define kBoltzman 1.3806226e-23
//#define eCharge 1.60219e-19

  temp = x[3] + tnom;//T
  rTd = temp/tnom;

  KT = kBoltzman*temp;
  vTd = eCharge/KT;

  AD rcx_T, rci_T, rbx_T, rbi_T, re_T, rex_T;
  AD bf_T, ne_T, br_T, nc_T, nex_T;
  AD is_T, isa_T, isb_T, ise_T, isc_T, isex_T;

  rcx_T = rcx*pow(rTd, xrc);
  rci_T = rci*pow(rTd, xrc);
  rbx_T = rci*pow(rTd, xrb);
  rbi_T = rbi*pow(rTd, xrb);
  re_T = re*pow(rTd, xre);
  rex_T = rex*pow(rTd, xrex);

  bf_T = bf*pow(rTd, xtb);
  ne_T = ne + tne*x[3];
  br_T = br*pow(rTd, xtb);
  nc_T = nc + tnc*x[3];
  nex_T = nex + tnex*x[3];

  is_T = is*exp(eCharge*eg/KT*(rTd - 1.0))*pow(rTd, xti);
  isa_T = isa*exp(eCharge*(eg + eaa)/KT*(rTd - 1.0))*pow(rTd, xti);
  isb_T = isb*exp(eCharge*(eg + eab)/KT*(rTd - 1.0))*pow(rTd, xti);
  ise_T = ise*pow(rTd, - xtb + xti/ne_T)*exp(eCharge*(eg + eae)
      /ne_T/KT - eCharge*(eg + eae)/ne/(kBoltzman*tnom));
  isc_T = isc*pow(rTd, - xtb + xti/nc_T)*exp(eCharge*(eg + eac)
      /nc_T/KT - eCharge*(eg + eac)/nc/(kBoltzman*tnom));
  isex_T = isex*pow(rTd, -xtb + xti/nex_T)*exp(eCharge*(eg + eax)
      /nex_T/KT - eCharge*(eg + eax)/nex/(kBoltzman*tnom));

  /****************current & voltage elements*****************/
  //delta
  q1 = 1.0/(1.0 + x[1]/vaf - x[0]/var);
  q2 = is_T/ik*(exp(eCharge*x[0]/nf/KT) - 1.0);
  qb = q1*(1.0 + pow((1.0 + 4.0*q2), 0.5))/2.0;
  delta = qb + is_T*exp(eCharge*x[0]/na/KT)/isa_T + is_T*exp(eCharge*x[1]/nb/KT)/isb_T;

  //Icf
  Icf = is_T*(exp(eCharge*x[0]/nf/KT) - 1.0)/delta;

  //Icr
  Icr = is_T*(exp(eCharge*x[1]/nr/KT) - 1.0)/delta;

  //Icc current from C to E
  Icc = Icf - Icr;

  //Ibei Current from Bi to E
  Ibei = Icf/bf_T + ise_T*(exp(eCharge*x[0]/ne_T/KT) - 1.0);

  //Ibci Current from Bi to C
  Ibci = Icr/br_T + isc_T*(exp(eCharge*x[1]/nc_T/KT) - 1.0);

  //Vbbi Voltage between B and Bi
  Vbbi = (Ibci + Ibei)*rbi_T;

  //Vbe & Vbc Voltage between B and E, C
  Vbe = Vbbi + x[0];
  Vbc = Vbbi + x[1];

  //Ibex Current from B to E
  Ibex = isex_T*(exp(eCharge*Vbe/nex_T/KT) - 1.0);

  //Ibcx Current from B to C
  Ibcx = iscx*(exp(eCharge*Vbc/ncx/KT) - 1.0);

  //Ibk Current from C node to Bi node
  if (bkdn == 0)
  {
    Ibk = 0.0;
  }
  else
  {
    if ((vTd < - x[1]) && (vTd < fa*bvc))
    {
      Mf = 1.0/(1.0 - pow( - x[1]/bvc, nbc));
    }
    else if (vTd > - x[1])
    {
      Mf = 1.0;
    }
    else if ( - x[1] > fa*bvc)
    {
      Mf1 = 1.0/(1.0 - pow(fa,nbc));
      g1 = Mf1*(Mf1 - 1.0)*nbc/fa/bvc;
      Mf = Mf1 + g1*( - x[1] - fa*bvc);
    }
    else
    {
      Mf = 0.0;
    }

    Ibk = (Mf - 1.0)*Icf;
  }

  //Ics Current from S to C
  Ics = ics*(exp(eCharge*x[2]/ncs/KT) - 1.0);

  /****************Charge Storage Elements*******************/
  AD tkrk_T, vkrk_T, ikrk_T, itc_T, itc2_T;
  tkrk_T = tkrk*pow(rTd, xttkrk);
  vkrk_T = vkrk*pow(rTd, xtvkrk);
  ikrk_T = ikrk*pow(rTd, xtikrk);
  itc_T = itc*pow(rTd, xtitc);
  itc2_T = itc2*pow(rTd, xtitc2);

  //Qkrk
  Qkrk = tkrk_T*Icf*exp(x[1]/vkrk_T + Icf/ikrk_T);
  ftt = pow(rTd, xttf);

  //Qbej
  AD vmin1;
  vmin1 =  vje*(1.0 - pow(cje/cemin, 1.0/mje));

  if ((x[0] < fce*vje) && (x[0] < vmin1))
  {
    Qbej = cemin*(x[0] - vje) + cemin*vje*mje/(mje - 1.0)*pow(cje/cemin, 1.0/mje);
  }
  else if ((x[0] < fce*vje) && (x[0] > vmin1))
  {
    Qbej = - cje*vje*pow(1.0 - x[0]/vje , 1.0 - mje)/(1.0 - mje);
  }
  else if ((x[0] > (fce*vje)) && (cje > cemin*pow(1.0 - fce, mje)))
  {
    Qbej = - cje*vje/pow(1.0 - fce, mje)*((1.0 - fce)/(1.0 - mje) + fce - x[0]/vje -
        mje*(fce - x[0]/vje)*(fce - x[0]/vje)/2.0/(1.0 - fce));
  }
  else if ((x[0] > (fce*vje)) && (cje < cemin*pow(1.0 - fce, mje)))
  {
    Qbej = cemin*(x[0] - vje) + cemin*vje*mje/(mje - 1.0)*pow(cje/cemin, 1.0/mje) +
      cje*vje*(x[0]/vje - fce)*(x[0]/vje - fce)*mje/2.0/pow(1.0 - fce, mje + 1.0);
  }
  else
  {
    Qbej = 0.0;
  }

  //Qbcj
  AD vmin2;
  vmin2 =  vjc*(1.0 - pow(cjc/ccmin, 1.0/mjc));

  if ((x[1] < fc*vjc) && (x[1] < vmin2))
  {
    Qbcj = ccmin*(x[1] - vjc) + ccmin*vjc*mjc/(mjc - 1.0)*pow(cjc/ccmin, 1.0/mjc);
  }
  else if ((x[1] < fc*vjc) && (x[1] > vmin2))
  {
    Qbcj = - cjc*vjc*pow(1.0 - x[1]/vjc , 1.0 - mjc)/(1.0 - mjc);
  }
  else if ((x[1] > (fc*vjc)) && (cjc > ccmin*pow(1.0 - fc, mjc)))
  {
    Qbcj = - cjc*vjc/pow(1.0 - fc, mjc)*((1.0 - fc)/(1.0 - mjc) + fc - x[1]/vjc -
        mjc*(fc - x[1]/vjc)*(fc - x[1]/vjc)/2.0/(1.0 - fc));
  }
  else if ((x[1] > (fc*vjc)) && (cjc < ccmin*pow(1.0 - fc, mjc)))
  {
    Qbcj = ccmin*(x[1] - vjc) + ccmin*vjc*mjc/(mjc - 1.0)*pow(cjc/ccmin, 1.0/mjc) +
      cjc*vjc*(x[1]/vjc - fc)*(x[1]/vjc - fc)*mjc/2.0/pow(1.0 - fc, mjc + 1.0);
  }
  else
  {
    Qbcj = 0.0;
  }

  qcc1 = (1.0 + (Icf/itc_T)*(Icf/itc_T));
  qcc2 = (1.0 + pow(Icf/itc2_T, 3) + (vjci - x[1])/vtc);
  qcc = qcc1/qcc2;
  Icrit = icrit0*qcc/ftt;

  //CjcH
  if ((1.0 - Icf/Icrit)<0)
  {
    CjcH = cjc*(-1.0)*pow((Icf/Icrit - 1.0), mjc);
  }
  else
  {
    CjcH = cjc*pow((1.0 - Icf/Icrit), mjc);
  }

  //Qbcf
  if (CjcH < 0)
  {
    Qbcf = ccmin*(x[1] - vjc) - ccmin*vjc*mjc/(mjc - 1.0)*pow( - CjcH/ccmin, 1.0/mjc);
  }
  else
  {
    Vmin = vjc*(1.0 - pow(CjcH/ccmin, 1.0/mjc));
    if ((x[1] < fc*vjc) && (x[1] < Vmin))
    {
      Qbcf = ccmin*(x[1] - vjc) + ccmin*vjc*mjc/(mjc - 1.0)*pow(CjcH/ccmin, 1.0/mjc);
    }
    else if ((x[1] < fc*vjc) && (x[1] > Vmin))
    {
      Qbcf = - CjcH*vjc*pow(1.0 - x[1]/vjc, 1.0 - mjc)/(1.0 - mjc);
    }
    else if ((x[1] > fc*vjc) && (CjcH > ccmin*pow(1.0 - fc, mjc)))
    {
      Qbcf = - CjcH*vjc/pow(1.0 - fc, mjc)*((1.0 - fc)/(1.0 - mjc) + fc
             - x[1]/vjc - mjc*pow(fc - x[1]/vjc, 2.0)/2.0/(1.0 - fc));
    }
    else if ((x[1] > fc*vjc) && (CjcH < ccmin*pow(1.0 - fc, mjc)))
    {
      Qbcf = ccmin*(x[1] - vjc) + ccmin*vjc*mjc/(mjc - 1.0)*pow(CjcH/ccmin, 1.0/mjc) +
             CjcH*vjc*pow(x[1]/vjc - fc, 2.0)*mjc/2.0/pow(1.0 - fc, mjc + 1.0);
    }
    else
    {
      Qbcf = 0.0;
    }
  }

  //Qbcm
  Qbcm = Qbcf - Qbcj;

  //Tfbt
  Tfbt = tfb*(1.0 + x[1]/vaf + x[0]/var) + tbexs*exp( - eCharge*(x[0] - vje)/na/KT)
         + tbcxs*exp(eCharge*(x[1] - vjc)/nb/KT);

  //Qfdiff
  Qfdiff = Icf*ftt*(Tfbt + tfc0/qcc) + Qbcm + Qkrk;

  //Qbe Charge beween Bi( + ) and E( - )
  Qbe = Qbej + (1.0 - fex)*Qfdiff;

  //Qbci Charge between Bi( + ) and C( - )
  Qbci = Qbcj + tre*Icr + fex*Qfdiff;

  //Qbcx Charge between B and C
  AD vmin3, Qbcxo;
  vmin3 =  vjcx*(1.0 - pow(cjcx/cxmin, 1.0/mjcx));

  if ((Vbc < fc*vjcx) && (Vbc < vmin3))
  {
    Qbcxo = cxmin*(Vbc - vjcx) + cxmin*vjcx*mjcx/(mjcx - 1.0)*pow(cjcx/cxmin, 1.0/mjcx);
  }
  else if ((Vbc < fc*vjcx) && (Vbc > vmin3))
  {
    Qbcxo = - cjcx*vjcx*pow(1.0 - Vbc/vjcx , 1.0 - mjcx)/(1.0 - mjcx);
  }
  else if ((Vbc > (fc*vjcx)) && (cjcx > cxmin*pow(1.0 - fc, mjcx)))
  {
    Qbcxo = - cjcx*vjcx/pow(1.0 - fc, mjcx)*((1.0 - fc)/(1.0 - mjcx) + fc - Vbc/vjcx -
        mjcx*(fc - Vbc/vjcx)*(fc - Vbc/vjcx)/2.0/(1.0 - fc));
  }
  else if ((Vbc > (fc*vjcx)) && (cjcx < cxmin*pow(1.0 - fc, mjcx)))
  {
    Qbcxo = cxmin*(Vbc - vjcx) + cxmin*vjcx*mjcx/(mjcx - 1.0)*pow(cjcx/cxmin, 1.0/mjcx) +
      cjcx*vjcx*(Vbc/vjcx - fc)*(Vbc/vjcx - fc)*mjcx/2.0/pow(1.0 - fc, mjcx + 1.0);
  }
  else
  {
    Qbcxo = 0.0;
  }

  Qbcx = trx*Ibcx + xcjc*Qbcxo;

  //Qbcxx Charge between B and C
  AD vmin4, Qbcxxo;
  vmin4 =  vjcx*(1.0 - pow(cjcx/cxmin, 1.0/mjcx));

  if ((Vbc < fc*vjcx) && (Vbc < vmin4))
  {
    Qbcxxo = cxmin*(Vbc - vjcx) + cxmin*vjcx*mjcx/(mjcx - 1.0)*pow(cjcx/cxmin, 1.0/mjcx);
  }
  else if ((Vbc < fc*vjcx) && (Vbc > vmin4))
  {
    Qbcxxo = - cjcx*vjcx*pow(1.0 - Vbc/vjcx , 1.0 - mjcx)/(1.0 - mjcx);
  }
  else if ((Vbc > (fc*vjcx)) && (cjcx > cxmin*pow(1.0 - fc, mjcx)))
  {
    Qbcxxo = - cjcx*vjcx/pow(1.0 - fc, mjcx)*((1.0 - fc)/(1.0 - mjcx) + fc - Vbc/vjcx -
        mjcx*(fc - Vbc/vjcx)*(fc - Vbc/vjcx)/2.0/(1.0 - fc));
  }
  else if ((Vbc > (fc*vjcx)) && (cjcx < cxmin*pow(1.0 - fc, mjcx)))
  {
    Qbcxxo = cxmin*(Vbc - vjcx) + cxmin*vjcx*mjcx/(mjcx - 1.0)*pow(cjcx/cxmin, 1.0/mjcx) +
      cjcx*vjcx*(Vbc/vjcx - fc)*(Vbc/vjcx - fc)*mjcx/2.0/pow(1.0 - fc, mjcx + 1.0);
  }
  else
  {
    Qbcxxo = 0.0;
  }

  Qbcxx = (1.0 - xcjc)*Qbcxxo;

  //Qcs Charge between C and S
  if (x[2] > - fc*vjs)
  {
    Qcs = - cjs*vjs*pow(1.0 + x[2]/vjs, 1.0 - mjs)/(1.0 - mjs);
  }
  else
  {
    Qcs = - cjs*vjs/pow(1.0 - fc, mjs)*((1.0 - fc)/(1.0 - mjs) + fc
        + x[2]/vjs - mjs/2.0/(1.0 - fc)*pow(fc + x[2]/vjs, 2.0));
  }

  /*****************************************************/
  AD Qcollector = - Qbci - Qbcx - Qbcxx + Qcs;
  AD Qbase = Qbe + Qbci + Qbcx + Qbcxx;
  AD Qbi = Qbe + Qbci;

  // Calculate dynamic current contributions due to charge
  // iqc is dynamic contribution to collector current
  // iqb is dynamic contribution to base current
  // iqbi is dynamic contribution to Bi(+) current
  AD iqc, iqb, iqbi;
  iqc = Qcollector.fastAccessDx(0)*x[4] + Qcollector.fastAccessDx(1)*x[5]
        + Qcollector.fastAccessDx(2)*x[6];
  iqb = Qbase.fastAccessDx(0)*x[4] + Qbase.fastAccessDx(1)*x[5] + Qbase.fastAccessDx(2)*x[6];
  iqbi = Qbi.fastAccessDx(0)*x[4] + Qbi.fastAccessDx(1)*x[5] + Qbi.fastAccessDx(2)*x[6];

  //effort[0] : Vc flow[0] : Ic
  //effort[1] : Vb flow[1] : Ib
  //effort[2] : Ve flow[2] : Ie
  //effort[3] : T  flow[3] : pout

  flow[0] = Icc - Ibci - Ibcx - Ics + Ibk + iqc; //Ic
  flow[1] = Ibci + Ibei + Ibcx + Ibex + iqb; //Ib
  flow[2] = - Icc - Ibei - Ibex - (iqc + iqb); //Ie

  AD rtd;

  effort[3] = x[3] + tnom; // Temperature
  // Apply Kirchhoff transformation (if requested)
  if (kirchhoff)
    effort[3] = ts/(bm1)*(b - pow(ts/effort[3],bm1));

  rtd = effort[3]/tnom;
  rcx_T = rcx*pow(rtd, xrc);
  rbi_T = rbi*pow(rtd, xrb);
  re_T = re*pow(rtd, xre);

  effort[0] = x[2] + rcx_T*flow[0]; //Vc
  effort[1] = Vbc + x[2] + rbi_T*iqbi; //Vb
  effort[2] = x[1] - x[0] + x[2] + re_T*flow[2]; //Ve

  flow[3] = -(flow[0]*effort[0] + flow[1]*effort[1]); //Pdiss
}
