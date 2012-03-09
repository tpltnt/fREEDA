#include "Mosnbsim4.h"

// Static members
const unsigned Mosnbsim4::n_par = 144;

// Element information
ItemInfo Mosnbsim4::einfo =
{
  "mosnbsim4",
  "MOSFET using BSIM4 model",
  "Nikhil Kriplani",
  DEFAULT_ADDRESS"transistor>mosfet"
  "2002_04_10"
};

//parameters
ParmInfo Mosnbsim4::pinfo[] =
{
  {"TOXE", "Electrical gate equivalent oxide thickness", TR_DOUBLE, false},
  {"TOXP", "Physical gate equivalent oxide thickness", TR_DOUBLE, false},
  {"EPSROX", "Gate dielectric constant relative to vacuum", TR_DOUBLE, false},
  {"VFB", "Flat-band voltage", TR_DOUBLE, false},
  {"VTH0", "Long-channel threshold voltage", TR_DOUBLE, false},
  {"NGATE", "Poly Si gate doping concentration", TR_DOUBLE, false},
  {"XL", "Channel length offset due to mask/etch effect", TR_DOUBLE, false},
  {"XW", "Channel width offset due to mask/etch effect", TR_DOUBLE, false},
  {"NF", "Number of device fingers", TR_DOUBLE, false},
  {"W", "Width of the device", TR_DOUBLE, false},
  {"L", "Length of the device", TR_DOUBLE, false},
  {"DWG", "Coefficient of gate bias dependence of Weff", TR_DOUBLE, false},
  {"DWB", "Coefficient of body bias dependence of Weff", TR_DOUBLE, false},
  {"WINT", "Channel-width offset parameter", TR_DOUBLE, false},
  {"WLN", "Power of length dependence of width offset", TR_DOUBLE, false},
  {"WL", "Coefficient of length dependence for width offset", TR_DOUBLE, false},
  {"WWN", "Power of width dependence of width offset", TR_DOUBLE, false},
  {"WW", "Coefficient of width dependence for width offset", TR_DOUBLE, false},
  {"WWL", "Coefficient of length and width cross term dependence for width offset", TR_DOUBLE, false},
  {"LINT", "Channel-length offset parameter", TR_DOUBLE, false},
  {"LLN", "Power of length dependence for length offset", TR_DOUBLE, false},
  {"LL", "Coefficient of length dependence for length offset", TR_DOUBLE, false},
  {"LW", "Coefficient of width dependence for length offset", TR_DOUBLE, false},
  {"LWN", "Power of width dependence for length offset", TR_DOUBLE, false},
  {"LWL", "Coefficient of length and width cross term dependence for length offset", TR_DOUBLE, false},
  {"K1", "First-order body bias coefficient", TR_DOUBLE, false},
  {"K2", "Second-order body bias coefficient", TR_DOUBLE, false},
  {"LPEB", "Lateral non-uniform doping effect on K1", TR_DOUBLE, false},
  {"LPE0", "Lateral non-uniform doping parameter at Vbs=0", TR_DOUBLE, false},
  {"K3", "Narrow width coefficient", TR_DOUBLE, false},
  {"K3B", "Body effect coefficient of K3", TR_DOUBLE, false},
  {"W0", "Narrow width parameter", TR_DOUBLE, false},
  {"DVT0W", "First coefficient of narrow width effect on Vth for small channel length", TR_DOUBLE, false},
  {"DVT0", "First coefficient of short-channel effect on Vth", TR_DOUBLE, false},
  {"DVT1W", "Second coefficient of narrow width effect on Vth for small channel length", TR_DOUBLE, false},
  {"DVT1", "Second coefficient of short-channel effect on Vth", TR_DOUBLE, false},
  {"DSUB", "DIBL coefficient exponent in sub-threshold region", TR_DOUBLE, false},
  {"ETA0", "DIBL coefficient in sub-threshold region", TR_DOUBLE, false},
  {"ETAB", "Body-bias coefficient for the sub-threshold region", TR_DOUBLE, false},
  {"TOXM", "Tox at which parameters are extracted", TR_DOUBLE, false},
  {"T", "Temperature", TR_DOUBLE, false},
  {"NDEP", "Channel doping concentration at depletion edge for zero body bias", TR_DOUBLE, false},
  {"PHIN", "Non-uniform vertical doping effect on surface potential", TR_DOUBLE, false},
  {"VBM", "Maximum applied body bias in VTH0 calculation", TR_DOUBLE, false},
  {"NSUB", "Substrate doping concentration", TR_DOUBLE, false},
  {"DVT2W", "Body-bias coefficient of narrow width effect for small channel length", TR_DOUBLE, false},
  {"NSD", "Source/drain doping concentration fatal error if not positive", TR_DOUBLE, false},
  {"DVT2", "Body-bias coefficient of short-channel effect on Vth", TR_DOUBLE, false},
  {"MINV", "Vgsteff fitting parameter for moderate inversion condition", TR_DOUBLE, false},
  {"NFACTOR", "Subthreshold swing factor", TR_DOUBLE, false},
  {"CDSC", "Coupling capacitance between source/drain and channel", TR_DOUBLE, false},
  {"CDSCD", "Drain-bias sensitivity of CDSC", TR_DOUBLE, false},
  {"CDSCB", "Body-bias sensitivity of CDSC", TR_DOUBLE, false},
  {"CIT", "Interface trap capacitance", TR_DOUBLE, false},
  {"KETA", "Body-bias coefficient of bulk charge effect", TR_DOUBLE, false},
  {"B0", "Bulk charge effect coefficient for channel width", TR_DOUBLE, false},
  {"B1", "Bulk charge effect width offset", TR_DOUBLE, false},
  {"A0", "Coefficient of channel-length dependence bulk charge effect", TR_DOUBLE, false},
  {"AGS", "Coefficient of Vgs dependence of bulk charge effect", TR_DOUBLE, false},
  {"XJ", "S/D junction depth", TR_DOUBLE, false},
  {"U0", "Low-field mobility", TR_DOUBLE, false},
  {"UA", "Coefficient of first-order mobility degradation due to vertical field", TR_DOUBLE, false},
  {"UB", "Coefficient of second-order mobility degradation due to vertical field", TR_DOUBLE, false},
  {"UC", "Coefficient of mobility degradation due to body-bias effect", TR_DOUBLE, false},
  {"EU", "Exponent for mobility degradation", TR_DOUBLE, false},
  {"DELTA", "Parameter for DC Vdseff", TR_DOUBLE, false},
  {"PDITS", "Impact of drain-induced Vth shift on Rout", TR_DOUBLE, false},
  {"FPROUT", "Effect of pocket implant on Rout degradation", TR_DOUBLE, false},
  {"PDITSL", "Channel-length dependence of drain-induced Vth shift for Rout", TR_DOUBLE, false},
  {"PDITSD", "Vds dependence of drain-induced Vth shift for Rout", TR_DOUBLE, false},
  {"PSCBE2", "Second substrate current induced body-effect parameter", TR_DOUBLE, false},
  {"PSCBE1", "First substrate current induced body-effect parameter", TR_DOUBLE, false},
  {"PDIBLCB", "Body bias coefficient of DIBL effect on Rout", TR_DOUBLE, false},
  {"PVAG", "Gate-bias dependence of Early voltage", TR_DOUBLE, false},
  {"PDIBL1", "Parameter for DIBL effect on Rout", TR_DOUBLE, false},
  {"PDIBL2", "Parameter for DIBL effect on Rout", TR_DOUBLE, false},
  {"DROUT", "Channel-length dependence of DIBL effect on Rout", TR_DOUBLE, false},
  {"PCLM", "Channel length modulation parameter", TR_DOUBLE, false},
  {"A1", "First non-saturation effect parameter", TR_DOUBLE, false},
  {"A2", "Second non-saturation factor", TR_DOUBLE, false},
  {"RDSWMIN", "Lightly-doped drain resistance per unit width at high Vgs and zero Vbs", TR_DOUBLE, false},
  {"RDSW", "Zero bias lightly-doped drain resistance per unit width", TR_DOUBLE, false},
  {"PRWG", "Gate-bias dependence of LDD resistance", TR_DOUBLE, false},
  {"PRWB", "Body-bias dependence of LDD resistance", TR_DOUBLE, false},
  {"WR", "Channel-width dependence parameter of LDD resistance", TR_DOUBLE, false},
  {"WLC", "Coefficient of length dependence for CV channel width offset", TR_DOUBLE, false},
  {"WWC", "Coefficient of width dependence for CV channel width offset", TR_DOUBLE, false},
  {"WWLC", "Coefficient of length and width crossterm dependence for CV channel width offset", TR_DOUBLE, false},
  {"DWJ", "Offset of the S/D junction width", TR_DOUBLE, false},
  {"CLC", "Constant term for the short channel model", TR_DOUBLE, false},
  {"CLE", "Exponential term for the short channel model", TR_DOUBLE, false},
  {"NOFF", "CV parameter in VgsteffCV for weak to strong inversion", TR_DOUBLE, false},
  {"VOFFCV", "CV parameter in VgsteffCV for weak to strong inversion", TR_DOUBLE, false},
  {"CF", "Fringing field capacitance", TR_DOUBLE, false},
  {"CKAPPAD", "Coefficient of bias-dependent overlap capacitance for the drain side", TR_DOUBLE, false},
  {"CKAPPAS", "Coefficient of bias-dependent overlap capacitance for the source side", TR_DOUBLE, false},
  {"LLC", "Coefficient of length dependence on CV channel length offset", TR_DOUBLE, false},
  {"LWC", "Coefficient of width dependence on CV channel length offset", TR_DOUBLE, false},
  {"LWLC", "Coefficient of length and width cross term dependence on CV channel length offset", TR_DOUBLE, false},
  {"WWLC", "Coefficient of length and width cross term dependence on CV channel width offset", TR_DOUBLE, false},
  {"VOFF", "Offset voltage in the subthreshold region for large W and L", TR_DOUBLE, false},
  {"VOFFL", "Channel length dependence of VOFF", TR_DOUBLE, false},
  {"POXEDGE", "Factor for the gate oxide thickness in the S/D overlap regions", TR_DOUBLE, false},
  {"TOXREF", "Nominal gate oxide thickness for gate dielectric tunneling current model", TR_DOUBLE, false},
  {"NTOX", "Exponent for the gate oxide ratio", TR_DOUBLE, false},
  {"DLCIG", "Source/drain overlap length for Igs and Igd", TR_DOUBLE, false},
  {"AIGSD", "parameter for Igs and Igd", TR_DOUBLE, false},
  {"BIGSD", "parameter for Igs and Igd", TR_DOUBLE, false},
  {"CIGSD", "parameter for Igs and Igd", TR_DOUBLE, false},
  {"MOIN", "Coefficient for gate-bias dependent surface potential", TR_DOUBLE, false},
  {"VSAT", "Saturation velocity", TR_DOUBLE, false},
  {"PDITSD", "Vds dependence of drain-induced Vth shift for Rout", TR_DOUBLE, false},
  {"AIGC", "Parameter for Igcs and Igcd", TR_DOUBLE, false},
  {"BIGC", "Parameter for Igcs and Igcd", TR_DOUBLE, false},
  {"CIGC", "Parameter for Igcs and Igcd", TR_DOUBLE, false},
  {"NIGC", "Parameter for Igcs, Igcd, Igs and Igd", TR_DOUBLE, false},
  {"PIGCD", "Vds dependence of Igcs and Igcd", TR_DOUBLE, false},
  {"DVTP0", "First coefficient of drain induced Vth shift due to long channel pocket devices", TR_DOUBLE, false},
  {"DVTP1", "First coefficient of drain induced Vth shift due to long channel pocket devices", TR_DOUBLE, false},
  {"PRT", "Temperature coefficient for RDSW", TR_DOUBLE, false},
  {"AT", "Temperature cpefficient for saturation velocity", TR_DOUBLE, false},
  {"XT", "Doping Depth", TR_DOUBLE, false},
  {"ALPHA0", "First parameter of impact ionization current", TR_DOUBLE, false},
  {"ALPHA1", "Isub parameter for length scaling", TR_DOUBLE, false},
  {"BETA0", "Second parameter of impact ionization current", TR_DOUBLE, false},
  {"AGIDL", "Pre-exponential coefficient for GIDL", TR_DOUBLE, false},
  {"BGIDL", "Exponential coefficient for GIDL", TR_DOUBLE, false},
  {"CGIDL", "Parameter for body-bias effect on GIDL", TR_DOUBLE, false},
  {"EGIDL", "Fitting parameter for band-bending for GIDL", TR_DOUBLE, false},
  {"ACDE", "Exponential coefficient for charge thickness", TR_DOUBLE, false},
  {"DLC", "Channel length offset parameter", TR_DOUBLE, false},
  {"DWC", "Channel width offset parameter", TR_DOUBLE, false},
  {"AIGBACC", "Parameter for Igb in accumulation", TR_DOUBLE, false},
  {"BIGBACC", "Parameter for Igb in accumulation", TR_DOUBLE, false},
  {"CIGBACC", "Parameter for Igb in accumulation", TR_DOUBLE, false},
  {"NIGBACC", "Parameter for Igb in accumulation", TR_DOUBLE, false},
  {"AIGBINV", "Parameter for Igb in inversion", TR_DOUBLE, false},
  {"BIGBINV", "Parameter for Igb in inversion", TR_DOUBLE, false},
  {"CIGBINV", "Parameter for Igb in inversion", TR_DOUBLE, false},
  {"EIGBINV", "Parameter for Igb in inversion", TR_DOUBLE, false},
  {"NIGBINV", "Parameter for Igb in inversion", TR_DOUBLE, false},
  {"KT1", "Temperature coefficient for threshold voltage", TR_DOUBLE, false},
  {"KT1l", "Channel Length dependence of KT1", TR_DOUBLE, false},
  {"KT2", "Body-bias coefficient of Vth with temperature effects", TR_DOUBLE, false}
};


Mosnbsim4::Mosnbsim4(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(TOXE = 3.0e-9);
  paramvalue[1] = &(TOXP = TOXE);
  paramvalue[2] = &(EPSROX = 3.9);
  paramvalue[3] = &(VFB = -1.0);
  paramvalue[4] = &(VTH0 = 0.7);
  paramvalue[5] = &(NGATE = 0.0);
  paramvalue[6] = &(XL = 0.0);
  paramvalue[7] = &(XW = 0.0);
  paramvalue[8] = &(NF = 1.0);
  paramvalue[9] = &(W = 5.0e-6);
  paramvalue[10] = &(L = 5.0e-6);
  paramvalue[11] = &(DWG = 0.0);
  paramvalue[12] = &(DWB = 0.0);
  paramvalue[13] = &(WINT = 0.0);
  paramvalue[14] = &(WLN = 1.0);
  paramvalue[15] = &(WL = 0.0);
  paramvalue[16] = &(WWN = 1.0);
  paramvalue[17] = &(WW = 0.0);
  paramvalue[18] = &(WWL = 0.0);
  paramvalue[19] = &(LINT = 0.0);
  paramvalue[20] = &(LLN = 1.0);
  paramvalue[21] = &(LL = 0.0);
  paramvalue[22] = &(LW = 0.0);
  paramvalue[23] = &(LWN = 1.0);
  paramvalue[24] = &(LWL = 0.0);
  paramvalue[25] = &(K1 = 0.53);
  paramvalue[26] = &(K2 = -0.0186);
  paramvalue[27] = &(LPEB = 0.0);
  paramvalue[28] = &(LPE0 = 1.74e-7);
  paramvalue[29] = &(K3 = 80.0);
  paramvalue[30] = &(K3B = 0.0);
  paramvalue[31] = &(W0 = 2.5e-6);
  paramvalue[32] = &(DVT0W = 0.0);
  paramvalue[33] = &(DVT0 = 2.2);
  paramvalue[34] = &(DVT1W = 5.3e6);
  paramvalue[35] = &(DVT1 = 0.53);
  paramvalue[36] = &(DSUB = 0.56);
  paramvalue[37] = &(ETA0 = 0.08);
  paramvalue[38] = &(ETAB = -0.07);
  paramvalue[39] = &(TOXM = TOXE);
  paramvalue[40] = &(T = 300.0);
  paramvalue[41] = &(NDEP = 1.7e17);
  paramvalue[42] = &(PHIN = 0.0);
  paramvalue[43] = &(VBM = -3.0);
  paramvalue[44] = &(NSUB = 6.0e16);
  paramvalue[45] = &(DVT2W = -0.032);
  paramvalue[46] = &(NSD = 1.0e20);
  paramvalue[47] = &(DVT2 = -0.032);
  paramvalue[48] = &(MINV = 0.0);
  paramvalue[49] = &(NFACTOR = 1.0);
  paramvalue[50] = &(CDSC = 2.4e-4);
  paramvalue[51] = &(CDSCD = 0.0);
  paramvalue[52] = &(CDSCB = 0.0);
  paramvalue[53] = &(CIT = 0.0);
  paramvalue[54] = &(KETA = -0.047);
  paramvalue[55] = &(B0 = 0.0);
  paramvalue[56] = &(B1 = 0.0);
  paramvalue[57] = &(A0 = 1.0);
  paramvalue[58] = &(AGS = 0.0);
  paramvalue[59] = &(XJ = 1.5e-7);
  paramvalue[60] = &(U0 = 0.067);
  paramvalue[61] = &(UA = 1.0e-15);
  paramvalue[62] = &(UB = 1.0e-19);
  paramvalue[63] = &(UC = -0.0465e-9);
  paramvalue[64] = &(EU = 1.67);
  paramvalue[65] = &(DELTA = 0.01);
  paramvalue[66] = &(PDITS = 0.0);
  paramvalue[67] = &(FPROUT = 0.0);
  paramvalue[68] = &(PDITSL = 0.0);
  paramvalue[69] = &(PDITSD = 0.0);
  paramvalue[70] = &(PSCBE2 = 1.0e-5);
  paramvalue[71] = &(PSCBE1 = 4.24e8);
  paramvalue[72] = &(PDIBLCB = 0.0);
  paramvalue[73] = &(PVAG = 0.0);
  paramvalue[74] = &(PDIBL1 = 0.0);
  paramvalue[75] = &(PDIBL2 = 0.0);
  paramvalue[76] = &(DROUT = 0.56);
  paramvalue[77] = &(PCLM = 1.3);
  paramvalue[78] = &(A1 = 0.0);
  paramvalue[79] = &(A2 = 1.0);
  paramvalue[80] = &(RDSWMIN = 0.0);
  paramvalue[81] = &(RDSW = 200.0);
  paramvalue[82] = &(PRWG = 1.0);
  paramvalue[83] = &(PRWB = 0.0);
  paramvalue[84] = &(WR = 1.0);
  paramvalue[85] = &(WLC = WL);
  paramvalue[86] = &(WWC = WW);
  paramvalue[87] = &(WWLC = WWL);
  paramvalue[88] = &(DWJ = WINT);
  paramvalue[89] = &(CLC = 1.0e-7);
  paramvalue[90] = &(CLE = 0.6);
  paramvalue[91] = &(NOFF = 1.0);
  paramvalue[92] = &(VOFFCV = 0.0);
  paramvalue[93] = &(CF = 1.08e-10);
  paramvalue[94] = &(CKAPPAD = 0.6);
  paramvalue[95] = &(CKAPPAS = 0.6);
  paramvalue[96] = &(LLC = 0.0);
  paramvalue[97] = &(LWC = 0.0);
  paramvalue[98] = &(LWLC = 0.0);
  paramvalue[99] = &(WWLC = 0.0);
  paramvalue[100] = &(VOFF = -0.08);
  paramvalue[101] = &(VOFFL = 0.0);
  paramvalue[102] = &(POXEDGE = 1.0);
  paramvalue[103] = &(TOXREF = TOXE);
  paramvalue[104] = &(NTOX = 1.0);
  paramvalue[105] = &(DLCIG = LINT);
  paramvalue[106] = &(AIGSD = 0.43);
  paramvalue[107] = &(BIGSD = 0.054);
  paramvalue[108] = &(CIGSD = 0.075);
  paramvalue[109] = &(MOIN = 15.0);
  paramvalue[110] = &(VSAT = 8.0e4);
  paramvalue[111] = &(PDITSD = 0.0);
  paramvalue[112] = &(AIGC = 0.43);
  paramvalue[113] = &(BIGC = 0.054);
  paramvalue[114] = &(CIGC = 0.075);
  paramvalue[115] = &(NIGC = 1.0);
  paramvalue[116] = &(PIGCD = 1.0);
  paramvalue[117] = &(DVTP0 = 0.0);
  paramvalue[118] = &(DVTP1 = 0.0);
  paramvalue[119] = &(PRT = 0.0);
  paramvalue[120] = &(AT = 3.3e4);
  paramvalue[121] = &(XT = 1.55e-7);
  paramvalue[122] = &(ALPHA0 = 0.0);
  paramvalue[123] = &(ALPHA1 = 0.0);
  paramvalue[124] = &(BETA0 = 30.0);
  paramvalue[125] = &(AGIDL = 0.0);
  paramvalue[126] = &(BGIDL = 2.3e9);
  paramvalue[127] = &(CGIDL = 0.5);
  paramvalue[128] = &(EGIDL = 0.8);
  paramvalue[129] = &(ACDE = 1.0);
  paramvalue[130] = &(DLC = LINT);
  paramvalue[131] = &(DWC = WINT);
  paramvalue[132] = &(AIGBACC = 0.43);
  paramvalue[133] = &(BIGBACC = 0.054);
  paramvalue[134] = &(CIGBACC = 0.075);
  paramvalue[135] = &(NIGBACC = 1.0);
  paramvalue[136] = &(AIGBINV = 0.35);
  paramvalue[137] = &(BIGBINV = 0.03);
  paramvalue[138] = &(CIGBINV = 0.006);
  paramvalue[139] = &(EIGBINV = 1.1);
  paramvalue[140] = &(NIGBINV = 3.0);
  paramvalue[141] = &(KT1 = -0.11);
  paramvalue[142] = &(KT1l = 0.0);
  paramvalue[143] = &(KT2 = 0.022);

  // Set the number of terminals
  setNumTerms(4);

  //Set number of state variables
  setNumberOfStates(3);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void Mosnbsim4::init() throw(string&)
{
  DenseIntVector var(3);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  initializeAD(var, var);
}

void Mosnbsim4::eval(AD * x, AD * effort, AD * flow)
{
  //state variables
  //vds x[0], vgs x[1], vbs x[2]

  AD lt, Vbseff, Cdsc_term, n, Theta0, Vbc, Phis, Weff;
  AD ltw, Vth, lambda, ueff, Esat, Xdep, sqrtPhis;
  AD Vgse, Vgsteff, Vb, Abulk, F_doping, Igidl;
  AD Rds, Vdsat, Vdseff, Vasat, Ids0, Ids, theta_rout;
  AD Vadits, Vadibl, Cclm, Vaclm, Va, Vascbe;
  AD T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, Isub;
  double t0, t1, t2, t3;

	const double MAX_EXP = 5.834617425e14;
	const double MIN_EXP = 1.713908431e-15;
	const double EXP_THRESHOLD = 34.0;

  double k = kBoltzman;
  double q = eCharge;

  double vt = k*T/q;

  double Eg0 = 1.16 - (7.02e-4 * T * T / (T + 1108.0));

  double ni = 1.45e10 * (T/300.15) * sqrt(T/300.15) * exp(21.5565981 - Eg0/(2.0 * vt));

  double e0 = epsilon0;
  double esi = 11.7 * e0;

  double coxe = EPSROX * (e0/TOXE);
  double coxp = EPSROX * (e0/TOXP);

  double Lnew = L + XL;
  double Wnew = W/NF + XW;
  t0 = pow(Lnew, LLN);
  t1 = pow(Wnew, LWN);
  double DL = LINT + LL/t0 + LW/t1 + LWL/(t0*t1);
  DLC += LLC/t0 + LWC/t1 + LWLC/(t0*t1);
  DLCIG += LLC/t0 + LWC/t1 + LWLC/(t0*t1);

  t2 = pow(Lnew, WLN);
  t3 = pow(Wnew,WWN);
  double DW = WINT + WL/t2 + WW/t3 + WWL/(t2*t3);
  DWC += WLC/t2 + WWC/t3 + WWLC/(t2*t3);
  DWJ += WLC/t2 + WWC/t3 + WWLC/(t2*t3);

  double Leff = Lnew - 2.0 * DL;

  Weff = Wnew - 2.0 * DW;

  double Weffcj = Wnew - 2.0 * DWJ;

  double phi_s = 0.4 + (k*T/q)*log(NDEP/ni) + PHIN;

  double m = 0.5 + atan(MINV)/pi;

  double Voff_prime = VOFF + VOFFL/Leff;

  double Vbi = (k*T/q)*log(NDEP*NSD/pow(ni,2.0));

  double Xdep0 = sqrt(2.0 * esi / (eCharge * NDEP * 1.0e6)) * sqrt(phi_s);

  double lt0 = sqrt(esi * TOXE * Xdep0 / (e0*EPSROX));

  double litl = sqrt(3.0 * XJ * TOXE);

  double gamma1 = 5.753e-12 * sqrt(NDEP)/coxe;
  double gamma2 = 5.753e-12 * sqrt(NSUB)/coxe;

  double VBX = phi_s - 7.7348e-4 * NDEP * XT * XT;
  if (VBX > 0.0)
    VBX = -VBX;

  t0 = gamma1 - gamma2;
  t1 = sqrt(phi_s - VBX) - sqrt(phi_s);
  t2 = sqrt(phi_s * (phi_s - VBM)) - phi_s;

  K2 = t0 * t1 / (2.0 * t2 + VBM);
  K1 = gamma2 - 2.0 * K2 * sqrt(phi_s - VBM);

  VTH0 = VFB + phi_s + K1 * sqrt(phi_s);

  double K1ox = K1 * TOXE / TOXM;
  double K2ox = K2 * TOXE / TOXM;

  //---------Calculation for Vbc--------
  if (K2 < 0.0)
	{
		Vbc = 0.9 * (phi_s - pow(0.5 * K1/K2,2.0));
		if (Vbc > -3.0)
			Vbc = -3.0;
		else if (Vbc < -30.0)
			Vbc = -30.0;
	}
  else
    Vbc = -30.0;
  //--------end Vbc-----------------

  //------Calculation for Vbseff-----------
  T0 = x[2] - Vbc - 0.001;
  T1 = sqrt(T0*T0 - 0.004 * Vbc);

  if (T0 > 0.0)
    Vbseff = Vbc + 0.5 * (T0 + T1);
  else
    Vbseff = Vbc * (1.0 + (-0.002 / (T1 - T0)));

  if (-Vbseff + x[2] > 0.0)
    Vbseff = x[2];
  //------end Vbseff----------------------

  //-----Calculation for Phis---------------
  if (Vbseff > 0.0)
    Phis = phi_s * phi_s / (phi_s + Vbseff);
  else
    Phis = phi_s - Vbseff;

  if (Vbseff > 0.0)
    sqrtPhis = sqrt(phi_s) * phi_s / (phi_s + 0.5*Vbseff);
  else
    sqrtPhis = sqrt(Phis);
  //------end Phis---------------------------

  //-----Calculation for Xdep---------------
  Xdep = Xdep0 * sqrtPhis / sqrt(phi_s);
  //-----end Xdep--------------------------

  //-------Calculation for Vth------------
  if (DVT2W*Vbseff + 0.5 > 0.0)
    T1 = 1.0 + DVT2W * Vbseff;
  else
    T1 = (1.0 + 3.0 * DVT2W*Vbseff) / (3.0 + 8.0 * DVT2W * Vbseff);

  ltw = sqrt(Xdep) * T1 * sqrt(esi/(EPSROX * e0) * TOXE);

  if (DVT2*Vbseff + 0.5 > 0.0)
    T1 = 1.0 + DVT2*Vbseff;
  else
    T1 = (1.0 + 3.0 * DVT2*Vbseff) / (3.0 + 8.0 * DVT2*Vbseff);

  lt = sqrt(Xdep) * T1 * sqrt(esi/(EPSROX * e0) * TOXE);

  if (-(DVT1*Leff/lt) + EXP_THRESHOLD > 0.0)
    Theta0 = exp(DVT1 * Leff/lt)/((exp(DVT1*Leff/lt) - 1.0) * (exp(DVT1*Leff/lt) - 1.0) + 2.0*exp(DVT1*Leff/lt)*MIN_EXP);
  else
    Theta0 = 1.0/(MAX_EXP - 2.0);

  if (-(DVT1W*Leff*Weff/ltw) + EXP_THRESHOLD > 0.0)
    T5 = exp(DVT1W*Leff*Weff/ltw)/((exp(DVT1W*Leff*Weff/ltw) - 1.0) * (exp(DVT1W*Leff*Weff/ltw) - 1.0) + 2.0*exp(DVT1W*Leff*Weff/ltw)*MIN_EXP);
  else
    T5 = 1.0/(MAX_EXP - 2.0);

  T0 = DVT0W * T5;
  T2 = (Vbi - phi_s) * T0;

  T0 = sqrt(1.0 + LPE0/Leff);
  T1 = K1ox * (T1 - 1.0) * sqrtPhis + (KT1 + KT1l/Leff + KT2*Vbseff) * (T/300.0 - 1.0);

  T8 = TOXE * phi_s /(Weff + W0);

  T6 = sqrt(esi / (EPSROX * e0) * TOXE * Xdep0);
  T0 = DSUB * Leff / T6;

  if (-T0+EXP_THRESHOLD > 0.0)
    T5 = exp(T0)/((exp(T0)-1.0)*(exp(T0)-1.0) + 2.0*exp(T0)*MIN_EXP);
  else
    T5 = 1.0/(MAX_EXP - 2.0);

  T3 = ETA0 + ETAB * Vbseff;

  if (-T3+1.0e-4 > 0.0)
    T3 = (2.0e-4 - T3) / (3.0 - 2.0e4 * T3);
  else
    T3 = ETA0 + ETAB * Vbseff;

  T7 = T3 * T5 * x[0];

  Vth = VTH0
	+ (K1ox * sqrtPhis - K1 * sqrtPhis) * sqrt(1.0 + LPEB/Leff)
	- K2ox * Vbseff
	- Theta0 * DVT0 * (Vbi - phi_s)
	- T2
	+ (K3 + K3B * Vbseff) * TOXE * phi_s /(Weff + W0)
	+ T1
	- T7;
  //--------------end Vth-------------------------


  //------Calculation for n-------------
  Cdsc_term = (CDSC + CDSCD * x[0] + CDSCB * Vbseff) * Theta0;
  T1 = (NFACTOR*esi/Xdep + Cdsc_term + CIT) / coxe;

  if (T1+0.5 > 0.0)
    n = 1.0 + T1;
  else
    n = (1.0 + 3.0 * T1) / (3.0 + 8.0 * T1);
  //------end n-------------------------

  //----Vth correction for pocket implant
  if (DVTP0 > 0.0)
	{
		T0 = -DVTP1 * x[0];

		if (T0 < -EXP_THRESHOLD)
			T2 = MIN_EXP;
		else
			T2 = exp(T0);

		T3 = Leff + DVTP0 * (1.0 + T2);
		T4 = vt * log(Leff/T3);
		Vth -= n * T4;
	}


  //----------Calculation for Vgse-------------
  T0 = VFB + phi_s;

  if((NGATE > 1.0e18) && (NGATE < 1.0e25) && (x[1] > T0))
	{
		T1 = 1.0e6 * eCharge * esi * NGATE / (coxe * coxe);
		T8 = x[1] - T0;
		T4 = sqrt(1.0 + 2.0 * T8 / T1);
		T2 = 2.0 * T8 / (T4 + 1.0);
		T3 = 0.5 * T2 * T2 / T1;
		T7 = 1.12 - T3 - 0.05;
		T6 = sqrt(T7 * T7 + 0.224);
		T5 = 1.12 - 0.5 * (T7 + T6);
	}
  else
    Vgse = x[1];
  //----------end Vgse-------------------------

  //---------Calculation for Vgsteff------------
  T2 = m * (Vgse - Vth) /(n * vt);

  if (T2 > EXP_THRESHOLD)
    T0 = m * (Vgse - Vth);
  else
    T0 = n * vt * log(1.0 + exp(T2));

  if (-T2 > EXP_THRESHOLD)
    T0 = vt * log(1.0 + MIN_EXP) * n;
  else
    T0 = n * vt * log(1.0 + exp(T2));

  T2 = (Voff_prime - (1.0 - m) * (Vgse - Vth)) / (n * vt);

  if (-T2 > EXP_THRESHOLD)
    T1 = m + n * coxe * MIN_EXP / sqrt((phi_s * q * NDEP * esi * 1.0e6)/2.0);
  else
    T1 = m + n * coxe * exp(T2) / sqrt((phi_s * q * NDEP * esi * 1.0e6)/2.0);

  if (T2 > EXP_THRESHOLD)
    T1 = m + n * coxe * MAX_EXP / sqrt((phi_s * q * NDEP * esi * 1.0e6)/2.0);
  else
    T1 = m + n * coxe * exp(T2) / sqrt((phi_s * q * NDEP * esi * 1.0e6)/2.0);

  Vgsteff = T0/T1;
  //--------end Vgsteff-------------------------

  //---Calculation for Weff------------------
  T9 = sqrtPhis - sqrt(phi_s);
  Weff += -2.0 * (DWJ * Vgsteff + DWB * T9);

  if (Weff < 2.0e-8)
	{
		T0 = 1.0 / (6.0e-8 - 2.0*Weff);
		Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
	}
  //-----end Weff---------------------------


	//--------Calculation for Rds----------------
  T0 = 1.0 + PRWG * Vgsteff;
  T1 = PRWB * (sqrtPhis - sqrt(phi_s));
  T2 = 1.0 / T0 + T1;
  T3 = T2 + sqrt(T2 * T2 + 0.01);
  T5 = RDSW + PRT * (T/300.0 - 1.0) * NF / (pow(Weffcj * 1.0e6,WR) * NF);
  T4 = 0.5 * T5;
  RDSWMIN += PRT * (T/300.0 - 1.0) * NF / (pow(Weffcj * 1.0e6,WR) * NF);
  Rds = RDSWMIN + T3 * T4;
	//------end Rds------------------------------

	//---------Calculation for lambda----
  if (A1 == 0.0)
    lambda = A2;
  else if (A1 > 0.0)
	{
		T0 = 1.0 - A2;
		T1 = T0 - A1 * Vgsteff - 0.0001;
		T2 = sqrt(T1 * T1 + 0.0004 * T0);
		lambda = A2 + T0 - 0.5 * (T1 + T2);
	}
  else
	{
		T1 = A2 + A1 * Vgsteff - 0.0001;
		T2 = sqrt(T1 * T1 + 0.0004 * A2);
		lambda = 0.5 * (T1 + T2);
	}
  //---------end lambda----------------


  //-----------Calculation for F_doping--------------
  F_doping = 0.5 * sqrt(1.0 + LPEB/Leff) * K1ox / (sqrt(phi_s - Vbseff))
	+ K2ox - K3B * TOXE * phi_s / (Weff + W0);
  //--------end of F-doping----------

  //--------Calculation for Abulk----------
  T0 = Leff/(Leff + 2.0 * sqrt(XJ * Xdep));
  T1 = A0 * T0 + B0/(Weff + B1);
  T2 = 1.0 + F_doping * T1;

  AD Abulk0 = T2;

  Abulk = Abulk0 - Vgsteff * F_doping * AGS * A0 * pow(T0,3.0);

  if (-Abulk + 0.1 > 0.0)
    Abulk = (0.2 - Abulk) / (3.0 - 20.0 * Abulk);

  if (-Abulk0 + 0.1 > 0.0)
    Abulk0 = (0.2 - Abulk0) / (3.0 - 20.0 * Abulk0);

  if (-Abulk + 0.1 > 0.0)
    Abulk = (0.2 - Abulk) / (3.0 - 20.0 * Abulk);

  T2 = KETA * Vbseff;

  if (T2 + 0.9 > 0.0)
    T0 = 1.0 / (1.0 + T2);
  else
    T0 = (17.0 + 20.0 * T2) / (0.8 + T2);

  Abulk *= T0;
  Abulk0 *= T0;
  //-------end Abulk-------------------

  //--------Calculation for Vb----------
  Vb = (Vgsteff + 2.0*vt)/Abulk;
  //-----------end Vb-------------------

  //----------Calculation for ueff--------------
  T1 = (Vgsteff + 2.0 * (VTH0 - VFB - phi_s))/TOXE;
  T2 = exp(EU * log(T1)) * (UA + UC * Vbseff);

  if (T2 + 0.8 > 0.0)
    T3 = 1.0 + T2;
  else
    T3 = (0.6 + T2) / (7.0 + 10.0 * T2);

  if (U0 > 1.0)
    U0 = U0 / 1.0e4;

  ueff = U0 / T3;
  //---------end ueff---------------------------

  //---------Calculation for Esat----------------
  Esat = 2.0 * (VSAT - AT * (T/300.0 - 1))/ueff;
  //----------end of Esat------------------------

  //---------Calculation for Vdsat-------------
  if((Rds == 0.0) && (lambda == 1.0))
	{
		T0 = 1.0 / (Abulk + Esat * Leff + Vgsteff + 2.0 * vt);
		T3 = Esat * Leff * (Vgsteff + 2.0 * vt);
		Vdsat = T3 * T0;
	}
  else
	{
		T6 = (Vgsteff + 2.0 * vt) * Weff * coxe * Rds * (VSAT - AT * (T/300.0 - 1.0));
		T9 = Abulk * Weff * coxe * Rds * (VSAT - AT * (T/300.0 - 1.0));
		T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0/lambda);
		T1 = (Vgsteff + 2.0 * vt) * (2.0/lambda - 1.0) + Abulk * Esat * Leff + 3.0 * (Vgsteff + 2.0 * vt) * Abulk * Rds * Weff * coxe * (VSAT - AT * (T/300.0 - 1.0));
		T2 = (Vgsteff + 2.0 * vt) * (Esat * Leff + 2.0 * T6);
		T3 = sqrt(T1*T1 - 2.0 * T0 * T2);
		Vdsat = (T1 - T3)/T0;
	}
  //---------end Vdsat-------------------------

  //----------Calculation for Vdseff---------
  T1 = Vdsat - x[0] - DELTA;
  T2 = sqrt(T1 * T1 + 4.0 * DELTA * Vdsat);
  T0 = T1/T2;
  T3 = 2.0 * DELTA / T2;

  if (T1 > 0.0)
    Vdseff = Vdsat - 0.5 * (T1 + T2);
  else
    Vdseff = Vdsat * (1.0 - 2.0 * DELTA / (T2 - T1));

  if (x[0] == 0.0)
    Vdseff = 0.0;

  if (Vdseff > x[0])
    Vdseff = x[0] - zero;

  //---------end Vdseff----------------------

  //------------Calculation for Vasat-------
  T0 = Esat*Leff + Vdsat + 2.0 * (Vgsteff * Rds * Weff * coxe * (VSAT - AT * (T/300.0 - 1.0))) * (1.0 - 0.5 * Abulk * Vdsat / (Vgsteff + 2.0 * vt));

  T1 = 2.0/lambda - 1.0 + Rds * Vgsteff * Weff * coxe * (VSAT - AT * (T/300.0 - 1.0));

  Vasat = T0 / T1;
  //-----------end Vasat--------------------

  //----------Calculation for Ids-----------
  T1 = 4.0 * (VTH0 - VFB - phi_s);
  T2 = (Vgsteff + T1) / (2.0e8 * TOXP);
  T3 = 1.9e-9 / (1.0 + exp(0.7 * log(T2)));

  T4 = esi * coxp /(esi + coxp * T3);
  T5 = ueff * T4 * Weff/Leff;
  T6 = Vgsteff * (1.0 - 0.5 * Vdseff * Abulk/(Vgsteff + 2.0 * vt));
  T7 = T5 * T6 / (1.0 + Vdseff/(Esat * Leff));

  Ids0 = T7 / (1.0 + T7 * Rds);

  AD F;

  if (FPROUT > 0.0)
    F = 1.0 / (1.0 + FPROUT * (sqrt(Leff)/(Vgsteff + 2.0 * vt)));
  else
    F = 1.0;

  //------Calculation of Cclm and Vaclm---------
  if (Vgsteff * PVAG/(Esat * L) + 0.9 > 0.0)
    T7 = 1.0 + Vgsteff * PVAG/(Esat * L);
  else
    T7 = (0.8 + (Vgsteff * PVAG/(Esat * L))) / (17.0 + 20.0 * (Vgsteff * PVAG/(Esat * L)));

  if (x[0] - Vdseff > 1.0e-10)
    Cclm = F * T7 * (1.0 + Rds * Ids0) * (Leff + Vdsat/Esat) / (PCLM * litl);
  else
    Cclm = MAX_EXP;

  if (x[0] - Vdseff > 1.0e-10)
    Vaclm = Cclm * (x[0] - Vdseff);
  else
    Vaclm = MAX_EXP;

  //--------end of Cclm and Vaclm----------

  Va = Vasat + Vaclm;

  //----------Calculation for Vadibl-------
  if (-(DROUT*Leff/lt0)+EXP_THRESHOLD > 0.0)
    T1 = exp(DROUT * Leff/lt0) / (2.0 * exp(DROUT * Leff/lt0) * MIN_EXP + pow(exp(DROUT * Leff/lt0)-1.0 ,2.0));
  else
    T1 = 1.0 / (MAX_EXP - 2.0);

  theta_rout = PDIBL1 * T1 + PDIBL2;

  if (theta_rout > 0.0)
	{
		T2 = (Vgsteff + 2.0 * vt) * Abulk * Vdsat;
		T3 = Vgsteff + 2.0 * vt + Abulk * Vdsat;
		Vadibl = (Vgsteff + 2.0 * vt - T2/T3) / theta_rout;

    if (PDIBLCB*Vbseff + 0.9 > 0.0)
      Vadibl = Vadibl * (1.0 / (1.0 + (PDIBLCB * Vbseff)));
    else
      Vadibl = Vadibl * (1.0 / (1.0 + (PDIBLCB * Vbseff)));

		Vadibl *= T7;
	}
  else
    Vadibl = MAX_EXP;
  //-----------end Vadibl---------------

  //-----Calculation for Vascbe----------
  if (PSCBE2 > 0.0)
	{
    if (x[0] - Vdseff > PSCBE1*litl/EXP_THRESHOLD)
      Vascbe = Leff * exp(PSCBE1 * litl/(x[0] - Vdseff)) / PSCBE2;
    else
      Vascbe = MAX_EXP * Leff / PSCBE2;
	}
  else
    Vascbe = MAX_EXP;
  //----------end Vascbe----------------

  //---------Calculation for Vadits------
  if (PDITSD * x[0] > EXP_THRESHOLD)
    T2 = MAX_EXP;
  else
    T2 = exp(PDITSD * x[0]);

  if (PDITS > 0.0)
    Vadits = F * (1/PDITS) * (1.0 + (1.0 + PDITSL * Leff) * T2);
  else
    Vadits = MAX_EXP;

  //--------end Vadits-------------------

  Ids = Ids0 * (1.0 + (x[0] - Vdseff)/Vadibl) * (1.0 + (x[0] - Vdseff)/Vadits) * (1.0 + log(Va/Vasat)/Cclm);

  //--------Calculation for Subthreshold current------
  T1 = ALPHA0 + ALPHA1 * Leff;
  if ((T1 <= 0.0) || (BETA0 <= 0.0))
    Isub = 0.0;
  else
	{   T2 = T1 / Leff;
    if ((x[0] - Vdseff) > BETA0 / EXP_THRESHOLD)
		{   T0 = -BETA0 / (x[0] - Vdseff);
      T1 = T2 * (x[0] - Vdseff) * exp(T0);
      T3 = T1 / (x[0] - Vdseff) * (T0 - 1.0);
		}
    else
		{   T3 = T2 * MIN_EXP;
      T1 = T3 * (x[0] - Vdseff);
		}
    T4 = Ids * Vdseff;
    Isub = T1 * T4;
	}
  //------------end Isub------------------------------

  //-------Calculation for Igidl---------------------
  T0 = 3.0 * TOXE;
  T1 = (x[0] - Vgse - EGIDL) / T0;

  if ((AGIDL <= 0.0) || (BGIDL <= 0.0) || (T1 < 0.0) || (CGIDL <= 0.0) || ((x[0] - x[2]) < 0.0))
    Igidl = 0.0;
  else
	{T2 = BGIDL / T1;
    if (T2 < 100.0)
		{   Igidl = AGIDL * Weffcj * T1 * exp(-T2);
      T3 = Igidl * (1.0 + T2) / T1;
		}
    else
      Igidl = T1 * AGIDL * Weffcj * 3.720075976e-44;

    T4 = pow((x[0] - x[2]),2.0);
    T5 = (x[0] - x[2]) * T4;
    T6 = CGIDL + T5;
    T7 = T5 / T6;
    Igidl *= T7;
	}
  //-------end Igidl---------------------------------

  //----and finally.....
  Ids =  NF * Vdseff * Ids * (1.0 + (x[0] - Vdseff)/Vascbe);
  //--------end Ids--------------------------

  //Equations for gate current

  double A = 4.97232e-7;
  double B = 7.45669e11;

  //-----------Calculation for Voxacc--------------------
  AD Vfbzb, Voxacc;

  T0 = DVT1W * Weff * Leff/(sqrt(esi*TOXE*Xdep0/e0));

  if (-T0+EXP_THRESHOLD > 0.0)
    T8 = exp(T0) / ((exp(T0)-1.0)*(exp(T0)-1.0) + 2.0*exp(T0)*MIN_EXP);
  else
    T8 = 1.0/(MAX_EXP - 2.0);

  T0 = DVT0W * T8;
  T8 = T0 * (Vbi - phi_s);

  T0 = DVT1 * Leff /(sqrt(esi*TOXE*Xdep0/e0));

  if (-T0+EXP_THRESHOLD > 0.0)
    T9 = exp(T0) / ((exp(T0)-1.0)*(exp(T0)-1.0) + 2.0*exp(T0)*MIN_EXP);
  else
    T9 = 1.0/(MAX_EXP - 2.0);

  double KT1 = -0.11;
  double KT1L = 0.0;

  T9 = DVT0 * T9 * (Vbi - phi_s);
  T4 = TOXE * phi_s / (Weff + W0);
  T0 = sqrt(1.0 + LPE0/Leff);
  T5 = K1ox * (T0 - 1.0) * sqrtPhis + (KT1 + KT1L / Leff) * (T/300.0 - 1.0);
  T6 = VTH0 - T8 - T9 + K3 * T4 + T5;

  Vfbzb = T6 - phi_s - K1 * sqrtPhis;

  T3 = Vfbzb - Vgse + Vbseff - 0.02;

  if (Vfbzb <= 0.0)
    T0 = sqrt(pow(T3,2.0) + 4*Vfbzb*0.02);
  else
    T0 = sqrt(pow(T3,2.0) - 4*Vfbzb*0.02);

  Voxacc = Vfbzb - (Vfbzb - 0.5*(T3 + T0));

  if (Voxacc < 0.0)
    Voxacc = 0.0;
  //---------end Voxacc--------------------------------

  //--------Calculation for Voxdepinv------------------
  AD Voxdepinv;

  T1 = Vgse - (Vfbzb - 0.5*(T3 + T0)) - Vbseff - Vgsteff;

  if (K1ox == 0.0)
    Voxdepinv = 0.0;

  if (T1 <= 0.0)
    Voxdepinv = -T1;
  else
    Voxdepinv = K1ox * (sqrt(T1 + 0.25*K1ox*K1ox) - 0.5*K1ox);

  Voxdepinv += Vgsteff;
  //-------end Voxdepinv------------------------------

  //-------Calculation for Igc---------------------
  AD Igc, Vaux;

  if ((Vgse-VTH0)/(vt*NIGC) > EXP_THRESHOLD)
    Vaux = Vgse - VTH0;
  else
    Vaux = vt*NIGC*log(1.0 + exp((Vgse - VTH0)/(vt*NIGC)));

  if (-(Vgse-VTH0)/(vt*NIGC) > EXP_THRESHOLD)
    Vaux = vt*NIGC*log(1.0 + MIN_EXP);
  else
    Vaux = vt*NIGC*log(1.0 + exp((Vgse - VTH0)/(vt*NIGC)));

  T3 = AIGC * CIGC - BIGC;
  T4 = BIGC * CIGC;
  T5 = (-B*TOXE) * (AIGC + T3 * Voxdepinv - T4 * Voxdepinv * Voxdepinv);

  if (T5 > EXP_THRESHOLD)
    T6 = MAX_EXP;
  else
    T6 = exp(T5);

  if (-T5 > EXP_THRESHOLD)
    T6 = MIN_EXP;
  else
    T6 = exp(T5);

  Igc = A* Weff* Leff* exp(NTOX * log(TOXREF/TOXE)) * Vgse * Vaux * T6;
  //--------end Igc---------------------------------

  //----Calculation for Igcs and Igcd--------------
  AD Igcs, Igcd;

  if ((-PIGCD*x[0]) > EXP_THRESHOLD)
    T1 = MAX_EXP;
  else if ((-PIGCD*x[0]) < -EXP_THRESHOLD)
    T1 = MIN_EXP;
  else
    T1 = exp(-PIGCD*x[0]);

  T8 = (-PIGCD*x[0]) * (-PIGCD*x[0]) + 2.0e-4;
  T0 = T8 * T8;
  T2 = T1 - 1.0 + 1.0e-4;
  T10 = (T2 - (-PIGCD*x[0])) / T8;

  Igcs = Igc * T10;

  T2 = T1 - 1.0 - 1.0e-4;
  T10 = ((-PIGCD * x[0]) * T1 - T2) / T8;

  Igcd = Igc * T10;
  //-----end Igcs and Igcd-------------------------

  //----Calculation for Igs and Igd---------------
  AD Igs, Igd;
  double Vfbsd;

  if (NGATE > 0.0)
    Vfbsd = (k*T/q)*log(NGATE/NSD);
  else
    Vfbsd = 0.0;

  T0 = x[1] - Vfbsd;
  Vgse = sqrt(T0 * T0 + 1.0e-4);

  T2 = x[1] * Vgse;
  T3 = AIGSD * CIGSD - BIGSD;
  T4 = BIGSD * CIGSD;
  T5 = (-B * TOXE * POXEDGE) * (AIGSD + T3 * Vgse - T4 * Vgse * Vgse);

  double  Toxratioedge = exp(NTOX * log(TOXREF / (TOXE * POXEDGE))) / TOXE / TOXE / POXEDGE / POXEDGE;

  if (T5 > EXP_THRESHOLD)
    T6 = MAX_EXP;
  else
    T6 = exp(T5);

  if (-T5 > EXP_THRESHOLD)
    T6 = MIN_EXP;
  else
    T6 = exp(T5);

  Igs = A * Weff * Toxratioedge * DLCIG * T6 * T2;

  T0 = (x[1] - x[0]) - Vfbsd;

  AD Vgde;
  Vgde = sqrt(T0 * T0 + 1.0e-4);

  T2 = Vgde * (x[1] - x[0]);
  T5 = (-B * TOXE * POXEDGE) * (AIGSD + T3 * Vgde - T4 * Vgde * Vgde);

  if (T5 > EXP_THRESHOLD)
    T6 = MAX_EXP;
  else if (T5 < -EXP_THRESHOLD)
    T6 = MIN_EXP;
  else
    T6 = exp(T5);

  Igd = A * Weff * Toxratioedge * DLCIG * T6 * T2;
  //------end Igs and Igd-----------------------------

  //---------Calculations for Igb--------------------
  T0 = vt * NIGBACC;
  T1 = -Vgse + Vbseff + Vfbzb;

  if (T1/T0 > EXP_THRESHOLD)
    Vaux = T1;
  else
    Vaux = T0 * log(1.0 + exp(T1/T0));

  if (-T1/T0 > EXP_THRESHOLD)
    Vaux = T0 * log(1.0 + MIN_EXP);
  else
    Vaux = T0 * log(1.0 + exp(T1/T0));

  T2 = (Vgse - Vbseff) * Vaux;
  T11 = 4.97232e-7 * Weff * Leff * exp(NTOX * log(TOXREF/TOXE));
  T12 = -7.45669e-11 * TOXE;
  T3 = AIGBACC * CIGBACC - BIGBACC;
  T4 = BIGBACC * CIGBACC;
  T5 = T12 * (AIGBACC + T3 * Voxacc - T4 * Voxacc * Voxacc);

  if (T5 > EXP_THRESHOLD)
    T6 = MAX_EXP;
  else
    T6 = exp(T5);

  if (-T5 > EXP_THRESHOLD)
    T6 = MIN_EXP;
  else
    T6 = exp(T5);

  AD Igbacc;
  Igbacc = T11 * T2 * T6;

  T0 = vt * NIGBINV;
  T1 = Voxdepinv - EIGBINV;

  if (T1/T0 > EXP_THRESHOLD)
    Vaux = T1;
  else
    Vaux = T0 * log(1.0 + exp(T1/T0));

  if (-T1/T0 > EXP_THRESHOLD)
    Vaux = T0 * log(1.0 + MIN_EXP);
  else
    Vaux = T0 * log(1.0 + exp(T1/T0));

  T2 = (Vgse - Vbseff) * Vaux;
  T11 *= 0.75610;
  T12 *= 1.31724;
  T3 = AIGBINV * CIGBINV - BIGBINV;
  T4 = BIGBINV * CIGBINV;
  T5 = T12 * (AIGBINV + T3 * Voxdepinv - T4 * Voxdepinv * Voxdepinv);

  if (T5 > EXP_THRESHOLD)
    T6 = MAX_EXP;
  else
    T6 = exp(T5);

  if (-T5 > EXP_THRESHOLD)
    T6 = MIN_EXP;
  else
    T6 = exp(T5);

  AD Igbinv;
  Igbinv = T11 * T2 * T6;

  AD Igb;
  Igb = Igbinv + Igbacc;

  //Equations for capacitance. CAPMOD=2 with a 40/60 charge partiton between the source and drain.
  //The equations for charge at the repsective transistor nodes begins here.
  //The derivatives of charge with respect to time are evaluated in eval2.
  //That corresponds to the current contribution at each node.

  //------Calculation for LeffCV and WeffCV----------
  Leff = Lnew - 2.0 * DLC;
  Weff = Wnew - 2.0 * DWC;
  //------end LeffCV and WeffCV---------------------

  //-----Calculation for VbseffCV---------------------
  AD VbseffCV;

  if (Vbseff > 0.0)
    VbseffCV = phi_s - Phis;
  else
    VbseffCV = Vbseff;
  //------end VbseffCV--------------------------------

  //------Calculation for VgsteffCV------------------
  T0 = n * NOFF * k * T / q;
  T1 = (Vgse - Vth) / T0;

  if (T1 > EXP_THRESHOLD)
    Vgsteff = Vgse - Vth - VOFFCV;
  else
    Vgsteff = T0 * log(1.0 + exp(T1));

  if (-T1 > EXP_THRESHOLD)
    Vgsteff = T0 * log(1.0 + MIN_EXP);
  else
    Vgsteff = T0 * log(1.0 + exp(T1));
  //-------End VgsteffCV------------------------------

  //---------------Calculation for Vfbeff-------------
  AD V3;
  V3 = Vfbzb - Vgse + VbseffCV - 0.02;

  if (Vfbzb > 0.0)
    T0 = sqrt(V3 * V3 + 4.0 * 0.02 * Vfbzb);
  else
    T0 = sqrt(V3 * V3 - 4.0 * 0.02 * Vfbzb);

  AD Vfbeff;
  Vfbeff = Vfbzb - 0.5 * (V3 + T0);
  //-----------end Vfbeff-----------------------------

  T0 = (Vgse - VbseffCV - Vfbzb) / TOXP;

  //-------------Calculation for Tcen----------------
  double LDEB = sqrt(esi * vt/(q * NDEP * 1.0e6)) / 3.0;

  AD Tcen;
  T1 = T0 * ACDE;

  if (EXP_THRESHOLD + T1 > 0.0)
    Tcen = LDEB * exp(T1);
  else
    Tcen = LDEB * MAX_EXP;

  if (-T1 + EXP_THRESHOLD > 0.0)
    Tcen = LDEB * exp(T1);
  else
    Tcen = LDEB * MAX_EXP;

  if (-T1 > EXP_THRESHOLD)
    Tcen = LDEB * MIN_EXP;
  else
    Tcen = LDEB * MAX_EXP;

  V3 = LDEB - Tcen - 1.0e-3 * TOXP;

  AD V4;
  V4 = sqrt(V3 * V3 + 4.0 * 1.0e-3 * TOXP * LDEB);

  Tcen = LDEB - 0.5 * (V3 + V4);
  //-------------end Tcen----------------------------

  AD Ccen;
  Ccen = esi / Tcen;

  AD Coxeff;
  Coxeff = Ccen * coxp / (Ccen + coxp);

  //--------Calculation for QoverlapCox-------------
  AD QoverlapCox, Qac0, CoxWLcen, Qsub0;

  CoxWLcen = coxp * Weff * Leff * NF * Coxeff / coxe;
  Qac0 = CoxWLcen * (Vfbeff - Vfbzb);
  QoverlapCox = Qac0 / Coxeff;

  T0 = 0.5 * K1ox;
  T3 = Vgse - Vfbeff - VbseffCV - Vgsteff;
  if (K1ox == 0.0)
	{ T1 = 0.0;
    T2 = 0.0;
	}
  else if (T3 < 0.0)
	{ T1 = T0 + T3 / K1ox;
    T2 = CoxWLcen;
	}
  else
	{ T1 = sqrt(T0 * T0 + T3);
    T2 = CoxWLcen * T0 / T1;
	}

  Qsub0 = CoxWLcen * K1ox * (T1 - T0);
  QoverlapCox = Qsub0 / Coxeff;

  //--------Calculation for Delta_phis------------
  if (K1ox > 0.0)
    T2 = MOIN * vt * K1ox * K1ox;
  else
    T2 = 0.25 * MOIN * vt;

  if (K1ox > 0.0)
    T0 = K1ox * sqrt(phi_s);
  else
    T0 = 0.5 * sqrt(phi_s);

  AD Delta_phis;

  T1 = 2.0 * T0 + Vgsteff;

  Delta_phis = vt * log(1.0 + T1 * Vgsteff / T2);
  //----------end Delta_phis-------------------

  //The calculation for Tcen must be done once more
  T0 = (Vgsteff + 4.0*(VTH0 - VFB - phi_s))/ (2.0 * TOXP);
  T1 = 1.0 + exp(0.7 * log(T0));
  T2 = 0.7 * exp(0.7 * log(T0)) / (T0 * 2.0 * TOXP);
  Tcen = 1.9e-9 / T1;

  Ccen = esi / Tcen;
  Coxeff = Ccen * coxp / (Ccen + coxp);
  CoxWLcen = coxp * Weff * Leff * Coxeff / coxe;

  AD AbulkCV;
  AbulkCV = Abulk0 * (1.0 + pow((CLC/Leff),CLE));

  AD VdsatCV;
  VdsatCV = (Vgsteff - Delta_phis) / AbulkCV;

  T0 = VdsatCV - x[0] - 0.02;
  T1 = sqrt(T0 * T0 + 4.0 * 0.02 * VdsatCV);

  AD VdseffCV;

  if (T0 > 0.0)
    VdseffCV = VdsatCV - 0.5 * (T0 + T1);
  else
    VdseffCV = VdsatCV * (1.0 - 0.04/(T1-T0));

  if (x[0] == 0.0)
    VdseffCV = 0.0;

  T0 = AbulkCV * VdseffCV;
  T1 = Vgsteff - Delta_phis;
  T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
  T3 = T0 / T2;
  T4 = 1.0 - 12.0 * T3 * T3;
  T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
  T6 = T5 * VdseffCV / AbulkCV;

  AD qgate;
  qgate = CoxWLcen * (T1 - T0 * (0.5 - T3));

  AD qbulk;
  qbulk = CoxWLcen * (1.0 - AbulkCV) * (0.5*VdseffCV - T0*VdseffCV/T2);

  QoverlapCox = qbulk / Coxeff;

  T2 = T2 / 12.0;
  T3 = 0.5 * CoxWLcen / (T2 * T2);
  T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0 * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;

  AD qsrc;
  qsrc = -T3 * T4;

  qgate += Qac0 + Qsub0 - qbulk;
  qbulk -= Qac0 + Qsub0;
  AD qdrn;
  qdrn = -(qbulk + qgate + qsrc);
  qsrc = -(qbulk + qgate + qdrn);

  // Calculate dynamic current contributions due to charge
  // iqd is dynamic contribution to drain current
  // iqg is dynamic contribution to gate current
  // iqs is dynamic contribution to source current
  AD iqd, iqg, iqs;
  iqd = qdrn.fastAccessDx(0)*x[3] + qdrn.fastAccessDx(1)*x[4] + qdrn.fastAccessDx(2)*x[5];
  iqg = qgate.fastAccessDx(0)*x[3] + qgate.fastAccessDx(1)*x[4] + qgate.fastAccessDx(2)*x[5];
  iqs = qsrc.fastAccessDx(0)*x[3] + qsrc.fastAccessDx(1)*x[4] + qsrc.fastAccessDx(2)*x[5];

  flow[0] = Ids + Isub + Igidl - Igcd - Igd + iqd; //Drain current
  flow[1] = Igs + Igcd + Igd + Igcs + Igb + iqg; //Gate current
  flow[2] = -Ids - Igs - Igcs - iqs; //Source current

  effort[0] = x[0] - x[2]; //Vdb
  effort[1] = x[1] - x[2]; //Vgb
  effort[2] = -x[2];//Vsb
}
