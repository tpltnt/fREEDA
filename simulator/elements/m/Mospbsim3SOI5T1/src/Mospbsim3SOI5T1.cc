#include "Mospbsim3SOI5T1.h"

const unsigned Mospbsim3SOI5T1::n_par = 116;

ItemInfo Mospbsim3SOI5T1::einfo =
{
  "mospbsim3soi5t1",
  "Mosfet model using Mospbsim3SOI5T1 level, version 3.2.4",
  "Ramya Mohan",
  DEFAULT_ADDRESS"transistor>mosfet",
  "2003_09_09"
};

ParmInfo Mospbsim3SOI5T1::pinfo[] =
{
  {"l", "Length (m)", TR_DOUBLE, false},
  {"w", "Width (m)", TR_DOUBLE, false},
  {"tsi", "Silicon film thickness (m)", TR_DOUBLE, false},
  {"tbox", "Buried Oxide thickness (m)", TR_DOUBLE, false},
  {"tox", "Gate oxide thickness (m)", TR_DOUBLE, false},
  {"toxqm", "Effective oxide thickness considering quantum effect",TR_DOUBLE, false},
  {"xj", "S/DJunction depth (m)", TR_DOUBLE, false},
  {"nch", "Channel doping concentration (1/cm^3)", TR_DOUBLE, false},
  {"nsub", "Substrate doping concentration (1/cm^3)", TR_DOUBLE, false},
  {"ngate", "Poly-gate doping concentration (1/cm^3)", TR_DOUBLE, false},
  {"vth0", "Threshold voltage @Vbs=0 for long and wide device", TR_DOUBLE, false},
  {"k1", "First order body effect coefficient (V^0.5)", TR_DOUBLE, false},
  {"k1w1", "First body effect width dependent parameter (m)", TR_DOUBLE, false},
  {"k1w2", "Second body effect width dependent parameter (m)", TR_DOUBLE, false},
  {"k2", "Second order body effect coefficient", TR_DOUBLE, false},
  {"k3", "Narrow width effect coefficient", TR_DOUBLE, false},
  {"k3b", "Body effect coefficient of k3 (1/V)", TR_DOUBLE, false},
  {"kb1", "Backgate body charge coefficient", TR_DOUBLE, false},
  {"w0", "Narrow width effect parameter (m)", TR_DOUBLE, false},
  {"nlx", "Lateral non-uniform doping parameter (m)", TR_DOUBLE, false},
  {"dvt0", "First coefficient of short-channel effect on Vth", TR_DOUBLE, false},
  {"dvt1", "Second coefficient of short-channel effect on Vth", TR_DOUBLE, false},
  {"dvt2", "Body-bias coefficient of short-channel effect on Vth (1/V)", TR_DOUBLE, false},
  {"dvt0w", "First coefficient of narrow width effect on Vth for small channel length", TR_DOUBLE, false},
  {"dvt1w", "Second coefficient of narrow width effect on Vth for small channel length", TR_DOUBLE, false},
  {"dvt2w", "Body-bias coefficient of narrow width effect on Vth for small channel length (1/V)", TR_DOUBLE, false},
  {"u0", "Mobility at Temp=Tnom (cm^2/V-sec)", TR_DOUBLE, false},
  {"ua", "First-oreder mobility degradation coefficient (m/V)", TR_DOUBLE, false},
  {"ub", "Second-order mobility degradation coefficient (m/V)^2", TR_DOUBLE, false},
  {"uc", "Body-effect of mobility degradation coefficient (1/V)", TR_DOUBLE, false},
  {"vsat", "Saturation velocity at Temp=Tnom (m/sec)", TR_DOUBLE, false},
  {"a0", "Bulk charge effect coefficient for channel length", TR_DOUBLE, false},
  {"ags", "Gate bias coefficient of Abulk (1/V)", TR_DOUBLE, false},
  {"b0", "Bulk charge effect coefficient for channel width (m)", TR_DOUBLE, false},
  {"b1", "Bulk charge effect width offset (m)", TR_DOUBLE, false},
  {"keta", "Body-bias coefficient of bulk charge effect (1/V)", TR_DOUBLE, false},
  {"ketas", "Surface potential adjustment for bulk charge effect (V)", TR_DOUBLE, false},
  {"a1", "First non-saturation effect parameter (1/V)", TR_DOUBLE, false},
  {"a2", "Second non-saturation effect parameter", TR_DOUBLE, false},
  {"rdsw", "Parasitic resistance per unit width (ohm-um)", TR_DOUBLE, false},
  {"prwb", "Body effect coefficient of rdsw (1/V)", TR_DOUBLE, false},
  {"prwg", "Gate-bias effect coefficient of rdsw", TR_DOUBLE, false},
  {"wr", "Width offset from Weff for Rds calculation", TR_DOUBLE, false},
  {"nfactor", "Subthreshold swing factor", TR_DOUBLE, false},
  {"wint", "Width offset fitting parameter from I-V without bias (m)", TR_DOUBLE, false},
  {"lint", "Length offset fitting parameter from I-V without bias (m)", TR_DOUBLE, false},
  {"dwg", "Coefficient of Weff's gate dependence (m/V)", TR_DOUBLE, false},
  {"dwb", "Coefficient of Weff's substrate body bias dependence (m/V)^0.5", TR_DOUBLE, false},
  {"dwbc", "Width offset for body contact isolation edge (m)", TR_DOUBLE, false},
  {"voff", "Offset voltage in the threshold region for large W and L (V)", TR_DOUBLE, false},
  {"eta0", "DIBL coefficeint subthreshold region", TR_DOUBLE, false},
  {"etab", "Body bias coefficeint for the subthreshold DIBL effect (1/V)", TR_DOUBLE, false},
  {"dsub", "DIBL coefficient in the subthreshold region", TR_DOUBLE, false},
  {"cit", "Interface trap capacitance (F/m^2)", TR_DOUBLE, false},
  {"cdsc", "Drain/Source to channel coupling capacitance (F/m^2)", TR_DOUBLE, false},
  {"cdscb", "Body-bias sensitivity of cdsc (F/m^2)", TR_DOUBLE, false},
  {"cdscd", "Drain-bias sensitivity of cdsc (F/m^2)", TR_DOUBLE, false},
  {"pclm", "Channel length modulation parameter", TR_DOUBLE, false},
  {"pdibl1", "First output resistance DIBL effect correction parameter", TR_DOUBLE, false},
  {"pdibl2", "Second output resistance DIBL effect correction parameter", TR_DOUBLE, false},
  {"pvag", "Gate dependence of Early voltage", TR_DOUBLE, false},
  {"delta", "Effective Vds parameter", TR_DOUBLE, false},
  {"alpha0", "The first parameter of impact ionization current (m/V)", TR_DOUBLE, false},
  {"beta0", "First Vds dependent parameter of impact ionization current (1/V)", TR_DOUBLE, false},
  {"beta1", "Second Vds dependent parameter of impact ionization current", TR_DOUBLE, false},
  {"beta2", "Third Vds dependent parameter of impact ionization current (V)", TR_DOUBLE, false},
  {"vdsatii0", "Nominal drain saturation voltage at threshold for impact ionization current (V)", TR_DOUBLE, false},
  {"cgeo", "Gate substrate overlap capacitance per unit channel length (F/m)", TR_DOUBLE, false},
  {"cjswg", "Source/Drain (gate side) sidewall junction capacitance per unit width (normalized to 100nm tsi) (F/m^2)",
  TR_DOUBLE, false},
  {"pbswg", "Source/Drain (gate side) sidewall junction capacitance built in potential (V)", TR_DOUBLE, false},
  {"mjswg", "Source/Drain (gate side) sidewall junction capacitance grading coefficient (V)", TR_DOUBLE, false},
  {"tt", "Diffusion capacitance transit time coefficient (sec)", TR_DOUBLE, false},
  {"csdesw", "S/D sidewall fringing capacitance per unit length (F/m)", TR_DOUBLE, false},
  {"cgsl", "Light doped source-gate region overlap capacitance (F/m)", TR_DOUBLE, false},
  {"cgdl", "Light doped drain-gate region overlap capacitance (F/m)", TR_DOUBLE, false},
  {"ckappa", "Coefficient for lightly doped region overlap capacitance fringing field capacitance (F/m)", TR_DOUBLE, false},
  {"clc", "Constant term for the short channel model (m)", TR_DOUBLE, false},
  {"cle", "Exponential term for the short channel model", TR_DOUBLE, false},
  {"dlc", "Length offset fitting parameter for gate charge (m)", TR_DOUBLE, false},
  {"dlcb", "Length offset fitting parameter for body charge (m)", TR_DOUBLE, false},
  {"dlbg", "Length offset fitting parameter for backgate charge (m)", TR_DOUBLE, false},
  {"dwc", "Width offset fitting parameter from C-V (m)", TR_DOUBLE, false},
  {"delvt", "Threshold voltage adjust for C-V (V)", TR_DOUBLE, false},
  {"fbody", "Scaling factor for body charge", TR_DOUBLE, false},
  {"moin", "Coefficient for the gate-bias dependent surface potential V^0.5", TR_DOUBLE, false},
  {"tnom", "Parameter measurement temperature (K)", TR_DOUBLE, false},
  {"ute", "Temperature coefficient of mobility", TR_DOUBLE, false},
  {"kt1", "Temperature coefficient of Vth (V)", TR_DOUBLE, false},
  {"kt1l", "Channel length dependence of the temperature coefficient of Vth (V*m)", TR_DOUBLE, false},
  {"kt2", "Body-bias coefficient of the Vth temperature effect", TR_DOUBLE, false},
  {"ua1", "Temperature coefficient for ua (m/V)", TR_DOUBLE, false},
  {"ub1", "Temperature coefficient for ub ((m/V)^2)", TR_DOUBLE, false},
  {"uc1", "Temperature coefficient for uc (1/V)", TR_DOUBLE, false},
  {"at", "Temperature coefficient of vsat (m/sec)", TR_DOUBLE, false},
  {"prt", "Temperature coefficient of rdsw (ohm-um)", TR_DOUBLE, false},
  {"vbm", "Maximum body voltage", TR_DOUBLE, false},
  {"xt1", "Doping depth", TR_DOUBLE, false},
  {"pdiblb", "Body-effect on drain induced barrier lowering", TR_DOUBLE, false},
  {"ll", "Length reduction parameter", TR_DOUBLE, false},
  {"llc", "Length reduction parameter for CV", TR_DOUBLE, false},
  {"lln", "Length reduction parameter", TR_DOUBLE, false},
  {"lw", "Length reduction parameter", TR_DOUBLE, false},
  {"lwc", "Length reduction parameter for CV", TR_DOUBLE, false},
  {"lwn", "Length reduction parameter", TR_DOUBLE, false},
  {"lwl", "Length reduction parameter", TR_DOUBLE, false},
  {"lwlc", "Length reduction parameter for CV", TR_DOUBLE, false},
  {"wl", "Width reduction parameter", TR_DOUBLE, false},
  {"wlc", "Width reduction parameter for CV", TR_DOUBLE, false},
  {"wln", "Width reduction parameter", TR_DOUBLE, false},
  {"ww", "Width reduction parameter", TR_DOUBLE, false},
  {"wwc", "Width reduction parameter for CV", TR_DOUBLE, false},
  {"wwn", "Width reduction parameter", TR_DOUBLE, false},
  {"wwl", "Width reduction parameter", TR_DOUBLE, false},
  {"wwlc", "Width reduction parameter for CV", TR_DOUBLE, false},
  {"temp", "Circuit temperature", TR_DOUBLE, false},
  {"acde", "Exponential coefficient for charge in m/V", TR_DOUBLE, false}
  //{"pmos", "True if PMOS", TR_BOOLEAN, false}
};

Mospbsim3SOI5T1::Mospbsim3SOI5T1(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  paramvalue[0] = &(l = 5e-6);
  paramvalue[1] = &(w = 5e-6);
  paramvalue[2] = &(tsi = 1e-7);
  paramvalue[3] = &(tbox = 3.0e-7);
  paramvalue[4] = &(tox = 1e-8);
  paramvalue[5] = &(toxqm = tox);
  paramvalue[6] = &(xj = tsi);
  paramvalue[7] = &(nch = 1.7e17);
  paramvalue[8] = &(nsub = 6e16);
  paramvalue[9] = &(ngate = 0.0);
  paramvalue[10] = &(vth0 = 0.7);
  paramvalue[11] = &(k1 = 0.6);
  paramvalue[12] = &(k1w1 = 0.0);
  paramvalue[13] = &(k1w2 = 0.0);
  paramvalue[14] = &(k2 = 0.0);
  paramvalue[15] = &(k3 = 0.0);
  paramvalue[16] = &(k3b = 0.0);
  paramvalue[17] = &(kb1 = 1);
  paramvalue[18] = &(w0 = 0);
  paramvalue[19] = &(nlx = 1.74e-7);
  paramvalue[20] = &(dvt0 = 2.2);
  paramvalue[21] = &(dvt1 = 0.53);
  paramvalue[22] = &(dvt2 = -0.032);
  paramvalue[23] = &(dvt0w = 0.0);
  paramvalue[24] = &(dvt1w = 5.3e6);
  paramvalue[25] = &(dvt2w = -0.032);
  paramvalue[26] = &(u0 = 0.025);
  paramvalue[27] = &(ua = 2.25e-9);
  paramvalue[28] = &(ub = 5.9e-19);
  paramvalue[29] = &(uc = -0.0465);
  paramvalue[30] = &(vsat = 8e4);
  paramvalue[31] = &(a0 = 1.0);
  paramvalue[32] = &(ags = 0.0);
  paramvalue[33] = &(b0 = 0.0);
  paramvalue[34] = &(b1 = 0.0);
  paramvalue[35] = &(keta = 0.0);
  paramvalue[36] = &(ketas = 0.0);
  paramvalue[37] = &(a1 = 0.0);
  paramvalue[38] = &(a2 = 1.0);
  paramvalue[39] = &(rdsw = 100);
  paramvalue[40] = &(prwb = 0.0);
  paramvalue[41] = &(prwg = 0.0);
  paramvalue[42] = &(wr = 1.0);
  paramvalue[43] = &(nfactor = 1.0);
  paramvalue[44] = &(wint = 0.0);
  paramvalue[45] = &(lint = 0.0);
  paramvalue[46] = &(dwg = 0.0);
  paramvalue[47] = &(dwb = 0.0);
  paramvalue[48] = &(dwbc = 0.0);
  paramvalue[49] = &(voff = -0.08);
  paramvalue[50] = &(eta0 = 0.08);
  paramvalue[51] = &(etab = -0.07);
  paramvalue[52] = &(dsub = 0.56);
  paramvalue[53] = &(cit = 0.0);
  paramvalue[54] = &(cdsc = 2.4e-4);
  paramvalue[55] = &(cdscb = 0.0);
  paramvalue[56] = &(cdscd = 0.0);
  paramvalue[57] = &(pclm = 1.3);
  paramvalue[58] = &(pdibl1 = 0.39);
  paramvalue[59] = &(pdibl2 = 0.0086);
  paramvalue[60] = &(pvag = 0.0);
  paramvalue[61] = &(delta = 0.01);
  paramvalue[62] = &(alpha0 = 0.0);
  paramvalue[63] = &(beta0 = 0.0);
  paramvalue[64] = &(beta1 = 0.0);
  paramvalue[65] = &(beta2 = 0.1);
  paramvalue[66] = &(vdsatii0 = 0.9);
  paramvalue[67] = &(cgeo = 0.0);
  paramvalue[68] = &(cjswg = 1e-10);
  paramvalue[69] = &(pbswg = 0.7);
  paramvalue[70] = &(mjswg = 0.5);
  paramvalue[71] = &(tt = 1e-12);
  paramvalue[72] = &(csdesw = 0.0);
  paramvalue[73] = &(cgsl = 0.0);
  paramvalue[74] = &(cgdl = 0.0);
  paramvalue[75] = &(ckappa = 0.6);
  paramvalue[76] = &(clc = 0.1e-7);
  paramvalue[77] = &(cle = 0.0);
  paramvalue[78] = &(dlc = lint);
  paramvalue[79] = &(dlcb = 0.0);
  paramvalue[80] = &(dlbg = 0.0);
  paramvalue[81] = &(dwc = wint);
  paramvalue[82] = &(delvt = 0.0);
  paramvalue[83] = &(fbody = 1.0);
  paramvalue[84] = &(moin = 15.0);
  paramvalue[85] = &(tnom = 300.15);
  paramvalue[86] = &(ute = -1.5);
  paramvalue[87] = &(kt1 = -0.11);
  paramvalue[88] = &(kt1l = 0.0);
  paramvalue[89] = &(kt2 = 0.022);
  paramvalue[90] = &(ua1 = 4.31e-9);
  paramvalue[91] = &(ub1 = -7.61e-18);
  paramvalue[92] = &(uc1 = -0.056);
  paramvalue[93] = &(at = 3.3e4);
  paramvalue[94] = &(prt = 0.0);
  paramvalue[95] = &(vbm = -3.0);
  paramvalue[96] = &(xt1 = 1.55e-7);
  paramvalue[97] = &(pdiblb = 0.0);
  paramvalue[98] = &(ll = 0.0);
  paramvalue[99] = &(llc = 0.0);
  paramvalue[100] = &(lln = 1.0);
  paramvalue[101] = &(lw = 0.0);
  paramvalue[102] = &(lwc = 0.0);
  paramvalue[103] = &(lwn = 1.0);
  paramvalue[104] = &(lwl = 0.0);
  paramvalue[105] = &(lwlc = 0.0);
  paramvalue[106] = &(wl = 0.0);
  paramvalue[107] = &(wlc = 0.0);
  paramvalue[108] = &(wln = 1.0);
  paramvalue[109] = &(ww = 0.0);
  paramvalue[110] = &(wwc = 0.0);
  paramvalue[111] = &(wwn = 1.0);
  paramvalue[112] = &(wwl = 0.0);
  paramvalue[113] = &(wwlc = 0.0);
  paramvalue[114] = &(temp = 300.15);
  paramvalue[115] = &(acde = 1.0);
  //paramvalue[98] = &(pmos=false);

  setNumTerms(5);

  setNumberOfStates(4);

  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void Mospbsim3SOI5T1::init() throw(string&)
{
  DenseIntVector var(4);
  var[0] = 0; //Vds
  var[1] = 1; //Vgs
  var[2] = 2; //Ves
  var[3] = 3; //Vbs
  DenseDoubleVector nodelay;
  initializeAD(var, var);
}

void Mospbsim3SOI5T1::eval(AD * x, AD * effort, AD * flow)
{
  AD Vgs_eff, V0, Vbsmos, Vbp, Vbsh, Vbseff, Phis, sqrtPhis, Xdep, ltl, ltw;
  AD Theta0, thetavth, Delt_vth, DeltVthw, DeltVthtemp, DIBL_Sft;
  AD sqrtPhisExt, Vth, n, Vgst, VgstNVt, ExpArg, Vgsteff, Vgst2Vtm, Weff, Rds;
  AD Abulk, Abulk0, Denomi, ueff, WVCox, WVCoxRds, Esat, EsatL, AbovVgst2Vtm;
  AD Vdsat, Vdseff, diffVds, Vasat, VACLM, VADIBL, Va, CoxWovL, beta, fgche1, fgche2;
  AD gche, Vfb, V3, Vfbeff, Qac0, Qsub0, AbulkCV, VdsatCV, V4, VdseffCV, qbulk;
  AD qinv, qsrc, qgate, qbody, qsub, qdrn;
  AD T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T13, T14;
  double t0, t1, t2, t3, tmp1, tmp2, vbsc;

	#define EPSOX 3.453133e-11
	#define EPSSI 1.03594e-10
	#define Charge_q 1.60219e-19
	#define KboQ 8.617087e-5  /*  Kb / q   */
	#define Eg300 1.115   /*  energy gap at 300K  */
	#define DELTA_1 0.02
	#define DELTA_2 0.02
	#define DELTA_3 0.02
	#define DELTA_3_SOI 0.02
	#define DELTA_4 0.02
	#define DELT_Vbseff  0.005
	#define DELTA_VFB  0.02
	#define OFF_Vbsitf 0.02
	#define CONST_2OV3 0.6666666666
	#define MAX_EXPL 2.688117142e+43
	#define MIN_EXPL 3.720075976e-44
	#define EXPL_THRESHOLD 100.0
	#define PI 3.141592654
	#define MAX_EXP 5.834617425e14
	#define MIN_EXP 1.713908431e-15
	#define EXP_THRESHOLD 34.0

	double TempRatio = temp / tnom;
	double factor1 = sqrt(EPSSI / EPSOX * tox);
	double Vtm0 = KboQ * tnom;
	double Vtm = KboQ * temp;

	double Eg = 1.16 - 7.02e-4 * temp * temp / (temp + 1108.0);
	/* ni is in cm^-3 */
	double ni = 1.45e10 * (temp / 300.15) * sqrt(temp / 300.15) * exp(21.5565981 - Eg / (2.0 * Vtm));

	double cox = 3.453133e-11 / tox;

	double ldrn = l;
	double wdrn = w;

	t0 = pow(ldrn, lln);
	t1 = pow(wdrn, lwn);
	tmp1 = ll / t0 + lw / t1 + lwl / (t0 * t1);
	double dl = lint + tmp1;
	tmp1 = llc / t0 + lwc / t1 + lwlc / (t0 * t1);
	double dlc = dlc + tmp1;

	t2 = pow(ldrn, wln);
	t3 = pow(wdrn, wwn);
	tmp2 = wl / t2 + ww / t3 + wwl / (t2 * t3);
	double dw = wint + tmp2;
	tmp2 = wlc / t2 + wwc / t3 + wwlc / (t2 * t3);
	double dwc = dwc + tmp2;

	double nbc = 0.0;
	double leff = l - 2.0 * dl;
	double weff = w - nbc * dwbc - (2.0 - nbc) * dw;

	double nseg = 1;
	double pdbcp = 0;
	double psbcp = 0;

	double leffCV = l - 2.0 * dlc;
	double weffCV = w - nbc * dwbc - (2.0 - nbc) * dwc;

	double wdiodCV = weffCV / nseg + pdbcp;
	double wdiosCV = weffCV / nseg + psbcp;
	double leffCVb = l - 2.0 * dlc - dlcb;
	double leffCVbg = leffCVb + 2 * dlbg;

	double T00 = (TempRatio - 1.0);

	//  double ua = ua + ua1 * T00;
	//  double ub = ub + ub1 * T00;
	//  double uc = uc + uc1 * T00;

	// if (u0 > 1.0)
	// u0 = u0 / 1.0e4;

	double u0temp = u0 * pow(TempRatio, ute);
	double vsattemp = vsat - at * T00;
	double rds0 = (rdsw + prt * T00) / pow(weff * 1e6, wr);

	double cf = 2.0 * EPSOX / PI * log(1.0 + 0.4e-6 / tox);

	double cgdo = (cgdo + cf) * wdiodCV;
	double cgso = (cgso + cf) * wdiosCV;
	double cgeo = cgeo * leffCV;
	double vtm = Vtm;
	double loggg = log(nch / ni);
	double phi = 2.0 * vtm * loggg;

	double sqrtPhi = sqrt(phi);
	double Xdep0 = sqrt(2.0 * EPSSI / (Charge_q * nch * 1.0e6)) * sqrtPhi;
	double litl = sqrt(3.0 * xj * tox);
	double vbi = vtm * log(1.0e20 * nch / (ni * ni));
	double cdep0 = sqrt(Charge_q * EPSSI * nch * 1.0e6 / 2.0 / phi);

	double vbx = phi - 7.7348e-4 * nch * xt1 * xt1;
	double gamma1 = 5.753e-12 * sqrt(nch) / cox;
	double gamma2 = 5.753e-12 * sqrt(nsub) / cox;

	t0 = gamma1 - gamma2;
	t1 = sqrt(phi - vbx) - sqrtPhi;
	t2 = sqrt(phi * (phi - vbm)) - phi;
	double k2 = t0 * t1 / (2.0 * t2 + vbm);
	double k1 = gamma2 - 2.0 * k2 * sqrt(phi - vbm);

  if (k2 < 0.0)
    vbsc = (0.9 * (phi - pow(0.5 * k1 / k2, 2.0)));
  else
    vbsc = -30.0;

  if (vbsc + 3.0 > 0.0)
    vbsc = -3.0;
  else
    vbsc = -30.0;

	if(vbsc > vbm)
		vbsc = vbm;

	if((t0 = weff + k1w2) < 1e-8)
		t0 = 1e-8;
	double k1eff = k1 * (1 + k1w1 / t0);

	double vth0 = 0.7;
	double vfb = -1.0 * vth0 - phi - k1eff * sqrtPhi;

	double drout = 0.56;
	t1 = sqrt(EPSSI / EPSOX * tox * Xdep0);
	t0 = exp(-0.5 * drout * leff / t1);
	t2 = (t0 + 2.0 * t0 * t0);
	double thetaRout = pdibl1 * t2 + pdibl2;

	t0 = exp(-0.5 * dsub * leff / t1);
	double theta0vb0 = (t0 + 2.0 * t0 * t0);
	double Tox = 1.0e8 * tox;

	/* Poly Gate Si DEpeletion Effect */

	T0 = vfb + phi;
	if((ngate > 1e18) && (ngate < 1e25) && (x[1] > T0))
	{
		T1 = 1.0e6 * Charge_q * EPSSI * ngate / (cox * cox);
		T4 = sqrt(1.0 + 2.0 * (x[1] - T0) / T1);
		T2 = T1 * (T4 -1.0);
		T3 = 0.5 * T2 * T2 / T1;
		T7 = 1.12 - T3 - 0.05;
		T6 = sqrt(T7 * T7 + 0.224);
		T5 = 1.12 - 0.5 * (T7 + T6);
		Vgs_eff = x[1] - T5;
	}
	else
	{
		Vgs_eff = x[1];
	}

	V0 = vbi - phi;

	/* B/S built-in potential lowering calculation */

	Vbsmos = x[2];

	T0 = Vbsmos + 5 - 0.001;
	T1 = sqrt(T0 * T0 - 0.004 * (-5.0));
	T2 = (-5.0) + 0.5 * (T0 + T1);

	/* Vbsh */

	T0 = 1.5;
	T1 = T0 - T2 - 0.002;
	T3 = sqrt(T1 * T1 + 0.008 * T0);
	Vbsh = T0 - 0.5 * (T1 + T3);

	/* Vbseff */

	T0 = 0.95 * phi;
	T1 = T0 - Vbsh - 0.002;
	T2 = sqrt(T1 * T1 + 0.008 * T0);
	Vbseff = T0 - 0.5 * (T1 + T2);

	Phis = phi - Vbseff;
	sqrtPhis = sqrt(Phis);
	Xdep = Xdep0 * sqrtPhis / sqrtPhi;

	/* Calculation of Threshold voltage-vth */

	T3 = sqrt(Xdep);
	T0 = dvt2 * Vbseff;
  if (T0 + 0.5 > 0.0)
    T1 = 1.0 + T0;
  else
    T1 = (1.0 + 3.0 * T0) / (3.0 + 8.0 * T0);
	ltl = factor1 * T3 * T1;

	T0 = dvt2w * Vbseff;
  if (T0 + 0.5 > 0.0)
    T1 = 1.0 + T0;
  else
    T1 = (1.0 + 3.0 * T0) / (3.0 + 8.0 * T0);
	ltw = factor1 * T3 * T1;

	T0 = -0.5 * dvt1 * leff / ltl;
  if (T0 + EXPL_THRESHOLD > 0.0)
    Theta0 = exp(T0) * (1.0 + 2.0 * exp(T0));
  else
    Theta0 = MIN_EXPL * (1.0 + 2.0 * MIN_EXPL);

	thetavth = dvt0 * Theta0;
	Delt_vth = thetavth * V0;

	T0 = -0.5 * dvt1w * weff * leff / ltw;
  if (T0 + EXPL_THRESHOLD > 0.0)
    T2 = exp(T0) * (1.0 + 2.0 * exp(T0));
  else
    T2 = MIN_EXPL * (1.0 + 2.0 * MIN_EXPL);

	T0 = dvt0w * T2;
	DeltVthw = T0 * V0; //DeltVthw is T2 in Bsim3 code

	double TempRatioMinus1 = temp / tnom - 1.0;

	T0 = sqrt(1.0 + nlx / leff);
	T1 = (kt1 + kt1l / leff + kt2 * Vbseff);
	DeltVthtemp = k1eff * (T0 - 1.0) * sqrtPhi + T1 * TempRatioMinus1;

	double TMP2 = tox * phi / (weff + w0);

	T3 = eta0 + etab * Vbseff;
	if (T3 < 1.0e-4)
	{
		T9 = 1.0 / (3.0 - 2.0e4 * T3);
		T3 = (2.0e-4 - T3) * T9;
		T4 = T9 * T9 * etab;
	}

	DIBL_Sft = T3 * theta0vb0 * x[0];

	T9 = 2.2361 / sqrtPhi;
	sqrtPhisExt = sqrtPhis - T9 * (Vbsh - Vbseff);

	Vth = -1.0 * vth0 + k1eff * (sqrtPhisExt - sqrtPhi) - k2 * Vbseff - Delt_vth - DeltVthw +
	(k3 + k3b * Vbseff) * TMP2 + DeltVthtemp - DIBL_Sft;

	/* Calculate n */

	T2 = nfactor * EPSSI / Xdep;
	T3 = cdsc + cdscb * Vbseff + cdscd * x[0];
	T4 = (T2 + T3 * Theta0 + cit) / cox;

  if (T4 + 0.5 > 0.0)
    n = 1.0 + T4;
  else
    n = (1.0 + 3.0 * T4) / (3.0 + 8.0 * T4);

	/* Effective Vgst (Vgsteff) Calculation */

	Vgst = Vgs_eff - Vth;
	T10 = 2.0 * n * Vtm;
	VgstNVt = Vgst / T10;
	ExpArg = (2.0 * voff - Vgst) / T10;
	AD ExpVgst, dT2_dVg;

	if (VgstNVt > EXPL_THRESHOLD)
	{
		Vgsteff = Vgst;
	}
	else if (ExpArg > EXPL_THRESHOLD)
	{
		T0 = (Vgst - voff) / (n * Vtm);
		ExpVgst = exp(T0);
		Vgsteff = Vtm * cdep0 / cox * ExpVgst;
	}
	else
	{
		ExpVgst = exp(VgstNVt);
		T1 = T10 * log(1.0 + ExpVgst);
		T3 = (1.0 / temp);
		dT2_dVg = cox / (Vtm * cdep0) * exp(ExpArg);
		T2 = 1.0 - T10 * dT2_dVg;
		Vgsteff = T1 / T2;
	}
	Vgst2Vtm = Vgsteff + 2.0 * Vtm;

	/* Calculate Effective Channel Geometry */

	T9 = sqrtPhis - sqrtPhi;
	Weff = weff - (2.0 - nbc) * (dwg * Vgsteff + dwb * T9);

	if (Weff < 2.0e-8)
	{
		T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
		Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
	}

	T0 = prwg * Vgsteff + prwb * T9;
  if (T0 + 0.9 > 0.0)
    Rds = rds0 * (1.0 + T0);
  else
    Rds = rds0 * (0.8 + T0) * (1.0 / (17.0 + 20.0 * T0));

	/* Calculate Abulk */
	if (a0 == 0.0)
		Abulk0 = Abulk = 1.0;
	else
	{
		T10 = keta * Vbsh;
    if (T10 + 0.9 > 0.0)
      T11 = 1.0 / (1.0 + T10);
    else
      T11 = (17.0 + 20.0 * T10)*(1.0/(0.8 + T10));
		T10 = phi + ketas;
		T13 = (Vbsh * T11) / T10;

		if(T13 < 0.96)
		{
			T14 = 1/ sqrt(1-T13);
			T10 = 0.5 * T14 / (1-T13);
		}
		else
		{
			T11 = 1.0 / (1.0 - 1.043406 * T13);
			T14 = (6.00167 - 6.26044 * T13) * T11;
			T10 = 0.001742 * T11 * T11;
		}

		T10 = 0.5 * k1eff / sqrt(phi + ketas);
		T1 = T10 * T14;
		T9 = sqrt(xj * Xdep);
		T5 = leff / (leff + 2.0 * T9);
		T2 = (a0 * T5) + (b0 / (weff + b1));
		T6 = T5 * T5;
		T7 = T5 * T6;

		Abulk0 = 1 + T1 * T2;
		T8 = ags * a0 * T7;
		Abulk = Abulk0 + (-T1 * T8) * Vgsteff;
	}

  if (-Abulk0 + 0.01 > 0.0)
    Abulk0 = ((0.02 - Abulk0) / (3.0 - 200.0 * Abulk0));

  if (-Abulk + 0.01 > 0.0)
    Abulk = ((0.02 - Abulk) / (3.0 - 200.0 * Abulk));

	/* Mobility calculation */
	ua = 2.25e-9 + 4.31e-9 * (temp/tnom - 1.0);
	ub = 5.87e-19 + (-7.61e-18) * (temp/tnom - 1.0);
	uc = -4.65e-11 + (-5.6e-11) * (temp/tnom - 1.0);

	T0 = Vgsteff + Vth + Vth;
	T2 = ua + uc * Vbseff;
	T3 = T0 / tox;
	T5 = T3 * (T2 + ub * T3);

  if (T5 + 0.8 > 0.0)
    Denomi = 1.0 + T5;
  else
    Denomi = (0.6 + T5) / (7.0 + 10.0 * T5);

	ueff = u0temp / Denomi;

	/* Saturation Drain Voltage Vdsat */

	WVCox = Weff * vsattemp * cox;
	WVCoxRds = WVCox * Rds;

	Esat = 2.0 * vsattemp / ueff;
	EsatL = Esat * leff;
	double Lambda;

	if (a1 == 0.0)
	{
		Lambda = a2;
	}
	AbovVgst2Vtm = Abulk / Vgst2Vtm;

	if ((Rds == 0.0) && (Lambda == 1.0))
	{
		T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
		T3 = EsatL * Vgst2Vtm;
		Vdsat = T3 * T0;
	}
	else
	{
		T9 = Abulk * WVCoxRds;
		T7 = Vgst2Vtm * T9;
		T6 = Vgst2Vtm * WVCoxRds;
		T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0 / Lambda);

		T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk * EsatL + 3.0 * T7;
		T2 = Vgst2Vtm * (EsatL + 2.0 * T6);
		T3 = sqrt(T1 * T1 - 2.0 * T0 * T2);
		Vdsat = (T1 -T3) / T0;
	}

	/* Effective Vds (Vdseff) Calculation */

	T1 = Vdsat - x[0] - delta;
	T2 = sqrt(T1 * T1 + 4.0 * delta * Vdsat);
	T0 = T1 / T2;
	T3 = 2.0 * delta / T2;

	Vdseff = Vdsat - 0.5 * (T1 + T2);

	if (Vdseff > x[0])
		Vdseff = x[0];

	diffVds = x[0] - Vdseff;

	/* Calculate VAsat */

	AD TMP4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;
	T9 = WVCoxRds * Vgsteff;
	T0 = EsatL + Vdsat + 2.0 * T9 * TMP4;

	T1 = 2.0 / Lambda - 1.0 + (WVCoxRds * Abulk);
	Vasat = T0 / T1;

	/* Calculate VACLM */

  if (diffVds > 1.0e-10)
    VACLM = ((1.0 / (pclm * Abulk * litl))*(leff * (Abulk + (Vgsteff / EsatL)))) * diffVds;
  else
    VACLM = MAX_EXPL;

	/* Calculate VADIBL */

	if (thetaRout > 0.0)
	{
		T8 = Abulk * Vdsat;
		T0 = Vgst2Vtm * T8;
		T1 = Vgst2Vtm + T8;
		T9 = T1 * T1;
		T2 = thetaRout;
		VADIBL = (Vgst2Vtm - T0 / T1) / T2;
		T7 = pdiblb * Vbseff;
		if (T7 >= -0.9){
			VADIBL = VADIBL * (1.0 / (1.0 + T7));
		}
		else{
			VADIBL = VADIBL * (17.0 + 20.0 * T7) * (1.0 / (0.8 + T7));
		}
	}
	else
		VADIBL = MAX_EXPL;

	/* Calculate Va */

	T8 = pvag / EsatL;
	T9 = T8 * Vgsteff;
  if (T9 + 0.9 > 0.0)
    T0 = 1.0 + T9;
  else
    T0 = (0.8 + T9) * (1.0 / (17.0 + 20.0 * T9));

	T1 = VACLM * VADIBL / (VACLM + VADIBL);
	Va = Vasat + T0 * T1;

	/* Calculate Ids */

	CoxWovL = cox * Weff / leff;
	beta = ueff * CoxWovL;

	T0 = 1.0 - 0.5 * Abulk * Vdseff / Vgst2Vtm;

	fgche1 = Vgsteff * T0;
	T9 = Vdseff / EsatL;
	fgche2 = 1.0 + T9;
	gche = beta * fgche1 / fgche2;

	T0 = 1.0 + gche * Rds;
	T9 = Vdseff / T0;
	AD Idl = gche * T9;

	T9 = diffVds / Va;
	T0 = 1.0 + T9;

	AD Ids = Idl * T0 / nseg;

	/***************************/
	/* C-V Model for capMod=3 */

	double dtoxcv = 0.0;
	double agbcp = 0.0;
	double CoxWL  = cox * (weffCV / nseg *leffCV + agbcp);
	CoxWL = CoxWL * (tox / (tox- dtoxcv));
	double fbody = 1.0;
	double CoxWLb = fbody * cox * (weffCV / nseg * leffCVb + agbcp);
	CoxWLb = CoxWLb * (tox / (tox - dtoxcv));

	/* vfbzb calculation for capMod 3 */

	double k1ox = k1;

	T0 = -0.5 * dvt1w * weff * leff / (factor1 * sqrt(Xdep0));

  if (T0 + EXPL_THRESHOLD > 0.0)
    T2 = exp(T0) * (1.0 + 2.0 * exp(T0));
  else
    T2 = MIN_EXPL * (1.0 + 2.0 * MIN_EXPL);

	T0 = dvt0w * T2;
	T2 = T0 * (vbi - phi);

	T0 = -0.5 * dvt1 * leff / (factor1 * sqrt(Xdep0));

  if (T0 + EXPL_THRESHOLD > 0.0)
    T3 = exp(T0) * (1.0 + 2.0 * exp(T0));
  else
    T3 = MIN_EXPL * (1.0 + 2.0 * MIN_EXPL);

	T3 = dvt0 * T3 * (vbi - phi);

	T4 = tox * phi / (weff + w0);

	T0 = sqrt(1.0 + nlx / leff);
	T5 = k1eff * (T0 - 1.0) * sqrtPhi + (kt1 + kt1l / leff) * (temp/tnom - 1.0);
	T6 = 1.0 * vth0 - T2 -T3 + k3 * T4 + T5;
	AD vfbzb = T6 - phi - k1eff * sqrtPhi;

	/* Calculation for Vfbeff */

	V3 = vfbzb - Vgs_eff + Vbseff - DELTA_3;

  if (vfbzb > 0.0)
    T0 = sqrt(V3 * V3 + 4.0 * DELTA_3 * vfbzb);
  else
    T0 = sqrt(V3 * V3 - 4.0 *	DELTA_3 * vfbzb);

	Vfbeff = vfbzb - 0.5 * (V3 + T0);

	T0 = (Vgs_eff - Vbseff - vfbzb) / Tox;

	/* Calculation for Tcen */

	double ldeb = sqrt(EPSSI * Vtm0/(Charge_q * nch * 1.0e6)) / 3.0;
	double acde = acde * pow((nch / 2.0e16), -0.25);

	AD Tcen;
	T1 = T0 * acde;

	if ((-EXPL_THRESHOLD < T1) && (T1 < EXPL_THRESHOLD))
		Tcen = ldeb * exp(T1);
	else if (T1 <= -EXPL_THRESHOLD)
		Tcen = ldeb * MIN_EXPL;
	else
		Tcen = ldeb * MAX_EXPL;

	double LINK = 1.0e-3 * (tox - 0.0);
	V3 = ldeb - Tcen - LINK;
	V4 = sqrt(V3 * V3 + 4.0 * LINK * ldeb);
	Tcen = ldeb - 0.5 * (V3 + V4);

	AD Ccen;
	Ccen = EPSSI / Tcen;
	T2 = cox / (cox + Ccen);
	AD Coxeff;
	Coxeff = T2 * Ccen;

	/* Calculation for QoverlapCox */

	AD CoxWLcenb = CoxWLb * Coxeff / cox;
	Qac0 = CoxWLcenb * (Vfbeff - vfbzb);
	AD QovCox = Qac0 / Coxeff;

	T0 = 0.5 * k1eff;
	T3 = Vgs_eff - Vfbeff - Vbseff - Vgsteff;
	if (k1eff == 0.0)
  {
		T1 = 0.0;
		T2 = 0.0;
  }
	else if (T3 < 0.0)
  {
		T1 = T0 + T3 / k1eff;
		T2 = CoxWLcenb;
  }
	else
  {
		T1 = sqrt(T0 * T0 + T3);
		T2 = CoxWLcenb * T0 / T1;
  }

	Qsub0 = CoxWLcenb * k1eff * (T1 - T0);
	QovCox = Qsub0 / Coxeff;

	/* Calculation for Delta_phis */

  if (k1eff > 0.0)
    Denomi = moin * Vtm * k1eff * k1eff;
  else
    Denomi = 0.25 * moin * Vtm;

  if (k1ox > 0.0)
    T0 = k1eff * sqrtPhi;
  else
    T0 = 0.5 * sqrtPhi;

	AD DeltaPhi;
	T1 = 2.0 * T0 + Vgsteff;
	DeltaPhi = Vtm * log(1.0 + T1 * Vgsteff / Denomi);

	T3 = 4.0 * (Vth - vfbzb - phi);
	T2 = sqrt(T3 * T3 + 0.0001);
	T5 = 0.5 * (1 + T3 / T2);
	T4 = 0.5 * (T3 + T2);

	Tox = Tox + Tox;
	T0 = (Vgsteff + T4) / Tox;
	AD tmp = exp(0.7 * log(T0));
	T1 = 1.0 + tmp;
	T2 = 0.7 * tmp / (T0 * Tox);
	Tcen = 1.9e-9 / T1;

	Ccen = EPSSI / Tcen;
	T0 = cox / (cox + Ccen);
	Coxeff = T0 * Ccen;
	AD CoxWLcen = CoxWL * Coxeff / cox;
	CoxWLcenb = CoxWLb * Coxeff / cox;

	double abulkCVfactor = 1.0 + pow((clc / leff), cle);
	AbulkCV = Abulk0 * abulkCVfactor;
	VdsatCV = (Vgsteff - DeltaPhi) / AbulkCV;
	V4 = VdsatCV - x[0] - DELTA_4;
	T0 = sqrt(V4 * V4 + 4.0 * DELTA_4 * VdsatCV);
	VdseffCV = VdsatCV - 0.5 * (V4 + T0);

	T0 = AbulkCV * VdseffCV;
	T1 = Vgsteff - DeltaPhi;
	T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
	T3 = T0 / T2;
	T4 = 1.0 - 12.0 * T3 * T3;
	T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
	T6 = T5 * VdseffCV / AbulkCV;

	qgate = CoxWLcen * (T1 - T0 * (0.5 - T3));
	QovCox = qgate / Coxeff;

	qbulk = CoxWLcenb * (1.0 - AbulkCV) * (0.5 * VdseffCV - T0 * VdseffCV / T2);
	QovCox = qbulk / Coxeff;

	/* xpartition 40/60 partition */

	T2 = T2 / 12.0;
	T3 = 0.5 * CoxWLcen / (T2 * T2);
	T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0 * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;

	qsrc = -T3 * T4;
	QovCox = qsrc / Coxeff;

	/* Backgate charge */

	double aebcp = 0.0;
	double tbox = 3e-7;
	double Cbox = 3.453133e-11 / tbox;
	double CboxWL = kb1 * fbody * Cbox * (weffCV / nseg * leffCVbg + aebcp);

	double vfbb = -1.0 * Vtm * log(nch / nsub);
	AD Vesfb = x[2] - vfbb;
	AD Qe1 = CboxWL * (Vesfb - x[3]);

	qgate = qgate + Qac0 + Qsub0 - qbulk;
	qbody = qbulk - Qac0 - Qsub0 - Qe1;
	qsub = Qe1;
	qdrn = -(qgate + qbody + qsub + qsrc);

  // Calculate dynamic current contributions due to charge
  // iqd is dynamic contribution to drain current
  // iqg is dynamic contribution to gate current
  // iqs is dynamic contribution to source current
  // iqsb is dynamic contribution to substrate current
  AD iqd, iqg, iqs, iqsb;
  iqd = qdrn.fastAccessDx(0)*x[4] + qdrn.fastAccessDx(1)*x[5]
        + qdrn.fastAccessDx(2)*x[6] + qdrn.fastAccessDx(3)*x[7];
  iqg = qgate.fastAccessDx(0)*x[4] + qgate.fastAccessDx(1)*x[5]
        + qgate.fastAccessDx(2)*x[6] + qgate.fastAccessDx(3)*x[7];
  iqs = qsrc.fastAccessDx(0)*x[4] + qsrc.fastAccessDx(1)*x[5]
        + qsrc.fastAccessDx(2)*x[6] + qsrc.fastAccessDx(3)*x[7];
  iqsb = qsub.fastAccessDx(0)*x[4] + qsub.fastAccessDx(1)*x[5]
        + qsub.fastAccessDx(2)*x[6] + qsub.fastAccessDx(3)*x[7];

	flow[0] = -Ids - iqd;
	flow[1] = -iqg;
	flow[2] = Ids + iqs;
	flow[3] = Ids + iqsb;

	effort[0] = -(x[0] - x[3]);
	effort[1] = -(x[1] - x[3]);
	effort[2] = -(x[2] - x[3]);
	effort[3] = x[3];
}
