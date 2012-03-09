#include "Mospbsim3.h"

//Static members
const unsigned Mospbsim3::n_par = 98;

//Element information
ItemInfo Mospbsim3::einfo =
{
  "mospbsim3",
  "Mosfet model using bsim3 level, version 3.2.4",
  "Ramya Mohan",
  DEFAULT_ADDRESS"transistor>mosfet",
  "2003_05_15"
};

//Parameters
ParmInfo Mospbsim3::pinfo[] =
{
  {"l", "Length", TR_DOUBLE, false},
  {"w", "Width", TR_DOUBLE, false},
  {"tox", "Gate oxide thickness (m)", TR_DOUBLE, false},
  {"toxm", "Gate oxide thickness used in extraction",TR_DOUBLE, false},
  {"cdsc", "Drain/Source and channel coupling capacitance", TR_DOUBLE, false},
  {"cdscb", "Body-bias dependence of cdsc", TR_DOUBLE, false},
  {"cdscd", "Drain-bias dependence of cdsc", TR_DOUBLE, false},
  {"cit", "Interface state capacitance", TR_DOUBLE, false},
  {"nfactor", "Subthreshold swing coefficient", TR_DOUBLE, false},
  {"xj", "Junction depth (m)", TR_DOUBLE, false},
  {"vsat", "Saturationvelocity at tnom", TR_DOUBLE, false},
  {"at", "Temperature coefficient of vsat", TR_DOUBLE, false},
  {"a0", "Non-uniform depletion width effect coefficient", TR_DOUBLE, false},
  {"ags", "Gate bias coefficient of Abulk", TR_DOUBLE, false},
  {"a1", "Non-saturation effect coefficient", TR_DOUBLE, false},
  {"a2", "Non-saturation effect coefficient", TR_DOUBLE, false},
  {"keta", "Body-bias coefficient of non-uniform depletion width effect", TR_DOUBLE, false},
  {"nsub", "Substrate doping concentration", TR_DOUBLE, false},
  {"nch", "Channel doping concentration", TR_DOUBLE, false},
  {"ngate", "Poly-gate doping concentration", TR_DOUBLE, false},
  {"vbm", "Maximum body voltage", TR_DOUBLE, false},
  {"xt1", "Doping depth", TR_DOUBLE, false},
  {"kt1", "Temperature coefficient of Vth", TR_DOUBLE, false},
  {"kt1l", "Temperature coefficient of Vth", TR_DOUBLE, false},
  {"kt2", "Body-coefficient of kt1", TR_DOUBLE, false},
  {"k3", "Narrow width effect coefficient", TR_DOUBLE, false},
  {"k3b", "Body effect coefficient of k3", TR_DOUBLE, false},
  {"w0", "Narrow width effect parameter", TR_DOUBLE, false},
  {"nlx", "Lateral non-uniform doping effect", TR_DOUBLE, false},
  {"dvt0", "Short channel effect coefficient 0", TR_DOUBLE, false},
  {"dvt1", "Short channel effect coefficient 1", TR_DOUBLE, false},
  {"dvt2", "Short channel effect coefficient 2", TR_DOUBLE, false},
  {"dvt0w", "Narrow width effect coefficient 0", TR_DOUBLE, false},
  {"dvt1w", "Narrow width effect coefficient 1", TR_DOUBLE, false},
  {"dvt2w", "Narrow width effect coefficient 2", TR_DOUBLE, false},
  {"drout", "DIBL coefficient of output resistance", TR_DOUBLE, false},
  {"dsub", "DIBL coefficient in the subthreshold region", TR_DOUBLE, false},
  {"ua", "Linear gate dependence of mobility", TR_DOUBLE, false},
  {"ub", "Quadratic gate dependence of mobility", TR_DOUBLE, false},
  {"uc", "Body-bias dependence of mobility", TR_DOUBLE, false},
  {"u0", "Low-field mobility at Tnom", TR_DOUBLE, false},
  {"voff", "Threshold voltage offset", TR_DOUBLE, false},
  {"tnom", "Parameter measurement temperature", TR_DOUBLE, false},
  {"elm", "Non-quasi-static Elmore Constant Parameter", TR_DOUBLE, false},
  {"delta", "Effective Vds parameter", TR_DOUBLE, false},
  {"rdsw", "Sorce-drain resistance per width", TR_DOUBLE, false},
  {"prwg", "Gate-bias effect on parasitic resistance", TR_DOUBLE, false},
  {"prwb", "Body-effect on parasitic resistance", TR_DOUBLE, false},
  {"prt", "Temperature coefficient of parasitic resistance", TR_DOUBLE, false},
  {"eta0", "Subthreshold region DIBL coefficeint", TR_DOUBLE, false},
  {"etab", "Subthreshold region DIBL coefficeint", TR_DOUBLE, false},
  {"pclm", "Channel length modulation coefficient", TR_DOUBLE, false},
  {"pdibl1", "Drain-induced barrier lowering oefficient", TR_DOUBLE, false},
  {"pdibl2", "Drain-induced barrier lowering oefficient", TR_DOUBLE, false},
  {"pdiblb", "Body-effect on drain induced barrier lowering", TR_DOUBLE, false},
  {"pscbe1", "Substrate current body-effect coeffiecient", TR_DOUBLE, false},
  {"pscbe2", "Substrate current body-effect coeffiecient", TR_DOUBLE, false},
  {"pvag", "Gate dependence of output resistance parameter", TR_DOUBLE, false},
  {"vfb", "Flat band voltage", TR_DOUBLE, false},
  {"acde", "Exponential coefficient for finite charge thickness", TR_DOUBLE, false},
  {"moin", "Coefficient for gate-bias dependent surface potential", TR_DOUBLE, false},
  {"noff", "C-V turn-on/off parameter", TR_DOUBLE, false},
  {"voffcv", "C-V lateral shift parameter", TR_DOUBLE, false},
  {"lint", "Length reduction parameter", TR_DOUBLE, false},
  {"ll", "Length reduction parameter", TR_DOUBLE, false},
  {"llc", "Length reduction parameter for CV", TR_DOUBLE, false},
  {"lln", "Length reduction parameter", TR_DOUBLE, false},
  {"lw", "Length reduction parameter", TR_DOUBLE, false},
  {"lwc", "Length reduction parameter for CV", TR_DOUBLE, false},
  {"lwn", "Length reduction parameter", TR_DOUBLE, false},
  {"lwl", "Length reduction parameter", TR_DOUBLE, false},
  {"lwlc", "Length reduction parameter for CV", TR_DOUBLE, false},
  {"wr", "Width dependence of rds", TR_DOUBLE, false},
  {"wint", "Width reduction parameter", TR_DOUBLE, false},
  {"dwg", "Width reduction parameter", TR_DOUBLE, false},
  {"dwb", "Width reduction parameter", TR_DOUBLE, false},
  {"wl", "Width reduction parameter", TR_DOUBLE, false},
  {"wlc", "Width reduction parameter for CV", TR_DOUBLE, false},
  {"wln", "Width reduction parameter", TR_DOUBLE, false},
  {"ww", "Width reduction parameter", TR_DOUBLE, false},
  {"wwc", "Width reduction parameter for CV", TR_DOUBLE, false},
  {"wwn", "Width reduction parameter", TR_DOUBLE, false},
  {"wwl", "Width reduction parameter", TR_DOUBLE, false},
  {"wwlc", "Width reduction parameter for CV", TR_DOUBLE, false},
  {"b0", "Abulk narrow width parameter", TR_DOUBLE, false},
  {"b1", "Abulk narrow width parameter", TR_DOUBLE, false},
  {"clc", "Vdsat paramater for C-V model", TR_DOUBLE, false},
  {"cle", "Vdsat paramater for C-V model", TR_DOUBLE, false},
  {"alpha0", "Substrate current model parameter", TR_DOUBLE, false},
  {"alpha1", "Substrate current model parameter", TR_DOUBLE, false},
  {"beta0", "Diode limiting current", TR_DOUBLE, false},
  {"ute", "Temperature coefficient of mobility", TR_DOUBLE, false},
  {"k1", "First order body effect coefficient", TR_DOUBLE, false},
  {"k2", "Second order body effect coefficient", TR_DOUBLE, false},
  {"temp", "Circuit temperature", TR_DOUBLE, false},
  {"ua1", "Temperature coefficient for ua in m/V", TR_DOUBLE, false},
  {"ub1", "Temperature coefficient for ub in (m/V)^2", TR_DOUBLE, false},
  {"uc1", "Temperature coefficient for uc in m/V^2", TR_DOUBLE, false}
  //{"pmos", "True if PMOS", TR_BOOLEAN, false}
};

Mospbsim3::Mospbsim3(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  //Set default parameter values
  paramvalue[0] = &(l = 1e-6);
  paramvalue[1] = &(w = 1e-6);
  paramvalue[2] = &(tox = 150.0e-10);
  paramvalue[3] = &(toxm = 150.0e-10);
  paramvalue[4] = &(cdsc = 2.4e-4);
  paramvalue[5] = &(cdscb = 0.0);
  paramvalue[6] = &(cdscd = 0.0);
  paramvalue[7] = &(cit = 0.0);
  paramvalue[8] = &(nfactor = 1.0);
  paramvalue[9] = &(xj = .15e-6);
  paramvalue[10] = &(vsat = 8.0e4);
  paramvalue[11] = &(at = 3.3e4);
  paramvalue[12] = &(a0 = 1.0);
  paramvalue[13] = &(ags = 0.0);
  paramvalue[14] = &(a1 = 0.0);
  paramvalue[15] = &(a2 = 1.0);
  paramvalue[16] = &(keta = -0.047);
  paramvalue[17] = &(nsub = 6.0e16);
  paramvalue[18] = &(nch = 1.7e17);
  paramvalue[19] = &(ngate = 0.0);
  paramvalue[20] = &(vbm = -3.0);
  paramvalue[21] = &(xt1 = 1.55e-7);
  paramvalue[22] = &(kt1 = -0.11);
  paramvalue[23] = &(kt1l = 0.0);
  paramvalue[24] = &(kt2 = 0.022);
  paramvalue[25] = &(k3 = 80.0);
  paramvalue[26] = &(k3b = 0.0);
  paramvalue[27] = &(w0 = 2.5e-6);
  paramvalue[28] = &(nlx = 1.74e-7);
  paramvalue[29] = &(dvt0 = 2.2);
  paramvalue[30] = &(dvt1 = 0.53);
  paramvalue[31] = &(dvt2 = -0.032);
  paramvalue[32] = &(dvt0w = 0.0);
  paramvalue[33] = &(dvt1w = 5.3e6);
  paramvalue[34] = &(dvt2w = -0.032);
  paramvalue[35] = &(drout = 0.56);
  paramvalue[36] = &(dsub = drout);
  paramvalue[37] = &(ua = 2.25e-9);
  paramvalue[38] = &(ub = 5.87e-19);
  paramvalue[39] = &(uc = -4.65e-11);
  paramvalue[40] = &(u0= 0.025);
  paramvalue[41] = &(voff = -0.08);
  paramvalue[42] = &(tnom = 300.15);
  paramvalue[43] = &(elm = 5.0);
  paramvalue[44] = &(delta = 0.01);
  paramvalue[45] = &(rdsw = 0.0);
  paramvalue[46] = &(prwg = 0.0);
  paramvalue[47] = &(prwb = 0.0);
  paramvalue[48] = &(prt = 0.0);
  paramvalue[49] = &(eta0 = 0.08);
  paramvalue[50] = &(etab = -0.07);
  paramvalue[51] = &(pclm = 1.3);
  paramvalue[52] = &(pdibl1 = .39);
  paramvalue[53] = &(pdibl2 = 0.0086);
  paramvalue[54] = &(pdiblb = 0.0);
  paramvalue[55] = &(pscbe1 = 4.24e8);
  paramvalue[56] = &(pscbe2 = 1.0e-5);
  paramvalue[57] = &(pvag = 0.0);
  paramvalue[58] = &(vfb = -1.0);
  paramvalue[59] = &(acde = 1.0);
  paramvalue[60] = &(moin = 15.0);
  paramvalue[61] = &(noff = 1.0);
  paramvalue[62] = &(voffcv = 0.0);
  paramvalue[63] = &(lint = 0.0);
  paramvalue[64] = &(ll = 0.0);
  paramvalue[65] = &(llc = ll);
  paramvalue[66] = &(lln = 1.0);
  paramvalue[67] = &(lw = 0.0);
  paramvalue[68] = &(lwc = lw);
  paramvalue[69] = &(lwn = 1.0);
  paramvalue[70] = &(lwl = 0.0);
  paramvalue[71] = &(lwlc = lwl);
  paramvalue[72] = &(wr = 1.0);
  paramvalue[73] = &(wint = 0.0);
  paramvalue[74] = &(dwg = 0.0);
  paramvalue[75] = &(dwb = 0.0);
  paramvalue[76] = &(wl = 0.0);
  paramvalue[77] = &(wlc = wl);
  paramvalue[78] = &(wln = 1.0);
  paramvalue[79] = &(ww = 0.0);
  paramvalue[80] = &(wwc = ww);
  paramvalue[81] = &(wwn = 1.0);
  paramvalue[82] = &(wwl = 0.0);
  paramvalue[83] = &(wwlc = wwl);
  paramvalue[84] = &(b0 = 0.0);
  paramvalue[85] = &(b1 = 0.0);
  paramvalue[86] = &(clc = 0.1e-6);
  paramvalue[87] = &(cle = 0.6);
  paramvalue[88] = &(alpha0 = 0.0);
  paramvalue[89] = &(alpha1 = 0.0);
  paramvalue[90] = &(beta0 = 30.0);
  paramvalue[91] = &(ute = -1.5);
  paramvalue[92] = &(k1 = 0.53);
  paramvalue[93] = &(k2 = -0.0186);
  paramvalue[94] = &(temp = 373.15);
  paramvalue[95] = &(ua1 = 4.31e-9);
  paramvalue[96] = &(ub1 = -7.61e-18);
  paramvalue[97] = &(uc1 = -5.6e-11);
  //paramvalue[98] = &(pmos=false);

  //Set number of terminals
  setNumTerms(4);

  //Set number of state variables
  setNumberOfStates(3);

  //Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void Mospbsim3::init() throw(string&)
{
  DenseIntVector var(3);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  initializeAD(var, var);
}

void Mospbsim3::eval(AD * x, AD * effort, AD * flow)
{
  AD vbsc, V0, Vbseff, Phis, sqrtPhis, Xdep, lt1, ltw, Theta0, thetavth;
  AD Delt_vth, dDIBL_Sft_dVd, DIBL_Sft, temp_tmp4, Vth, n, Vgs_eff, Vgst;
  AD VgstNVt, ExpArg, Vgsteff, Weff, Rds, rds0, Abulk0, Abulk, Denomi, ueff;
  AD Esat, Vgst2Vtm, EsatL, WVCox, WVCoxRds, Lambda, Vdsat, Vdseff;
  AD diffVds, Vasat, VACLM, VADIBL, Va, VASCBE, CoxWovL, beta, fgche1;
  AD fgche2, gche, Idl, Idsa, Isub, Ids, vfbzb;
  AD T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10;
  double t0, t1, t2, t3, t4, t11, tmp1, tmp2, tmp3, tmp4, TMP2;

	#define Kb 1.3806226e-23
	#define KboQ 8.617087e-5  /* Kb / q  where q = 1.60219e-19 */
	#define EPSOX 3.453133e-11
	#define EPSSI 1.03594e-10 //this is in meters --> 1.03594e-8 F/cm
	#define PI 3.141592654
	#define MAX_EXP 5.834617425e14
	#define MIN_EXP 1.713908431e-15
	#define EXP_THRESHOLD 34.0
	#define Charge_q 1.60219e-19

	double factor1 = sqrt(EPSSI / EPSOX * tox);
	double Vtm0 = KboQ * tnom;
	double Eg0 = 1.16 - (7.02e-4 * tnom * tnom) / (tnom + 1108.0);
	double ni = 1.45e10 * (tnom / 300.15) * sqrt(tnom / 300.15) * exp(21.5565981 - Eg0 / (2.0 * Vtm0));

	double e0 = epsilon0;
	double esi = 11.7 * e0;
	double vt = KboQ * tnom;

	double Vtm = KboQ * temp;//Temp=tnom

	double ldrn = l;
	double wdrn = w;

	t0 = pow(ldrn, lln);
	t1 = pow(wdrn, lwn);
	tmp1 = ll / t0 + lw / t1 + lwl / (t0 * t1);
	double dl = lint + tmp1;
	tmp2 = llc / t0 + lwc / t1 + lwlc / (t0 * t1);
	double dlc = dlc + tmp2;

	t2 = pow(ldrn, wln);
	t3 = pow(wdrn, wwn);
	tmp3 = wl / t2 + ww / t3 + wwl / (t2 * t3);
	double dw = wint + tmp3;
	tmp4 = wlc / t2 + wwc / t3 + wwlc / (t2 * t3);
	double dwc = dwc + tmp4;

	double leff = l - 2.0 * dl;
	double weff = w - 2.0 * dw;

	double leffCV = l - 2.0 * dlc;

	t4 = (temp/tnom - 1.0);//t4 =TempRatio below

	double vsattemp = vsat - at * t4;
	rds0 = (rdsw + prt * t4) / pow(weff * 1e6, wr);
	t11 = leffCV * leffCV;

	double cox = 3.453133e-11 / tox;

	double phi = 2.0 * Vtm0 * log(nch / ni);
	double sqrtPhi = sqrt(phi);
	double phis3 = sqrtPhi * phi;

	double Xdep0 = sqrt(2.0 * EPSSI / (Charge_q * nch * 1.0e6)) * sqrtPhi;
	double litl = sqrt(3.0 * xj * tox);
	double vbi = Vtm0 * log(1.0e20 * nch / (ni * ni));
	double cdep0 = sqrt(Charge_q * EPSSI * nch * 1.0e6 / 2.0 / phi);

	double vbx = phi - 7.7348e-4  * nch * xt1 * xt1;
	if (vbx > 0.0)
		vbx = -vbx;

	double gamma1 = 5.753e-12 * sqrt(nch) / cox;
	double gamma2 = 5.753e-12 * sqrt(nsub) / cox;

	t0 = gamma1 - gamma2;
	t1 = sqrt(phi - vbx) - sqrtPhi;
	t2 = sqrt(phi * (phi - vbm)) - phi;
	//k2 = t0 * t1 / (2.0 * t2 + vbm);
	//k1 = gamma2 - 2.0 * k2 * sqrt(phi - vbm);
	double vth0 = vfb + phi + k1 * sqrtPhi;
	//cout << "vth0 = " << vth0 << endl;
	double k1ox = k1 * tox / toxm;
	double k2ox = k2 * tox / toxm;

	t1 = sqrt(EPSSI / EPSOX * tox * Xdep0);
	t0 = exp(-0.5 * dsub * leff / t1);
	double theta0vb0 = (t0 + 2.0 * t0 * t0);
	double Tox = 1.0e8 * tox;

	//Calculation of vbsc(Vbc) and Vbseff
  if (-k2 > 0.0)
    vbsc = (0.9 * (phi - pow(0.5 * k1 / k2, 2.0)));
  else
    vbsc = -30.0;

  if (vbsc + 3.0 > 0.0)
    vbsc = -3.0;
  else
    vbsc = -30.0;

	T0 = x[2] - vbsc - 0.001; //Vbs = x[2]
	T1 = sqrt(T0 * T0 - 0.004 * vbsc);
	Vbseff = vbsc + 0.5 * (T0 + T1);

  if (-Vbseff + x[2] > 0.0)
    Vbseff = x[2] + zero;

	//Calculation of Phis, sqrtPhis and Xdep
  if (Vbseff > 0.0)
    Phis = phi * phi / (phi + Vbseff);
  else
    Phis = phi - Vbseff;

  if (Vbseff > 0.0)
    sqrtPhis = phis3 / (phi + 0.5 * Vbseff);
  else
    sqrtPhis = sqrt(Phis);

	Xdep = Xdep0 * sqrtPhis / sqrtPhi;

	//Calculation of Threshold voltage-vth
	T3 = sqrt(Xdep);
	V0 = vbi - phi;

  if (dvt2 * Vbseff + 0.5 > 0.0)
    T1 = 1.0 + dvt2 * Vbseff;
  else
    T1 = (1.0 + 3.0 * dvt2 * Vbseff) / (3.0 + 8.0 * dvt2 * Vbseff);

	lt1 = factor1 * T3 * T1;

  if (dvt2w * Vbseff + 0.5 > 0.0)
    T1 = 1.0 + dvt2w * Vbseff;
  else
    T1 = (1.0 + 3.0 * dvt2w * Vbseff) / (3.0 + 8.0 * dvt2w * Vbseff);

	ltw = factor1 * T3 * T1;

  if (-0.5 * dvt1 * leff / lt1 + EXP_THRESHOLD > 0.0)
    Theta0 = exp(-0.5 * dvt1 * leff / lt1) * (1.0 + 2.0 * exp(-0.5 * dvt1 * leff / lt1));
  else
    Theta0 = MIN_EXP * (1.0 + 2.0 * MIN_EXP);

	thetavth = dvt0 * Theta0;
	Delt_vth = thetavth * V0;

  if (-0.5 * dvt1w * weff * leff / ltw + EXP_THRESHOLD > 0.0)
    T2 = exp(-0.5 * dvt1w * weff * leff / ltw) * (1.0 + 2.0 * exp(-0.5 * dvt1w * weff * leff / ltw));

	T0 = dvt0w * T2;
	T2 = T0 * V0;
	double TempRatio = temp / tnom - 1.0;

	T0 = sqrt(1.0 + nlx / leff);
	T1 = k1ox * (T0 - 1.0) * sqrtPhi + (kt1 + kt1l / leff + kt2 * Vbseff) * TempRatio;
	TMP2 = tox * phi / (weff + w0);

	T3 = eta0 + etab * Vbseff;
  if (-T3 + 1.0e-4 > 0.0)
    T4 = 1.0 / (3.0 - 2.0e4 * eta0 + etab * Vbseff);
  else
    T4 = 1.0;

	dDIBL_Sft_dVd = T3 * theta0vb0;
	DIBL_Sft = dDIBL_Sft_dVd * x[0];

	Vth = vth0 - k1 * sqrtPhi + k1ox * sqrtPhis
	- k2ox * Vbseff - Delt_vth - T2 + (k3 + k3b * Vbseff) * TMP2 + T1 - DIBL_Sft;

	//Calculate n
	AD temp_tmp2 = nfactor * EPSSI / Xdep;
	AD temp_tmp3 = cdsc + cdscb * Vbseff + cdscd * x[0];
	temp_tmp4 = (temp_tmp2 + temp_tmp3 * Theta0 + cit) / cox;

  if (temp_tmp4 + 0.5 > 0.0)
    n = 1.0 + temp_tmp4;
  else
    n = (1.0 + 3.0 * temp_tmp4) * (1.0 / (3.0 + 8.0 * temp_tmp4));

	//Poly Gate Si Depletion Effect
	Vgs_eff = x[1];
	Vgst = Vgs_eff - Vth;//not in Nikhil's code

	//Effective Vgst (Vgsteff) Calculation
	T10 = 2.0 * n * Vtm;
	VgstNVt = Vgst / T10;
	ExpArg = (2.0 * (-0.08) - Vgst) / T10;

	AD ExpVgst, dT2_dVg;

	T1 = T10 * log(1.0 + (exp(VgstNVt)));
	T2 = 1.0 - T10 * (-cox / (Vtm * cdep0) * exp(ExpArg));

  if (VgstNVt > EXP_THRESHOLD)
    Vgsteff = Vgst;
  else
    Vgsteff = (T1/T2);

  if (ExpArg > EXP_THRESHOLD)
    Vgsteff = ((Vtm * cdep0) / (cox * (exp((Vgst - (-0.08)) / (n * Vtm)))));
  else
    Vgsteff = (T1/T2);

	T3 = T2 * T2;

	//Calculate Effective Channel Geometry
	T9 = sqrtPhis - sqrtPhi;
	Weff = weff - 2.0 * (dwg * Vgsteff + dwb * T9);

  if (-Weff + 2.0e-8 > 0.0)
    Weff = (2.0e-8 * (4.0e-8 - Weff) * T0);

	T0 = prwg * Vgsteff + prwb * (sqrtPhis - sqrtPhi);

  if (T0 + 0.9 > 0.0)
    Rds = rds0 * (1.0 + T0);
  else
    Rds = rds0 * (0.8 +T0) * (1.0 / (17.0 + 20.0 * T0));

	//Calculate Abulk
	T1 = 0.5 * k1ox / sqrtPhis;

	T9 = sqrt(xj * Xdep);
	T5 = leff / (leff + 2.0 * T9);
	T2 = (a0 * T5) + (b0 / (weff + b1));
	T6 = T5 * T5;
	T7 = T5 * T6;

	Abulk0 = 1.0 + T1 * T2;

	T8 = ags * a0 * T7;
	Abulk = Abulk0 + (-T1 * T8) * Vgsteff;

  if (-Abulk0 + 0.1 > 0.0)
    Abulk0 = ((0.2 - Abulk0) / (3.0 - 20.0 * Abulk0));

  if (-Abulk + 0.1 > 0.0)
    Abulk = ((0.2 - Abulk) / (3.0 - 20.0 * Abulk));
	T2 = keta * Vbseff;

  if (T2 + 0.9 > 0.0)
    T0 = 1.0 / (1.0 + T2);
  else
    T0 = (17.0 + 20.0 * T2) / (0.8 + T2);
	Abulk *= T0;
	Abulk0 *= T0;

	//Mobility calculation
	ua = 2.25e-9 + 4.31e-9 * (temp/tnom - 1.0);
	ub = 5.87e-19 + (-7.61e-18) * (temp/tnom - 1.0);
	uc = -4.65e-11 + (-5.6e-11) * (temp/tnom - 1.0);

	// ua = 6.47e-9 + 3.31e-10 * (temp/tnom - 1.0);
	// ub = 4.23e-18 + 2.61e-19 * (temp/tnom - 1.0);
	// uc = -4.706281e-11 + -3.42e-10 * (temp/tnom - 1.0);

	T0 = Vgsteff + Vth + Vth;
	T2 = ua + uc * Vbseff;
	T3 = T0 / tox;
	T5 = T3 * (T2 + ub * T3);

  if (T5 + 0.8 > 0.0)
    Denomi = 1.0 + T5;
  else
    Denomi = (0.6 + T5) / (7.0 + 10.0 * T5);

  if (u0 > 1.0)
		u0 = u0 / 1.0e4;

	double u0temp = u0 * pow((temp/tnom), ute); //TRatio = temp/tnom

	ueff = u0temp / Denomi;

	Esat = 2.0 * vsattemp / ueff;

	//Saturation Drain Voltage Vdsat

	WVCox = Weff * vsattemp * cox;
	WVCoxRds = WVCox * Rds;

	EsatL = Esat * leff;

	//calculation of Lambda
  Lambda = a2;

	Vgst2Vtm = Vgsteff + 2.0 * Vtm;

	T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
	T3 = EsatL * Vgst2Vtm;
	Vdsat = T3 * T0;

	//Effective Vds(Vdseff) Calculation
	T1 = Vdsat - x[0] - delta;
	T2 = sqrt(T1 * T1 + 4.0 * delta * Vdsat);
	T0 = T1 / T2;
	T3 = 2.0 * delta / T2;

	Vdseff = Vdsat - 0.5 * (T1 + T2);

	//Added to eliminate non-zero Vdseff at Vds=0.0
	if (x[0] == 0.0)
		Vdseff = 0.0;

	//Calculate Vasat
	T6 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;//T6=tmp4
	T9 = WVCoxRds * Vgsteff;//expanded
	T0 = EsatL + Vdsat + 2.0 * T9 * T6;

	T9 = WVCoxRds * Abulk;
	T1 = 2.0 / Lambda - 1.0 + T9;

	Vasat = T0 / T1;

  if (Vdseff > x[0])
    Vdseff = x[0] - zero;
	diffVds = x[0] - Vdseff;

	//Calculate VACLM
  if (diffVds > 1.0e-10)
    VACLM = ((1.0 / (pclm * Abulk * litl))*(leff * (Abulk + (Vgsteff / EsatL)))) * diffVds;
  else
    VACLM = MAX_EXP;

	//Calculate VADIBL
	T1 = sqrt(EPSSI / EPSOX * tox * Xdep0);
	T0 = exp(-0.5 * drout * leff / T1);
	T2 = (T0 + 2.0 * T0 * T0);
	AD thetaRout = pdibl1 * T2 + pdibl2;//drout, pdibl1, pdibl2 are given

  if (thetaRout > 0.0)
    VADIBL = ((Vgst2Vtm - (Vgst2Vtm * Abulk * Vdsat) / (Vgst2Vtm + (Abulk * Vdsat))) / (thetaRout));
  else
    VADIBL = MAX_EXP;

  if ((pdiblb * Vbseff) + 0.9 > 0.0)
    VADIBL = (VADIBL * (1.0 / (1.0 + (pdiblb * Vbseff))));
  else
    VADIBL = (VADIBL * ((17.0 + 20.0 * (pdiblb * Vbseff)) * (1.0 / (0.8 + (pdiblb * Vbseff)))));

	//Calculate Va
	T8 = pvag / EsatL;
	T9 = T8 * Vgsteff;

  if (T9 + 0.9 > 0.0)
    T0 = 1.0 + T9;
  else
    T0 = (0.8 + T9) * (1.0 / (17.0 + 20.0 * T9));

	T3 = VACLM + VADIBL;//tmp3 = T3
	T1 = VACLM * VADIBL / T3;
	Va = Vasat + T0 * T1;

	//Calculate VASCBE
  if (diffVds > (pscbe1 * litl/EXP_THRESHOLD))
    VASCBE = leff * exp(pscbe1 * litl/diffVds) / pscbe2;
  else
    VASCBE = MAX_EXP * leff / pscbe2;

	//Calculate Ids
	CoxWovL = cox * Weff / leff;
	beta = ueff * CoxWovL;

	T0 = 1.0 - 0.5 * Abulk * Vdseff / Vgst2Vtm;

	fgche1 = Vgsteff * T0;

	T9 = Vdseff / EsatL;
	fgche2 = 1.0 + T9;

	gche = beta * fgche1 / fgche2;

	T0 = 1.0 + gche * Rds;
	T9 = Vdseff / T0;
	Idl = gche * T9;

	T9 = diffVds / Va;
	T0 = 1.0 + T9;
	Idsa = Idl * T0;

	T9 = diffVds / VASCBE;
	T0 = 1.0 + T9;
	Ids = Idsa * T0;

	//Substrate current begins
	T1 = alpha0 + alpha1 * leff;

  if (diffVds > (beta0/EXP_THRESHOLD))
    T1 = (((alpha0 + alpha1 * leff)/leff) * diffVds * exp(-beta0 / diffVds));
  else
    T1 = ((((alpha0 + alpha1 * leff)/leff) * MIN_EXP) * diffVds);

  if ((T1 <= 0.0) || (beta0 <= 0.0))
    Isub = 0.0;
  else
    Isub = (T2 * Idsa);

	//calculation of vfbzb
	T0 = -0.5 * dvt1w * weff * leff / (factor1 * sqrt(Xdep0));

  if (T0 + EXP_THRESHOLD > 0.0)
    T2 = exp(T0) * (1.0 + 2.0 * exp(t0));
	else
    T2 = MIN_EXP * (1.0 + 2.0 * MIN_EXP);

	T0 = dvt0w * T2;
	T2 = T0 * (vbi - phi);

	T0 = -0.5 * dvt1 * leff / (factor1 * sqrt(Xdep0));

  if (T0 + EXP_THRESHOLD > 0.0)
    T3 = exp(T0) * (1.0 + 2.0 * exp(T0));
  else
    T3 = MIN_EXP * (1.0 + 2.0 * MIN_EXP);

	T3 = dvt0 * T3 * (vbi - phi);
	T4 = tox * phi / (weff + w0);
	T0 = sqrt(1.0 + nlx / leff);
	T5 = k1ox * (T0 - 1.0) * sqrtPhi + (kt1 + kt1l / leff) * (temp/tnom - 1.0);
	T6 = vth0 - T2 -T3 + k3 * T4 + T5;
	vfbzb = T6 - phi - k1 * sqrtPhi;

	//Calculation for VbseffCV
	AD VbseffCV;
  if (Vbseff > 0.0)
    VbseffCV = phi - Phis;
  else
    VbseffCV = Vbseff;

	//Calculation for VgsteffCV
	T0 = n * noff * vt;
	T1 = (Vgs_eff - Vth) / T0;

  if (T1 > EXP_THRESHOLD)
    Vgsteff = Vgs_eff - Vth - voffcv;
  else
    Vgsteff = T0 * log(1.0 + exp(T1));

  if (-T1 > EXP_THRESHOLD)
    Vgsteff = T0 * log(1.0 + MIN_EXP);
  else
    Vgsteff = T0 * log(1.0 + exp(T1));

	//Calculation for Vfbeff
	AD V3;
	V3 = vfbzb - Vgs_eff + VbseffCV - 0.02;

  if (vfbzb > 0.0)
    T0 = sqrt(V3 * V3 + 4.0 * 0.02 * vfbzb);
  else
    T0 = sqrt(V3 * V3 - 4.0 * 0.02 * vfbzb);

	AD Vfbeff;
	Vfbeff = vfbzb - 0.5 * (V3 + T0);

	Tox = 1.0e8 * tox;
	T0 = (Vgs_eff - VbseffCV - vfbzb) / Tox;

	//Calculation for Tcen

	double ldeb = sqrt(esi * vt/(Charge_q * nch * 1.0e6)) / 3.0;

	AD Tcen;
	T1 = T0 * acde;

  if (EXP_THRESHOLD + T1 > 0.0)
    Tcen = ldeb * exp(T1);
  else
    Tcen = ldeb * MAX_EXP;

  if (EXP_THRESHOLD > T1)
    Tcen = ldeb * exp(T1);
  else
    Tcen = ldeb * MAX_EXP;

  if (-T1 > EXP_THRESHOLD)
    Tcen = ldeb * MIN_EXP;
  else
    Tcen = ldeb * MAX_EXP;

	V3 = ldeb - Tcen - 1.0e-3 * Tox;

	AD V4;
	V4 = sqrt(V3 * V3 + 4.0 * 1.0e-3 * Tox * ldeb);

	Tcen = ldeb - 0.5 * (V3 + V4);

	AD Ccen;
	Ccen = esi / Tcen;

	AD Coxeff;
	Coxeff = Ccen * cox / (Ccen + cox);

	//Calculation for QoverlapCox
	AD QovCox, Qac0, CoxWLcen, Qsub0;

	CoxWLcen = cox * Weff * leff * Coxeff / cox;
	Qac0 = CoxWLcen * (Vfbeff - vfbzb);
	QovCox = Qac0 / Coxeff;

	T0 = 0.5 * k1ox;
	T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
	T1 = sqrt(T0 * T0 + T3);
	T2 = CoxWLcen * T0 / T1;

  if (-T3 > 0.0)
    T1 = (T0 + T3 / k1ox);
  else
    T1 = (sqrt(T0 * T0 + T3));

  if (-T3 > 0.0)
    T2 = CoxWLcen;
  else
    T2 = (CoxWLcen * T0 / T1);

	Qsub0 = CoxWLcen * k1ox * (T1 - T0);
	QovCox = Qsub0 / Coxeff;

	//Calculation for Delta_phis
  if (k1ox > 0.0)
  {
    T2 = moin * vt * k1ox * k1ox;
    T0 = k1ox * sqrt(phi);
  }
  else
  {
    T2 = 0.25 * moin * vt;
    T0 = 0.5 * sqrt(phi);
  }

	AD DeltaPhi;

	T1 = 2.0 * T0 + Vgsteff;

	DeltaPhi = vt * log(1.0 + T1 * Vgsteff / T2);

	//The calculation for Tcen must be done once more

	T0 = (Vgsteff + 4.0*(vth0 - vfb - phi))/ (2.0 * Tox);
	T1 = 1.0 + exp(0.7 * log(T0));
	T2 = 0.7 * exp(0.7 * log(T0)) / (T0 * 2.0 * Tox);
	Tcen = 1.9e-9 / T1;

	Ccen = esi / Tcen;
	Coxeff = Ccen * cox / (Ccen + cox);
	CoxWLcen = cox * Weff * leff * Coxeff / cox;

	AD AbulkCV;
	AbulkCV = Abulk0 * (1.0 + pow((clc/leff),cle));

	AD VdsatCV;
	VdsatCV = (Vgsteff - DeltaPhi) / AbulkCV;

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
	T1 = Vgsteff - DeltaPhi;
	T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
	T3 = T0 / T2;
	T4 = 1.0 - 12.0 * T3 * T3;
	T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
	T6 = T5 * VdseffCV / AbulkCV;

	AD qgate;
	qgate = CoxWLcen * (T1 - T0 * (0.5 - T3));

	AD qbulk;
	qbulk = CoxWLcen * (1.0 - AbulkCV) * (0.5*VdseffCV - T0*VdseffCV/T2);

	QovCox = qbulk / Coxeff;

	T2 = T2 / 12.0;
	T3 = 0.5 * CoxWLcen / (T2 * T2);
	T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0 * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;

	AD qsrc;
	qsrc = -T3 * T4;

	qgate += Qac0 + Qsub0 - qbulk;
	qbulk -= (Qac0 + Qsub0);
	AD qdrn;
	qdrn = -(qbulk + qgate + qsrc);

  // Calculate dynamic current contributions due to charge
  // iqd is dynamic contribution to drain current
  // iqg is dynamic contribution to gate current
  // iqs is dynamic contribution to source current
  AD iqd, iqg, iqs;
  iqd = qdrn.fastAccessDx(0)*x[3] + qdrn.fastAccessDx(1)*x[4] + qdrn.fastAccessDx(2)*x[5];
  iqg = qgate.fastAccessDx(0)*x[3] + qgate.fastAccessDx(1)*x[4] + qgate.fastAccessDx(2)*x[5];
  iqs = qsrc.fastAccessDx(0)*x[3] + qsrc.fastAccessDx(1)*x[4] + qsrc.fastAccessDx(2)*x[5];

	flow[0] = -(Ids + Isub + iqd);
	flow[1] = -iqg;
	flow[2] = Ids - iqs;

	effort[0] = -(x[0] - x[2]);
	effort[1] = -(x[1] - x[2]);
	effort[2] = x[2];
}
