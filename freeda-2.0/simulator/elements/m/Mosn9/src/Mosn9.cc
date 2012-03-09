#include "Mosn9.h"

//Static members
const unsigned Mosn9::n_par = 125;

//element information
ItemInfo Mosn9::einfo =
{
  "mosn9",
  "Philips MOS9 model",
  "Yogesh Kuldip Ajit Xuemin Dapeng Xin",
  DEFAULT_ADDRESS"transistor>mosfet",
  "2003_05_15"
};

//Parameter information
ParmInfo Mosn9::pinfo[]=
{
  //Parameters of the electrical model
  {"LEVEL","Model Level",TR_DOUBLE,false},
  {"VT0","Threshold voltage at zero back-bias for the actual transistor at the actual temperature(volts) ",TR_DOUBLE,false},
  {"K0","low back-bias body factor for the actual transistor(V^1/2)",TR_DOUBLE,false},
  {"K","high back-bias body factor for the actual transistor(V^1/2)",TR_DOUBLE,false},
  {"PHIB","Surface potential at strong inversion for the actual transistor at the actual temperature(volts)",TR_DOUBLE,false},
  {"VSBX","Transition voltage for the dual-k-factor model for the actual transistor(volts)",TR_DOUBLE,false},
  {"BET","gain factor for the actual transistor at the actual temperature(A/V^2)",TR_DOUBLE,false},
  {"THE1","coefficient of the mobility reduction due to the gate induced field for the actual transistor at the actual temperature(1/V)",TR_DOUBLE,false},
  {"THE2","coefficient of the mobility reduction due to the back-biased for the actual transistor at the actual temperature(1/V^0.5)",TR_DOUBLE,false},
  {"THE3","coefficient of the mobility reduction due to the lateral field for the actual transistor at the actual temperature(1/V)",TR_DOUBLE,false},
  {"GAM1","coefficient for the drain induced threshold shift for large gate drive for the actual transistor (V^(1-ETADS))",TR_DOUBLE,false},
  {"ETADS","exponent of the vds dependence of GAM1 for the actual transistor",TR_DOUBLE,false},
  {"ALP","factor of the channel length modulation for the actual transistor",TR_DOUBLE,false},
  {"VP","characteristic voltage of the channel length modulation for the actual transistor (V)",TR_DOUBLE,false},
  {"GAM00","coefficient for the drain induced threshold shift at zero gate drive for the actual transistor",TR_DOUBLE,false},
  {"ETAGAM","exponent of the back-bias dependence of gam0 for the actual transistor",TR_DOUBLE,false},
  {"M0","factor of the subthreshold slope for the actual transistor at the actual temperature",TR_DOUBLE,false},
  {"ETAM","exponent of the back-bias dependence of m for the actual transistor",TR_DOUBLE,false},
  {"PHIT","thermal voltage at the actual temperature(volts)",TR_DOUBLE,false},
  {"ZET1","weak inversion correction factor for the actual transistor",TR_DOUBLE,false},
  {"VSBT","limiting voltage of vsb dependence of m and gam0 for the actual transistor(volts)",TR_DOUBLE,false},
  {"A1","factor of the weak avalanche current for the actual transistor",TR_DOUBLE,false},
  {"A2","exponent of the weak avalanche current for the actual transistor(volts)",TR_DOUBLE,false},
  {"A3","factor of the drain source voltage above which weak avalanche occurs for the actual transistor",TR_DOUBLE,false},
  {"COX","gate to channel capacitance for the actual transistor(farads)",TR_DOUBLE,false},
  {"CGDO","gate drain overlap capacitance for the actual transistor(farads)",TR_DOUBLE,false},
  {"CGSO","gate source overlap capacitance for the actual transistor(farads)",TR_DOUBLE,false},
  {"NT","coefficient of thermal noise for the actual transistor(J)",TR_DOUBLE,false},
  {"NFMOD","switch that selects either new or old flicker noise model for the actual transistor",TR_DOUBLE,false},
  {"NF","flicker noise coefficient for the actual transistor(for NFMOD=0)(V^2)",TR_DOUBLE,false},
  {"NFA","first coefficient of flicker noise coefficient for the actual transistor(for NFMOD=1)(1/V)(1/m^4))",TR_DOUBLE,false},
  {"NFB","second coefficient of flicker noise coefficient for the actual transistor(for NFMOD=1)(1/V)(1/m^2)",TR_DOUBLE,false},
  {"NFC","third coefficient of flicker noise coefficient for the actual transistor(for NFMOD=1)(1/V)",TR_DOUBLE,false},
  {"TOX","thickness of gate oxide layer (meters)",TR_DOUBLE,false},
  {"MULT","number of devices operating in parallel",TR_DOUBLE,false},
  //Parameters of the geometrical model
  {"LER", "Effective channel length of the reference transistor (m)", TR_DOUBLE, false},
  {"WER", "Effective channel width of the reference transistor (m)",TR_DOUBLE, false},
  {"LVAR", "Difference between the actual and the programmed poly-silicon gate length (m)", TR_DOUBLE, false},
  {"LAP", "Effective channel length reduction per side due to the lateral diffusion of the source/drain dopant ions (m)", TR_DOUBLE, false},
  {"WVAR", "Difference between the actual and the programmed field-oxide opening (m)", TR_DOUBLE, false},
  {"WOT", " Effective reduction of the channel width per side due to the lateral diffusion of the channel-stop dopant ions (m)", TR_DOUBLE, false},
  {"TR", "Temperature at which the parameters for the reference transistor have been determined (0C)", TR_DOUBLE, false},
  {"VT0R", "Threshold voltage at zero back-bias for the reference transistor at the reference temperature (V)", TR_DOUBLE, false},
  {"STVT0", "Coefficient of the temperature dependence VT0 (VK-1)", TR_DOUBLE, false},
  {"SLVT0", "Coefficient of the length dependence of VT0 (Vm)", TR_DOUBLE, false},
  {"SL2VT0", "Second coefficient of the length dependence of VT0 (Vm2)", TR_DOUBLE, false},
  {"SL3VT0", "Third coefficient of the length dependence of VT0 (V)", TR_DOUBLE, false},
  {"SWVT0", "Coefficient of the width dependence of VT0 (Vm)", TR_DOUBLE, false},
  {"K0R", "Low-backbias body factor for the reference transistor (V1/2)",TR_DOUBLE, false},
  {"SLK0", "Coefficient of the length dependence of K0 (V1/2m)",TR_DOUBLE, false},
  {"SL2K0", "Second coefficient of the length dependence of K0 (V1/2m2)", TR_DOUBLE, false},
  {"SWK0", "Coefficient of the width dependence of K0 (V1/2m)", TR_DOUBLE, false},
  {"KR", "High-backbias body factor for the reference transistor (V1/2)", TR_DOUBLE, false},
  {"SLK", "Coefficient of the length dependence of K (V1/2m)", TR_DOUBLE, false},
  {"SL2K", "Second coefficient of the length dependence of K (V1/2m2)", TR_DOUBLE, false},
  {"SWK", "Coefficient of the width dependence of K (V1/2m)", TR_DOUBLE, false},
  {"PHIBR", "Surface potential at strong inversion for the reference transistor at the reference temperature (V)", TR_DOUBLE, false},
  {"VSBXR", "Transition voltage for the dual-k-factor model for the reference transistor (V)", TR_DOUBLE, false},
  {"SLVSBX", "Coefficient of the length dependence of VSBX (Vm)", TR_DOUBLE, false},
  {"SWVSBX", "Coefficient of the width dependence VSBX (Vm)", TR_DOUBLE, false},
  {"BETSQ", "Gain factor for an infinite square transistor at the reference temperature (AV-2)", TR_DOUBLE, false},
  {"ETABET", "Exponent of the temperature dependence of the gain factor (-)", TR_DOUBLE, false},
  {"LP1", "Characteristic length of first profile (m)", TR_DOUBLE, false},
  {"FBET1", "Relative mobility decrease due to first profile (-)", TR_DOUBLE, false},
  {"LP2", "Characteristic length of second profile (m)", TR_DOUBLE, false},
  {"FBET2", "Relative mobility decrease due to second profile (-)", TR_DOUBLE, false},
  {"THE1R", "Coefficient of the mobility reduction due to the gate-induced field for the reference transistor at the reference temperature (V-1)", TR_DOUBLE, false},
  {"STTHE1R", "Coefficient of the temperature dependence of THE1 for the reference transistor (V-1K-1)", TR_DOUBLE, false},
  {"SLTHE1R", "Coefficient of the length dependence of THE1 at the reference temperature (V-1m)", TR_DOUBLE, false},
  {"STLTHE1", "Coefficient of the temperature dependence of the length dependence of THE 1 (V-1mK-1)", TR_DOUBLE, false},
  {"GTHE1", "Parameter that selects either the old ( ) or the new ( ) scaling gTHE0 = gTHE1 = rule of THE 1 (-)", TR_DOUBLE, false},
  {"SWTHE1", "Coefficient of the width dependence of THE1 (V-1m)", TR_DOUBLE, false},
  {"WDOG", "Characteristic drawn gate width, below which dogboning appears (m)", TR_DOUBLE, false},
  {"FTHE1", "Coefficient describing the geometry dependence of THE1 for W < WDOG (-)", TR_DOUBLE, false},
  {"THE2R", "Coefficient of the mobility reduction due to the back-bias for the reference transistor at the reference temperature (V-1/2)", TR_DOUBLE, false},
  {"STTHE2R", "Coefficient of the temperature dependence of THE2 for the reference transistor (V-1/2K-1)", TR_DOUBLE, false},
  {"SLTHE2R", "Coefficient of the length dependence of THE2 at the reference temperature (V-1/2m)", TR_DOUBLE, false},
  {"STLTHE2", "Coefficient of the temperature dependence of the length dependence of THE2 (V-1/2mK-1)", TR_DOUBLE, false},
  {"SWTHE2", "Coefficient of the width dependence of THE2 (V-1/2m)", TR_DOUBLE, false},
  {"THE3R", "Coefficient of the mobility reduction due to the lateral field for the reference transistor at the reference temperature (V-1)", TR_DOUBLE, false},
  {"STTHE3R", "Coefficient of the temperature dependence of THE3 for the reference temperature (V-1K-1)", TR_DOUBLE, false},
  {"SLTHE3R", "Coefficient of the length dependence of THE3 at the reference temperature (V-1m)", TR_DOUBLE, false},
  {"STLTHE3", "Coefficient of the temperature dependence of the length dependence of THE3 (V-1mK-1)", TR_DOUBLE, false},
  {"SWTHE3", "Coefficient of the width dependence of THE3 (V-1m)", TR_DOUBLE, false},
  {"GAM1R", "Coefficient for the drain induced threshold shift for large gate drive for the reference transistor (V(1-ETADS))", TR_DOUBLE, false},
  {"SLGAM1", "Coefficient of the length dependence of GAM1 (V(1-ETADS)m)", TR_DOUBLE, false},
  {"SWGAM1", "Coefficient of the width dependence of GAM1 (V(1-ETADS)m)", TR_DOUBLE, false},
  {"ETADSR", "Exponent of the VDS dependence of GAM1 for the reference transistor (-)", TR_DOUBLE, false},
  {"ALPR", "Factor of the channel-length modulation for the reference transistor (-)", TR_DOUBLE, false},
  {"ETAALP", "Exponent of the length dependence of ALP (-)", TR_DOUBLE, false},
  {"SLALP", "Coefficient of the length dependence of ALP (mETAALP)", TR_DOUBLE, false},
  {"SWALP", "Coefficient of the width dependence of ALP (m)", TR_DOUBLE, false},
  {"VPR", "Characteristic voltage of the channel length modulation for the reference transistor (V)", TR_DOUBLE, false},
  {"GAM00R", "Coefficient of the drain induced threshold shift at zero gate drive for the reference transistor (-)", TR_DOUBLE, false},
  {"SLGAM00", "Coefficient of the length dependence of GAM00 (m2)", TR_DOUBLE, false},
  {"SL2GAM00", "Second coefficient of the length dependence of GAM00 (-)", TR_DOUBLE, false},
  {"ETAGAMR", "Exponent of the back-bias dependence of GAM0 for the reference transistor (-)", TR_DOUBLE, false},
  {"M0R", "Factor of the sub threshold slope for the reference transistor at the reference temperature (-)", TR_DOUBLE, false},
  {"STM0", "Coefficient of the temperature dependence of m0 (K-1)", TR_DOUBLE, false},
  {"SLM0", "Coefficient of the length dependence of m0 (m1/2)", TR_DOUBLE, false},
  {"ETAMR", "Exponent of the back-bias dependence of m for the reference transistor (-)", TR_DOUBLE, false},
  {"ZET1R", "Weak-inversion correction factor for the reference transistor (-)", TR_DOUBLE, false},
  {"ETAZET", "Exponent of the length dependence of ZET1 (-)", TR_DOUBLE, false},
  {"SLZET1", "Coefficient of the length dependence of ZET1 (mETAZET)", TR_DOUBLE, false},
  {"VSBTR", "Limiting voltage of the VSB dependence of m and GAM0 for the reference transistor (V)", TR_DOUBLE, false},
  {"SLVSBT", "Coefficient of the length dependence of VSBT (Vm)", TR_DOUBLE, false},
  {"A1R", "Factor of the weak-avalanche current for the reference transistor at the reference temperature (-)", TR_DOUBLE, false},
  {"STA1", "Coefficient of the temperature dependence of a1 (K-1)", TR_DOUBLE, false},
  {"SLA1", "Coefficient of the length dependence of a1 (m)", TR_DOUBLE, false},
  {"SWA1", "Coefficient of the width dependence of a1 (m)", TR_DOUBLE, false},
  {"A2R", "Exponent of the weak-avalanche current for the reference transistor (V)", TR_DOUBLE, false},
  {"SLA2", "Coefficient of the length dependence of a2 (Vm)", TR_DOUBLE, false},
  {"SWA2", "Coefficient of the width dependence of a2 (Vm)", TR_DOUBLE, false},
  {"A3R", "Factor of the drain-source voltage above which weak-avalanche occurs, for the reference transistor (-)", TR_DOUBLE, false},
  {"SLA3", "Coefficient of the length dependence of a3 (m)", TR_DOUBLE, false},
  {"SWA3", "Coefficient of the width dependence of a3 (m)", TR_DOUBLE, false},
  // {"TOX", "Thickness of the gate-oxide layer. tox is used for calculation
  // of 1/f noise and Cox, not for BET!!! (m)", TR_DOUBLE, false},
  {"COL", "Gate overlap capacitance per unit channel width (Fm-1)", TR_DOUBLE, false},
  {"NTR", "Coefficient of the thermal noise for the reference transistor (J)", TR_DOUBLE, false},
  // {"NFMOD", "Switch that selects either old or new flicker noise model (-)", TR_DOUBLE, false},
  {"NFR", "Flicker noise coefficient of the reference transistor (for NFMOD = 0) (V2)", TR_DOUBLE, false},
  {"NFAR", "First coefficient of the flicker noise of the reference transistor (for NFMOD = 1) (V-1m-4)", TR_DOUBLE, false},
  {"NFBR", "Second coefficient of the flicker noise of the reference transistor (for NFMOD = 1) (V-1m-2)", TR_DOUBLE, false},
  {"NFCR", "Third coefficient of the flicker noise of the reference transistor (for NFMOD = 1) () V-1", TR_DOUBLE, false},
  {"l", "Drawn channel length in the lay-out of the actual transistor (m)", TR_DOUBLE, false},
  {"w", "Drawn channel width in the lay-out of the actual transistor (m)", TR_DOUBLE, false},
  {"DTA", "Temperature offset of the device with respect to TA (0C)", TR_DOUBLE, false}
};


Mosn9::Mosn9(const string& iname) : ADInterface(&einfo,pinfo,n_par,iname)
{
  //Set default parameter values (electrical model) n-channel
  paramvalue[0] = &(LEVEL = 903);
  paramvalue[1] = &(VT0 = 7.099154e-1);
  paramvalue[2] = &(K0 = 6.478116e-1);
  paramvalue[3] = &(K = 4.280174e-1);
  paramvalue[4] = &(PHIB = 6.225999e-1);
  paramvalue[5] = &(VSBX = 6.599578e-1);
  paramvalue[6] = &(BET = 1.418789e-3);
  paramvalue[7] = &(THE1 = 1.923533e-1);
  paramvalue[8] = &(THE2 = 1.144632e-2);
  paramvalue[9] = &(THE3 = 1.381597e-1);
  paramvalue[10] = &(GAM1 = 1.476930e-1);
  paramvalue[11] = &(ETADS = 0.6);
  paramvalue[12] = &(ALP = 2.878165e-3);
  paramvalue[13] = &(VP = 3.338182e-1);
  paramvalue[14] = &(GAM00 = 1.861785e-2);
  paramvalue[15] = &(ETAGAM = 2);
  paramvalue[16] = &(M0 = 5.024606e-1);
  paramvalue[17] = &(ETAM = 2);
  paramvalue[18] = &(PHIT = 2.662680e-2);
  paramvalue[19] = &(ZET1 = 4.074464e-1);
  paramvalue[20] = &(VSBT = 2.025926);
  paramvalue[21] = &(A1 = 6.022073);
  paramvalue[22] = &(A2 = 3.801696e+1);
  paramvalue[23] = &(A3 = 6.407407e-1);
  paramvalue[24] = &(COX = 2.979787e-14);
  paramvalue[25] = &(CGDO = 6.392e-15);
  paramvalue[26] = &(CGSO = 6.392e-15);
  paramvalue[27] = &(NT = 2.563182e-20);
  paramvalue[28] = &(NFMOD = 0);
  paramvalue[29] = &(NF);
  paramvalue[30] = &(NFA = 7.15e+22);
  paramvalue[31] = &(NFB = 2.16e+7);
  paramvalue[32] = &(NFC = 0);
  paramvalue[33] = &(TOX = 25e-9);
  paramvalue[34] = &(MULT = 1);
  //Geometrical model. n-channel
  paramvalue[35] = &(LER = 1.10e-6);
  paramvalue[36] = &(WER = 20e-6);
  paramvalue[37] = &(LVAR = -0.220e-6);
  paramvalue[38] = &(LAP = 0.100e-6);
  paramvalue[39] = &(WVAR = -0.025e-6);
  paramvalue[40] = &(WOT = 0.000e-6);
  paramvalue[41] = &(TR = 21.00);
  paramvalue[42] = &(VT0R = 0.730);
  paramvalue[43] = &(STVT0 = -1.20e-3);
  paramvalue[44] = &(SLVT0 = -0.135e-6);
  paramvalue[45] = &(SL2VT0 = 0.0);
  paramvalue[46] = &(SL3VT0 = 0.0);
  paramvalue[47] = &(SWVT0 = 0.130e-6);
  paramvalue[48] = &(K0R = 0.650);
  paramvalue[49] = &(SLK0 = -0.130e-6);
  paramvalue[50] = &(SL2K0 = 0.0);
  paramvalue[51] = &(SWK0 = 0.002e-6);
  paramvalue[52] = &(KR = 0.110);
  paramvalue[53] = &(SLK = -0.280e-6);
  paramvalue[54] = &(SL2K = 0.0);
  paramvalue[55] = &(SWK = 0.275e-6);
  paramvalue[56] = &(PHIBR = 0.650);
  paramvalue[57] = &(VSBXR = 0.660);
  paramvalue[58] = &(SLVSBX = 0.000e-6);
  paramvalue[59] = &(SWVSBX = -0.675e-6);
  paramvalue[60] = &(BETSQ = 83.00e-6);
  paramvalue[61] = &(ETABET = 1.600);
  paramvalue[62] = &(LP1 = 1.0e-6);
  paramvalue[63] = &(FBET1 = 0.0);
  paramvalue[64] = &(LP2 = 1.0e-8);
  paramvalue[65] = &(FBET2 = 0.0);
  paramvalue[66] = &(THE1R = 0.190);
  paramvalue[67] = &(STTHE1R = 0.000e-3);
  paramvalue[68] = &(SLTHE1R = 0.140e-6);
  paramvalue[69] = &(STLTHE1 = 0.000e-3);
  paramvalue[70] = &(GTHE1 = 0.0);
  paramvalue[71] = &(SWTHE1 = -0.058e-6);
  paramvalue[72] = &(WDOG = 0.0);
  paramvalue[73] = &(FTHE1 = 0.0);
  paramvalue[74] = &(THE2R = 0.012);
  paramvalue[75] = &(STTHE2R = 0.000e-9);
  paramvalue[76] = &(SLTHE2R = -0.033e-6);
  paramvalue[77] = &(STLTHE2 = 0.000e-3);
  paramvalue[78] = &(SWTHE2 = 0.030e-6);
  paramvalue[79] = &(THE3R = 0.145);
  paramvalue[80] = &(STTHE3R = -0.660e-3);
  paramvalue[81] = &(SLTHE3R = 0.185e-6);
  paramvalue[82] = &(STLTHE3 = -0.620e-9);
  paramvalue[83] = &(SWTHE3 = 0.020e-6);
  paramvalue[84] = &(GAM1R = 0.145);
  paramvalue[85] = &(SLGAM1 = 0.160e-6);
  paramvalue[86] = &(SWGAM1 = -0.010e-6);
  paramvalue[87] = &(ETADSR = 0.600);
  paramvalue[88] = &(ALPR = 0.003);
  paramvalue[89] = &(ETAALP = 0.15);
  paramvalue[90] = &(SLALP = -5.65e-3);
  paramvalue[91] = &(SWALP = 1.67e-9);
  paramvalue[92] = &(VPR = 0.340);
  paramvalue[93] = &(GAM00R = 0.018);
  paramvalue[94] = &(SLGAM00 = 20.00e-15);
  paramvalue[95] = &(SL2GAM00 = 0.0);
  paramvalue[96] = &(ETAGAMR = 2.0);
  paramvalue[97] = &(M0R = 0.500);
  paramvalue[98] = &(STM0 = 0.000);
  paramvalue[99] = &(SLM0 = 0.280e-3);
  paramvalue[100] = &(ETAMR = 2.0);
  paramvalue[101] = &(ZET1R = 0.420);
  paramvalue[102] = &(ETAZET = 0.17);
  paramvalue[103] = &(SLZET1 = -0.390);
  paramvalue[104] = &(VSBTR = 2.10);
  paramvalue[105] = &(SLVSBT = -4.40e-6);
  paramvalue[106] = &(A1R = 6.00);
  paramvalue[107] = &(STA1 = 0.000);
  paramvalue[108] = &(SLA1 = 1.30e-6);
  paramvalue[109] = &(SWA1 = 3.00e-6);
  paramvalue[110] = &(A2R = 38.0);
  paramvalue[111] = &(SLA2 = 1.00e-6);
  paramvalue[112] = &(SWA2 = 2.00e-6);
  paramvalue[113] = &(A3R = 0.650);
  paramvalue[114] = &(SLA3 = -0.550e-6);
  paramvalue[115] = &(SWA3 = 0.000);
  //paramvalue[82] = &(TOX = 25.0e-9);
  paramvalue[116] = &(COL = 0.320e-9);
  paramvalue[117] = &(NTR = 0.244e-19);
  //paramvalue[118] = &(NFMOD = 0);
  paramvalue[118] = &(NFR = 0);
  paramvalue[119] = &(NFAR = 7.15e22);
  paramvalue[120] = &(NFBR = 2.16e7);
  paramvalue[121] = &(NFCR = 0);
  paramvalue[122] = &(l = 1.5e-6);
  paramvalue[123] = &(w = 20.0e-6);
  paramvalue[124] = &(DTA = 0.0);

  //Set the number of terminals
  setNumTerms(4);

  //Set number of state variables
  // vgs,vds,vsb
  setNumberOfStates(3);

  //Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
};

void Mosn9::init() throw (string&)
{
  DenseIntVector var(3);
  var[0]=0;
  var[1]=1;
  var[2]=2;
  initializeAD(var,var);
}

void Mosn9:: eval(AD * x, AD * effort, AD * flow)
{
  //state variables
  AD vds,vgs,vsb;

  vds=x[0]+1e-15;
  vgs=x[1];
  vsb=x[2];

  //Set constants values
  double T0 = 273.15; //Offset for conversion from Celsius to Kelvin temperature scale
  double BOLTZMAN = 1.3806226e-23; //Boltzmann constant(JK-1)
  double Q = 1.6021918e-19; //Elementary unit charge(C)
  double EPIOX = 3.453143800e-11; //Absolute permittivity of the oxide layer(Fm-1)
  double eps2 = 0.1;

  //=========Geometrical scaling and temperature scaling=========
  double LE, WE ;
  double TKR,TKD;
  double TA = 25.0;

  //---------Calculation of Transistor Geometry---------
  LE = l + LVAR - 2.0 * LAP;
  WE = w + WVAR - 2.0 * WOT;

  // clip levels
  if (LE < 1e-9)
    LE = 1e-9;
  if (WE < 1e-9)
    WE = 1e-9;
  //---------end Transistor Geometry---------

  //---------Calculation of Transistor Temperature---------
  TKR = T0 + TR;
  TKD = T0 + TA + DTA;
  //---------end Transistor Temperature---------

  //---------Calculation of Threshold-Voltage Parameters---------
  double VT0TILDA, GPE, GPER;
  double STPHIB;

  VT0TILDA = VT0R + (TKD-TKR)*STVT0;
  GPE = FBET1 * LP1/LE * (1.0 - exp(-LE/LP1)) + FBET2 * LP2/LE * (1.0 - exp(-LE/LP2));
  GPER = FBET1 * LP1/LER * (1.0 - exp(-LER/LP1)) + FBET2 * LP2/LER
         * (1.0 - exp(-LER/LP2));

  // clip levels
  if (1.0 + GPE < 1e-15)
    GPE = 1e-15 - 1.0;
  if (1.0 + GPER < 1e-15)
    GPER = 1e-15 - 1.0;

  VT0 = VT0TILDA + (1.0/LE - 1.0/LER) * SLVT0 + (1.0/pow(LE,2) - 1.0/pow(LER,2))
        * SL2VT0 + (GPE-GPER) * SL3VT0 + (1.0/WE - 1.0/WER) * SWVT0;
  K0 = K0R + (1.0/LE - 1.0/LER) * SLK0 + (1.0/pow(LE,2) - 1.0/pow(LER,2))
       * SL2K0 + (1.0/WE - 1.0/WER) * SWK0;
  K = KR + (1.0/LE - 1.0/LER) * SLK + (1.0/pow(LE,2) - 1.0/pow(LER,2))
      * SL2K + (1.0/WE - 1.0/WER) * SWK;

  STPHIB = (PHIBR - 1.13 - 2.5e-4 * TKR)/300.0;
  PHIB = PHIBR + (TKD - TKR) * STPHIB;
  VSBX = VSBXR + (1.0/LE - 1.0/LER) * SLVSBX + (1.0/WE - 1.0/WER) * SWVSBX;

  AD tmp = K0 * sqrt(hyp1(-VSBX,eps2))/sqrt(VSBX + PHIB);
  if (K < tmp.val())
    K = tmp.val();
  //---------end Threshold-Voltage Parameters---------

  //---------Calculation of Channel-Current Parameters---------
  double BETTILDA, THE1TILDA;
  double SLTHE1, WEDOG;

  BETTILDA = BETSQ * pow((TKR/TKD),ETABET);
  BET = BETTILDA/(1.0 + GPE) * WE/LE;

  THE1TILDA = THE1R + (TKD-TKR) * STTHE1R;
  SLTHE1 = SLTHE1R + (TKD-TKR) * STLTHE1;
  WEDOG = WDOG + WVAR - 2.0 * WOT;

  if (w > WDOG)
  {
    THE1 = THE1TILDA+(1.0/(LE * (1.0 + GTHE1*GPE))
           - 1.0/(LER*(1.0 + GTHE1*GPER)))*SLTHE1 + (1.0/WE - 1.0/WER)*SWTHE1;
  }
  else
  {
    THE1 = THE1TILDA + (1.0/(LE*(1.0 + GTHE1*GPE))
           - 1.0/(LER * (1.0+GTHE1*GPER))) * SLTHE1
           + (1.0/WE - 1.0/WER) * SWTHE1 + (WE/WEDOG - 1.0)
           * FTHE1/(LE * (1.0 + GTHE1*GPE)) * SLTHE1;
  }

  //Calculation of Drain-Feedback Parameters
  GAM1 = GAM1R + (1.0 / LE - 1.0 / LER) * SLGAM1 + (1.0 / WE - 1.0 / WER) * SWGAM1;
  ETADS = ETADSR;
  ALP = ALPR + (1.0 / pow(LE,ETAALP) - 1.0 / pow(LER,ETAALP) * SLALP
        + (1.0 / WE - 1.0 / WER)) * SWALP;
  VP = VPR * (LE / LER);
  //---------end Drain-Feedback Parameters---------

  //---------Calculation of Sub-Threshold Parameters---------
  double M0TILDA;

  GAM00 = GAM00R + (1.0/pow(LE,2) - 1.0/pow(LER,2)) * SLGAM00 + (GPE - GPER) * SL2GAM00;

  ETAGAM = ETAGAMR;
  M0TILDA = M0R + (TKD - TKR) * STM0;
  M0 = M0TILDA + (1.0/sqrt(LE) - 1.0/sqrt(LER)) * SLM0;
  ETAM = ETAMR;
  PHIT = BOLTZMAN * TKD / Q;

  ZET1 = ZET1R + (1.0/pow(LE,ETAZET) - 1.0/pow(LER,ETAZET)) * SLZET1;
  VSBT = VSBTR + (1.0/LE - 1.0/LER) * SLVSBT;

  //---------end Sub-Threshold Parameters---------

  //---------Calculation of Week-Avalanche Parameters---------
  double A1TILDA;

  A1TILDA = A1R + (TKD - TKR) * STA1;
  A1 = A1TILDA + (1.0/LE - 1.0/LER) * SLA1 + (1.0/WE - 1.0/WER) * SWA1;
  A2 = A2R + (1.0/LE - 1.0/LER) * SLA2 + (1.0/WE - 1.0/WER) * SWA2;
  A3 = A3R + (1.0/LE - 1.0/LER) * SLA3 + (1.0/WE - 1.0/WER) * SWA3;
	//---------end Sub-Threshold Parameters---------

  //---------Calculation of Charge Parameters---------
  COX = EPIOX * WE * LE / TOX;
  CGDO = WE * COL;
  CGSO = WE * COL;

  if (THE1 < 0 )
    THE1 = 0.0;
  if (THE2 < 0 )
    THE2 = 0.0;
  if (GAM00 < 0)
    GAM00 = 0.0;
  if (THE3 < 0 )
    THE3 = 0.0;
  if (GAM1 < 0)
    GAM1 = 0.0;
  //---------end Charge Parameters---------
  //=========end Geometrical scaling and temperature scaling=========

  //=========Model Equations=========
  //---------DC current model---------
  double lambda10 = 0.9;
  double eps1 = 1e-2;
  double us0, ust, usx;
  AD h1, us;
  h1 = hyp1(vsb + lambda10 * PHIB,eps1) + (1.0 - lambda10) * PHIB;
  us = sqrt(h1);
  us0 = sqrt(PHIB);
  ust = sqrt(VSBT + PHIB);
  usx = sqrt(VSBX + PHIB);

  AD deltavt0, vt1;
  deltavt0 = K * (sqrt(hyp4(vsb,VSBX,eps2) + pow((K/K0),2) * pow(usx,2))
             - (K/K0) * usx) + K0 * (sqrt(h1 - hyp4(vsb,VSBX,eps2)) - us0);
  vt1 = VT0+deltavt0;

  double eps3 = 0.01;
  AD us1,gam0;
  us1 = hyp2(us,ust,eps3);
  gam0 = GAM00 * pow((us1/us0),ETAGAM);

  double eps4 = 5e-4;
  double lambda1 = 0.1;
  double lambda2 = 1e-4;
  AD vgt1, deltavt1, vgtx, vt2;
  vgt1 = hyp1(vgs-vt1,eps4);
  vgtx = 1.0/sqrt(2.0);
  deltavt1 = (-gam0 - (GAM1 * pow(vds + lambda2,ETADS-1) - gam0)
             * (pow(vgt1,2) / (pow(vgtx,2) + pow(vgt1,2))))
             *(pow(vds,2) / (vds+lambda1));
  vt2 = vt1 + deltavt1;

  double lambda3 = 1e-8;
  double lambda7 = 37.0;
  AD m, vgt2, vgta, g1, vgt3;
  m = 1.0 + M0 * pow((us0/us1),ETAM);
  vgt2 = vgs - vt2;
  vgta = 2.0 * m * PHIT * lambda7;
  g1 = exp((vgt2/(2.0 * m * PHIT)));
  if (vgt2 > vgta)
    vgt3 = vgt2 + lambda3;
  else
    vgt3 = 2.0 * m * PHIT * log(1.0 + g1) + lambda3;
  double lambda4 = 0.3, lambda5 = 0.1;
  AD del1;
  del1 = (lambda4/us) * (K + ((K0 - K)*(pow(VSBX,2))/(pow(VSBX,2)
         + pow(lambda5 * (vgt1 + vsb),2))));

  double lambda9 =0.1;
  double eps8= 0.001;
  AD vdss1;
  vdss1 = (vgt3/(1.0 + del1)) * (2.0/(1.0+sqrt(lambda9
          + hyp1(1.0 - lambda9 + (2.0 * THE3 * vgt3/(1.0 + del1)),eps8))));

  double lambda6 = 0.3;
  double vdssx = 1.0;
  AD vds1;
  AD eps5;
  eps5=lambda6 * (vdss1/(vdss1 + vdssx));
  vds1=hyp5(vds,vdss1,eps5);

  AD g2;
  g2=1.0 + ALP * log(1.0 + (vds-vds1)/VP);

  AD g3;
  if (vgt2 > vgta)
    g3 = g2;
  else
    g3 = (ZET1*(1.0 - exp(-(vds/PHIT))) + (g1*g2))/(g1 + (1.0/ZET1));

  AD ids;
  ids = BET * g3 * ((vgt3 * vds1 - ((1.0 + del1)/2)
        *pow(vds1,2))/((1.0 + THE1 * vgt1 + THE2 * (us-us0))
        * (lambda9 + hyp1(1.0 - lambda9 + THE3 * vds1,eps8))));

  AD vdsa;
  vdsa=A3 * vdss1;

  AD iavl;
  if (vds > vdsa)
    iavl = ids * A1 * exp(-A2/(vds - vdsa));
  else
    iavl = 0.0;
  //---------end DC current model---------

  //---------Charge model---------
  AD vdb;
  vdb=vds + vsb;

  AD h2;
  h2=hyp1(vdb + lambda10 * PHIB,eps1) + (1.0 - lambda10) * PHIB;

  AD deltavt0d, vt1d;
  deltavt0d = K * (sqrt(hyp4(vdb,VSBX,eps2) + pow((K/K0),2) * pow(usx,2))
              - (K/K0) * usx) + K0 * (sqrt(h2-hyp4(vdb,VSBX,eps2)) - us0);
  vt1d = VT0 + deltavt0d;

  AD diff_delta_vt0_vsb, diff_delta_vt0_vgs, diff_delta_vt0_vds;
  AD diff_delta_vt1_vsb, diff_delta_vt1_vgs, diff_delta_vt1_vds;
  AD diff_vt2_vsb, diff_vt2_vgs, diff_vt2_vds;

  diff_delta_vt0_vsb = 0.5 * K * (diff_hyp4(vsb,VSBX,eps2)/sqrt(hyp4(vsb,VSBX,eps2)
                      + pow((K/K0)*usx,2))) +	0.5*K0*(1.0/(sqrt(h1 - hyp4(vsb,VSBX,eps2))))*
                      (diff_hyp1(vsb + lambda10 * PHIB,eps1) - diff_hyp4(vsb,VSBX,eps2));

  AD diff_vgt1_vsb;
  AD diff_gam0_vsb;

  diff_vgt1_vsb = -diff_hyp1(vsb-vt1,eps4) * diff_delta_vt0_vsb;
  diff_gam0_vsb = GAM00 * (ETAGAM - 1.0) * (1.0/us0) * pow(us1,ETAGAM - 1.0)
                  * diff_hyp2(us,ust,eps3) * (0.5/sqrt(h1))
                  * diff_hyp1(vsb + lambda10 * PHIB,eps1);

  diff_delta_vt1_vsb = (pow(vds,2)/(vds+lambda1))
                       * (-diff_gam0_vsb - (GAM1*pow(vds + lambda2,ETADS - 1.0) - gam0)
                       * 2.0 * vgt1 * diff_vgt1_vsb * pow(vgtx,2)/(pow(pow(vgtx,2) + pow(vgt1,2),2))
                       + diff_gam0_vsb * pow(vgt1,2)/(pow(vgtx,2) + pow(vgt1,2)));

  diff_vt2_vsb = diff_delta_vt0_vsb + diff_delta_vt1_vsb;

  diff_vt2_vgs = -((2.0 * vgt1 * pow(vgtx,2))/(pow(pow(vgtx,2) + pow(vgt1,2),2)))
                 * diff_hyp1(vgs - vt1,eps4) * (pow(vds,2)/(vds + lambda1))
                 * (GAM1 * pow(vds + lambda2,ETADS - 1.0) - gam0);

  diff_vt2_vds = (-gam0 - (GAM1 * pow(vds + lambda2,ETADS-1) - gam0)
                 * (pow(vgt1,2)/(pow(vgtx,2) + pow(vgt1,2))))
                 * ((pow(vds,2) + 2.0 * lambda1 * vds)/pow(vds + lambda1,2))
                 + (-GAM1 * (ETADS-1.0) * (pow(vds + lambda2,ETADS-2)
                 * pow(vgt1,2)/(pow(vgtx,2) + pow(vgt1,2))
                 * (pow(vds,2)/(vds + lambda1))));

  AD del2;
  del2 = diff_vt2_vsb - diff_vt2_vgs - diff_vt2_vds;

  AD temp=2.0 * m * PHIT;
  AD diff_vgt3_vsb, diff_vgt3_vgs, diff_vgt3_vds;
  AD diff_m_vsb, diff_g1_vsb;

  diff_m_vsb = -M0 * ETAM * pow(us0/us1,ETAM-1) * us0 * pow(us1,-2)
               * diff_hyp2(us,ust,eps3) * (0.5/sqrt(h1))
               * diff_hyp1(vsb + lambda10 * PHIB,eps1);

  diff_g1_vsb = exp(vgt2/temp) * (temp * diff_vt2_vsb - diff_m_vsb * vgt2
                * 2.0 * PHIT)/pow(temp,2);

  if (vgt2 > vgta)
  {
    diff_vgt3_vsb = -diff_vt2_vsb;
    diff_vgt3_vgs = 1.0 - diff_vt2_vgs;
    diff_vgt3_vds = -diff_vt2_vds;
  }
  else
  {
    diff_vgt3_vsb = (2.0 * diff_m_vsb * PHIT * log(1.0 + g1))
                    + (temp*diff_g1_vsb/(1.0+g1));
    diff_vgt3_vgs = exp(vgt2/temp) * (1.0 - diff_vt2_vgs)/(1.0 + g1);
    diff_vgt3_vds = -exp(vgt2/temp) * diff_vt2_vds/(1.0 + g1);
  }

  AD delta2;
  delta2 = diff_vgt3_vsb + diff_vgt3_vgs + diff_vgt3_vds;

  AD vdss2;
  vdss2 = (vgt3/(1.0 + del2))*(2.0/(1.0 + sqrt(lambda9 + hyp1(1.0 - lambda9
          + (2.0 * THE3 * vgt3/(1.0 + del2)),eps8))));

  double lambda8 = 0.1;
  AD eps7; // clarify data type
  eps7=lambda8 * (vdss2/(vdssx + vdss2));

  AD vds2;
  vds2 = hyp5(vds,vdss2,eps7);

  AD fj;
  fj=((1.0 + del2) * (lambda9 + hyp1(1.0 - lambda9 + THE3 * vds2,eps8)) * vds2)
     /(2.0 * vgt3 - (1.0 + del2) * vds2);

  AD qd, qs;
  qd = -COX * (0.5 * vgt3 + delta2 * vds2 * (fj/12.0 + (pow(fj,2)/60.0) - 1.0/3.0));
  qs = -COX * (0.5 * vgt3 + delta2 * vds2 * (fj/12.0 - (pow(fj,2)/60.0) - 1.0/6.0));

  double eps6 = 0.03;
  AD vgb;
  vgb = vgs + vsb;

  double vfb;
  vfb = VT0 - PHIB - K0 * sqrt(PHIB);

  AD qbs, qbd;
  if (vgb > vfb)
  {
    qbs = -COX * K0 * (-K0/2.0 + (sqrt(pow(K0/2.0,2) + hyp3(vgb-vfb,vsb+vt1-vfb,eps6))));
    qbd = -COX * K0 * (-K0/2.0 + (sqrt(pow(K0/2.0,2) + hyp3(vgb-vfb,vds2+vsb+vt1d-vfb,eps6))));
  }
  else
  {
    qbs = -COX * hyp3(vgb - vfb,vsb + vt1 - vfb,eps6);
    qbd = -COX * hyp3(vgb - vfb,vds2 + vsb + vt1d - vfb,eps6);
  }

  AD qb, qg;
  qb = 0.5 * (qbs + qbd);
  qg = -(qd + qs + qb);
  //---------end DC current model---------
  //=========end Model Equations=========

  // Calculate dynamic current contributions due to charge
  // iqd is dynamic contribution to drain current
  // iqg is dynamic contribution to gate current
  // iqs is dynamic contribution to source current
  AD iqd, iqg, iqs;
  iqd = qd.fastAccessDx(0)*x[3] + qd.fastAccessDx(1)*x[4] + qd.fastAccessDx(2)*x[5];
  iqg = qg.fastAccessDx(0)*x[3] + qg.fastAccessDx(1)*x[4] + qg.fastAccessDx(2)*x[5];
  iqs = qs.fastAccessDx(0)*x[3] + qs.fastAccessDx(1)*x[4] + qs.fastAccessDx(2)*x[5];

  flow[0] = ids+iavl + iqd; //Drain current
  flow[1] = iqg; //Gate current
  flow[2] = -ids - iqs; //Source current

  effort[0] = x[0]+x[2]; //Vdb
  effort[1] = x[1]+x[2]; //Vgb
  effort[2] = x[2]; //Vsb
}

AD Mosn9:: hyp1( AD x, double EPI)
{
	return (0.5*(x+sqrt(x * x + 4 * EPI * EPI)));
}

AD Mosn9:: hyp1( AD x, AD EPI)
{
	return (0.5*(x+sqrt(x * x + 4 * EPI * EPI)));
}

AD Mosn9::diff_hyp1( AD x, double EPI)
{
	return (0.5*(1+ x/sqrt(x * x + 4 * EPI * EPI)));
}

AD Mosn9::hyp2( AD x, AD x0, double EPI)
{
	return(x - EPI - hyp1(x-x0,EPI));
}

AD Mosn9::diff_hyp2( AD x, AD x0, double EPI)
{
	return(1- diff_hyp1(x-x0,EPI));
}

AD Mosn9::hyp3( AD x, AD x0, double EPI)
{
	return(hyp2(x,x0,EPI) - hyp2(0,x0,EPI));
}

AD Mosn9::diff_hyp3( AD x, AD x0, double EPI)
{
	return(diff_hyp2(x,x0,EPI));
}

AD Mosn9::hyp4( AD x, double x0, double EPI)
{
	return(hyp1(x-x0, EPI) - hyp1(-x0, EPI));
}

AD Mosn9::diff_hyp4( AD x, double x0, double EPI)
{
	return(diff_hyp1(x-x0,EPI));
}

AD Mosn9::hyp5( AD x, AD x0, AD EPI)
{
	return(x0 - hyp1(x0 - x - EPI * EPI/x0 ,EPI));
}

AD Mosn9::diff_hyp5( AD x, AD x0, double EPI)
{
	return(-diff_hyp1(x0 - x - EPI * EPI/x0 ,EPI));
}
