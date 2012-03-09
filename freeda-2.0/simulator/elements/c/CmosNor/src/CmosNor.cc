#include "CmosNor.h"

// Define the number of parameters
const unsigned CmosNor :: n_par = 12;

// Element Information
ItemInfo CmosNor::einfo =
{
  "cmosnor",
  "CMOS NOR Gate",
  "Mazen M Kharbutli",
  DEFAULT_ADDRESS"category:cmos,digital",
  "2003_05_15"
};

// Parameter Information
ParmInfo CmosNor::pinfo[] =
{
  {"vtn","NMOS Threshold Voltage (V)", TR_DOUBLE, false},
  {"vtp","PMOS Threshold Voltage (V)", TR_DOUBLE, false},
  {"un","Effective Mobility of Electrons in NMOS ((cm^2)/(V-sec))", TR_DOUBLE, false},
  {"up","Effective Mobility of Holes in PMOS ((cm^2)/(V-sec))", TR_DOUBLE, false},
  {"en","Permittivity of the Gate Insultor in NMOS (F/cm)", TR_DOUBLE, false},
  {"ep","Permittivity of the Gate Insultor in PMOS (F/cm)", TR_DOUBLE, false},
  {"tox","Thickness of the Gate Insulator (cm)", TR_DOUBLE, false},
  {"wn","Channel Width of NMOS (cm)", TR_DOUBLE, false},
  {"ln","Channel Length of NMOS (cm)", TR_DOUBLE, false},
  {"wp","Channel Width of PMOS (cm)", TR_DOUBLE, false},
  {"lp","Channel Length of PMOS (cm)", TR_DOUBLE, false},
  {"td","Response Delay Time (sec)", TR_DOUBLE, false},
};

// Constructor
CmosNor::CmosNor(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(vtn = 1);
  paramvalue[1] = &(vtp = -1);
  paramvalue[2] = &(un = 500);
  paramvalue[3] = &(up = 200);
  paramvalue[4] = &(en = 34.515e-14);
  paramvalue[5] = &(ep = 34.515e-14);
  paramvalue[6] = &(tox = 2e-6);
  paramvalue[7] = &(ln = 2.0e-6);
  paramvalue[8] = &(wn = 50e-6);
  paramvalue[9] = &(lp = 2.0e-6);
  paramvalue[10] = &(wp = 50e-6);
  paramvalue[11] = &(td = 0);

  // Set Number of Terminals
  setNumTerms(5);

  // Set Flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set Number of States
  setNumberOfStates(4);
}

// Initialization Function
void CmosNor::init() throw(string&)
{
  // Find Beta for the NMOS and PMOS
  betan = (un * en * wn) / (tox * ln);
  betap = (up * ep * wp) / (tox * lp);

  DenseIntVector var(4);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  var[3] = 3;
  DenseIntVector dvar(4);
  dvar[0] = 0;
  dvar[1] = 1;
  dvar[2] = 2;
  dvar[3] = 3;
  DenseIntVector none(4);
  DenseIntVector tvar(4);
  tvar[0] = 0;
  tvar[1] = 1;
  tvar[2] = 2;
  tvar[3] = 3;
  DenseDoubleVector delay(4);
  delay[0] = 0;
  delay[1] = td;
  delay[2] = td;
  delay[3] = 0;
  initializeAD(var, dvar, none, tvar, delay);
}

// Evaluate Function
void CmosNor::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]  : Vdd    x[4] : d(Vdd)/dt
  // x[1]  : Vin1   x[5] : d(Vin1)/dt
  // x[2]  : Vin2   x[6] : d(Vin2)/dt
  // x[3]  : Vout   x[7] : d(Vout)/dt
  // effort[0] : Vdd    flow[0] : Isdp
  // effort[1] : Vin1   flow[1] : Ig1
  // effort[2] : Vin2   flow[2] : Ig2
  // effort[3] : Vout   flow[3] : Iout

  AD Idsn, Idsp, Vout, VoutTime;
  AD Idsn1, Idsp1, Vout1;
  AD IdsnA1, IdspA1, VoutA1;
  AD IdsnB1, IdspB1, VoutB1, tempB1;
  AD IdsnC1, IdspC1, VoutC1, tempC1;
  AD IdsnD1, IdspD1, VoutD1, tempD1;
  AD IdsnE1, IdspE1, VoutE1;
  AD Idsn2, Idsp2, Vout2;
  AD IdsnA2, IdspA2, VoutA2;
  AD IdsnB2, IdspB2, VoutB2, tempB2;
  AD IdsnC2, IdspC2, VoutC2, tempC2;
  AD IdsnD2, IdspD2, VoutD2, tempD2;
  AD IdsnE2, IdspE2, VoutE2;

   // Region A1: 0 <= Vin1 <= Vtn
  Vout1 = x[0];
  Idsn1 = 0;
  Idsp1 = 0;

  // Region B1: Vtn <= Vin1 <= Vdd/2
  IdsnB1 = (betan/2) * (x[13] - vtn) * (x[13] - vtn);
  IdspB1 = -betap * (((x[13] - x[0] - vtp) * (x[3] - x[0]))
                        - ((x[3] - x[0]) * (x[3] - x[0]) / 2));
  tempB1 = (x[13] - vtp)*(x[13] - vtp);
  tempB1 = tempB1 - (2 * x[0] * (x[13] - (x[0]/2) - vtp));
  tempB1 = tempB1 - ((betan/betap) * (x[13] - vtn) * (x[13] - vtn));
  VoutB1 = (x[13] - vtp) + sqrt(tempB1);

  if (x[13] > vtn)
  {
    Vout1 = VoutB1;
    Idsn1 = IdsnB1;
    Idsp1 = IdspB1;
  }

  // Region C1: Vin1 = Vdd/2
  IdsnC1 = betan * (x[13] - vtn) * (x[13] - vtn) / 2;
  IdspC1 = -IdsnC1;
  VoutC1 = x[3];

  // Region D1: Vdd/2 <= Vin1 <= Vdd + Vtp
  IdsnD1 = betan * (((x[13] - vtn) * x[3]) - (x[3] * x[3] / 2));
  IdspD1 = - (betap/2) * (x[13] - x[0] - vtp) * (x[13] - x[0] - vtp);
  tempD1 = (x[13] - vtn) * (x[13] - vtn);
  tempD1 = tempD1 - ((betap/betan) * (x[13] - x[0] - vtp) * (x[13] - x[0] - vtp));
  VoutD1 = (x[13] - vtn) - sqrt(tempD1);

  if (x[13] > (x[0]/2))
  {
    Vout1 = VoutD1;
    Idsn1 = IdsnD1;
    Idsp1 = IdspD1;
  }

  // Region E1: Vin1 >= Vdd + Vtp
  if (x[13]-x[0]-vtp > 0.0)
  {
    Vout1 = 0.0;
    Idsn1 = 0.0;
    Idsp1 = 0.0;
  }

  // Region A2: 0 <= Vin2 <= Vtn
  Vout2 = x[0];
  Idsn2 = 0;
  Idsp2 = 0;

  // Region B2: Vtn <= Vin2 <= Vdd/2
  IdsnB2 = (betan/2) * (x[14] - vtn) * (x[14] - vtn);
  IdspB2 = -betap * (((x[14] - x[0] - vtp) * (x[3] - x[0]))
                        - ((x[3] - x[0]) * (x[3] - x[0]) / 2));
  tempB2 = (x[14] - vtp)*(x[14] - vtp);
  tempB2 = tempB2 - (2 * x[0] * (x[14] - (x[0]/2) - vtp));
  tempB2 = tempB2 - ((betan/betap) * (x[14] - vtn) * (x[14] - vtn));
  VoutB2 = (x[14] - vtp) + sqrt(tempB2);

  if (x[14] > vtn)
  {
    Vout2 = VoutB2;
    Idsn2 = IdsnB2;
    Idsp2 = IdspB2;
  }

  // Region C2: Vin2 = Vdd/2
  IdsnC2 = betan * (x[14] - vtn) * (x[14] - vtn) / 2;
  IdspC2 = -IdsnC2;
  VoutC2 = x[3];

  // Region D2: Vdd/2 <= Vin2 <= Vdd + Vtp
  IdsnD2 = betan * (((x[14] - vtn) * x[3]) - (x[3] * x[3] / 2));
  IdspD2 = - (betap/2) * (x[14] - x[0] - vtp) * (x[14] - x[0] - vtp);
  tempD2 = (x[14] - vtn) * (x[14] - vtn);
  tempD2 = tempD2 - ((betap/betan) * (x[14] - x[0] - vtp) * (x[14] - x[0] - vtp));
  VoutD2 = (x[14] - vtn) - sqrt(tempD2);

  if (x[14] > (x[0]/2))
  {
    Vout2 = VoutD2;
    Idsn2 = IdsnD2;
    Idsp2 = IdspD2;
  }

  // Region E2: Vin2 >= Vdd + Vtp
  if (x[14]-x[0]-vtp > 0.0)
  {
    Vout2 = 0.0;
    Idsn2 = 0.0;
    Idsp2 = 0.0;
  }

  // Assign output voltage
  effort[0] = x[0];
  effort[1] = x[1];
  effort[2] = x[2];
  if (x[14] > x[13])
    Vout = Vout2;
  else
    Vout = Vout1;
  effort[3] = Vout;

  // Assign output current
  if (x[14] > x[13])
  {
    Idsp = Idsp2;
    Idsn = Idsn2;
  }
  else
  {
    Idsp = Idsp1;
    Idsn = Idsn1;
  }
  flow[0] = -Idsp;
  flow[1] = 0;
  flow[2] = 0;
  flow[3] = Idsn1 + Idsn2 + Idsp;
}

