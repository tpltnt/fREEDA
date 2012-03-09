// Thermal CMOS NOR Gate
// By Travis Lentz and Tony Mulder

#include "CmosNorT.h"

// Define the number of parameters
const unsigned CmosNorT :: n_par = 23;

// Element Information
ItemInfo CmosNorT::einfo = {
  "cmosnort",
  "CMOS NOR Gate, Electro-thermal",
  "Travis Lentz, Tony Mulder",
  DEFAULT_ADDRESS"category:cmos,electrothermal,digital",
  "2005_04_19"
};

// Parameter Information
ParmInfo CmosNorT::pinfo[] = {
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
  {"thermal", "Thermal element flag", TR_BOOLEAN, false},
  {"ts", "Kirchhoff transformation temperature (K)", TR_DOUBLE, false},
  {"tnom", "Ambient temperature (K)", TR_DOUBLE, false},
  {"zt", "Thermal impedence summation (K/W)", TR_DOUBLE, false},
  {"c1", "P channel GS capacitance (F)", TR_DOUBLE, false},
  {"c2", "N channel GS capacitance (F)", TR_DOUBLE, false},
  {"c3", "output GS (Miller) capacitance (F)", TR_DOUBLE, false},
  {"c4", "parasitic diode capacitance from out to Vdd (F)", TR_DOUBLE, false},
  {"c5", "parasitic diode capacitance from out to Gnd(F)", TR_DOUBLE, false},
  {"freq", "operating frequency (Hz)", TR_DOUBLE, false},
  {"lk", "leakage current (A)", TR_DOUBLE, false}
};

// Constructor
CmosNorT::CmosNorT(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
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
  paramvalue[10] = &(wp = 100e-6);
  paramvalue[11] = &(td = 0);
  paramvalue[12] = &(thermal = false);
  paramvalue[13] = &(ts = 300.);
  paramvalue[14] = &(tnom = 300);
  paramvalue[15] = &(zt = 310);
  paramvalue[16] = &(c1 = 1.0e-12);
  paramvalue[17] = &(c2 = 1.0e-12);
  paramvalue[18] = &(c3 = 1.0e-12);
  paramvalue[19] = &(c4 = 1.0e-12);
  paramvalue[20] = &(c5 = 1.0e-12);
  paramvalue[21] = &(freq = 1.0e6);
  paramvalue[22] = &(lk = 8e-6);

  // Set Number of Terminals
  setNumTerms(7);

  // Set Flags
  setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);

  // Set Number of States
  setNumberOfStates(5);
}

// Initialization Function
void CmosNorT::init() throw(string&)
{
  // Find Beta for the NMOS and PMOS
  betan = (un * en * wn) / (tox * ln);
  betap = (up * ep * wp) / (tox * lp);
  power = ( freq * (c1+c2+c3+c4+c5));

  DenseIntVector var(5);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  var[3] = 3;
  var[4] = 4;
  DenseIntVector dvar(5);
  dvar[0] = 0;
  dvar[1] = 1;
  dvar[2] = 2;
  dvar[3] = 3;
  dvar[4] = 4;
  DenseIntVector none(5);
  DenseIntVector tvar(5);
  tvar[0] = 0;
  tvar[1] = 1;
  tvar[2] = 2;
  tvar[3] = 3;
  tvar[4] = 4;
  DenseDoubleVector delay(5);
  delay[0] = 0;
  delay[1] = td;
  delay[2] = td;
  delay[3] = 0;
  delay[4] = 0;
  initializeAD(var, dvar, none, tvar, delay);
}

void CmosNorT::getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1));
  term_list.push_back(getTerminal(2));
  term_list.push_back(getTerminal(3));
  term_list.push_back(getTerminal(4)); // Local reference terminal
  term_list.push_back(getTerminal(5));
  term_list.push_back(getTerminal(6)); // Local reference terminal

  local_ref_vec.push_back(4); // Local reference index
  local_ref_vec.push_back(6); // Local reference index
}

// Evaluate Function
void CmosNorT::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]  : Vdd    x[5] : d(Vdd)/dt
  // x[1]  : Vin1   x[6] : d(Vin1)/dt
  // x[2]  : Vin2   x[7] : d(Vin2)/dt
  // x[3]  : Vout   x[8] : d(Vout)/dt
  // x[4]  : temp
  // effort[0] : Vdd    flow[0] : Isdp
  // effort[1] : Vin1   flow[1] : Ig1
  // effort[2] : Vin2   flow[2] : Ig2
  // effort[3] : Vout   flow[3] : Iout
  // effort[4] : temp   flow[4] : pout

  AD Idsn, Idsp, Vout, VoutTime, tj, pwr;
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
  Idsn1 = lk;
  Idsp1 = -lk;

  // Region B1: Vtn <= Vin1 <= Vdd/2
  IdsnB1 = ((betan/2) * (x[16] - vtn) * (x[16] - vtn))+(power*x[0])+lk;
  IdspB1 = (-betap * (((x[16] - x[0] - vtp) * (x[3] - x[0]))
    - ((x[3] - x[0]) * (x[3] - x[0]) / 2)))-(power*x[0])-lk;
  tempB1 = (x[16] - vtp)*(x[16] - vtp);
  tempB1 = tempB1 - (2 * x[0] * (x[16] - (x[0]/2) - vtp));
  tempB1 = tempB1 - ((betan/betap) * (x[16] - vtn) * (x[16] - vtn));
  VoutB1 = (x[16] - vtp) + sqrt(tempB1);

  if (x[16] > vtn)
  {
    Vout1 = VoutB1;
    Idsn1 = IdsnB1;
    Idsp1 = IdspB1;
  }

  // Region C1: Vin1 = Vdd/2
  IdsnC1 = (betan * (x[16] - vtn) * (x[16] - vtn) / 2)+(power*x[0])+lk;
  IdspC1 = -IdsnC1;
  VoutC1 = x[3];

  // Region D1: Vdd/2 <= Vin1 <= Vdd + Vtp
  IdsnD1 = (betan * (((x[16] - vtn) * x[3]) - (x[3] * x[3] / 2))) +
   (power*x[0])+lk;
  IdspD1 = (- (betap/2) * (x[16] - x[0] - vtp) * (x[16] - x[0] - vtp))
       - (power*x[0])-lk;
  tempD1 = (x[16] - vtn) * (x[16] - vtn);
  tempD1 = tempD1 - ((betap/betan) * (x[16] - x[0] - vtp) * (x[16] - x[0] - vtp));
  VoutD1 = (x[16] - vtn) - sqrt(tempD1);

  if (x[16] > (x[0]/2))
  {
    Vout1 = VoutD1;
    Idsn1 = IdsnD1;
    Idsp1 = IdspD1;
  }

  // Region E1: Vin1 >= Vdd + Vtp
  if (x[16]-x[0]-vtp > 0.0)
  {
    Vout1 = 0.0;
    Idsn1 = lk;
    Idsp1 = -lk;
  }

  // Region A2: 0 <= Vin2 <= Vtn
  Vout2 = x[0];
  Idsn2 = lk;
  Idsp2 = -lk;

  // Region B2: Vtn <= Vin2 <= Vdd/2
  IdsnB2 = ((betan/2) * (x[17] - vtn) * (x[17] - vtn))+(power*x[0])+lk;
  IdspB2 = (-betap * (((x[17] - x[0] - vtp) * (x[3] - x[0]))
     - ((x[3] - x[0]) * (x[3] - x[0]) / 2)))-(power*x[0])-lk;
  tempB2 = (x[17] - vtp)*(x[17] - vtp);
  tempB2 = tempB2 - (2 * x[0] * (x[17] - (x[0]/2) - vtp));
  tempB2 = tempB2 - ((betan/betap) * (x[17] - vtn) * (x[17] - vtn));
  VoutB2 = (x[17] - vtp) + sqrt(tempB2);

  if (x[17] > vtn)
  {
    Vout2 = VoutB2;
    Idsn2 = IdsnB2;
    Idsp2 = IdspB2;
  }

  // Region C2: Vin2 = Vdd/2
  IdsnC2 = (betan * (x[17] - vtn) * (x[17] - vtn) / 2)+(power*x[0])+lk;
  IdspC2 = -IdsnC2;
  VoutC2 = x[3];

  // Region D2: Vdd/2 <= Vin2 <= Vdd + Vtp
  IdsnD2 = (betan * (((x[17] - vtn) * x[3]) - (x[3] * x[3] / 2))) +
   (power*x[0])+lk ;
  IdspD2 = (- (betap/2) * (x[17] - x[0] - vtp) * (x[17] - x[0] - vtp))
    - (power*x[0])-lk ;
  tempD2 = (x[17] - vtn) * (x[17] - vtn);
  tempD2 = tempD2 - ((betap/betan) * (x[17] - x[0] - vtp) * (x[17] - x[0] - vtp));
  VoutD2 = (x[17] - vtn) - sqrt(tempD2);

  if (x[17] > (x[0]/2))
  {
    Vout2 = VoutD2;
    Idsn2 = IdsnD2;
    Idsp2 = IdspD2;
  }

  // Region E2: Vin2 >= Vdd + Vtp
  if (x[17]-x[0]-vtp > 0.0)
  {
    Vout2 = 0.0;
    Idsn2 = lk;
    Idsp2 = -lk;
  }

  // Assign output voltage
  effort[0] = x[0];
  effort[1] = x[1];
  effort[2] = x[2];
  if (x[17] > x[16])
    Vout = Vout2;
  else
    Vout = Vout1;
  effort[3] = Vout;
  effort[4] = tnom + x[4]*zt;

  // Assign output current
  if (x[17] > x[16])
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
  flow[4] = ((effort[0]-effort[3])*Idsp)-(effort[3]*(Idsn1+Idsn2))-(power*x[0])-lk;
}

