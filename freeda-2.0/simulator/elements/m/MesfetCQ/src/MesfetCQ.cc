#include "MesfetCQ.h"

// Static members
const unsigned MesfetCQ::n_par = 29;

// Element information
ItemInfo MesfetCQ::einfo =
{
  "mesfetcq",
  "Intrinsic MESFET using Curtice-Ettemberg cubic model",
  "Senthil Velu",
  DEFAULT_ADDRESS"transistor>mesfet",
  "2002_12_10"
};

// Parameter information
ParmInfo MesfetCQ::pinfo[] =
{
  {"a0", "Drain saturation current for Vgs=0 (A)", TR_DOUBLE, false},
  {"a1", "Coefficient for V1 (A/V)", TR_DOUBLE, false},
  {"a2", "Coefficient for V1^2 (A/V^2)", TR_DOUBLE, false},
  {"a3", "Coefficient for V1^3 (A/V^3)", TR_DOUBLE, false},
  {"beta", "V1 dependance on Vds (1/V)", TR_DOUBLE, false},
  {"vds0", "Vds at which BETA was measured (V)", TR_DOUBLE, false},
  {"gama", "Slope of drain characteristic in the linear region (1/V)", TR_DOUBLE, false},
  {"vt0", "Voltage at which the channel current is forced to be zero for Vgs<=Vto (V)", TR_DOUBLE, false},
  {"cgs0", "Gate-source Schottky barrier capacitance for Vgs=0 (F)", TR_DOUBLE, false},
  {"cgd0", "Gate-drain Schottky barrier capacitance for Vgd=0 (F)", TR_DOUBLE, false},
  {"is", "Diode saturation current (A)", TR_DOUBLE, false},
  {"n", "Diode ideality factor", TR_DOUBLE, false},
  {"ib0", "Breakdown current parameter (A)", TR_DOUBLE, false},
  {"nr", "Breakdown ideality factor", TR_DOUBLE, false},
  {"t", "Channel transit time (s)", TR_DOUBLE, false},
  {"vbi", "Built-in potential of the Schottky junctions (V)", TR_DOUBLE, false},
  {"fcc", "Forward-bias depletion capacitance coefficient (V)", TR_DOUBLE, false},
  {"vbd", "Breakdown voltage (V)", TR_DOUBLE, false},
  {"tnom", "Reference Temperature (K)", TR_DOUBLE, false},
  {"avt0", "Pinch-off voltage (VP0 or VT0) linear temp. coefficient (1/K)",	TR_DOUBLE, false},
  {"bvt0", "Pinch-off voltage (VP0 or VT0) quadratic temp. coefficient (1/K^2)", TR_DOUBLE, false},
  {"tbet", "BETA power law temperature coefficient (1/K)", TR_DOUBLE, false},
  {"tm", "Ids linear temp. coeff. (1/K)", TR_DOUBLE, false},
  {"tme", "Ids power law temp. coeff. (1/K^2)", TR_DOUBLE, false},
  {"eg", "Barrier height at 0.K (eV)", TR_DOUBLE, false},
  {"m", "Grading coefficient", TR_DOUBLE, false},
  {"xti", "Diode saturation current temperature exponent", TR_DOUBLE, false},
  {"tj", "Junction Temperature (K)", TR_DOUBLE, false},
  {"area", "Area multiplier", TR_DOUBLE, false}
};


MesfetCQ::MesfetCQ(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(a0 = .1);
  paramvalue[1] = &(a1 = .05);
  paramvalue[2] = &(a2 = zero);
  paramvalue[3] = &(a3 = zero);
  paramvalue[4] = &(beta = zero);
  paramvalue[5] = &(vds0 = 4.);
  paramvalue[6] = &(gama = 1.5);
  paramvalue[7] = &(vt0 = -1e10);
  paramvalue[8] = &(cgs0 = zero);
  paramvalue[9] = &(cgd0 = zero);
  paramvalue[10] = &(is = zero);
  paramvalue[11] = &(n = one);
  paramvalue[12] = &(ib0 = zero);
  paramvalue[13] = &(nr = 10.);
  paramvalue[14] = &(t = zero);
  paramvalue[15] = &(vbi = .8);
  paramvalue[16] = &(fcc = .5);
  paramvalue[17] = &(vbd = 1e10);
  paramvalue[18] = &(tnom = 293.);
  paramvalue[19] = &(avt0 = zero);
  paramvalue[20] = &(bvt0 = zero);
  paramvalue[21] = &(tbet = zero);
  paramvalue[22] = &(tm = zero);
  paramvalue[23] = &(tme = zero);
  paramvalue[24] = &(eg = .8);
  paramvalue[25] = &(m = .5);
  paramvalue[26] = &(xti = 2.);
  paramvalue[27] = &(tj = 293.);
  paramvalue[28] = &(area = one);

  // Set the number of terminals
  setNumTerms(3);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(2);
}


void MesfetCQ::init() throw(string&)
{
  cgs0 = cgs0 * .1e13;
  cgd0 = cgd0 * .1e13;

  k2 = cgs0 / sqrt(one - fcc);
  k3 = cgd0 / sqrt(one - fcc);

  delta_T = tj - tnom;
  tn = tj / tnom;
  Vt = kBoltzman * tj / eCharge;
  k5 = n  * Vt;
  k6 = nr * Vt;
  Vt0 = vt0 * (one + (delta_T * avt0) + (delta_T * delta_T * bvt0));
  if (tbet)
    Beta = beta * pow(1.01, (delta_T * tbet));
  else
    Beta = beta;
  Ebarr = eg -.000702 * tj * tj / (tj + 1108.);
  EbarrN = eg -.000702 * tnom * tnom / (tnom + 1108.);
  Nn = eCharge / 38.696 / kBoltzman / tj;
  Is = is * exp((tn - one) * Ebarr / Nn / Vt);
  if (xti)
    Is *= pow(tn, xti / Nn);
  Vbi = vbi * tn -3. * Vt * log(tn) + tn * EbarrN - Ebarr;
  k1 = fcc * Vbi;
  k4 = 2. * Vbi * (one - fcc);
  M1 = 2.*(k4 - k1);
  M2 = 2.*k4/k2;

  DenseIntVector var(2);
  var[0] = 0;
  var[1] = 1;
  DenseIntVector dvar(2);
  dvar[0] = 0;
  dvar[1] = 1;
  DenseIntVector d2var;
  DenseIntVector tvar(1);
  DenseDoubleVector delay(1);
  delay[0] = t;
  initializeAD(var, dvar, d2var, tvar, delay);
}

void MesfetCQ::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: qgs
  // x[1]: qgd
  // x[2]: dqgs/dt
  // x[3]: dqgd/dt
  // x[4]: qgs(t-tau)
  // effort[0]: ugs , flow[0]: ig
  // effort[1]: uds , flow[1]: id

  AD vgs,vgst,vgd;

  //calculate vgs
  if ((k1*k1+k1*M1)/M2 > x[0])
    vgs = vbi*(one - pow(one - x[0] / 2. / cgs0 / vbi,2.));
  else
    vgs = (-M1 + sqrt(M1*M1 + 4. * M2 * x[0])) / 2.;

  //calculate vgst
  if ((k1*k1+k1*M1)/M2 > x[0])
    vgst = vbi*(one - pow(one - x[4] / 2. / cgs0 / vbi,2.));
  else
    vgst = (-M1 + sqrt(M1*M1 + 4. * M2 * x[4])) / 2.;

  //calculate vgd
  if ((k1*k1+k1*M1)/M2 > x[1])
    vgd = vbi*(one - pow(one - x[1] / 2. / cgd0 / vbi,2.));
  else
    vgd = (-M1 + sqrt(M1*M1 + 4. * M2 * x[1])) / 2.;

  // Assign known output voltages
  effort[0] = vgs;
  effort[1] = vgs - vgd;

  AD ids, igd, igs, itmp, vx, cgs, cgd;

  // static igs current
  igs = Is *(exp(vgs / k5) - one) - ib0 * exp(-(vgs +vbd) / k6);

  //cgs *= (one + m * (0.0004 * delta_T + one - Vbi / vbi));

	// Calculate the total current igs = static + dq_dt
  igs += x[3];

  // static igd current
  igd = Is * (exp(vgd / k5) - one) - ib0 * exp(-(vgd + vbd) / k6);

  //cgd *= (one + m * (0.0004 * delta_T + one - Vbi / vbi));

  // Calculate the total current igd = static + dq_dt
  igd += x[4];

  // Calculate ids. Include temperature effects.
  vx = vgst * (one + Beta * (vds0 - effort[1]));

  itmp = (a0 + vx*(a1 + vx*(a2 + vx  * a3)))* tanh(gama * effort[1]);
  if ((itmp * effort[1]) * (vgs - Vt0) > 0.0)
    ids = itmp;
  else
    ids = zero;

  if (tme && tm)
    ids *= pow((1 + delta_T * tm), tme);

  // Calculate the output currents
  flow[0]  = (igd + igs) * area / .1e13;
  flow[1]  = (ids - igd) * area / .1e13;
}

