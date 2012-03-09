#include "MesfetCT.h"

// Static members
const unsigned MesfetCT::n_par = 31;

// Element information
ItemInfo MesfetCT::einfo =
{
  "mesfetct",
  "Intrinsic MESFET using Curtice-Ettemberg cubic model with thermal effects",
  "Carlos E. Christoffersen, Hector Gutierrez",
  DEFAULT_ADDRESS"transistor>mesfet,electrothermal",
  "2000_07_20"
};

// Parameter information
ParmInfo MesfetCT::pinfo[] =
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
  {"avt0", "Pinch-off voltage (VP0 or VT0) linear temp. coefficient (1/K)", TR_DOUBLE, false},
  {"bvt0", "Pinch-off voltage (VP0 or VT0) quadratic temp. coefficient (1/K^2)", TR_DOUBLE, false},
  {"tbet", "BETA power law temperature coefficient (1/K)", TR_DOUBLE, false},
  {"tm", "Ids linear temp. coeff. (1/K)", TR_DOUBLE, false},
  {"tme", "Ids power law temp. coeff. (1/K^2)", TR_DOUBLE, false},
  {"eg", "Barrier height at 0.K (eV)", TR_DOUBLE, false},
  {"m", "Grading coefficient", TR_DOUBLE, false},
  {"xti", "Diode saturation current temperature exponent", TR_DOUBLE, false},
  {"kirchhoff", "Use Kirchhoff transformation flag", TR_BOOLEAN, false},
  {"b", "Thermal conductivity temperature exponent", TR_DOUBLE, false},
  {"ts", "Kirchhoff transformation temperature (K)", TR_DOUBLE, false},
  {"area", "Area multiplier", TR_DOUBLE, false}
};


MesfetCT::MesfetCT(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
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
  paramvalue[27] = &(kirchhoff = false);
  paramvalue[28] = &(b = 1.22);
  paramvalue[29] = &(ts = 300.0);
  paramvalue[30] = &(area = one);

  // Set the number of terminals
  setNumTerms(5);

  // Set flags
  setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(3);
}


void MesfetCT::init() throw(string&)
{
  k2 = cgs0 / sqrt(one - fcc);
  k3 = cgd0 / sqrt(one - fcc);

  bm1 = b - one;

  DenseIntVector var(3);
  var[0] = 0;
  var[1] = 1;
  var[2] = 2;
  DenseIntVector dvar(2);
  dvar[0] = 0;
  dvar[1] = 1;
  DenseIntVector d2var;
  DenseIntVector tvar(1);
  DenseDoubleVector delay(1);
  delay[0] = t;
  initializeAD(var, dvar, d2var, tvar, delay);
}

void MesfetCT::getLocalRefIdx(UnsignedVector& local_ref_vec,
TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1));
  term_list.push_back(getTerminal(2)); // Local reference terminal
  term_list.push_back(getTerminal(3));
  term_list.push_back(getTerminal(4)); // Local reference terminal

  local_ref_vec.push_back(2); // Local reference index
  local_ref_vec.push_back(4); // Local reference index
}


void MesfetCT::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: vgs
  // x[1]: vgd
  // x[2]: temperature - tnom
  // x[3]: dvgs/dt
  // x[4]: dvgd/dt
  // x[5]: vgs(t-tau)
  // effort[0]: ugs , flow[0]: ig
  // effort[1]: uds , flow[1]: id
  // effort[2]: temperature (theta if kirchhoff is true), flow[2]: pout

  // Assign known output voltages
  effort[0] = x[0];
  effort[1] = x[0] - x[1];
  // Shift input temperature by 293K
  effort[2] = x[2] + tnom;

  AD delta_T, tn, Vt, k1, k4, k5, k6, Vt0, Beta, Ebarr, EbarrN, Nn;
  AD ids, igd, igs, itmp, vx, cgs, cgd, Is, Vbi;

  delta_T = effort[2] - tnom;
  tn = effort[2] / tnom;
  Vt = kBoltzman * effort[2] / eCharge;
  k5 = n  * Vt;
  k6 = nr * Vt;
  Vt0 = vt0 * (one + (delta_T * avt0) + (delta_T * delta_T * bvt0));
  if (tbet)
    Beta = beta * pow(1.01, (delta_T * tbet));
  else
    Beta = beta;
  Ebarr = eg -.000702 * effort[2] * effort[2] / (effort[2] + 1108.);
  EbarrN = eg -.000702 * tnom * tnom / (tnom + 1108.);
  Nn = eCharge / 38.696 / kBoltzman / effort[2];
  Is = is * exp((tn - one) * Ebarr / Nn / Vt);
  if (xti)
    Is *= pow(tn, xti / Nn);
  Vbi = vbi * tn -3. * Vt * log(tn) + tn * EbarrN - Ebarr;
  k1 = fcc * Vbi;
  k4 = 2. * Vbi * (one - fcc);

  // static igs current
  igs = Is *(exp(x[0] / k5) - one) - ib0 * exp(-(x[0] +vbd) / k6);

  // Power dissipated on this junction
  flow[2] = -igs * x[0];

  // Calculate cgs, including temperature effect.
  if (k1 > x[0])
    cgs = cgs0 / sqrt(one - x[0] / Vbi);
  else
    cgs = k2 * (one + (x[0] - k1) / k4);
  cgs *= (one + m * (0.0004 * delta_T + one - Vbi / vbi));

  // Calculate the total current igs = static + dq_dt
  igs += cgs * x[3];

  // static igd current
  igd = Is * (exp(x[1] / k5) - one) - ib0 * exp(-(x[1] + vbd) / k6);

  // Power dissipated on this junction
  flow[2] -= igd * x[1];

  // Calculate cgd, including temperature effect.
  if (k1 > x[1])
    cgd = cgd0 / sqrt(one - x[1] / Vbi);
  else
    cgd = k3 * (one + (x[1] - k1) / k4);
  cgd *= (one + m * (0.0004 * delta_T + one - Vbi / vbi));

	// Calculate the total current igd = static + dq_dt
  igd += cgd * x[4];

  // Calculate ids. Include temperature effects.
  vx = x[5] * (one + Beta * (vds0 - effort[1]));

  itmp = (a0 + vx*(a1 + vx*(a2 + vx  * a3)))* tanh(gama * effort[1]);
  if ((itmp * effort[1]) * (x[0] - Vt0) > 0.0)
    ids = itmp;
  else
    ids = zero;

  if (tme && tm)
    ids *= pow((1 + delta_T * tm), tme);

  // Apply Kirchhoff transformation (if requested)
  if (kirchhoff)
    effort[2] = ts/(bm1) * (b - pow(ts/effort[2], bm1));

  // Calculate the output currents
  flow[0]  = (igd + igs) * area;
  flow[1]  = (ids - igd) * area;
  flow[2] -= effort[1] * ids;
  flow[2] *= area;
}
