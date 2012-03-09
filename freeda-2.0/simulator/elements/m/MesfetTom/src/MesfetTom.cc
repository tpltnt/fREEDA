#include "MesfetTom.h"

// Static members
const unsigned MesfetTom::n_par = 29;

// Element information
ItemInfo MesfetTom::einfo =
{
  "mesfettom",
  "Triquint TOM model",
  "Shuping Zhang",
  DEFAULT_ADDRESS"transistor>mesfet"
};

// Parameter information
ParmInfo MesfetTom::pinfo[] =
{
  {"vgsi", "Intrinsic gate-source voltage(V)", TR_DOUBLE, false},
  {"vgdi", "Intrinsic gate-drain voltage(V)", TR_DOUBLE, false},
  {"vdsi", "Intrinsic drain-source voltage(V)", TR_DOUBLE, false},
  {"vmax", "max voltage(V)", TR_DOUBLE, false},
  {"k", "Boltzmann's constant", TR_DOUBLE, false},
  {"beta", "Transconductance coefficient (A/V^2)", TR_DOUBLE, false},
  {"vt0", "Pinch-off voltage (V)", TR_DOUBLE, false},
  {"gama", "Slope parameter of pinch-off voltage", TR_DOUBLE, false},
  {"q", "Power law parameter", TR_DOUBLE, false},
  {"delt", "Slope of drain characteristic in the saturated region (/A V)", TR_DOUBLE, false},
  {"alfa", "Slope of drain characteristic in the linear region (/V)", TR_DOUBLE, false},
  {"t", "Channel transit-time delay (sec)", TR_DOUBLE, false},
  {"cgs0", "Gate-source schottly barrier capacitance at Vgs=0 (farad)", TR_DOUBLE, false},
  {"cgd0", "Gate-drain schottly barrier capacitance at Vgd=0 (farad)", TR_DOUBLE, false},
  {"vbi", "Built-in barrier potential (V)", TR_DOUBLE, false},
  {"is", "Saturation current of diodes (A)", TR_DOUBLE, false},
  {"n", "Ideality  factor", TR_DOUBLE, false},
  {"ib0", "Reverse breakdown saturation current (A)", TR_DOUBLE, false},
  {"nr", "Reverse breakdown ideality factor", TR_DOUBLE, false},
  {"vbd", "Reverse breakdown voltage (V)", TR_DOUBLE, false},
  {"tj", "Junction temperature (K)", TR_INT, false},
  {"t1", "Scan time (ms)", TR_DOUBLE, false},
  {"tnom", "Reference Temperature (K)", TR_DOUBLE, false},
  {"tbet", "BETA power law temperature coefficient (1/K)", TR_DOUBLE, false},
  {"xti", "Diode saturation current temperature exponent", TR_DOUBLE, false},
  {"avt0", "Pinch-off voltage(vt0)linear temperature cofficient(/K)", TR_DOUBLE, false},
  {"bvt0", "Pinch-off voltage(vt0) quadratic tempature cofficent(/K^2)", TR_DOUBLE, false},
  {"eg", "Brrier height at 0K", TR_DOUBLE, false},
  {"tm", "Ids linear temp. coeff. (1/K)", TR_DOUBLE, false},
  {"tme", "Ids power law temp. coeff. (1/K^2)", TR_DOUBLE, false},
  {"m", "Grading coefficient", TR_DOUBLE, false},
  {"area", "Area multglier", TR_DOUBLE, false},
  {"a0", "Drain saturation current for Vgs=0 (A)", TR_DOUBLE, false},
  {"a1", "Coefficient for V1 (A/V)", TR_DOUBLE, false},
  {"a2", "Coefficient for V1^2 (A/V^2)", TR_DOUBLE, false},
  {"a3", "Coefficient for V1^3 (A/V^3)", TR_DOUBLE, false},
  {"fcc", "Forward-bias depletion capacitance coefficient (V)", TR_DOUBLE, false},
  {"vds0", "Vds at which BETA was measured (V)", TR_DOUBLE, false}
};


MesfetTom::MesfetTom(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Value of cap is required
  paramvalue[0] = &(vgsi =- 0.5);
  paramvalue[1] = &(vgdi =- 3.5);
  paramvalue[2] = &(vdsi = 3.);
  paramvalue[3] = &(vmax = 0.5);
  paramvalue[4] = &(k = 1.38e-23);
  paramvalue[5] = &(beta = 0.1);
  paramvalue[6] = &(vt0 = -2.0);
  paramvalue[7] = &(gama = 1.5);
  paramvalue[8] = &(q = 2.0);
  paramvalue[9] = &(delt = 0.2);
  paramvalue[10] = &(alfa = 2.0);
  paramvalue[11] = &(t = zero);
  paramvalue[12] = &(cgs0 = zero);
  paramvalue[13] = &(cgd0 = zero);
  paramvalue[14] = &(vbi = 0.8);
  paramvalue[15] = &(is = zero);
  paramvalue[16] = &(n = one);
  paramvalue[17] = &(ib0 = zero);
  paramvalue[14] = &(vbi = 0.8);
  paramvalue[15] = &(is = zero);
  paramvalue[16] = &(n = one);
  paramvalue[17] = &(ib0 = zero);
  paramvalue[18] = &(nr = 10.0);
  paramvalue[19] = &(vbd = 1e10);
  paramvalue[20] = &(tj = 293);
  paramvalue[21] = &(t1 = 100.0);
  paramvalue[22] = &(tnom = 293);
  paramvalue[23] = &(tbet = zero);
  paramvalue[24] = &(xti = 2.0);
  paramvalue[25] = &(avt0 = zero);
  paramvalue[26] = &(bvt0 = zero);
  paramvalue[27] = &(eg = 0.8);
  paramvalue[28] = &(tm = zero);
  paramvalue[29] = &(tme = zero);
  paramvalue[30] = &(m = 0.5);
  paramvalue[31] = &(area = one);
  paramvalue[32] = &(a0 = .016542);
  paramvalue[33] = &(a1 = .0500214);
  paramvalue[34] = &(a2 = .02012);
  paramvalue[35] = &(a3 = -.00806592);
  paramvalue[36] = &(fcc = 0.5);
  paramvalue[37] = &(vds0 = 4.0);

  // Set the number of terminals
  setNumTerms(3);
  // Set number of states
  setNumberOfStates(2);
  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void MesfetTom::init() throw(string&)
{
  k2 = cgs0 / sqrt(one - fcc);
  k3 = cgd0 / sqrt(one - fcc);
  Veff = 1/2 *(vgsi + vgdi + (sqrt((vgsi - vgdi)*(vgsi-vgdi)+(1/alfa)*(1/alfa))));
  VT = vt0 - gama * vdsi;
  A1 = 1/2 * (Veff + VT + (sqrt((Veff -VT)*(Veff-VT) +( delt * delt))));
  if(A1 < vmax)
    Vnew = A1;
  else
    Vnew = vmax;

  F1 = 1/2 * (1+((Veff - VT)/(sqrt((Veff-VT)*(Veff-VT) +( delt * delt)))));
  F2 = 1/2 * (1+((vgsi - vgdi)/(sqrt((vgsi-vgdi)*(vgsi-vgdi) + (1/alfa)*(1/alfa)))));
  F3 = 1/2 * (1-((vgsi - vgdi)/(sqrt((vgsi-vgdi)*(vgsi-vgdi) + (1/alfa)*(1/alfa)))));

  delta = tj - tnom;
  tn = tj / tnom;
  Vt = k * tj / eCharge;
  k5 = n  * Vt;
  k6 = nr * Vt;
  Vt0 = vt0 * (one + (delta * avt0) + (delta * delta * bvt0));

  if (tbet)
    Beta = beta * pow(1.01, (delta * tbet));
  else
    Beta = beta;

  Ebarr = eg -.000702 * tj * tj / (tj + 1108.);
  EbarrN = eg -.000702 * tnom * tnom / (tnom + 1108.);
  Nn = eCharge / 38.696 / k / tj;

  if (xti)
    is1 *= pow(tn, xti / Nn);
    Vbi = vbi * tn -3. * Vt * log(tn) + tn * EbarrN - Ebarr;

  if((vdsi >= 0) && (vdsi <= (3/alfa)))
      KTanh = 1 - (1 - ((alfa * vdsi)/3)) * (1-((alfa * vdsi)/3)) * (1-((alfa * vdsi)/3));
  else
      KTanh = 1;

  Vp = vt0 - gama * vdsi;
  if((vgsi*(t1-t)) >= Vp)
	Ids0 = Beta * pow((vgsi * (t1-t)- Vp), q) * KTanh;
  else
        Ids0 = 0;
  k1 = fcc * Vbi;
  k4 = 2. * Vbi * (one - fcc);

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

void MesfetTom::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: vgs
  // x[1]: vgd
  // x[2]: dvgs/dt
  // x[3]: dvgd/dt
  // x[4]: vgs(t-tau)
  // effort[0]: ugs , flow[0]: ig
  // effort[1]: uds , flow[1]: id

  // Assign known output voltages
  effort[0] = x[0];
  effort[1] = x[0] - x[1];

  AD cgs, cgd, igs, igd, ids, itmp, vx;

  //get capacity value
  cgs = (AD)(cgs0 * ((F1 * F2)/(sqrt(1-Vnew/vbi))) + cgd0 * F3);
  cgd = (AD)(cgs0 * ((F1 * F3)/(sqrt(1-Vnew/vbi))) + cgd0 * F2);

  //static igs current
  igs = (AD)((is1 * (exp(vgsi/(n*Vt))-1))-(ib0 * (exp(-((vgsi + vbd)/(nr * Vt))))));
  igd = (AD)((is1 * (exp(vgdi/(n*Vt))-1))-(ib0 * (exp(-((vgdi + vbd)/(nr * Vt))))));

  ids = (AD)(Ids0 / (1+ (delt *Ids0 * vdsi)));
  // Calculate cgs, including temperature effect.
  if (k1 > x[0])
    cgs = cgs0 / sqrt(one - x[0] / Vbi);
  else
    cgs = k2 * (one + (x[0] - k1) / k4);
  cgs *= (one + m * (0.0004 * delta + one - Vbi / vbi));

  // Calculate the total current igs = static + dq_dt
  igs += cgs * x[2];

  // Calculate cgd, including temperature effect.
  if (k1 > x[1])
    cgd = cgd0 / sqrt(one - x[1] / Vbi);
  else
    cgd = k3 * (one + (x[1] - k1) / k4);
  cgd *= (one + m * (0.0004 * delta + one - Vbi / vbi));

  // Calculate the total current igd = static + dq_dt
  igd += cgd * x[3];

  // Calculate ids. Include temperature effects.
  vx = x[4] * (one + Beta * (vds0 - effort[1]));

  itmp = (a0 + vx*(a1 + vx*(a2 + vx  * a3)))* tanh(gama * effort[1]);
  if ((itmp * effort[1]) * (x[0] - Vt0) > 0.0)
    ids = itmp;
  else
    ids = zero;

  if (tme && tm)
    ids *= pow((1 + delta * tm), tme);
  // Calculate the output currents
  flow[0]  = (igd + igs) * area;
  flow[1]  = (ids - igd) * area;
}

