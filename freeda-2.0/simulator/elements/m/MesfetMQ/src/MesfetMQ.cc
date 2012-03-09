#include "MesfetMQ.h"

// Static members
const unsigned MesfetMQ::n_par = 22;

// Element information
ItemInfo MesfetMQ::einfo =
{
  "mesfetmq",
  "Intrinsic MESFET using the Materka-Kacprzac model",
  "Senthil Velu",
  DEFAULT_ADDRESS"transistor>mesfet",
  "2002_12_10"
};

// Parameter information
ParmInfo MesfetMQ::pinfo[] =
{
  {"idss", "Drain saturation current for Vgs=0 (A)", TR_DOUBLE, false},
  {"vp0", "Pinch-off voltage for Vds=0 (V)", TR_DOUBLE, false},
  {"gama", "Voltage slope parameter of pinch-off voltage (1/V)", TR_DOUBLE, false},
  {"e", "Constant part of power law parameter", TR_DOUBLE, false},
  {"ke", "Dependence of power law on Vgs (1/V)", TR_DOUBLE, false},
  {"sl", "Slope of the Vgs=0 drain characteristic in the saturated region (S)",	TR_DOUBLE, false},
  {"kg", "Drain dependence on Vgs in the linear region (1/V)", TR_DOUBLE,	false},
  {"ss", "Slope of the drain characteristic in the saturated region (S)",	TR_DOUBLE, false},
  {"t", "Channel transit-time delay (s)", TR_DOUBLE, false},
  {"ig0", "Saturation current of gate-source Schottky barrier (A)",	TR_DOUBLE, false},
  {"afag", "Slope factor of gate conduction current (1/V)", TR_DOUBLE, false},
  {"ib0", "Current parameter of gate-drain breakdown source (A)",	TR_DOUBLE, false},
  {"afab", "Slope factor of breakdown current (1/V)", TR_DOUBLE, false},
  {"vbc", "Breakdown voltage (V)", TR_DOUBLE, false},
  {"r10", "Intrinsic channel resistance for Vgs = 0 (Ohm)", TR_DOUBLE, false},
  {"kr", "Slope factor of intrinsic channel resistance (1/V)", TR_DOUBLE,	false},
  {"c10", "Gate-source Schottky barrier capacitance for Vgs = 0 (F)",	TR_DOUBLE, false},
  {"k1", "Slope parameter of gate-source capacitance (1/V)", TR_DOUBLE, false},
  {"c1s", "Constant parasitic component of gate-source capacitance (F)", TR_DOUBLE, false},
  {"cf0", "Gate-drain feedback capacitance for Vgd = 0 (F)", TR_DOUBLE, false},
  {"kf", "Slope parameter of gate-drain feedback capacitance (1/V)", TR_DOUBLE, false},
  {"area", "Area multiplier", TR_DOUBLE, false},
};


MesfetMQ::MesfetMQ(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(idss = .1);
  paramvalue[1] = &(vp0 = -2.);
  paramvalue[2] = &(gama = zero);
  paramvalue[3] = &(e = 2.);
  paramvalue[4] = &(ke = zero);
  paramvalue[5] = &(sl = .15);
  paramvalue[6] = &(kg = zero);
  paramvalue[7] = &(ss = zero);
  paramvalue[8] = &(t = zero);
  paramvalue[9] = &(ig0 = zero);
  paramvalue[10] = &(afag = 38.696);
  paramvalue[11] = &(ib0 = zero);
  paramvalue[12] = &(afab = zero);
  paramvalue[13] = &(vbc = 1e10);
  paramvalue[14] = &(r10 = zero);
  paramvalue[15] = &(kr = zero);
  paramvalue[16] = &(c10 = zero);
  paramvalue[17] = &(k1 = 1.25);
  paramvalue[18] = &(c1s = zero);
  paramvalue[19] = &(cf0 = zero);
  paramvalue[20] = &(kf = 1.25);
  paramvalue[21] = &(area = one);

  // Set the number of terminals
  setNumTerms(3);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(2);
}


void MesfetMQ::init() throw(string&)
{
  // calculate some constants.
  k01 = r10 * kr;
  k02 = c10 * sqrt(5.) + c1s;
  k03 = cf0 * sqrt(5.);

  DenseIntVector var(2);
  var[0] = 0;
  var[1] = 1;
  DenseIntVector dvar(2);
  dvar[0] = 0;
  dvar[1] = 1;
  DenseIntVector novar;
  DenseIntVector tvar(1);
  DenseDoubleVector delay(1);
  delay[0] = t;
  initializeAD(var, var, novar, tvar, delay);
}

void MesfetMQ::eval(AD * x, AD * effort, AD * flow)
{
  // x[0]: qgs
  // x[1]: qgd
  // x[2]: dqgs/dt
  // x[3]: dqgd/dt
  // x[4]: qgs(t-tau)
  // effort[0]: ugs , flow[0]: ig
  // effort[1]: uds , flow[1]: id

  AD tmp1,vgs,vgd,vgst,tmp2;
  AD ids, igd, igs;
  AD cgs, cgd, Ri, base;

	//calculate temp_vgs & temp_vgst.
	//temp_vgs = c1s * (k1 * x[0] / c10 - 2. - c1s);
	//temp_vgst = c1s * (k1 * x[4] / c10 - 2. - c1s);

	// calculate Vgsi
  if (x[0] > (1.1*c10/k1))
    vgs = (one - pow(one - (k1*x[0]/2./c10),2.)) / k1;
  else
    vgs = x[0] / c10 / sqrt(5.);

	// calculate Vgst
  if (x[0] > (1.1*c10/k1))
    vgst = (one - pow(one - (k1*x[4]/2./c10),2.)) / k1;
  else
    vgst = x[4] / c10 /sqrt(5.);

	// static igs current
	igs = ig0 * (exp(afag * vgs) - one);
	igs += x[2];

  // Calculate Ri
  if (kr * vgs < one)
    Ri = r10 - k01 * vgs;
  else
    Ri = zero;

  // Calculate Vgdi
  if (x[1] > (1.1*cf0/kf))
    vgd = (one - pow(one - kf * x[1] / 2. / cf0,2.)) / kf;
  else
    vgd = x[1] / cf0 / sqrt(5.);

	// calculate uds = - x[1] + ugs
	effort[0] = vgs + igs * Ri;
	effort[1] = effort[0] - vgd;

	// calculate ids
	tmp1 = one - vgst / (vp0 + gama * effort[1]);
	tmp2 = tanh(sl*effort[1]/(idss*(one - kg * vgst)));
  if (tmp1 > 0.0)
    ids = (idss + ss * effort[1]) * pow(tmp1,e + ke * vgst) *	tmp2;
  else
    ids = zero;

  // static igd current
  igd = - ib0 * (exp(afab * ( -vgd - vbc)) - one);
  igd += x[3];

  // Now calculate the output currents.
  flow[0] = area * (igd + igs);
  flow[1] = area * (ids - igd);
}
