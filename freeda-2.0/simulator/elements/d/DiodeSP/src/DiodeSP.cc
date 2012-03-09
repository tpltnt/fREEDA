#include "DiodeSP.h"

// Static members
const unsigned DiodeSP::n_par = 12;

// Element information
ItemInfo DiodeSP::einfo =
{
  "diodesp",
  "Spice diode model (conserves charge)",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:diode",
  "2000_07_20"
};

// Parameter information
ParmInfo DiodeSP::pinfo[] =
{
  {"is", "Saturation current (A)", TR_DOUBLE, false},
  {"n", "Emission coefficient", TR_DOUBLE, false},
  {"ibv", "Current magnitude at the reverse breakdown voltage (A)",
	TR_DOUBLE, false},
  {"bv", "Breakdown voltage (V)", TR_DOUBLE, false},
  {"fc", "Coefficient for forward-bias depletion capacitance",
	TR_DOUBLE, false},
  {"cj0", "Zero-bias depletion capacitance (F)", TR_DOUBLE, false},
  {"vj", "Built-in junction potential (V)", TR_DOUBLE, false},
  {"m", "PN junction grading coefficient", TR_DOUBLE, false},
  {"tt", "Transit time (s)", TR_DOUBLE, false},
  {"area", "Area multiplier", TR_DOUBLE, false},
  {"charge", "Use charge-conserving model", TR_BOOLEAN, false},
  {"rs", "Series resistance (ohms)", TR_DOUBLE, false}
};


DiodeSP::DiodeSP(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(is = 1e-14);
  paramvalue[1] = &(n = one);
  paramvalue[2] = &(ibv = 1e-10);
  paramvalue[3] = &(bv = 40.0);
  paramvalue[4] = &(fc = .5);
  paramvalue[5] = &(cj0 = zero);
  paramvalue[6] = &(vj = one);
  paramvalue[7] = &(m = .5);
  paramvalue[8] = &(tt = zero);
  paramvalue[9] = &(area = one);
  paramvalue[10] = &(charge = true);
  paramvalue[11] = &(rs = zero);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void DiodeSP::init() throw(string&)
{
  if (charge)
	{
    // Add one terminal
    Circuit* cir = getCircuit();
    unsigned tref_id = getTerminal(1)->getId();
    unsigned term_id1 = cir->addTerminal(getInstanceName()
		+ ":extra");
    cir->connect(getId(), term_id1);
    // Connect an external 1K resistor (not really needed if using the
    // augmentation network)
    unsigned newelem_id =
		cir->addElement("resistor", getInstanceName() + ":res");
    cir->connect(newelem_id, tref_id);
    cir->connect(newelem_id, term_id1);
    Element* elem = cir->getElement(newelem_id);
    double res = 1e3;
    elem->setParam("r", &res, TR_DOUBLE);
    elem->init();

    // Set the number of terminals
    setNumTerms(3);
    // Set number of states
    setNumberOfStates(2);

    DenseIntVector var(2);
    var[0] = 0;
    var[1] = 1;
    DenseIntVector dvar(1);
    dvar[0] = 1;
    initializeAD(var, dvar);
  }
  else
	{
    // Set the number of terminals
    setNumTerms(2);
    // Set number of states
    setNumberOfStates(1);

    DenseIntVector var(1);
    initializeAD(var, var);
  }
}

void DiodeSP::eval(AD * x, AD * effort, AD * flow)
{
  double alfa = eCharge / n / kBoltzman / 300.; // tnom = 300K
  double v1 = log(5e8 / alfa) / alfa; // normal is .5e9
  double k3 = exp(alfa * v1);

  // x[0]: x as in Rizolli's equations
  AD vd,id,ibreakdown;
  // Diode voltage
  if (v1 > x[0])
    vd = x[0] + zero;
  else
    vd = v1 + log(one + alfa*(x[0]-v1))/alfa;

  // Static current
  if (v1 > x[0])
    id = is * (exp(alfa * x[0]) - one);
  else
    id = is * k3 * (one + alfa * (x[0] - v1)) - is;

  // breakdown current (taken from Antognetti's SPICE book)
  if (vd == -bv)
    ibreakdown = -ibv;
  else if (vd < -bv)
    ibreakdown = -is * exp(-alfa*(vd + bv) - 1.0 + bv/alfa);
  else
    ibreakdown = 0.0;

  // include the breakdown current contribution
  flow[0] = id + ibreakdown;

  if (charge)
	{
    // x[1]: q
    // x[2]: dq/dt
    // Form the additional error function
    AD qvj;
    double km;
    if (isSet(&cj0))
		{
      if (fc * vj > vd)
        qvj = vj * cj0 / (one - m) * (one - pow(one - vd / vj, one - m));
      else
        qvj = cj0 * pow(one - fc, - m - one) *
              (((one - fc * (one + m)) * vd + .5 * m * vd * vd / vj) -
              ((one - fc * (one + m)) * vj * fc + .5 * m * vj * fc * fc)) +
              vj * cj0 / (one - m) * (one - pow(one - fc, one - m));
      km = vj * cj0 / (one - m) * 1e1;
    }
    else
		{
      qvj = zero;
      km = 1e-12;
    }
    if (isSet(&tt))
		{
      qvj += tt * flow[0];
      km += tt * 1e-2;
    }
    // Add capacitor current
    flow[0] += x[2] * km;
    flow[1] = - flow[0];
    effort[1] = qvj / km - x[1];
    effort[0] = vd + flow[0] * rs + effort[1];

    // scale the current according to area.
    flow[0] *= area;
    flow[1] *= area;
  }
  else
	{
    AD dvd_dx;
    if (v1 > x[0])
      dvd_dx = one;
    else
      dvd_dx = one / (one + alfa*(x[0]-v1));
    AD cd;
    if (isSet(&cj0))
    {
      if (fc * vj > vd)
        cd = cj0 * pow(one - vd / vj, -m);
      else
        cd = cj0 * pow(one - fc, - m - one) * ((one - fc * (one + m)) + m * vd / vj);
    }
    else
      cd = zero;

    if (isSet(&tt))
      cd += alfa * tt * flow[0];
    // Add capacitor current
    flow[0] += cd * dvd_dx * x[1];
    effort[0] = vd + flow[0] * rs;
    // scale the current according to area.
    flow[0] *= area;
  }
}
