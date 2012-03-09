#include "Jfetp.h"

//Static members
const unsigned Jfetp::n_par = 20;

//Information about this element
ItemInfo Jfetp::einfo=
{
  "jfetp",
  "The n-channel JFET based on the SPICE model",
  "Nikhil Kriplani",
  DEFAULT_ADDRESS"category:transistor>mesfet",
  "2001_04_15"
};

ParmInfo Jfetp::pinfo[] =
{
  {"af", "flicker noise exponent", TR_DOUBLE, false},
  {"area", "The Area", TR_DOUBLE, false},
  {"beta", "Transconductance parameter", TR_DOUBLE, false},
  {"cgs", "zero bias gate-source junction capacitance", TR_DOUBLE, false},
  {"cgd", "zero bias gate-drain junction capacitance", TR_DOUBLE, false},
  {"eg", "Barrier height at 0.K (eV)", TR_DOUBLE, false},
  {"fc", "coefficient for forward bias depletion capacitance", TR_DOUBLE, false},
  {"is", "gate junction saturation current", TR_DOUBLE, false},
  {"kf", "flicker noise coefficient", TR_DOUBLE, false},
  {"lambda", "Channel length modulation parameter", TR_DOUBLE, false},
  {"pb", "gate junction potential", TR_DOUBLE, false},
  {"rd", "drain ohmic resistance", TR_DOUBLE, false},
  {"rs", "source ohmic resistance", TR_DOUBLE, false},
  {"vt0", "Threshold voltage", TR_DOUBLE, false},
  {"betatc", "Temperature coefficient of transconductance parameter beta",
	TR_DOUBLE, false},
  {"m", "Gate p-n Grading coefficient", TR_DOUBLE, false},
  {"vt0tc", "Temperature corfficient for vt0", TR_DOUBLE, false},
  {"tnom", "temperature measurement parameter", TR_DOUBLE, false},
  {"b", "doping tail parameter", TR_DOUBLE, false},
  {"T", "circuit temperature", TR_DOUBLE, false}
};


Jfetp::Jfetp(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  paramvalue[0] = & (af = one);
  paramvalue[1] = & (area = one);
  paramvalue[2] = & (beta = 1.0e-4);
  paramvalue[3] = & (cgs0 = 0.0);
  paramvalue[4] = & (cgd0 = 0.0);
  paramvalue[5] = & (eg = 0.8);
  paramvalue[6] = & (fc = 0.5);
  paramvalue[7] = & (is = 1.0e-14);
  paramvalue[8] = & (kf = 0.0);
  paramvalue[9] = & (lambda = 0.0);
  paramvalue[10]= & (pb = 1.0);
  paramvalue[11]= & (rd = 0.0);
  paramvalue[12]= & (rs = 0.0);
  paramvalue[13]= & (vt0 = -2.0);
  paramvalue[14]= & (betatc = 0.0);
  paramvalue[15]= & (m = 0.5);
  paramvalue[16]= & (vt0tc = 0.0);
  paramvalue[17]= & (tnom = 298.0);
  paramvalue[18]= & (b = 1.0);
  paramvalue[19]= & (T = 300);

  setNumTerms(3);
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
  setNumberOfStates(2);
}

void Jfetp::init() throw(string&)
{
  DenseIntVector var(2);
  var[1] = 1;
  DenseIntVector novar;
  DenseDoubleVector nodelay;
  initializeAD(var,var);
}

void Jfetp::eval(AD * x, AD * effort, AD * flow)
{
  // x[0] : vds
  // x[1] : vgs
  // x[2] : dvds/dt
  // x[3] : dvgs/dt
  // effort[0]: vds, flow[0] : id
  // effort[1]: vgs, flow[1] : ig

	//Assign the known output voltages
  effort[0] = x[0];
  effort[1] = x[1];

  double vt = (kBoltzman * T)/eCharge;

  AD cgs, cgd, igd, igs, idrain;
  AD ggs, csat, gmin, evgs, cg;

  csat = is * exp((T/tnom -1) * 1.11/vt);
  gmin = 1.0e-12;

  if (-x[1] <= -5 * vt)
	{
		ggs = csat/x[1] + gmin;
		cg = ggs * (-x[1]);
	}
  else
	{
		evgs = exp((-x[1]) / vt);
		ggs = csat * evgs/vt + gmin;
		cg = csat * (evgs - 1) + gmin * (-x[1]);
	}

  AD ggd, evgd;

  if (-(x[1]-x[0]) <= -5 * vt)
	{
		ggd = -csat/(-(x[1]-x[0])) + gmin;
		cgd = ggd * (-(x[1]-x[0]));
	}
  else
	{
		evgd = exp(-(x[1] - x[0]) / vt);
		ggd = csat * evgd/vt + gmin;
		cgd = csat * (evgd - 1) + gmin * (-(x[1]-x[0]));
	}

  cg = cg + cgd;

	//Forward operation
  AD vgst, vgdt, bfac, cdrain, gm, gds, betap, apart, cpart;

  bfac = 1.0 - b;

  if (-x[0] >= 0.0)
	{
		vgst = -x[1] - vt0;
		if (vgst <= 0)
		{
			cdrain = 0;
			gm = 0.0;
			gds = 0.0;
		}
		else
		{
			betap = beta * (1 + lambda * -x[0]);
			if (vgst >= -x[0])
	    {
	      apart = 2 * b + 3 * bfac*(vgst + x[0]);
	      cpart = -x[0] * (-x[0]*(bfac*(-x[0]) - b) + vgst*apart);
	      cdrain = betap * cpart;
	      gm = betap * (-x[0]) * (apart + 3*bfac*vgst);
	      gds = betap * (vgst + x[0])*apart + beta*lambda*cpart;
	    } else {
	      bfac = vgst * bfac;
	      gm = betap * vgst * (2 * b + 3 * bfac);
	      cpart = vgst * vgst * (b + bfac);
	      cdrain = betap * cpart;
	      gds = lambda * beta * cpart;
	    }
		}
	}
  else //reverse operation
	{
    vgdt = -(x[1]-x[0]) - vt0;
    if (vgdt <= 0)
		{
			cdrain = 0;
			gm = 0.0;
			gds = 0.0;
		}
    else
		{
			betap = beta * (1 + lambda * (-x[0]));
			if (vgdt-x[0] >= 0)
			{
				apart = 2 * b + 3 * bfac*(vgdt - x[0]);
				cpart = -x[0] * (x[0]*(bfac*x[0] - b) + vgdt*apart);
				cdrain = betap*cpart;
				gm = betap * (-x[0]) * (apart + 3*bfac*vgdt);
				gds = betap * (vgst - x[0])*apart - beta*lambda*cpart - gm;
			} else {
				bfac = vgdt * bfac;
				gm = -betap * vgdt * (2 * b + 3 * bfac);
				cpart = vgdt * vgdt * (b + bfac);
				cdrain = -betap * cpart;
				gds = lambda * beta * cpart - gm;
			}
		}
  }

  idrain = cdrain - cgd;

  if (-x[1] < fc * pb)
    cgs = cgs0 / sqrt(1.0 + x[1]/pb);
  else
    cgs = cgs0 / pow(1.0 - fc,1.5) * (1.0 - 1.5*fc - x[1]/(2.0 * pb));

  if (-(x[1]-x[0]) < fc * pb)
    cgd = cgd0 / sqrt(1.0 + (x[1]-x[0])/pb);
  else
    cgd = cgd0 / pow(1.0 - fc,1.5) * (1.0 - 1.5*fc - (x[1]-x[0])/(2.0 * pb));

  igs = area * is * (exp(-x[1]/vt - one));
  igd = area * is * (exp(-(x[1] - x[0])/vt - one));

  flow[0] = -idrain - cgd * x[2] + igd;
  flow[1] = -area * cgs * x[3] - cgd * x[3] + igs;
}

