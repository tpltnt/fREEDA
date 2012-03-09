#include "DiodeQk.h"

// Static members
const unsigned DiodeQk::n_par = 4;

// Element information
ItemInfo DiodeQk::einfo =
{
	"diodeqk",
	"Microwave Diode",
	"James Murray",
	"category:diode",
	"2008_05_26"
};

// Parameter information
ParmInfo DiodeQk::pinfo[] =
{
	{"js", "Saturation current (A)", TR_DOUBLE, false},
	{"alfa", "Slope factor of conduction current (1/V)", TR_DOUBLE, false},
	{"r0", "Series resistance (Ohms)", TR_DOUBLE, false},
	{"area", "Area multiplier", TR_DOUBLE, false},
};


DiodeQk::DiodeQk(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
	// Set default parameter values
	paramvalue[0] = &(js = 1e-12);
	paramvalue[1] = &(alfa = 38.696);
	paramvalue[2] = &(r0 = 2);
	paramvalue[3] = &(area = one);

	// Set the number of terminals
	setNumTerms(2);

	// Set flags
	setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

	// Set number of states
	setNumberOfStates(1);
}


void DiodeQk::init() throw(string&)
{
	Vo = log(5e8 / alfa) / alfa; // normal is .612
	DenseIntVector var(1,0);
	initializeAD(var, var);
}

void DiodeQk::eval(AD * x, AD * effort, AD * flow)
{
	// x[0]: state variable
	// x[1]: time derivative of x[0]
	AD Vj, itmp;

	Vj=x[0];

	if(Vj>Vo)
		flow[0]=js * (exp(alfa*(Vj-Vo)) - 1);
	else
		flow[0]=0;
	effort[0] = Vj + flow[0] * r0;

	flow[0] *= area;
}
