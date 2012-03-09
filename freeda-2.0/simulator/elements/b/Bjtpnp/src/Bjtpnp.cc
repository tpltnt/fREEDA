#include "Bjtpnp.h"

const unsigned Bjtpnp :: n_par = 50;

ItemInfo Bjtpnp::einfo =
{
	"bjtpnp",
	"Gummel Poon model, PNP",
	"Senthil Velu",
	DEFAULT_ADDRESS"category:transistor>bjt",
	"2002_05_15"
};

// Parameter Information
ParmInfo Bjtpnp :: pinfo[] =
{
	{"is","Transport saturation current (A)", TR_DOUBLE, false},
	{"bf","Ideal maximum forward beta", TR_DOUBLE, false},
	{"nf","Forward current emission coefficient", TR_DOUBLE, false},
	{"vaf","Forward early voltage (V)", TR_DOUBLE, false},
	{"ikf","Forward-beta high current roll-off knee current (A)", TR_DOUBLE, false},
	{"ise","Base-emitter leakage saturation current (A)", TR_DOUBLE, false},
	{"ne","Base-emitter leakage emission coefficient", TR_DOUBLE, false},
	{"br","Ideal maximum reverse beta", TR_DOUBLE, false},
	{"nr","Reverse current emission coefficient", TR_DOUBLE, false},
	{"var","revers early voltage (V)", TR_DOUBLE, false},
	{"ikr","Corner for reverse-beta high current roll off (A)", TR_DOUBLE, false},
	{"isc","Base collector leakage saturation current (A)", TR_DOUBLE, false},
	{"nc","Base-collector leakage emission coefficient", TR_DOUBLE, false},
	{"re","Emitter ohmic resistance (W)", TR_DOUBLE, false},
	{"rb","Zero bias base resistance (W)", TR_DOUBLE, false},
	{"rbm","Minimum base resistance (W)", TR_DOUBLE, false},
	{"irb","Current at which rb falls to half of rbm (A)", TR_DOUBLE, false},
	{"rc","Collector ohmic resistance (W)", TR_DOUBLE, false},
	{"eg","Badgap voltage (eV)", TR_DOUBLE, false},
	{"cje","Base emitter zero bias p-n capacitance (F)", TR_DOUBLE, false},
	{"vje","Base emitter built in potential (V)", TR_DOUBLE, false},
	{"mje","Base emitter p-n grading factor", TR_DOUBLE, false},
	{"cjc","Base collector zero bias p-n capacitance (F)", TR_DOUBLE, false},
	{"vjc","Base collector built in potential (V)", TR_DOUBLE, false},
	{"mjc","Base collector p-n grading factor", TR_DOUBLE, false},
	{"xcjc","Fraction of cbc connected internal to rb", TR_DOUBLE, false},
	{"fc","Forward bias depletion capacitor coefficient", TR_DOUBLE, false},
	{"tf","Ideal forward transit time (S)", TR_DOUBLE, false},
	{"xtf","Transit time bias dependence coefficient", TR_DOUBLE, false},
	{"vtf","Transit time dependency on vbc (V)", TR_DOUBLE, false},
	{"itf","Transit time dependency on ic (A)", TR_DOUBLE, false},
	{"tr","Ideal reverse transit time (S)", TR_DOUBLE, false},
	{"xtb","Forward and reverse beta temperature coefficient", TR_DOUBLE, false},
	{"xti","IS temperature effect exponent", TR_DOUBLE, false},
	{"tre1","RE temperature coefficient (linear)", TR_DOUBLE, false},
	{"tre2","RE temperature coefficient (quadratic)", TR_DOUBLE, false},
	{"trb1","RB temperature coefficient (linear)", TR_DOUBLE, false},
	{"trb2","RB temperature coefficient  (quadratic)", TR_DOUBLE, false},
	{"trm1","RBM temperature coefficient (linear)", TR_DOUBLE, false},
	{"trm2","RBM temperature coefficient (quadratic)", TR_DOUBLE, false},
	{"trc1","RC temperature coefficient (linear)", TR_DOUBLE, false},
	{"trc2","RC temperature coefficient (quadratic)", TR_DOUBLE, false},
	{"tnom","Nominal temperature (K)", TR_DOUBLE, false},
	{"t","temperature (K)", TR_DOUBLE, false},
	{"cjs", "Collector substrate capacitance", TR_DOUBLE, false},
	{"mjs", "substrate junction exponential factor", TR_DOUBLE, false},
	{"vjs", "substrate junction built in potential", TR_DOUBLE, false},
	{"area","Current multiplier", TR_DOUBLE, false},
	{"ns","substrate p-n coefficient",TR_DOUBLE,false},
	{"iss","Substrate saturation current",TR_DOUBLE,false}
};

Bjtpnp::Bjtpnp(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
	// Set default parameter values
	paramvalue[0] = &(is = 1e-16);
	paramvalue[1] = &(bf = 100.);
	paramvalue[2] = &(nf = one);
	paramvalue[3] = &(vaf = zero);
	paramvalue[4] = &(ikf = zero);
	paramvalue[5] = &(ise = zero);
	paramvalue[6] = &(ne = 1.5);
	paramvalue[7] = &(br = one);
	paramvalue[8] = &(nr = one);
	paramvalue[9] = &(var = zero);
	paramvalue[10] = &(ikr = zero);
	paramvalue[11] = &(isc = zero);
	paramvalue[12] = &(nc = 2.);
	paramvalue[13] = &(re = zero);
	paramvalue[14] = &(rb = zero);
	paramvalue[15] = &(rbm = zero);
	paramvalue[16] = &(irb = zero);
	paramvalue[17] = &(rc = zero);
	paramvalue[18] = &(eg = 1.11);
	paramvalue[19] = &(cje = zero);
	paramvalue[20] = &(vje = 0.75);
	paramvalue[21] = &(mje = 0.33);
	paramvalue[22] = &(cjc = zero);
	paramvalue[23] = &(vjc = 0.75);
	paramvalue[24] = &(mjc = 0.33);
	paramvalue[25] = &(xcjc = one);
	paramvalue[26] = &(fc = 0.5);
	paramvalue[27] = &(tf = zero);
	paramvalue[28] = &(xtf = zero);
	paramvalue[29] = &(vtf = zero);
	paramvalue[30] = &(itf = zero);
	paramvalue[31] = &(tr = zero);
	paramvalue[32] = &(xtb = zero);
	paramvalue[33] = &(xti = 3.);
	paramvalue[34] = &(tre1 = zero);
	paramvalue[35] = &(tre2 = zero);
	paramvalue[36] = &(trb1 = zero);
	paramvalue[37] = &(trb2 = zero);
	paramvalue[38] = &(trm1 = zero);
	paramvalue[39] = &(trm2 = zero);
	paramvalue[40] = &(trc1 = zero);
	paramvalue[41] = &(trc2 = zero);
	paramvalue[42] = &(tnom = 300.);
	paramvalue[43] = &(t = 300.);
	paramvalue[44] = &(cjs = zero);
	paramvalue[45] = &(mjs = zero);
	paramvalue[46] = &(vjs = 0.75);
	paramvalue[47] = &(area = one);
	paramvalue[48] = &(ns = one);
	paramvalue[49] = &(iss = zero);

	setNumTerms(4);
	setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
	setNumberOfStates(3);
}

void Bjtpnp::init() throw(string&)
{
	DenseIntVector var2(3);
	var2[0] = 0;
	var2[1] = 1;
	var2[2] = 2;
	DenseIntVector novar;
	DenseDoubleVector nodelay;
	initializeAD(var2, var2, var2, novar, nodelay);
}

//evaluate function
void Bjtpnp :: eval(AD * x, AD * effort, AD * flow)
{
	//x[vbe,vbc,vcjs,dvbe,dvbc,dvcjs,d2vbe,d2vbc,d2vcjs]
	//x[0] : Veb
	//x[1] : Vcb
	//x[2] : Vcjs
	//x[3] : dveb/dt
	//x[4] : dvcb/dt
	//x[5] : dvcjs/dt
	//x[6] : d2veb/dt
	//x[7] : d2vcb/dt
	//x[8] : d2vcjs/dt
	//effort[0] : Vcjs flow[0] : Ic
	//effort[1] : Vbjs flow[1] : Ib
	//effort[2] : Vejs flow[2] : Ie

	AD Ibe,Ibc,Ice,Ibf,Ile,Ibr,Ilc,kqb,kqbtemp;
	AD cbej,cbet,cbcj,cbct,ibe,ibc;
	AD vth = (kBoltzman * tnom) / eCharge;

	if (!isSet(&rbm))
    rbm = rb;
	//dc current equations

	Ibf = is * (exp(x[0]/nf/vth) - one);
	Ile = ise * (exp(x[0]/ne/vth) - one);
	Ibr = is * (exp(x[1]/nr/vth) - one);
	Ilc = isc * (exp(x[1]/nc/vth) - one);
	Ibe = Ibf / bf + Ile;
	Ibc = Ibr / br + Ilc;
	AD Ibf1 = (Ibe - Ile) * bf;
	AD Ibr1 = (Ibc - Ilc) * br;
	AD kqbtem = zero;
	if (ikf)
		kqbtem =  4.*(Ibf1/ikf);
	if (ikr)
		kqbtem += 4.*(Ibr1/ikr);
	kqbtemp = sqrt(one + kqbtem);
	AD tempvaf = zero;
	if (vaf)
		tempvaf = x[1] / vaf;
	AD tempvar = zero;
	if (var)
		tempvar = x[0] / var;
	kqb = 0.5 * (one / (one - tempvaf - tempvar)) * (one + kqbtemp);
	Ice = (Ibf1 - Ibr1) / kqb;

	// Charge base-collector
  if (x[1] > fc*vjc)
    cbcj = cjc*pow(one-fc,-one-mjc)*(one-fc*(one+mjc)+mjc*x[1]/vjc);
  else
    cbcj = cjc*pow(one-x[1]/vjc,-mjc);
	cbct=tr*is*exp(x[1]/nr/vth)/(nr*vth);
	ibc=area*(cbct+xcjc*cbcj)*x[4];
	//end

	//Current Base-Emitter
  if (x[0] > fc*vje)
    cbej = cjc*pow(one-fc,-one-mjc)*(one-fc*(one+mje)+mje*x[0]/vje);
  else
    cbej = cje*pow(one-x[0]/vje,-mje);

	AD tZ = one;
	if (vtf)
		tZ = exp(.69*x[1]/vtf);

	cbet=is*exp(x[0]/nf/vth)*tf*(one+xtf*tZ)/(nf*vth);

	ibe=area*(cbet+cbej)*x[3];

	//end

	AD Tibe,Tibc;
	Tibe = Ibe + ibe;
	Tibc = Ibc + ibc;
	AD rb1;
	AD vbx,cbx,ibx;

	rb1 = (rbm + (rb - rbm) / kqb) / area;
	vbx = x[1] + (Tibe + Tibc) * rb;

  if (vbx > fc*vjc)
    cbx = (one-xcjc)*cjc*pow(one-fc,-one-mjc)*(one-fc*(one+mjc)+mjc*vbx/vjc);
  else
    cbx = (one-xcjc)*cjc*pow(one-vbx/vjc,-mjc);

	AD dibc,dibe,dIbc,dIbe;
	dibc = area*(cbct+xcjc*cbcj) * x[7];
	dibe = area*(cbet+cbej) * x[6];
	dIbe = is/bf * exp(x[0]/nf/vth) * (x[3]/nf/vth)
	+ ise * exp(x[0]/ne/vth) * (x[3]/ne/vth);
	dIbc = is/br * exp(x[1]/nr/vth) * (x[4]/nr/vth)
	+ isc * exp(x[1]/nc/vth) * (x[4]/nc/vth);

	ibx = cbx * (x[3] + rb*(dibc+dibe+dIbc+dIbe));

	AD Ijs;
	AD  cjjs, ijs;
	Ijs = area * iss * (exp(x[2]/ns/vth) - one);

  if (x[2] > zero)
    cjjs = cjs*(one+mjs*x[2]/vjs);
  else
    cjjs = cjs*pow(one-x[2]/vjs,-mjs);
	ijs = Ijs + area*(cjjs * x[5]);

	flow[0] = Tibc + ibx + ijs - Ice;  // collector current
	flow[1] = -ibx - Tibc - Tibe;     // base current
	flow[2] = (Ice + Tibe - ijs);    // emitter current

	effort[0] = -x[2] + rc*flow[0];                // Collector Substrate Voltage
	effort[1] = -x[2] - x[1] + flow[1]*rb;        // Base Substrate Voltage
	effort[2] = x[0] - x[1] - x[2] + flow[2]*re; // Emitter Substrate Voltage
}

