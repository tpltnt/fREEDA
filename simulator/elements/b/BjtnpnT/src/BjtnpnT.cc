#include "BjtnpnT.h"

const unsigned BjtnpnT :: n_par = 52;

ItemInfo BjtnpnT::einfo =
{
	"bjtnpnt",
	"Gummel Poon model, NPN, Electro-thermal",
	"Senthil Velu",
	DEFAULT_ADDRESS"category:transistor>bjt,electrothermal",
  "2002_05_15"
};

// Parameter Information
ParmInfo BjtnpnT :: pinfo[] =
{
	{"is","Transport saturation current (A)", TR_DOUBLE, false},
	{"bf","Ideal maximum forward beta", TR_DOUBLE, false},
	{"nf","Forward current emission coefficient", TR_DOUBLE, false},
	{"vaf","Forward early voltage (V)", TR_DOUBLE, false},
	{"ikf","Forward-beta high current roll-off knee current (A)", TR_DOUBLE,
	false},
	{"ise","Base-emitter leakage saturation current (A)", TR_DOUBLE, false},
	{"ne","Base-emitter leakage emission coefficient", TR_DOUBLE, false},
	{"br","Ideal maximum reverse beta", TR_DOUBLE, false},
	{"nr","Reverse current emission coefficient", TR_DOUBLE, false},
	{"var","revers early voltage (V)", TR_DOUBLE, false},
	{"ikr","Corner for reverse-beta high current roll off (A)", TR_DOUBLE,
	false},
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
	{"kirchhoff", "Use Kirchhoff transformation flag", TR_BOOLEAN, false},
	{"b", "Thermal conductivity temperature exponent", TR_DOUBLE, false},
	{"ts", "Kirchhoff transformation temperature (K)", TR_DOUBLE, false},
	{"cjs", "Collector substrate capacitance", TR_DOUBLE, false},
	{"mjs", "substrate junction exponential factor", TR_DOUBLE, false},
	{"vjs", "substrate junction built in potential", TR_DOUBLE, false},
	{"area","Current multiplier", TR_DOUBLE, false},
	{"ns","substrate p-n coefficient",TR_DOUBLE,false},
	{"iss","Substrate saturation current",TR_DOUBLE,false}
};

BjtnpnT::BjtnpnT(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
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
	paramvalue[43] = &(kirchhoff = false);
	paramvalue[44] = &(b = 1.22);
	paramvalue[45] = &(ts = 300.);
	paramvalue[46] = &(cjs = zero);
	paramvalue[47] = &(mjs = zero);
	paramvalue[48] = &(vjs = 0.75);
	paramvalue[49] = &(area = one);
	paramvalue[50] = &(ns = one);
	paramvalue[51] = &(iss = zero);

	setNumTerms(6);
	setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);
	setNumberOfStates(4);
}

void BjtnpnT::init() throw(string&)
{
	bm1 = b - one;
	DenseIntVector var2(4);
	var2[0] = 0;
	var2[1] = 1;
	var2[2] = 2;
	var2[3] = 3;
	DenseIntVector novar;
	DenseDoubleVector nodelay;
	initializeAD(var2, var2, var2, novar, nodelay);
}

void BjtnpnT::getLocalRefIdx(UnsignedVector& local_ref_vec,
TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1));
  term_list.push_back(getTerminal(2));
  term_list.push_back(getTerminal(3)); // Local reference terminal
  term_list.push_back(getTerminal(4));
  term_list.push_back(getTerminal(5)); // Local reference terminal

  local_ref_vec.push_back(3); // Local reference index
  local_ref_vec.push_back(5); // Local reference index
}

//evaluate function
void BjtnpnT :: eval(AD * x, AD * effort, AD * flow)
{
	//x[vbe,vbc,vcjs,dvbe,dvbc,dvcjs,d2vbe,d2vbc,d2vcjs]
	//x[0] : Vbe
	//x[1] : Vbc
	//x[2] : Vcjs
	//x[3] : temperature - tnom
	//x[4] : dvbe/dt
	//x[5] : dvbc/dt
	//x[6] : dvcjs/dt
	//x[7] : d2vbe/dt
	//x[8] : d2vbc/dt
	//x[9] : d2vcjs/dt
	//effort[0] : Vcjs flow[0] : Ic
	//effort[1] : Vbjs flow[1] : Ib
	//effort[2] : Vejs flow[2] : Ie
	//effort[3] :temperature(theta if kirchhoff is true), flow[2] : pout

	AD Ibe,Ibc,Ice,Ibf,Ile,Ibr,Ilc,kqb,kqbtemp;
	AD cbej,cbet,cbcj,cbct,ibe,ibc;
	AD egt,ist,iset,isct,isst,vjet,vjct,vjst,cjct,cjet,cjst,bft,brt;
	AD rbt,rbmt,rct,ret;
	egt = eg;
	ist = is;
	iset = ise;
	isct = isc;
	isst = iss;
	vjet = vje;
	vjct = vjc;
	vjst = vjs;
	cjct = cjc;
	cjet = cje;
	cjst = cjs;
	rbt = rb;
	rbmt = rbm;
	rct = rc;
	ret = re;
	bft = bf;
	brt = br;

	//**********Shift input temperature
	effort[3] = x[3] + tnom;
	AD delta_T = effort[3] / tnom;
	AD vth = (kBoltzman * effort[3]) / eCharge;
	egt = egt - (.000702*effort[3]*effort[3] / (1108*effort[3]));
	ist = ist*(eg*delta_T - eg) / vth + pow(delta_T,xti/nf);
	iset = iset*(eg*delta_T - eg) / vth + pow(delta_T,xti/ne);
	isct = isct*(eg*delta_T - eg) / vth + pow(delta_T,xti/nc);
	isst = isst*(eg*delta_T - eg) / vth + pow(delta_T,xti/ns);
	vjet = vje*(effort[3]-tnom) - 3.*vth*log(delta_T)*eg*tnom*delta_T - egt;
	vjct = vjc*(effort[3]-tnom) - 3.*vth*log(delta_T)*eg*tnom*delta_T - egt;
	vjst = vjs*(effort[3]-tnom) - 3.*vth*log(delta_T)*eg*tnom*delta_T - egt;
	cjct = cjct*(one+mjc*(.0004*(effort[3]-tnom)+(one-vjct/vjc)));
	cjet = cjet*(one+mjc*(.0004*(effort[3]-tnom)+(one-vjet/vje)));
	cjst = cjst*(one+mjc*(.0004*(effort[3]-tnom)+(one-vjst/vjs)));
	bft = pow(bf,xtb);
	brt = pow(br,xtb);
	rbt = rbt*(one+trb1*(effort[3]-tnom)+trb2*pow(effort[3]-tnom,2));
	rbmt = rbm*(one+trm1*(effort[3]-tnom)+trm2*pow(effort[3]-tnom,2));
	rct = rc*(one+trc1*(effort[3]-tnom)+trc2*pow(effort[3]-tnom,2));
	ret = re*(one+tre1*(effort[3]-tnom)+tre2*pow(effort[3]-tnom,2));

	if (!isSet(&rbm))
    rbmt = rbt;
	//dc current equations

	Ibf = ist * (exp(x[0]/nf/vth) - one);
	Ile = iset * (exp(x[0]/ne/vth) - one);
	Ibr = ist * (exp(x[1]/nr/vth) - one);
	Ilc = isct * (exp(x[1]/nc/vth) - one);
	Ibe = Ibf / bft + Ile;
	Ibc = Ibr / brt + Ilc;
	AD Ibf1 = (Ibe - Ile) * bft;
	AD Ibr1 = (Ibc - Ilc) * brt;
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
  if (x[1] > fc*vjct)
    cbcj = cjct*pow(one-fc,-one-mjc)*(one-fc*(one+mjc)+mjc*x[1]/vjct);
  else
    cbcj = cjct*pow(one-x[1]/vjct,-mjc);

	cbct=tr*is*exp(x[1]/nr/vth)/(nr*vth);

	ibc=area*(cbct+xcjc*cbcj)*x[5];
	//end

	//Current Base-Emitter
  if (x[0] > fc*vjet)
    cbej = cjct*pow(one-fc,-one-mjc)*(one-fc*(one+mje)+mje*x[0]/vjet);
  else
    cbej = cjet*pow(one-x[0]/vjet,-mje);

	AD tZ = one;
	if (vtf)
		tZ = exp(.69*x[1]/vtf);

	cbet=is*exp(x[0]/nf/vth)*tf*(one+xtf*tZ)/(nf*vth);

	ibe=area*(cbet+cbej)*x[4];

	//end

	AD Tibe,Tibc;
	Tibe = Ibe + ibe;
	Tibc = Ibc + ibc;
	AD rb1;
	AD vbx,cbx,ibx;

	rb1 = (rbmt + (rbt - rbmt) / kqb) / area;
	vbx = x[1] + (Tibe + Tibc) * rbt;

  if (vbx > fc*vjc)
    cbx = (one-xcjc)*cjc*pow(one-fc,-one-mjc)*(one-fc*(one+mjc)+mjc*vbx/vjct);
  else
    cbx = (one-xcjc)*cjct*pow(one-vbx/vjct,-mjc);

	AD dibc,dibe,dIbc,dIbe;
	dibc = area*(cbct+xcjc*cbcj) * x[8];
	dibe = area*(cbet+cbej) * x[7];
	dIbe = ist/bft * exp(x[0]/nf/vth) * (x[4]/nf/vth)
	+ iset * exp(x[0]/ne/vth) * (x[4]/ne/vth);
	dIbc = ist/brt * exp(x[1]/nr/vth) * (x[5]/nr/vth)
	+ isct * exp(x[1]/nc/vth) * (x[5]/nc/vth);

	ibx = cbx * (x[4] + rbt*(dibc+dibe+dIbc+dIbe));

	AD Ijs;
	AD  cjjs, ijs;
	Ijs = area * isst * (exp(x[2]/ns/vth) - one);

  if (x[2] > zero)
    cjjs = cjst*(one+mjs*x[2]/vjst);
  else
    cjjs = cjst*pow(one-x[2]/vjst,-mjs);
	ijs = Ijs + area*(cjjs * x[6]);

	// Apply Kirchhoff transformation (if requested)
	if (kirchhoff)
		effort[3] = ts/(bm1) *(b - pow(ts/effort[3],bm1));

	flow[0] = Ice - Tibc - ibx - ijs; // collector current
	flow[1] = ibx + Tibc + Tibe;      // base current
	flow[2] = -(Ice + Tibe - ijs);    // emitter current
	flow[3] -= effort[0]*Tibc;

	effort[0] = x[2] + rc*flow[0];                // Collector Substrate Voltage
	effort[1] = x[2] + x[1] + flow[1]*rb;        // Base Substrate Voltage
	effort[2] = x[1] - x[0] + x[2] + flow[2]*re; // Emitter Substrate Voltage
}

