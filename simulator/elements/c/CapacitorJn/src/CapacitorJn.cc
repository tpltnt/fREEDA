#include "CapacitorJn.h"

// Static members
const unsigned CapacitorJn::n_par = 26;

// Element information
ItemInfo CapacitorJn::einfo =
{
  "capacitorjn",
  "Philips MOS9 -- CapacitorJn model",
  "Yogesh Ramdoss, Kuldip Gothi, Xuemin Yang, Ajit Rajagopalan, Dapeng ding, Xin Cai",
  DEFAULT_ADDRESS"category:lumped,semiconductor,diode",
  "2003_05_15"
};

// Parameter information
ParmInfo CapacitorJn::pinfo[] =
{
	{"level", "level of this model. Must be set to 1", TR_DOUBLE, false},
	{"ab", "Diffusion area (m2)", TR_DOUBLE, false},
	{"ls", "Length of the side-wall of the diffusion area ab which is not under the gate(m)",TR_DOUBLE, false},
	{"lg", "Length of the side-wall of the diffusion area ab which is under the gate (m)", TR_DOUBLE, false},
	{"dta", "Temperature offset of the JUNCAP element with respect to TA (oC)", TR_DOUBLE, false},
	{"tr", "Temperature at which the parameters have been determined (oC)", TR_DOUBLE, false},
	{"vr", " Voltage at which parameters have been determined(V)", TR_DOUBLE, false},
	{"jsgbr", "Bottom saturation-current density due to electron-hole generation at V = vr (Am-2)", TR_DOUBLE, false},
	{"jsdbr", "Bottom saturation-current density due to diffusion from back contact (Am-2)", TR_DOUBLE, false},
	{"jsgsr", "Sidewall saturation-current density due to electron-hole generation at V = vr (Am-1)", TR_DOUBLE, false},
	{"jsdsr", "Sidewall saturation-current density due to diffusion from back contact (Am-1)", TR_DOUBLE, false},
	{"jsggr", "Gate edge saturation-current density due to electron-hole generation at V = vr (Am-1)", TR_DOUBLE, false},
	{"jsdgr", "Gate edge saturation-current density due to diffusion from back contact (Am-1)", TR_DOUBLE, false},
	{"nb", "Emission coefficient of the bottom forward current", TR_DOUBLE, false},
	{"ns", "Emission coefficient of the sidewall forward current",TR_DOUBLE, false},
	{"ng", "Emission coefficient of the gate edge forward current",TR_DOUBLE, false},
	{"vb", "Reverse breakdown voltage (V)", TR_DOUBLE, false},
	{"cjbr", "Bottom junction capacitance at V = vr (Fm-2)", TR_DOUBLE, false},
	{"cjsr", "Sidewall junction capacitance at V = vr (Fm-1)", TR_DOUBLE, false},
	{"cjgr", "Gate edge junction capacitance at V = vr (Fm-1)", TR_DOUBLE, false},
	{"vdbr", "Diffusion voltage of the bottom junction at T = tr(V)", TR_DOUBLE, false},
	{"vdsr", "Diffusion voltage of the sidewall junction at T = tr(V)", TR_DOUBLE, false},
	{"vdgr", "Diffusion voltage of the gate edge junction at T = tr(V)", TR_DOUBLE, false},
	{"pb", "Bottom-junction grading coefficient", TR_DOUBLE, false},
	{"ps", "Sidewall-junction grading coefficient", TR_DOUBLE, false},
	{"pg", "Gate-edge-junction grading coefficient", TR_DOUBLE, false}
};


CapacitorJn::CapacitorJn(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
	// Set default parameter values
	paramvalue[0] = &(level = 1.00);
	paramvalue[1] = &(ab = 1e-12);
	paramvalue[2] = &(ls = 1e-6);
	paramvalue[3] = &(lg = 1e-6);
	paramvalue[4] = &(dta = 0.00);
	paramvalue[5] = &(tr = 25.00);
	paramvalue[6] = &(vr = 0.00);
	paramvalue[7] = &(jsgbr = 1e-3);
	paramvalue[8] = &(jsdbr = 1e-3);
	paramvalue[9] = &(jsgsr = 1e-3);
	paramvalue[10] = &(jsdsr = 1e-3 );
	paramvalue[11] = &(jsggr = 1e-3);
	paramvalue[12] = &(jsdgr = 1e-3);
	paramvalue[13] = &(nb = 1.00);
	paramvalue[14] = &(ns = 1.00);
	paramvalue[15] = &(ng = 1.00);
	paramvalue[16] = &(vb = 0.90); //the value for VB is not used
	paramvalue[17] = &(cjbr = 1e-12);
	paramvalue[18] = &(cjsr = 1e-12);
	paramvalue[19] = &(cjgr = 1e-12);
	paramvalue[20] = &(vdbr = 1.00);
	paramvalue[21] = &(vdsr = 1.00);
	paramvalue[22] = &(vdgr = 1.00);
	paramvalue[23] = &(pb = 0.40);
	paramvalue[24] = &(ps = 0.40);
	paramvalue[25] = &(pg = 0.40);

	// Set the number of terminals
	setNumTerms(2);

	// Set number of states
	setNumberOfStates(1);

	// Set flags
	setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);
}

void CapacitorJn::init() throw(string&)
{
	DenseIntVector var(1);
	var[0] = 0;
	DenseIntVector novar;
	DenseDoubleVector nodelay;
	initializeAD(var, var, var, novar, nodelay);
}

void CapacitorJn :: eval(AD * x, AD * effort, AD * flow)
{
	// x[0]:    V(V) Potential applied to the anode-cathode
	// x[1]:    dV/dt
	// effort[0]:    V
	// flow[0]:    I  DC current into the anode + dQ/dt

	//All the active variables "AD" are initialiazed here

	//All the Internal Variables and Parameters
	AD vdb;   //VDB (V) Diffusion voltage of bottom area ab
	AD vds;   //VDS (V) Diffusion voltage of Locos-edge ls
	AD vdg;   //VDG (V) Diffusion voltage of gate-edge lg
	AD cjb;   //CJB (F) Capacitance of bottom area ab
	AD cjs;   //CJS (F) Capacitance of Locos-edge ls
	AD cjg;   //CJG (F) Capacitance of gate-edge lg
	AD isdb;  //ISDB (A) Diffusion saturation-current of bottom area ab
	AD isds;  //ISDS (A) Diffusion saturation-current of Locos-edge ls
	AD isdg;  //ISDG (A) Diffusion saturation-current of gate-edge lg
	AD isgb;  //ISGB (A) Generation saturation-current of bottom area ab
	AD isgs;  //ISGS (A) Generation saturation-current of Locos-edge ls
	AD isgg;  //ISGG (A) Generation saturation-current of gate-edge lg
	AD ta;    //TA (0C) Ambient circuit temperature
	AD tkd;   //TKD (K) absolute temperature of the junction/device
	AD v;     //V(V) Diode bias voltage (V = VA - VK)
	AD i;     //I(A) Total DC current from anode to cathode (I = IA =- IK)
	AD Q ;    //Q(c) Total junction charge (Q = QA = - QK)

	AD tkr,phitr,phitd,phigr,phigd,ftd,c;

	AD qlb,qjbv,qls,qjsv,qlg,qjgv,qjds,qjdg,qjdb;
	AD cjbv,cjsv,cjgv,clb,cls,clg;
	AD idb,igb,ids,igs,idg,igg;
	AD vlb,vls,vlg,vab,vas,vag;
	AD fcb,fcs,fcg,fsb,fss,fsg;
	AD dqjbvdt,dqjsvdt,dqjgvdt,dvdt,dQdt;

	//constant
	double t0 = 273.15; //Offset for conversion from Celsius to Kelvin temperature scale
	double k = 1.3806226e-23; //Boltzmann constant(JK-1)
	double q = 1.6021918e-19; //Elementary unit charge(C)
	//Absolute permittivity of the oxide layer(Fm-1)
	// double EPIox = 3.453143800e-11;
	double b = 2 ;

	//v = va - vk
	v = x[0];
	// dvdt = dvadt - dvkdt
	dvdt = x[1];

	//general scaling rules, which apply to all three components of the JUNCAP model
	tkr = t0 + tr;
	tkd = t0 + ta + dta;
	phitr = k * tkr / q;
	phitd = k * tkd / q;
	phigr = 1.16 - (7.02e-4 * tkr * tkr/ (1108.0 + tkr));
	phigd = 1.16 - (7.02e-4 * tkd * tkd/ (1108.0 + tkd));
	ftd = pow(tkd / tkr,1.5) * exp(phigr / (2 * phitr) - phigd / (2 * phitd));

	//The internal reference parameters for the bottom component
	vdb = vdbr * tkd / tkr - 2 * phitd * log(ftd);
	cjb = cjbr * ab * pow((vdbr - vr) / vdb,pb);
	isgb = jsgbr * ftd * ab * pow(vdb / (vdbr - vr),pb);
	isdb = jsdbr * pow(ftd,2) * ab;

	//The internal reference parameters for the locos-edge component
	vds = vdsr * tkd / tkr - 2 * phitd * log(ftd);
	cjs = cjsr * ls * pow((vdsr - vr) / vds,ps);
	isgs = jsgsr * ftd * ls * pow(vds / (vdsr - vr),ps);
	isds = jsdsr * pow(ftd,2) * ls;

	//The internal reference parameters for the gate-edge component
	vdg = vdgr * tkd / tkr - 2 * phitd * log(ftd);
	cjg = cjgr * lg * pow((vdgr - vr) / vdg,pg);
	isgg = jsggr * ftd * ls * pow(vdg / (vdgr - vr),pg);
	isdg = jsdgr * pow(ftd,2) * lg;

	//Capacitor and Leakage Current Model for the bottom component

	// charge
	qjdb = cjb * vdb / (1 - pb);
	fcb = 1 - pow((1 + pb)/3,1 / pb);
	vlb = fcb * vdb;
	clb = cjb * pow(1 - fcb, - pb);
	qlb = qjdb * (1 - pow(1 - fcb, 1 - pb));
  if (v > vlb)
  {
    qjbv = qlb + clb * (v - vlb) * (1 + pb * (v - vlb)/ (2 * vdb * (1 - fcb)));
    dqjbvdt = clb * dvdt + clb * pb * (v - vlb) * dvdt / ( vdb * (1 - fcb));
    cjbv = clb + clb * pb * (v - vlb) / ( vdb * (1 - fcb));
  }
  else
  {
    qjbv = qjdb * (1 - pow(1 - v/vdb, 1 - pb));
    dqjbvdt = qjdb * (1 - pb) * pow (1 - v / vdb, -pb) *(-1/vdb) * dvdt;
    cjbv = cjb / pow(1 - v / vb, pb);
  }

	// bulk to source or bulk to drain diode current
	vab = b * vdb / pb;
	fsb = isgb * pow (vab, b);
	idb = isdb * (exp (v / (nb * phitd)) - 1);
  if (v > 0.0)
    igb = isgb * pow (vab/(v + vab), b) * (1 - exp ( -v / (nb * phitd)));
  else
    igb = isgb * pow ((vdb - v)/vdb,pb) * (exp (v / (nb * phitd)) - 1);

	//Capacitor and Leakage Current Model for the locos-edge component
	// charge
	qjds = cjs * vds / (1 - ps);
	fcs = 1 - pow((1 + ps)/3,1 / ps);
	vls = fcs * vds;
	cls = cjs * pow(1 - fcs, - ps);
	qls = qjds * (1 - pow(1 - fcs, 1 - ps));
  if (v > vls)
  {
    qjsv = qls + cls * (v - vls) * (1 + ps * (v - vls)/ (2 * vds * (1 - fcs)));
    dqjsvdt = cls * dvdt + cls * ps * (v - vls) * dvdt/(vds * (1 - fcs));
    cjsv = cls + cls * ps * (v - vls) / (vds * (1 - fcs));
  }
  else
  {
    qjsv = qjds * (1 - pow(1 - v/vds, 1 - ps));
    dqjsvdt = qjds * (1-ps) * pow (1-v/vds, -ps)*(-1/vds) * dvdt;
    cjsv = cjs / pow(1 - v / vds, ps);
  }

	// bulk to source or bulk to drain diode current
	vas = b * vds / ps;
	fss = isgs * pow (vas, b);
	ids = isds * (exp (v / (ns * phitd)) - 1);
  if (v > 0.0)
    igs = isgs * pow (vas/(v + vas), b) * (1 - exp ( -v / (ns * phitd)));
  else
    igs = isgs * pow ((vds - v)/vds,ps) * (exp (v / (ns * phitd)) - 1);

	//Capacitor and Leakage Current Model for the gate-edge component
	// charge
	qjdg = cjg * vdg / (1 - pg);
	fcg = 1 - pow((1 + pg)/3,1 / pg);
	vlg = fcg * vdg;
	clg = cjg * pow(1 - fcg, - pg);
	qlg = qjdg * (1 - pow(1 - fcg, 1 - pg));
  if (v > vlg)
  {
    qjgv = qlg + clg * (v - vlg) * (1 + pg * (v - vlg)/ (2 * vdg * (1 - fcg)));
    dqjgvdt = clg * dvdt + clg * pg * (v - vlg) * dvdt / ( vdg * (1 - fcg));
    cjgv = clg + clg * pg * (v - vlg) / ( vdg * (1 - fcg));
  }
  else
  {
    qjgv = qjdg * (1 - pow(1 - v/vdg, 1 - pg));
    dqjgvdt = qjdg * (1 - pg) * pow (1 - v / vdg, -pg) *(-1/vdg) * dvdt;
    cjgv = cjg / pow(1 - v / vdg, pg);
  }

	// bulk to source or bulk to drain diode current
	vag = b * vdg / pg;
	fsg = isgg * pow (vag, b);
	idg = isdg * (exp (v / (ng * phitd)) - 1);
  if (v > 0.0)
    igg = isgg * pow (vag/(v + vag), b) * (1 - exp (-v / (ng * phitd)));
  else
    igg = isgg * pow ((vdg - v)/vdg,pg) * (exp (v / (ng * phitd)) - 1);

	// total charge
	Q = qjbv + qjsv + qjgv;
	// total capacitance
	c = cjbv + cjsv + cjgv;
	// total junction current
	i = idb + igb + ids + igs + idg + igg;

	dQdt = dqjbvdt + dqjsvdt + dqjgvdt;

	//current = DC current + dQ/dt
	flow[0] = i + dQdt;
	//voltage = x[0]
	effort[0] = x[0];
}
