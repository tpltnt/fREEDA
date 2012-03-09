#include "MosnldMet.h"

const unsigned MosnldMet :: n_par = 60;

ItemInfo MosnldMet::einfo =
{
	"mosnldmet",
	"MET LDMOS",
	"Jiankai Chang/Jason Thurston",
	DEFAULT_ADDRESS"transistor>mosfet",
	"2003_05_15"
};

// Parameter Information
ParmInfo MosnldMet :: pinfo[] =
{
	{"rg_0","Gate Resistance Evaluated at Tnom", TR_DOUBLE, false},
	{"rg_1","Gate Resistance Coefficient  ", TR_DOUBLE, false},
	{"rs_0","Source Resistance Evaluated at Tnom", TR_DOUBLE, false},
	{"rs_1","Source Resistance Coefficient  ", TR_DOUBLE, false},
	{"rd_0","Drain Resistance Evaluated at Tnom", TR_DOUBLE, false},
	{"rd_1","Drain Resistance Coefficient  ", TR_DOUBLE, false},
	{"vtO_0","Forward Threshold Voltage Evaluated at Tnom", TR_DOUBLE, false},
	{"vtO_1","Forward Threshold Voltage Coefficient", TR_DOUBLE, false},
	{"gamma","IDS Equation Coefficient  ", TR_DOUBLE, false},
	{"vst","Sub-Threshold Slope Coefficient  ", TR_DOUBLE, false},
	{"beta_0","IDS Equation Coefficient. BETA Evaluated at Tnom", TR_DOUBLE, false},
	{"beta_1","IDS Equation Coefficient  ", TR_DOUBLE, false},
	{"lambda","IDS Equation Coefficient  ", TR_DOUBLE, false},
	{"vgexp","IDS Equation Coefficient  ", TR_DOUBLE, false},
	{"alpha","IDS Equation Coefficient  ", TR_DOUBLE, false},
	{"vk","IDS Equation Coefficient  ", TR_DOUBLE, false},
	{"delta","IDS Equation Coefficient  ", TR_DOUBLE, false},
	{"vbr_0","Breakdown Voltage Evaluated at Tnom", TR_DOUBLE, false},
	{"vbr_1","Breakdown Coefficient @ Vgs=0V", TR_DOUBLE, false},
	{"k1","Breakdown Parameter   ", TR_DOUBLE, false},
	{"k2","Breakdown Parameter   ", TR_DOUBLE, false},
	{"m1","Breakdown Parameter   ", TR_DOUBLE, false},
	{"m2","Breakdown Parameter   ", TR_DOUBLE, false},
	{"m3","Breakdown Parameter   ", TR_DOUBLE, false},
	{"br","Reverse IDS Equation Coefficient", TR_DOUBLE, false},
	{"rdiode_0","Reverse Diode Series Resistance Evaluated at Tnom", TR_DOUBLE, false},
	{"rdiode_1","Reverse Diode Series Resistance Coefficient", TR_DOUBLE, false},
	{"isr","Reverse Diode Leakage Current", TR_DOUBLE, false},
	{"nr","Reverse Diode Ideality Factor", TR_DOUBLE, false},
	{"vtO_r","Reverse Threshold Voltage Coefficient", TR_DOUBLE, false},
	{"rth","Thermal Resistance Coefficient  ", TR_DOUBLE, false},
	{"ggs","Gate To Source Conductance", TR_DOUBLE, false},
	{"ggd","Gate to Drain Conductance", TR_DOUBLE, false},
	{"tau","Transit Time Under Gate", TR_DOUBLE, false},
	{"tnom","Temperature at Which Model Parameters are Extracted", TR_DOUBLE, false},
	{"tsnk","Heat Sink Temp.  ", TR_DOUBLE, false},
	{"cgst","Cgs Temperature Coefficient  ", TR_DOUBLE, false},
	{"cdst","Cds Temperature Coefficient  ", TR_DOUBLE, false},
	{"cgdt","Cgd Temperature Coefficient  ", TR_DOUBLE, false},
	{"cth","Thermal Capacitance   ", TR_DOUBLE, false},
	{"kf","Flicker Noise Coefficient  ", TR_DOUBLE, false},
	{"af","Flicker Noise Exponent  ", TR_DOUBLE, false},
	{"ffe","Flicker Noise Frequency Exponent", TR_DOUBLE, false},
	{"n","Forward Diode Ideality Factor", TR_DOUBLE, false},
	{"iss","Forward Diode Leakage Current", TR_DOUBLE, false},
	{"cgs1","Cgs Equation Coefficient  ", TR_DOUBLE, false},
	{"cgs2","Cgs Equation Coefficient  ", TR_DOUBLE, false},
	{"cgs3","Cgs Equation Coefficient  ", TR_DOUBLE, false},
	{"cgs4","Cgs Equation Coefficient  ", TR_DOUBLE, false},
	{"cgs5","Cgs Equation Coefficient  ", TR_DOUBLE, false},
	{"cgs6","Cgs Equation Coefficient  ", TR_DOUBLE, false},
	{"cgd1","Cgd Equation Coefficient  ", TR_DOUBLE, false},
	{"cgd2","Cgd Equation Coefficient  ", TR_DOUBLE, false},
	{"cgd3","Cgd Equation Coefficient  ", TR_DOUBLE, false},
	{"cgd4","Cgd Equation Coefficient  ", TR_DOUBLE, false},
	{"cds1","Cds Equation Coefficient  ", TR_DOUBLE, false},
	{"cds2","Cds Equation Coefficient  ", TR_DOUBLE, false},
	{"cds3","Cds Equation Coefficient  ", TR_DOUBLE, false},
	{"area","Gate Periphery Scaling Parameter", TR_DOUBLE, false},
	{"n_fing","Gate Finger Scaling Parameter", TR_DOUBLE, false}
};

MosnldMet::MosnldMet(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
	// Set default parameter values
	paramvalue[0]= &(RG_0 = 1);            //in units of ohms
	paramvalue[1]= &(RG_1 = 0.001);        //in units of ohms/K
	paramvalue[2]= &(RS_0 = 0.1);          //in units of ohms
	paramvalue[3]= &(RS_1 = 0.0001);       //in units of ohms/K
	paramvalue[4]= &(RD_0 = 1.5);          //in units of ohms
	paramvalue[5]= &(RD_1 = 0.0015);       //in units of ohms/K
	paramvalue[6]= &(VTO_0 = 3.5);         //in units of V
	paramvalue[7]= &(VTO_1 = -0.001);      //in units of V/K
	paramvalue[8]= &(GAMMA = -0.02);       //is a constant
	paramvalue[9]= &(VST = 0.15);          //in units of V
	paramvalue[10]= &(BETA_0 = 0.2);       //in units of 1/ohms
	paramvalue[11]= &(BETA_1 = -0.0002);   //in units of 1/(ohms*K)
	paramvalue[12]= &(LAMBDA = -0.0025);   //in units of 1/V
	paramvalue[13]= &(VGEXP = 1.1);        //is a constant
	paramvalue[14]= &(ALPHA = 1.5);        //is a constant
	paramvalue[15]= &(VK = 7);             //in units of V
	paramvalue[16]= &(DELTA = 0.9);        //in units of V
	paramvalue[17]= &(VBR_0 = 75);         //in units of V
	paramvalue[18]= &(VBR_1 = 0.01);       //in units of V/K
	paramvalue[19]= &(K1 = 1.5);           //is a constant
	paramvalue[20]= &(K2 = 1.15);          //in units of 1/V
	paramvalue[21]= &(M1 = 9.5);           //is a constant
	paramvalue[22]= &(M2 = 1.2);           //in units of 1/V
	paramvalue[23]= &(M3 = 0.001);         //is a constant
	paramvalue[24]= &(BR = 0.5);           //in units of 1/(V*ohms)
	paramvalue[25]= &(RDIODE_0 = 0.5);     //in units of ohms
	paramvalue[26]= &(RDIODE_1 = 0.001);   //in units of ohms/K
	paramvalue[27]= &(ISR = 1.00E-13);     //in units of A
	paramvalue[28]= &(NR = 1);             //is a constant
	paramvalue[29]= &(VTO_R = 3);          //in units of V
	paramvalue[30]= &(RTH = 10);           //in units of degrees C/Watts
	paramvalue[31]= &(GGS = 1.00E+05);     //in units of 1/ohms
	paramvalue[32]= &(GGD = 1.00E+05);     //in units of 1/ohms
	paramvalue[33]= &(TAU = 1.00E-12);     //in units of seconds
	paramvalue[34]= &(TNOM = 298);         //in units of K
	paramvalue[35]= &(TSNK = 25);          //in units of degrees C
	paramvalue[36]= &(CGST = 0.001);       //in units of 1/K
	paramvalue[37]= &(CDST = 0.001);       //in units of 1/K
	paramvalue[38]= &(CGDT = 0);           //in units of 1/K
	paramvalue[39]= &(CTH = 0);            //in units of J/degrees C
	paramvalue[40]= &(KF = 0);             //is a constant
	paramvalue[41]= &(AF = 1);             //is a constant
	paramvalue[42]= &(FFE = 1);            //is a constant
	paramvalue[43]= &(N = 1);              //is a constant
	paramvalue[44]= &(ISS = 1.00E-13);     //in units of A
	paramvalue[45]= &(CGS1 = 2.00E-12);    //in units of F
	paramvalue[46]= &(CGS2 = 1.00E-12);    //in units of F
	paramvalue[47]= &(CGS3 = -4);          //in units of V
	paramvalue[48]= &(CGS4 = 1.00E-12);    //in units of F
	paramvalue[49]= &(CGS5 = 0.25);        //in units of 1/V
	paramvalue[50]= &(CGS6 = 3.5);         //in units of 1/V
	paramvalue[51]= &(CGD1 = 4.00E-13);    //in units of F
	paramvalue[52]= &(CGD2 = 1.00E-13);    //in units of F
	paramvalue[53]= &(CGD3 = 0.1);         //in units of 1/V^2
	paramvalue[54]= &(CGD4 = 4);           //in units of V
	paramvalue[55]= &(CDS1 = 1.00E-12);    //in units of F
	paramvalue[56]= &(CDS2 = 1.50E-12);    //in units of F
	paramvalue[57]= &(CDS3 = 0.1);         //in units of 1/V^2
	paramvalue[58]= &(AREA = 1);           //is a constant
	paramvalue[59]= &(N_FING = 1);         //is a constant

	setNumTerms(5);

	//Set Flags
	setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);

	setNumberOfStates(3);
}

void MosnldMet::getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list)
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

void MosnldMet::init() throw(string&)
{
	//takes default parameters and scales it with the area of the device
	RD_0 = RD_0/AREA;
	RS_0 = RS_0/AREA;
	RG_0 = RG_0*AREA*N_FING*N_FING;
	RD_1 = RD_1/AREA;
	RS_1 = RS_1/AREA;
	RG_1 = RG_1*AREA*N_FING*N_FING;
	//RDS0 = RDS0/AREA;
	GGD  = GGD*AREA;
	GGS  = GGS*AREA;
	RTH  = RTH/AREA;
	CTH  = CTH*AREA;
	BETA_0 = BETA_0*AREA;
	BETA_1 = BETA_1*AREA;
	CGS1 = CGS1*AREA;
	CGS2 = CGS2*AREA;
	CGS4 = CGS4*AREA;
	CGD1 = CGD1*AREA;
	CGD2 = CGD2*AREA;
	CDS1 = CDS1*AREA;
	CDS2 = CDS2*AREA;
	ISS = ISS * AREA;
	ISR = ISR * AREA;
	BR  = BR * AREA;
	RDIODE_0 = RDIODE_0/AREA;
	RDIODE_1 = RDIODE_1/AREA;

	DenseIntVector var(3);
	var[0] = 0;
	var[1] = 1;
	var[2] = 2;
	DenseIntVector dvar(3);
	dvar[0] = 0;
	dvar[1] = 1;
	dvar[2] = 2;
	DenseIntVector d2var(2);
	d2var[0] = 0;
	d2var[1] = 1;

	DenseIntVector tvar(1);
	DenseDoubleVector delay(1);
	delay[0] = TAU;
	initializeAD(var, dvar, d2var, tvar, delay);
}

//evaluate function
void MosnldMet :: eval(AD * x, AD * effort, AD * flow)
{
	//x[vbe,vbc,vcjs,dvbe,dvbc,dvcjs,d2vbe,d2vbc,d2vcjs]
	//x[0] : Vgs
	//x[1] : Vgd
	//x[2] : temperature - tnom
	//x[3] : dvgs/dt
	//x[4] : dvgd/dt
	//x[5] : dvth_rise/dt
	//x[6] : dvgs/d2t
	//x[7] : dvgd/d2t
	//x[8] : vgs(t-tau)
	//effort[0] : Vgs flow[0] : Ig
	//effort[1] : Vds flow[1] : Id
	//effort[2] : temperature, flow[2]: pout

	AD Vth_rise, T, Rg, Rd, Rs;
	AD Vto_f, Beta, Vbr, Vgst2, Vgst1, Vgst, Vbreff,Vbreff1, Ids, Ids_r, Ids_f;
	AD Vt, Idiode_f, Vto_r, Ism, Idiode_r, Rdiode, Vgst2_r, Vgst1_r, Vgst_r;
	AD Cgs, Cgd, Cds, Igs, Igd, Idiode, Vt2;
	AD Vds, Vdst, Cgst, Cgdt, Vth_rise2;

	double Eg = 1.11;  // energy gap for Silicon
	double k = 1.381E-23;  // Boltzmann's constant
	double q = 1.602E-19;   // electron charge

	T = x[2] + TSNK + 273;

	Rg = RG_0 + RG_1 * ( T - TNOM);
	Rd = RD_0 + RD_1 * ( T - TNOM);
	Rs = RS_0 + RS_1 * ( T - TNOM);

	// gate to source capacitance
	Cgs = (CGS1+CGS2*(1+tanh(CGS6*(x[0]+CGS3))) + CGS4*(1-tanh(x[0]*CGS5)))*(1+CGST*(T-TNOM));
	Cgst = (CGS2*pow(1/cosh(CGS6*(x[0]+CGS3)),2)*CGS6*x[3] - CGS4*pow(1/cosh(x[0]*CGS5),2)*CGS5*x[3])*(1+CGST*(T-TNOM));
	// gate to drain capacitance
	Cgd = (CGD1+CGD2/(1+CGD3*pow((x[1] - CGD4), 2))) * (1+CGDT*(T-TNOM));
	Cgdt = -(1+CGDT*(T-TNOM))*CGD2*pow((1+CGD3*pow((x[1] - CGD4), 2)),-2) * (2*CGD3*(x[1]-CGD4))*x[4];

	Igs = x[3] * Cgs;
	Igd = x[4] * Cgd;
	Vds = x[0] + Igs/GGS - x[1] - Igd/GGD;
	Vdst = x[3] + (x[6]*Cgs + x[3]*Cgst)/GGS - x[4] - (x[7]*Cgd + x[4]*Cgdt)/GGD;

	// drain to source capacitance
	Cds = (CDS1+CDS2/(1+CDS3*Vds*Vds)) * (1+CDST*(T-TNOM));

	// forward bias drain to source current
	Vto_f = VTO_0 + VTO_1 * ( T - TNOM);
	Beta = BETA_0 + BETA_1 * ( T - TNOM);
	Vbr = VBR_0 + VBR_1 * ( T - TNOM);
	Vgst2 = x[8] - (Vto_f + (GAMMA * Vds));
	Vgst1 = Vgst2 - 1/2 * (Vgst2 + sqrt(pow((Vgst2-VK),2) + pow(DELTA,2)) - sqrt(pow(VK,2)+pow(DELTA,2)));
	Vgst = VST * log(exp(Vgst1/VST) + 1);
	Vbreff = Vbr/2 * (1 + tanh(M1 - Vgst * M2));
	Vbreff1 = 1/K2 * (Vds - Vbreff) + M3*(Vds/Vbreff);
	Ids_f = Beta * pow(Vgst, VGEXP) * (1+LAMBDA*Vds) * tanh(Vds * ALPHA/ Vgst) * (1+K1* exp(Vbreff1));

	// reverse bias drain to source current
	Vto_r = VTO_R + VTO_1 * (T-TNOM);
	Vgst2_r = x[8] - (Vto_r - (GAMMA*Vds));
	Vgst1_r = Vgst2_r - 1/2 * (Vgst2_r + sqrt( pow((Vgst2_r-VK),2) + pow(DELTA,2)) - sqrt(pow(VK,2)+pow(DELTA,2)));
	Vgst_r = VST * log(exp(Vgst1_r/VST) + 1);
	Ids_r = BR * Vds * Vgst_r;

	// drain to source current
  if (Vds > 0.0)
    Ids = Ids_f;
  else
    Ids = Ids_r;

	// forward bias drain to source diode
	Vt = k * T / q;
	Idiode_f = ISS * exp((Vds - Vbr)/ (N * Vt));

	// reverse bias drain to source diode
	Vt2 = k * T / q;
	Ism = ISR * pow(T/TNOM, 3/NR) * exp(-Eg/NR/Vt2*(1-T/TNOM));
	Idiode_r = Ism * (exp(Vds/(NR*Vt2)) - 1);          // Vdiode_r  -> Vds
	//Rdiode = RDIODE_0 + RDIODE_1 * (T - TNOM);

	// drain to source diode current
  if (Vds > 0.0)
    Idiode = Idiode_f;
  else
    Idiode = 0-Idiode_r;

	flow[0] = Igs + Igd;      // Ig
	flow[1] = Ids - Igd + Vdst*Cds + Idiode;      // Id
	flow[2] = flow[0]*flow[0]*Rg + Igs*Igs/GGS + Igd*Igd/GGD
            + (flow[0]+flow[1])*(flow[0]+flow[1])*Rs
            + Igs*x[0]+Igd*x[2]+ Vds*(Ids + Vdst*Cds + Idiode);

	//Vgs is calculated using the current throught the gate and the source and the resistance
	//and charge of the device
	effort[0] = x[0] + Rg*flow[0] + Igs/GGS + (flow[0]+flow[1])*Rs;      // Gate source Voltage

	//Vds is calculated using the current throught the drain and the source and the resistance
	//and charge of the device
	effort[1] = Vds + (flow[0]+flow[1])*Rs + flow[1]*Rd;        // Drain Source Voltage

	AD Vth_rise3, Vth_rise4;
	Vth_rise3 = x[2];
  if (Vth_rise3 > 250.0)
    Vth_rise2 = 250.0 + 50.0*tanh((Vth_rise3-250.0)/50.0);
  else
    Vth_rise2 = Vth_rise3;

  if (Vth_rise3 > 0.0)
    Vth_rise4 = Vth_rise2;
  else
    Vth_rise4 = Vth_rise3;

	//Thermal effects is calculated using the current throught the device
	//the temperature and the thermal resistance and capacitance.
	effort[2] = Vth_rise4 + TSNK + 273 - (flow[2] - x[5]*CTH) * RTH;
}

