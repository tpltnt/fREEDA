/* This version edited for Semi-Therm XVII form of Rth00(t). */


#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "ThermalShunt.h"
#include <cstdio>

#define Np 6			/* Number of terms used in Stehfest
                                   algorithm for numerical Laplace
                                   inversion. */
#define NR_END 0		/* For use in dynamic allocation of arrays. */

#define LN2 0.6931471806

// Static members
// Number of netlist paramenters (remember to update!)
const unsigned ThermalShunt::n_par = 17;

// Element information
ItemInfo ThermalShunt::einfo = {
  "thermalshunt",
  "2-port MMIC die with variable base temperature (for interface matching), and a single averaged surface heating element",
  "Bill Batty, Carlos Christoffersen",
  DEFAULT_ADDRESS"category:thermal",
  "2001_06_20"
 };

// Parameter information
// true means the paramenter is required and false means it 
// can be omitted in the netlist.
// See ../network/NetListItem.h for a list of types (INT, DOUBLE, etc.)
ParmInfo ThermalShunt::pinfo[] = {
  {"ntimesteps", "Number of time steps in transient simulation", 
   TR_INT, false},
  {"dt", "Length of timestep (s)", TR_DOUBLE, false},
  {"tambient", "Ambient temperature (K)", TR_DOUBLE, false},
  {"time_d", "Flag, if true, calculate in the time domain.", TR_BOOLEAN, 
   false},
  {"read_input", "Flag, read_input thermal resistance matrices from file.", 
   TR_BOOLEAN, false},
  {"l", "Substrate x-dimension in meters.", TR_DOUBLE, false},
  {"w", "Substrate y-dimension in meters.", TR_DOUBLE, false},
  {"d", "Substrate z-dimension in meters.", TR_DOUBLE, false},
  {"xl", "x-coordinate of left edge of heating element.", TR_DOUBLE, false},
  {"xr", "x-coordinate of right edge of heating element.", TR_DOUBLE, false},
  {"yu", "y-coordinate of upper edge of heating element.", TR_DOUBLE, false},
  {"yd", "y-coordinate of lower edge of heating element.", TR_DOUBLE, false},
  {"ks", "Thermal conductivity of die material (W/m.K).", TR_DOUBLE, false},
  {"rho", "Density of die material (kg.m-3).", TR_DOUBLE, false},
  {"c", "Specific heat of die material (J/kg.K).", TR_DOUBLE, false},
  {"nfingers", "Number of power transistor fingers", TR_INT, false},
  {"b","Exponent in power law temperature dependence of thermal conductivity",TR_DOUBLE,false}
};


ThermalShunt::ThermalShunt(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Setup parameters
  paramvalue[0] = &(Ntimesteps = 0);
  paramvalue[1] = &(dt = zero);
  paramvalue[2] = &(Tambient=300.);
  paramvalue[3] = &(time_d = false);
  paramvalue[4] = &(read_input = false);
  paramvalue[5] = &(L = 400.0e-6);
  paramvalue[6] = &(W = 400.0e-6);
  paramvalue[7] = &(D = 400.0e-6);
  paramvalue[8] = &(xL = 220.0e-6);
  paramvalue[9] = &(xR = 180.0e-6);
  paramvalue[10] = &(yU = 220.0e-6);
  paramvalue[11] = &(yD = 180.0e-6);
  paramvalue[12] = &(Ks = 46.0);
  paramvalue[13] = &(rho = 5320.0);
  paramvalue[14] = &(C = 350.0);
  paramvalue[15] = &(Nfingers = 1);
  paramvalue[16] = &(b = 1.22);

  // Set the number of terminals
  setNumTerms(3);
  // Set flags
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);
  // Set number of states
  setNumberOfStates(2);

  // Set unallocated pointers to 0
  Rth00 = Rth0D = RthD0 = RthDD = P0store = PDstore = NULL;
}


ThermalShunt::~ThermalShunt()
{
  if (PDstore)
    free_vector(P0store,1,Ntimesteps);
  if (P0store)
    free_vector(PDstore,1,Ntimesteps);
  if (RthDD)
    free_vector(RthDD,0,Ntimesteps);
  if (RthD0)
    free_vector(RthD0,0,Ntimesteps);
  if (Rth0D)
    free_vector(Rth0D,0,Ntimesteps);
  if (Rth00)
    free_vector(Rth00,0,Ntimesteps);
}

// This is the initialization routine.
void ThermalShunt::init() throw(string&)
{ 
  FILE *fp;
  double input1, input2, input3, input4;

  Pnorm0 = double(Nfingers) / ((xL-xR)*(yU-yD)); /* W.m-2 Nfingers scales
  the normalised power density to treat multi-finger power devices. */

  /* Areal power dissipation over grid array heating area (total power
     set at 1 Watt for normalisation). */
  PnormD = -1.0 / (L*W); /* NOTE sign: in Transim, flux is positive into element.  In Rth approach flux is positive along z-axis, i.e. INTO top of MMIC, z=0, and OUT of MMIC base, z=D. */
 
  k = Ks/(rho*C);	/* Diffusivity. */ 

  if (time_d) {
    if (!Ntimesteps || !dt)
      throw(getInstanceName() + 
	    ": Both Ntimesteps and dt must be specified if time_d=1");
    // Change element flags
    setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

    Rth00 = newvector(0,Ntimesteps); /* Dynamically allocate storage for
                                     thermal impedances Rth(t)
                                     evaluated at specified times
                                     t. */
    Rth0D = newvector(0,Ntimesteps);
    RthD0 = newvector(0,Ntimesteps);
    RthDD = newvector(0,Ntimesteps);
    
    Rth00[0] = 0.0; /* At t = 0 temperature rise is always 0. */
    Rth0D[0] = 0.0;
    RthD0[0] = 0.0;
    RthDD[0] = 0.0;
    
    P0store = newvector(1,Ntimesteps); /* Dynamically allocate storage
                                       for the self-consistently
                                       calculated power dissipation
                                       (Watts) at each timestep. */
    PDstore = newvector(1,Ntimesteps);


    //----------- Allow storage and re-use of precomputed Rth(t).   
    if (!(read_input)) fp = fopen("Rth_elements.dat","a");
    if (read_input) fp = fopen("Rth_elements.dat","r");

    for (int n = 1; n <= Ntimesteps; n++) {	
      double t = n*dt;
      /* Calculate temperature rise from
	 Tambient as a function of time t. */
      /* deltaT = Rth.P and power P is normalised to 1 W (over
         specified area of heating element). */
      if (!(read_input) ){
	if (n==1) 
	  printf("\nPrecomputation for the Rth(t) ...\n\n");	
//	  cout<<"print time step "<<t<<endl;
	Rth00[n] = TimeDomain00(t);
	Rth0D[n] = TimeDomain0D(t);
	RthD0[n] = TimeDomainD0(t);
	RthDD[n] = TimeDomainDD(t);	
  		printf("%f %f %f %f %f\n",
  		       t,Rth00[n],Rth0D[n],RthD0[n],RthDD[n]);
	fprintf(fp,"%f %f %f %f\n",
		Rth00[n],Rth0D[n],RthD0[n],RthDD[n]);
      }
      else 
	if (read_input) {
	  if (n==1) 
	    printf("\nRead file input for the Rth(t) ...\n\n");
	  fscanf(fp,"%lf %lf %lf %lf",&input1,&input2,&input3,&input4);
	  Rth00[n] = input1;
	  Rth0D[n] = input2;
	  RthD0[n] = input3;
	  RthDD[n] = input4;
	  // printf("%f %f %f %f %f\n",t,Rth00[n],Rth0D[n],RthD0[n],RthDD[n]);
	}
    }
    printf("\nPrecomputation/read for the Rth(t) complete.\n\n");

    fclose(fp);
  }
}
  
unsigned ThermalShunt::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  // Keep the equation number assigned to this element
  my_row1 = eqn_number;
  my_row2 = eqn_number + 1;
  // Add one extra RC
  return 2;
}

void ThermalShunt::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row1);
  first_eqn = my_row1;
  n_rows = 2;
}


void ThermalShunt::fillMNAM(FreqMNAM* mnam)
{
double SMALL=1.0e0; // Should be a safe `small' number, essentially zero. 
  double_complex rth00, rth0D, rthD0, rthDD;
  double_complex s(zero, twopi * mnam->getFreq());
  
  if (abs(s) < SMALL) {
    // The element behave like this at very low frequencies
    //
    //                  rth00
    //  P0,T0 o--------/\/\/\/\-------o PD,TD
    // 
    //        o-----------------------o
    // 
    mnam->setOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row1);
    rth00 = - TimeIndependent00();
    mnam->setElement(my_row1, my_row1, rth00);
    // Additional equation P0 = PD
    mnam->setElement(my_row2, my_row1, one);
    mnam->setElement(my_row2, my_row2, -one);
  }
  else if (abs(s) >= SMALL) {
    // Set structural ones
    mnam->setOnes(getTerminal(0)->getRC(), getTerminal(2)->getRC(), my_row1);
    mnam->setOnes(getTerminal(1)->getRC(), getTerminal(2)->getRC(), my_row2);

    // Set z elements 
    rth00 = - CLaplaceDomain00(s);
    mnam->setElement(my_row1, my_row1, rth00);
    rth0D = - CLaplaceDomain0D(s);
    mnam->setElement(my_row1, my_row2, rth0D);
    rthD0 = - CLaplaceDomainD0(s);
    mnam->setElement(my_row2, my_row1, rthD0);
    rthDD = - CLaplaceDomainDD(s);
    mnam->setElement(my_row2, my_row2, rthDD);
  }  	
}
  

// This routine is called at each time step
// (Actually, several times at each time step because the Jacobian
// is evaluated numerically).
void ThermalShunt::svTran(TimeDomainSV* tdsv)
{
  if(tdsv->DC()) {
    tdsv->i(0) = tdsv->getX(0);
    tdsv->u(0) = 0;
    tdsv->i(1) = tdsv->getX(1);
    tdsv->u(1) = 0;
  }
  else {
    int n = tdsv->getIndex();
    // Check consistency with svtr parameters
    assert(n <= Ntimesteps);
    assert(tdsv->getdt() == dt);
    // The state variable is equal to the input power
    tdsv->i(0) = tdsv->getX(0);
    tdsv->i(1) = tdsv->getX(1);
    /* Value of P and the timestep n would be set in Transim. */
    P0store[n] = tdsv->i(0); /* The self-consistently obtained value
                                of P at the given timestep would be
                                stored in Transim. */
    PDstore[n] = tdsv->i(1); /* Two values of P set, corresponding to
                                z = 0 (top) and z = D (base) of
                                MMIC. */
    /* P should be supplied as its numerical value in Watts. */
    if (tdsv->firstTime()) {
      storedDeltaT0 = zero; 
      storedDeltaTD = zero; 
      // This is the whole thermal impedance calculation (after precomputation)
      for (int m = n; m > 1; m--) {
	storedDeltaT0 += (Rth00[m] - Rth00[m-1])*P0store[n-m+1] 
	  + (Rth0D[m] - Rth0D[m-1])*PDstore[n-m+1];
	storedDeltaTD += (RthD0[m] - RthD0[m-1])*P0store[n-m+1] 
	  + (RthDD[m] - RthDD[m-1])*PDstore[n-m+1];
      }
    }
    double deltaT0 = storedDeltaT0 
      + (Rth00[1]-Rth00[0]) * P0store[n] + (Rth0D[1]-Rth0D[0]) * PDstore[n];
    double deltaTD = storedDeltaTD 
      + (RthD0[1]-RthD0[0]) * P0store[n] + (RthDD[1]-RthDD[0]) * PDstore[n];

    //Implement inverse Kirchhoff transformation
    deltaT0 = Tambient*(pow(1.0-(b-1.0)*deltaT0/Tambient,-1.0/(b-1.0))-1.0);
    deltaTD = Tambient*(pow(1.0-(b-1.0)*deltaTD/Tambient,-1.0/(b-1.0))-1.0);

    // Set temperature
    tdsv->u(0) = deltaT0;
    tdsv->u(1) = deltaTD;
  }
}

void ThermalShunt::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJi()(0,0) = one;
  tdsv->getJi()(1,1) = one;
  tdsv->getJi()(0,1) = zero;
  tdsv->getJi()(1,0) = zero;
  if(tdsv->DC()) {
    tdsv->getJu()(0,0) = zero;
    tdsv->getJu()(0,1) = zero;
    tdsv->getJu()(1,0) = zero;
    tdsv->getJu()(1,1) = zero;
  }
  else {
    assert(!tdsv->firstTime());

    double deltaT0 = storedDeltaT0 
      + (Rth00[1]-Rth00[0]) * tdsv->getX(0) 
      + (Rth0D[1]-Rth0D[0]) *tdsv->getX(1);
    double deltaTD = storedDeltaTD 
      + (RthD0[1]-RthD0[0]) * tdsv->getX(0) 
      + (RthDD[1]-RthDD[0]) * tdsv->getX(1);

    tdsv->getJu()(0,0) = pow(one-(b-one)/Tambient*deltaT0, b/(one-b)) 
      * (Rth00[1] - Rth00[0]);
    tdsv->getJu()(0,1) = pow(one-(b-one)/Tambient*deltaT0, b/(one-b)) 
      * (Rth0D[1] - Rth0D[0]);
    tdsv->getJu()(1,0) = pow(one-(b-one)/Tambient*deltaTD, b/(one-b)) 
      * (RthD0[1] - RthD0[0]);
    tdsv->getJu()(1,1) = pow(one-(b-one)/Tambient*deltaTD, b/(one-b)) 
      * (RthDD[1] - RthDD[0]);
  }
}

double ThermalShunt::TimeIndependent00()
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, lambda_m, mu_n, gamma_mn, Imn;

  mMAX = nMAX = 75;	/* Number of terms in the m and n summations of the double series for the thermal impedance. */	
	
  result = 0.0;
  for (m = 0; m <= mMAX; m++){
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++){
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n);

      /* Imn are area integrals over the heating element. */

      if (m != 0 && n != 0) Imn = (sin(lambda_m*xL) - sin(lambda_m*xR))*(sin(mu_n*yU) - sin(mu_n*yD))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imn = (xL - xR)*(sin(mu_n*yU) - sin(mu_n*yD))/mu_n;
      if (m != 0 && n == 0) Imn = (sin(lambda_m*xL) - sin(lambda_m*xR))*(yU - yD)/lambda_m;
      if (m == 0 && n == 0) Imn = (xL - xR)*(yU - yD);

      /* Insert into analytical expression for thermal impedance. */

      if (!(m == 0 && n == 0)) result += 4.0*Imn*Imn*tanh(gamma_mn*D)/(gamma_mn*(1.0 + delta_m0)*(1.0 + delta_n0)*(xL - xR)*(yU - yD));
				
    }
  }

  result += D*(xL - xR)*(yU - yD);    
  result = result*Pnorm0/(Ks*L*W);
	
  return result;
}


/* Function LaplaceDomain returns the Laplace s-space thermal impedance. */
/* This subroutine contains all the geometrical and material details required for analytical construction of the s-space thermal impedance. */
double ThermalShunt::LaplaceDomain00(double s)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, lambda_m, mu_n, gamma_mn, Imn;

  mMAX = nMAX = 75;	/* Number of terms in the m and n summations of the double series for the thermal impedance. */	
	
  result = 0.0;
  for (m = 0; m <= mMAX; m++){
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++){
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n + s/k);

      /* Imn are area integrals over the heating element. */

      if (m != 0 && n != 0) Imn = (sin(lambda_m*xL) - sin(lambda_m*xR))*(sin(mu_n*yU) - sin(mu_n*yD))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imn = (xL - xR)*(sin(mu_n*yU) - sin(mu_n*yD))/mu_n;
      if (m != 0 && n == 0) Imn = (sin(lambda_m*xL) - sin(lambda_m*xR))*(yU - yD)/lambda_m;
      if (m == 0 && n == 0) Imn = (xL - xR)*(yU - yD);

      /* Insert into analytical expression for thermal impedance. */

      result += 4.0*Imn*Imn*(1.0/tanh(gamma_mn*D))/(gamma_mn*s*(1.0 + delta_m0)*(1.0 + delta_n0)*(xL - xR)*(yU - yD));				
    }
  }

/****  result += 2.0*sqrt(k/s)*(xL - xR)*(yU - yD)/(s*sinh(2.0*sqrt(s/k)*D)); ****/    
  result = result*Pnorm0/(Ks*L*W);
	
  return result;
}



double ThermalShunt::LaplaceDomain0D(double s)
{
  double result;

  result = -(1.0/Ks)*sqrt(k/s)*PnormD/(s*sinh(sqrt(s/k)*D));

  return result;
}



double ThermalShunt::LaplaceDomainD0(double s)
{
  double result;

  result = (1.0/Ks)*sqrt(k/s)*(xL - xR)*(yU - yD)*Pnorm0/(s*L*W*sinh(sqrt(s/k)*D));

  return result;
}



double ThermalShunt::LaplaceDomainDD(double s)
{
  double result;
double r_cosh = 0;
double r_sinh=0;

	r_cosh = cosh(sqrt(s/k)*D);
	r_sinh = sinh(sqrt(s/k)*D);
	
//	cout<<"R_cosh "<<r_cosh<<" R_sinh "<<r_sinh<<endl;
	
  result = -(1.0/Ks)*sqrt(k/s)*(r_cosh/r_sinh)*PnormD/s;

  return result;
}



double_complex ThermalShunt::CLaplaceDomain00(double_complex s)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double lambda_m, mu_n, Imn;
  double_complex  gamma_mn, result;

  mMAX = nMAX = 75;	/* Number of terms in the m and n summations of the double series for the thermal impedance. */	
	
  result = zero;
  for (m = 0; m <= mMAX; m++){
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++){
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n + s/k);

      /* Imn are area integrals over the heating element. */

      if (m != 0 && n != 0) Imn = (sin(lambda_m*xL) - sin(lambda_m*xR))*(sin(mu_n*yU) - sin(mu_n*yD))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imn = (xL - xR)*(sin(mu_n*yU) - sin(mu_n*yD))/mu_n;
      if (m != 0 && n == 0) Imn = (sin(lambda_m*xL) - sin(lambda_m*xR))*(yU - yD)/lambda_m;
      if (m == 0 && n == 0) Imn = (xL - xR)*(yU - yD);

      /* Insert into analytical expression for thermal impedance. */


	    // Calculate exp(-gamma_mn*D)
      double_complex expx = -gamma_mn*D, imag_i(zero,one);
      expx = exp(expx.real())*(cos(expx.imag())+imag_i*sin(expx.imag()));
	
      result += 4.0*Imn*Imn*((one-expx*expx)/(one+expx*expx))/(gamma_mn*(1.0 + delta_m0)*(1.0 + delta_n0)*(xL - xR)*(yU - yD));
				
    }
  }

	// Calculate exp(-2.0*sqrt(s/k)*D)
	double_complex exps = -2.0*sqrt(s/k)*D, imag_i(zero,one);
	exps = exp(exps.real())*(cos(exps.imag())+imag_i*sin(exps.imag()));

  result += 2.0*sqrt(k/s)*(xL - xR)*(yU - yD)*2.0*exps/(1.0-exps*exps);    
  result = result*Pnorm0/(Ks*L*W);
	
  return result;
}



double_complex ThermalShunt::CLaplaceDomain0D(double_complex s)
{
  double_complex result;

	// Calculate exp(-sqrt(s/k)*D)
	double_complex exps = -sqrt(s/k)*D, imag_i(zero,one);
	exps = exp(exps.real())*(cos(exps.imag())+imag_i*sin(exps.imag()));

  result = -(1.0/Ks)*sqrt(k/s)*PnormD*2.0*exps/(1.0-exps*exps);

  return result;
}



double_complex ThermalShunt::CLaplaceDomainD0(double_complex s)
{
  double_complex result;

	// Calculate exp(-sqrt(s/k)*D)
	double_complex exps = -sqrt(s/k)*D, imag_i(zero,one);
	exps = exp(exps.real())*(cos(exps.imag())+imag_i*sin(exps.imag()));

  result = (1.0/Ks)*sqrt(k/s)*(xL - xR)*(yU - yD)*Pnorm0*2.0*exps/(L*W*(1.0-exps*exps));

  return result;
}



double_complex ThermalShunt::CLaplaceDomainDD(double_complex s)
{
  double_complex result;

// Calculate exp(-sqrt(s/k)*D)
	double_complex exps = -sqrt(s/k)*D, imag_i(zero,one);
	exps = exp(exps.real())*(cos(exps.imag())+imag_i*sin(exps.imag()));


  result = -(1.0/Ks)*sqrt(k/s)*((1.0+exps*exps)/(1.0-exps*exps))*PnormD;

  return result;
}




/* Function TimeDomain performs the inverse Laplace transform numerically based on Stehfest's algorithm. */
/* In fact the inversion is easy to do analytically, but the numerical approach is general and computationally cheap. */
double ThermalShunt::TimeDomain00(double t)
{
  int v;
  double result, increment;
	
  increment = result = 0;
	
  for (v=1;v<=Np;v++)
    {
      increment = CalculateWeight(v) * LaplaceDomain00(v*LN2/t);
      result = result + increment;
    }			
  result = result*LN2/t;
  return result;
}


double ThermalShunt::TimeDomain0D(double t)
{
  int v;
  double result, increment;
	
  increment = result = 0;
	
  for (v=1;v<=Np;v++)
    {
      increment = CalculateWeight(v) * LaplaceDomain0D(v*LN2/t);
      result = result + increment;
    }			
  result = result*LN2/t;
  return result;
}


double ThermalShunt::TimeDomainD0(double t)
{
  int v;
  double result, increment;
	
  increment = result = 0;
	
  for (v=1;v<=Np;v++)
    {
      increment = CalculateWeight(v) * LaplaceDomainD0(v*LN2/t);
      result = result + increment;
    }			
  result = result*LN2/t;
  return result;
}



double ThermalShunt::TimeDomainDD(double t)
{
  int v;
  double result, increment;
	
  increment = result = 0;

  for (v=1;v<=Np;v++)
    {
      increment = CalculateWeight(v) * LaplaceDomainDD(v*LN2/t);
      result = result + increment;

    }			
  result = result*LN2/t;
  return result;
}


/* Function CalculateWeight constructs the weights for the numerical
   Laplace inversion. */
double ThermalShunt::CalculateWeight(int v)
{
  int NumberOfTerms, k, start;
  double result, increment;
	
	
  if (v<(Np/2))
    NumberOfTerms = v;
  else 
    NumberOfTerms = Np/2;
	
  if (v==1) start = 1;
  if (v==2) start = 1;
  if (v==3) start = 2;
  if (v==4) start = 2;
  if (v==5) start = 3;
  if (v==6) start = 3;
  result = increment = 0;	
	
  for(k = start; k<=NumberOfTerms; k++)
    {
	
      increment = pow((double)k,Np/2) * fact(2*k);
      increment = increment / (fact(Np/2-k) * fact(k) * fact(k-1) * fact(v-k) * fact(2*k-v));
      result = result + increment;
		
    }
  result = result*pow(-1.0,Np/2.0+v);
  return result;
}


/* Function factorial returns x! */
double ThermalShunt::fact(int x)
{	
  int i;
  double cumul = 1;
  for(i=1;i<=x;i++) cumul = cumul *  i;
  if (x==0) cumul = 1.0;
  if (x < 0) {
    printf("Argument less than zero in factorial\n");
    exit(1);
  }
  return cumul;
}	


/* nrerror prints an error message to stderr on failure. */

void ThermalShunt::nrerror(const char *error_text)
{
  fprintf(stderr,"Run time error....\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"... now exiting to system ..\n");
  exit(1);
}

/* Dynamic memory allocation and deallocation routines. */

double* ThermalShunt::newvector(int nl,int nh)	
{ 
  double *v;
  v=(double *) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("Allocation failure in newvector()"); 
  return v-nl+NR_END;
}

void ThermalShunt::free_vector(double *v,int nl,int nh)	
{
  free((char*) (v+nl-NR_END));
}

