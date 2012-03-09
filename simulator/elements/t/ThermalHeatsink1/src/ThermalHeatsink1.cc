#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "ThermalHeatsink1.h"
#include <cstdio>

#define Np 6			
// Number of terms used in Stehfest 
// algorithm for numerical Laplace inversion.
#define NR_END 0		/* For use in dynamic allocation of arrays. */
#define LN2 0.6931471806
#define sigma 5.67051e-8 	/* Stefan-Boltzmann constant W.m-2.K-4. */


// Static members
// Number of netlist paramenters (remember to update!)
const unsigned ThermalHeatsink1::n_par = 19;

// Element information
ItemInfo ThermalHeatsink1::einfo =
{
  "thermalheatsink1",
  "Grid array substrate: 1-port thermal element with a single averaged surface heating element",
  "Bill Batty, Carlos Christoffersen",
  DEFAULT_ADDRESS"category:thermal",
  "2001_04_20"
};

// Parameter information
// true means the paramenter is required and false means it
// can be omitted in the netlist.
// See ../network/NetListItem.h for a list of types (TR_INT, TR_DOUBLE, etc.)
ParmInfo ThermalHeatsink1::pinfo[] =
{
  {"ntimesteps", "Number of time steps in transient simulation", TR_INT, false},
  {"dt", "Length of timestep (s)", TR_DOUBLE, false},
  {"tambient", "Ambient temperature (K)", TR_DOUBLE, false},
  {"time_d", "Flag, if true, calculate in the time domain.", TR_BOOLEAN, false},
  {"l", "Substrate x-dimension in meters.", TR_DOUBLE, false} ,
  {"w", "Substrate y-dimension in meters", TR_DOUBLE, false},
  {"d", "Substrate z-dimension in meters.", TR_DOUBLE, false},
  {"xl", "x-coordinate of left edge of heating element (m)", TR_DOUBLE, false},
  {"xr", "x-coordinate of right edge of heating element (m)", TR_DOUBLE, false},
  {"yu", "y-coordinate of upper edge of heating element (m)", TR_DOUBLE, false},
  {"yd", "y-coordinate of lower edge of heating element (m)", TR_DOUBLE, false},
  {"xi", "Adjustment for T^4 non linearity.", TR_DOUBLE, false},
  {"eta", "Adjustment for natural convection.", TR_DOUBLE, false},
  {"epsilon", "Emissivity.", TR_DOUBLE, false},
  {"ks", "Thermal conductivity W/m.K.", TR_DOUBLE, false},
  {"rho", "Density kg.m-3.", TR_DOUBLE, false},
  {"c", "Specific heat J/kg.K.", TR_DOUBLE, false},
  {"ndevices", "Number of heat dissipating devices.", TR_INT, false},
	{"b","Exponent in power law temperature dependence of thermal conductivity",TR_DOUBLE,false}
};


ThermalHeatsink1::ThermalHeatsink1(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Setup parameters
  paramvalue[0] = &(Ntimesteps = 0);
  paramvalue[1] = &(dt = zero);
  paramvalue[2] = &(Tambient=300.);
  paramvalue[3] = &(time_d=false);
  paramvalue[4] = &(L = 0.05);
  paramvalue[5] = &(W = 0.05);
  paramvalue[6] = &(D = 0.0016);
  paramvalue[7] = &(xL = 0.04);
  paramvalue[8] = &(xR = 0.01);
  paramvalue[9] = &(yU = 0.04);
  paramvalue[10] = &(yD = 0.01);
  paramvalue[11] = &(xi = 1.3);
  paramvalue[12] = &(eta = 3.0);
  paramvalue[13] = &(epsilon = 0.7);
  paramvalue[14] = &(Ks = 0.294);
  paramvalue[15] = &(rho = 1900.0);
  paramvalue[16] = &(C = 1150.0);
  paramvalue[17] = &(Ndevices = 1);
  paramvalue[18] = &(b = 0.0);

  // Set the number of terminals
  setNumTerms(2);
  // Set flags
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);
  // Set number of states
  setNumberOfStates(1);

  // Set unallocated pointers to 0
  Rth = Pstore = NULL;
}


ThermalHeatsink1::~ThermalHeatsink1()
{
  if (Pstore)
    free_vector(Pstore,1,Ntimesteps);
  if (Rth)
    free_vector(Rth,0,Ntimesteps);
}

// This is the initialization routine.
void ThermalHeatsink1::init() throw(string&)
{
  if (time_d) {
    if (!Ntimesteps || !dt)
      throw(getInstanceName() +
		"Both Ntimesteps and dt must be specified if time_d=1");
    // Change element flags
    setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

    Rth = newvector(0,Ntimesteps); /* Dynamically allocate storage for
		thermal impedance Rth(t) evaluated
		at specified times t. */
    Rth[0] = 0.0; /* At t = 0 temperature rise is always 0. */

    Pstore = newvector(1,Ntimesteps); /* Dynamically allocate storage
		for the self-consistently
		calculated power dissipation
		(Watts) at each timestep. */

    printf("\nPrecomputation for the Rth(t) ...\n\n");
    for (int n = 1; n <= Ntimesteps; n++) {
      double t = n*dt;
      double deltaT = TimeDomain(Tambient,t); /* Calculate temperature rise
			from Tambient as a function
			of time t. */
      Rth[n] = deltaT; /* deltaT = Rth.P and power P is normalised to
			1 W (over specified area of heating
			element). */

      printf("%f %f \n",t,Rth[n] + Tambient);
    }
    printf("\nPrecomputation for the Rth(t) complete.\n\n");
  }
}


void ThermalHeatsink1::fillMNAM(FreqMNAM* mnam)
{
  double SMALL=1.0e0;	/* Should be a safe `small' number, essentially zero. */
  double_complex rth;
  double_complex s(zero, twopi * mnam->getFreq());

  if (abs(s) < SMALL)
  	rth = TimeIndependent(Tambient);
  else if (abs(s) >= SMALL)
  	rth = CLaplaceDomain(Tambient, s);

  // For debugging purposes only
  cout << s << "\t" << rth << endl;
  double_complex gth = one / rth;
  // Ask my terminals the row numbers
  mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), gth);
}

// This routine is called at each time step
// (Actually, several times at each time step because the Jacobian
// is evaluated numerically).
void ThermalHeatsink1::svTran(TimeDomainSV* tdsv)
{
  if(tdsv->DC())
	{
    tdsv->i(0) = tdsv->getX(0);
    tdsv->u(0) = zero;
  }
  else
	{
    int n = tdsv->getIndex();
    // Check consistency with svtr parameters
    assert(n <= Ntimesteps);
    assert(tdsv->getdt() == dt);
    // The state variable is equal to the input power
    tdsv->i(0) = tdsv->getX(0);
    /* The self-consistently obtained value of P at the given timestep
		would be stored in Transim. */
    Pstore[n] = tdsv->i(0);
    // Avoid doing the convolution at each nonlinear iteration.
    if (tdsv->firstTime())
		{
      storedDeltaT = zero;
      // This is the whole thermal impedance calculation (after precomputation)
      for (int m = n; m > 1; m--)
				storedDeltaT += (Rth[m] - Rth[m-1])*Pstore[n-m+1];
    }
    double deltaT = storedDeltaT + (Rth[1] - Rth[0]) * Pstore[n];

    //Implement inverse Kirchhoff transformation
    deltaT = Tambient*(pow(one-(b-one)*deltaT/Tambient,-one/(b-one))-one);

    // Set temperature
    tdsv->u(0) = deltaT;
  }
}

void ThermalHeatsink1::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJi()(0,0) = one;
  if(tdsv->DC())
    tdsv->getJu()(0,0) = zero;
  else
	{
    assert(!tdsv->firstTime());
    // Derivative with no Kirchoff transformation
    //  tdsv->getJu()(0,0) = Rth[1] - Rth[0];
    // With transformation
    double deltaT = storedDeltaT + (Rth[1] - Rth[0]) * tdsv->getX(0);
    tdsv->getJu()(0,0) = pow(one-(b-one)/Tambient*deltaT, b/(one-b))
		* (Rth[1] - Rth[0]);
  }
}

double ThermalHeatsink1::TimeIndependent(double Tambient)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imn, I00;
  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body
	radiation and
	natural convection
	to order (T -
	T0). */

  k = Ks/(rho*C);	/* Diffusivity. */
	P = Ndevices*1.0/((xL-xR)*(yU-yD)); /* W.m-2 */ /* Areal power dissipation
	over grid array heating area (total
	power set at 1 Watt for
	normalisation). */

  mMAX = nMAX = 100;	/* Number of terms in the m and n summations
	of the double series for the thermal
	impedance. */
  I00 = (xL - xR)*(yU - yD);

  result = 0.0;
  for (m = 0; m <= mMAX; m++)
	{
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++)
		{
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n);

      /* Imn are area integrals over the heating element. */

      if (m != 0 && n != 0)
				Imn = (sin(lambda_m*xL) - sin(lambda_m*xR))
			*(sin(mu_n*yU) - sin(mu_n*yD))/(lambda_m*mu_n);
      if (m == 0 && n != 0)
				Imn = (xL - xR)*(sin(mu_n*yU) - sin(mu_n*yD))/mu_n;
      if (m != 0 && n == 0)
				Imn = (sin(lambda_m*xL) - sin(lambda_m*xR))*(yU - yD)/lambda_m;
      if (m == 0 && n == 0)
				Imn = (xL - xR)*(yU - yD);

      /* Insert into analytical expression for thermal impedance. */
			if (!(m == 0 && n == 0))
			{
				result += ( 4.0 / (I00*(1.0 + delta_m0)*(1.0 + delta_n0)) )
				*(H*tanh(gamma_mn*D) + Ks*gamma_mn)
				*Imn*Imn/( (H*H+gamma_mn*gamma_mn*Ks*Ks)*tanh(gamma_mn*D)+ 2.0*H*Ks*gamma_mn );
			}
    }
  }

  result = result + (I00/H)*(1.0 + (H*D/Ks))/(2.0 + (H*D/Ks));
  result = result*P/(L*W);

  return result;
}

/* Function LaplaceDomain returns the Laplace s-space thermal
impedance. */
/* This subroutine contains all the geometrical and material details
required for analytical construction of the s-space thermal
impedance. */
double ThermalHeatsink1::LaplaceDomain(double Tambient, double s)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imn;

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body radiation and natural convection to order (T - T0). */

  k = Ks/(rho*C);	/* Diffusivity. */

  P = Ndevices*1.0/((xL-xR)*(yU-yD)); /* W.m-2 */ /* Areal power dissipation
	over grid array heating area (total
	power set at 1 Watt for
	normalisation). */

  mMAX = nMAX = 100;	/* Number of terms in the m and n summations
	of the double series for the thermal
	impedance. */

  result = 0.0;
  for (m = 0; m <= mMAX; m++)
	{
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++)
		{
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n + s/k);

      /* Imn are area integrals over the heating element. */
      if (m != 0 && n != 0)
				Imn = (sin(lambda_m*xL) - sin(lambda_m*xR))
			*(sin(mu_n*yU) - sin(mu_n*yD))/(lambda_m*mu_n);
      if (m == 0 && n != 0)
				Imn = (xL - xR)*(sin(mu_n*yU) - sin(mu_n*yD))/mu_n;
      if (m != 0 && n == 0)
				Imn = (sin(lambda_m*xL) - sin(lambda_m*xR))*(yU - yD)/lambda_m;
      if (m == 0 && n == 0)
				Imn = (xL - xR)*(yU - yD);

      /* Insert into analytical expression for thermal impedance. */

      result +=
			(cos(m*pi/2.0)*cos(n*pi/2.0) / ((1.0 + delta_m0)*(1.0 + delta_n0)))
			*Imn*(gamma_mn/s) / ((H*H + gamma_mn*gamma_mn*Ks*Ks)*sinh(gamma_mn*D)
			+ 2.0*H*Ks*gamma_mn*cosh(gamma_mn*D));
    }
  }

  result = result*4.0*P*Ks/(L*W);

  return result;
}

double_complex ThermalHeatsink1::CLaplaceDomain(double Tambient, double_complex s)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double_complex result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imn, i(0.0,1.0);

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body radiation and natural convection to order (T - T0). */
  k = Ks/(rho*C);	/* Diffusivity. */

  P = Ndevices*1.0/((xL-xR)*(yU-yD)); /* W.m-2 */ /* Areal power dissipation
	over grid array heating area (total
	power set at 1 Watt for
	normalisation). */

  mMAX = nMAX = 100;	/* Number of terms in the m and n summations
	of the double series for the thermal
	impedance. */

  result = 0.0;
  for (m = 0; m <= mMAX; m++)
	{
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++)
		{
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n + s/k);

      /* Imn are area integrals over the heating element. */
      if (m != 0 && n != 0)
				Imn = (sin(lambda_m*xL) - sin(lambda_m*xR))
			*(sin(mu_n*yU) - sin(mu_n*yD))/(lambda_m*mu_n);
      if (m == 0 && n != 0)
				Imn = (xL - xR)*(sin(mu_n*yU) - sin(mu_n*yD))/mu_n;
      if (m != 0 && n == 0)
				Imn = (sin(lambda_m*xL) - sin(lambda_m*xR))*(yU - yD)/lambda_m;
      if (m == 0 && n == 0)
				Imn = (xL - xR)*(yU - yD);

      /* Insert into analytical expression for thermal impedance. */
      // Calculate exp(-gamma_mn*D)
      double_complex expx = -gamma_mn*D;
      expx = exp(expx.real())*(cos(expx.imag())+i*sin(expx.imag()));

      result +=
			2.0*(cos(m*pi/2.0)*cos(n*pi/2.0) / ((1.0 + delta_m0)*(1.0 + delta_n0)))
			*Imn*(gamma_mn/s)*expx / ((H*H + gamma_mn*gamma_mn*Ks*Ks)*(1.0 - expx*expx)
			+ 2.0*H*Ks*gamma_mn*(1.0 + expx*expx));
    }
  }
  result = result*4.0*P*Ks/(L*W);
  return result;
}

/* Function TimeDomain performs the inverse Laplace transform
numerically based on Stehfest's algorithm. */
/* In fact the inversion is easy to do analytically, but the numerical
approach is general and computationally cheap. */
double ThermalHeatsink1::TimeDomain(double Tambient, double t)
{
  int v;
  double result, increment;

  increment = result = 0;

  for (v=1;v<=Np;v++)
	{
		increment = CalculateWeight(v) * LaplaceDomain(Tambient, v*LN2/t);
		result = result + increment;
	}
  result = result*LN2/t;
  return result;
}


/* Function CalculateWeight constructs the weights for the numerical
Laplace inversion. */
double ThermalHeatsink1::CalculateWeight(int v)
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
  result = result*pow((double)-1,(double)(Np/2+v));
  return result;
}

/* Function factorial returns x! */
double ThermalHeatsink1::fact(int x)
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
void ThermalHeatsink1::nrerror(const char *error_text)
{
  fprintf(stderr,"Run time error....\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"... now exiting to system ..\n");
  exit(1);
}

/* Dynamic memory allocation and deallocation routines. */
double* ThermalHeatsink1::newvector(int nl,int nh)
{
  double *v;
  v=(double *) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("Allocation failure in newvector()");
  return v-nl+NR_END;
}

void ThermalHeatsink1::free_vector(double *v,int nl,int nh)
{
  free((char*) (v+nl-NR_END));
}

