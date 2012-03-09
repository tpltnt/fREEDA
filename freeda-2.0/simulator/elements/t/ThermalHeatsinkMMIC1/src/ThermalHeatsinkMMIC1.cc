#include "../../../../network/ElementManager.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "ThermalHeatsinkMMIC1.h"
#include <cstdio>

#define Np 6			/* Number of terms used in Stehfest
algorithm for numerical Laplace
inversion. */
#define NR_END 0		/* For use in dynamic allocation of arrays. */

#define LN2 0.6931471806

// Static members
// Number of netlist paramenters (remember to update!)
const unsigned ThermalHeatsinkMMIC1::n_par = 18;

// Element information
ItemInfo ThermalHeatsinkMMIC1::einfo =
{
  "thermalheatsinkmmic1",
  "Heatsink mounted MMIC: 1-port thermal element with a single averaged surface heating element",
  "Bill Batty, Carlos Christoffersen mailto:w.batty@elec-eng.leeds.ac.uk",
  DEFAULT_ADDRESS"category:thermal",
  "2001_04_20"
};

// Parameter information
// true means the paramenter is required and false means it
// can be omitted in the netlist.
// See ../network/NetListItem.h for a list of types (TR_INT, TR_DOUBLE, etc.)
ParmInfo ThermalHeatsinkMMIC1::pinfo[] =
{
  {"ntimesteps", "Number of time steps in transient simulation",
	TR_INT, false},
  {"dt", "Length of timestep (s)", TR_DOUBLE, false},
  {"tambient", "Constant heatsink mount temperature (K)", TR_DOUBLE, false},
  {"time_d", "Flag, if true, calculate in the time domain.", TR_BOOLEAN,
	false},
  {"l", "Die x-dimension in meters.", TR_DOUBLE, false},
  {"w", "Die y-dimension in meters.", TR_DOUBLE, false},
  {"d", "Die z-dimension in meters.", TR_DOUBLE, false},
  {"xl", "x-coordinate of left edge of heating element.", TR_DOUBLE, false},
  {"xr", "x-coordinate of right edge of heating element.", TR_DOUBLE, false},
  {"yu", "y-coordinate of upper edge of heating element.", TR_DOUBLE, false},
  {"yd", "y-coordinate of lower edge of heating element", TR_DOUBLE, false},
  {"ks", "Thermal conductivity (W/m.K).", TR_DOUBLE, false},
  {"rho", "Density (kg.m-3).", TR_DOUBLE, false},
  {"c", "Specific heat (J/kg.K).", TR_DOUBLE, false},
  {"nfingers", "Number of power transistor fingers", TR_INT, false},
  {"b","Exponent in power law temperature dependence of thermal conductivity",TR_DOUBLE,false},
  {"zero_sv", "Flag, if true element does not add state variables.",
	TR_BOOLEAN, false},
  {"kt","Flag, if false Kirchhoff transformation is not used with zero_sv=1.",
	TR_BOOLEAN, false}
};


ThermalHeatsinkMMIC1::ThermalHeatsinkMMIC1(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Setup parameters
  paramvalue[0] = &(Ntimesteps = 0);
  paramvalue[1] = &(dt = zero);
  paramvalue[2] = &(Tambient=300.);
  paramvalue[3] = &(time_d = false);
  paramvalue[4] = &(L = 400.0e-6);
  paramvalue[5] = &(W = 400.0e-6);
  paramvalue[6] = &(D = 400.0e-6);
  paramvalue[7] = &(xL = 220.0e-6);
  paramvalue[8] = &(xR = 180.0e-6);
  paramvalue[9] = &(yU = 220.0e-6);
  paramvalue[10] = &(yD = 180.0e-6);
  paramvalue[11] = &(Ks = 46.0);
  paramvalue[12] = &(rho = 5320.0);
  paramvalue[13] = &(C = 350.0);
  paramvalue[14] = &(Nfingers = 1);
  paramvalue[15] = &(b = 1.22);
  paramvalue[16] = &(zero_sv = false);
  paramvalue[17] = &(kt = true);

  // Set the number of terminals
  setNumTerms(2);

  // Set flags
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);

  // Set number of states
  setNumberOfStates(1);

  // Set unallocated pointers to 0
  Rth = Pstore = NULL;

  my_row = 0;
}


ThermalHeatsinkMMIC1::~ThermalHeatsinkMMIC1()
{
  if (Pstore)
    free_vector(Pstore,1,Ntimesteps);
  if (Rth)
    free_vector(Rth,0,Ntimesteps);
}

// This is the initialization routine.
void ThermalHeatsinkMMIC1::init() throw(string&)
{
  if (time_d) {
    if (!Ntimesteps || !dt)
      throw(getInstanceName() +
		": Both Ntimesteps and dt must be specified if time_d=1");
    // Change element flags
    if (zero_sv)
      setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN | SOURCE | CONVOLUTION);
    else
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

    // Initialize history variables
    nstep = 1;
    storedDeltaT = zero;
  }
}

unsigned ThermalHeatsinkMMIC1::getExtraRC(const unsigned& eqn_number,
const MNAMType& type)
{
  if (type == TIME_DOMAIN)
	{
    // Keep the equation number assigned to this element
    my_row = eqn_number;
    // Add one extra RC
    return 1;
  }
  else
    return 0;
}

void ThermalHeatsinkMMIC1::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row);
  first_eqn = my_row;
  n_rows = 1;
}

void ThermalHeatsinkMMIC1::fillMNAM(TimeMNAM* mnam)
{
  assert(my_row);
  // Ask my terminals the row numbers
  mnam->setMOnes(getTerminal(0)->getRC(), getTerminal(1)->getRC(), my_row);

  // Add thermal resistance here
  mnam->setMElement(my_row, my_row, Rth[0] - Rth[1]);
  return;
}

void ThermalHeatsinkMMIC1::setLastResult(DenseDoubleVector& res, const double& time)
{
  // Assume the time step is constant for now
  assert(abs(time-storedTime-dt) < dt*1e-15*nstep);
  // Remember mnam index starts from 1 (res vector starts from 0)
  Pstore[nstep] = res[my_row - 1];
  nstep++;
  assert(nstep <= Ntimesteps);
  storedTime = time;
  return;
}

void ThermalHeatsinkMMIC1::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime();

  // Add summation here
  //  int n = tdsv->getIndex();
  // Check consistency with svtr parameters
  assert(abs(ctime-storedTime-dt) < dt*1e-15*nstep);
  assert(nstep <= Ntimesteps);
  // Perform summation
  storedDeltaT = zero;
  // This is the whole thermal impedance calculation (after precomputation)
  for (int m = nstep; m > 1; m--)
    storedDeltaT += (Rth[m] - Rth[m-1]) * Pstore[nstep-m+1];

  mnam->setSource(my_row, storedDeltaT);
  return;
}

void ThermalHeatsinkMMIC1::fillMNAM(FreqMNAM* mnam)
{
  /* Should be a safe `small' number, essentially zero. */
  const double SMALL=1.0e0;
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
// (Actually, several times at each time step because of the nonlinear
// iterations)
void ThermalHeatsinkMMIC1::svTran(TimeDomainSV* tdsv)
{
  if(tdsv->DC())
	{
    tdsv->i(0) = tdsv->getX(0);
    tdsv->u(0) = 0;
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
				storedDeltaT += (Rth[m] - Rth[m-1]) * Pstore[n-m+1];
    }
    double deltaT = storedDeltaT + (Rth[1] - Rth[0]) * Pstore[n];
    if (kt)
      //Implement inverse Kirchhoff transformation
		deltaT = Tambient
		* (pow(one-(b-one)*deltaT/Tambient, one/(one-b)) - one);

    // Set temperature
    tdsv->u(0) = deltaT;
  }
}

void ThermalHeatsinkMMIC1::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJi()(0,0) = one;
  if(tdsv->DC())
    tdsv->getJu()(0,0) = zero;
  else
	{
    assert(!tdsv->firstTime());
    if (kt)
		{
      // With transformation
      double deltaT = storedDeltaT + (Rth[1] - Rth[0]) * tdsv->getX(0);
      tdsv->getJu()(0,0) = pow(one-(b-one)/Tambient*deltaT, b/(one-b))
			* (Rth[1] - Rth[0]);
    }
    else
      // Derivative with no Kirchoff transformation
		tdsv->getJu()(0,0) = Rth[1] - Rth[0];
  }
}

double ThermalHeatsinkMMIC1::TimeIndependent(double Tambient)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, k, Imn;
  k = Ks/(rho*C);	/* Diffusivity. */
	P = Nfingers*1.0/((xL-xR)*(yU-yD)); /* W.m-2 */ /* Areal power dissipation
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
			if (!(m == 0 && n == 0)){
				result +=
				4.0*(Imn*Imn*(1.0/((xL - xR)*(yU - yD))) / ((1.0 + delta_m0)*(1.0 + delta_n0)))
				*tanh(gamma_mn*D)*(Imn/((xL - xR)*(yU - yD)))/(gamma_mn);
			}
    }
  }

  result = result + D*(xL - xR)*(yU - yD);
  result = result*P/(Ks*L*W);

  return result;
}

/* Function LaplaceDomain returns the Laplace s-space thermal
impedance. */
/* This subroutine contains all the geometrical and material details
required for analytical construction of the s-space thermal
impedance. */
double ThermalHeatsinkMMIC1::LaplaceDomain(double Tambient, double s)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, k, Imn;

  k = Ks/(rho*C);	/* Diffusivity. */

  P = Nfingers*1.0/((xL-xR)*(yU-yD)); /* W.m-2 */ /* Areal power dissipation
	over grid array heating area (total
	power set at 1 Watt for
	normalisation). */

  mMAX = nMAX = 100;	/* Number of terms in the m and n summations
	of the double series for the thermal
	impedance. */

  result = 0.0;
  for (m = 0; m <= mMAX; m++){
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++){
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
			(Imn*Imn*(1.0/((xL - xR)*(yU - yD))) / ((1.0 + delta_m0)*(1.0 + delta_n0)))*tanh(gamma_mn*D)/(gamma_mn*s);
    }
  }

  result = result*4.0*P/(Ks*L*W);

  return result;
}

double_complex ThermalHeatsinkMMIC1::CLaplaceDomain(double Tambient, double_complex s)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double_complex result, gamma_mn, i(0.0,1.0);
  double P, lambda_m,mu_n, k, Imn;

  k = Ks/(rho*C);	/* Diffusivity. */

  P = Nfingers*1.0/((xL-xR)*(yU-yD)); /* W.m-2 */ /* Areal power dissipation
	over grid array heating area (total
	power set at 1 Watt for
	normalisation). */

  mMAX = nMAX = 100;	/* Number of terms in the m and n summations
	of the double series for the thermal
	impedance. */

  result = 0.0;
  for (m = 0; m <= mMAX; m++){
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++){
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

			result += (1.0 / ((1.0+delta_m0)*(1.0+delta_n0)))*
			Imn*Imn*(1.0/((xL - xR)*(yU - yD)))*((1.0-expx*expx)/(1.0+expx*expx))  /(gamma_mn);

    }
  }

  result = result*4.0*P/(Ks*L*W);
  return result;
}

/* Function TimeDomain performs the inverse Laplace transform
numerically based on Stehfest's algorithm. */
/* In fact the inversion is easy to do analytically, but the numerical
approach is general and computationally cheap. */
double ThermalHeatsinkMMIC1::TimeDomain(double Tambient, double t)
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
double ThermalHeatsinkMMIC1::CalculateWeight(int v)
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
double ThermalHeatsinkMMIC1::fact(int x)
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

void ThermalHeatsinkMMIC1::nrerror(const char *error_text)
{
  fprintf(stderr,"Run time error....\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"... now exiting to system ..\n");
  exit(1);
}

/* Dynamic memory allocation and deallocation routines. */

double* ThermalHeatsinkMMIC1::newvector(int nl,int nh)
{
  double *v;
  v=(double *) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("Allocation failure in newvector()");
  return v-nl+NR_END;
}

void ThermalHeatsinkMMIC1::free_vector(double *v,int nl,int nh)
{
  free((char*) (v+nl-NR_END));
}

