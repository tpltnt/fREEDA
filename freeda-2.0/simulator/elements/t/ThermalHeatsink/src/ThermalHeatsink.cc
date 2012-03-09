#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "ThermalHeatsink.h"
#include <cstdio>

#define Np 6 /* Number of terms used in Stehfest algorithm for
                numerical Laplace inversion. */
#define NR_END 0 /* For use in dynamic allocation of arrays. */

#define FREE_ARG char*

#define LN2 0.6931471806

// Static members
// Number of netlist paramenters (remember to update!)
const unsigned ThermalHeatsink::n_par = 21;
// Stefan-Boltzmann constant W.m-2.K-4. 
const double ThermalHeatsink::sigma = 5.67051e-8;

// Element information
ItemInfo ThermalHeatsink::einfo = {
  "thermalheatsink",
  "N-port grid array substrate with NxN surface heating elements",
  DEFAULT_ADDRESS"category:thermal",
  "2001_04_20"
};

// Parameter information
// true means the paramenter is required and false means it 
// can be omitted in the netlist.
// See ../network/NetListItem.h for a list of types (TR_INT, TR_DOUBLE, etc.)
ParmInfo ThermalHeatsink::pinfo[] = {
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
  {"yd", "y-coordinate of lower edge of heating element", TR_DOUBLE, false},
  {"ks", "Thermal conductivity (W/m.K).", TR_DOUBLE, false},
  {"rho", "Density  (kg.m-3).", TR_DOUBLE, false},
  {"c", "Specific heat (J/kg.K).", TR_DOUBLE, false},
  {"xi","Adjustment for T**4 non linearity", TR_DOUBLE, false},
  {"eta", "Adjustment for natural convection", TR_DOUBLE, false},
  {"epsilon", "Emissivity", TR_DOUBLE, false},
  {"narray", "Order of NxN grid array", TR_INT, false},
  {"ndevices", "Number of heat dissipating devices", TR_INT, false},
  {"b","Exponent in power law temperature dependence of thermal conductivity",TR_DOUBLE,false}
  
};


ThermalHeatsink::ThermalHeatsink(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Setup parameters
  paramvalue[0] = &(Ntimesteps = 0);
  paramvalue[1] = &(dt = zero);
  paramvalue[2] = &(Tambient=300.);
  paramvalue[3] = &(time_d = false);
  paramvalue[4] = &(read_input = false);
  paramvalue[5] = &(L = .05);
  paramvalue[6] = &(W = .05);
  paramvalue[7] = &(D = .0016);
  paramvalue[8] = &(xl = 220.0e-6);
  paramvalue[9] = &(xr = 180.0e-6);
  paramvalue[10] = &(yu = 220.0e-6);
  paramvalue[11] = &(yd = 180.0e-6);
  paramvalue[12] = &(Ks = .294);
  paramvalue[13] = &(rho = 1900.);
  paramvalue[14] = &(C = 1150.);
  paramvalue[15] = &(xi = 1.3);
  paramvalue[16] = &(eta = 3.);
  paramvalue[17] = &(epsilon = .7);
  paramvalue[18] = &(Narray=3);
  paramvalue[19] = &(Ndevices = 1);
  paramvalue[20] = &(b = 0.0);

  // Set flags
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);

  // Set some unallocated pointers to 0
  Rth = Pstore = 0;
  Rthtemp = Pstoretotal = storedDeltaT = 0;

  Nsquared = 0;
}


ThermalHeatsink::~ThermalHeatsink()
{
  if (Pstore)
    free_matrix(Pstore,1,Ntimesteps,1,Nsquared);
  if (Rth)
    free_matrix(Rth,0,Ntimesteps,1,Nsquared);
  if (Rthtemp)
    free_vector(Rthtemp,1,Nsquared);
  if (Pstoretotal)
    free_vector(Pstoretotal,1,Ntimesteps);
  if (storedDeltaT)
    delete [] storedDeltaT;
  if (Nsquared) {
    delete [] xL;
    delete [] xR;
    delete [] yU;
    delete [] yD;
  }
}

// This is the initialization routine.
void ThermalHeatsink::init() throw(string&)
{
  FILE *fp;
  double input;

  if (Narray < 1)
    throw(getInstanceName() + 
	  ": Narray must be greater than zero.");
  Nsquared = Narray * Narray;
  // Allocate memory for vectors
  xL = new double[Nsquared+1];
  xR = new double[Nsquared+1];
  yU = new double[Nsquared+1];
  yD = new double[Nsquared+1];

  // Set number of terminals and state variables
  setNumTerms(Nsquared + 1);
  // resize vector with MNAM equation row numbers.
  my_row.resize(Nsquared);
  // Set number of state variables
  setNumberOfStates(Nsquared);

  /* Element coordinates could easily be read from a data file, to
     allow easy input of different grid array layouts. */

  /* Use a single averaged heating element, to represent the effects
   of heat spreading by surface metallisation, between power
   dissipating elements. */



  for (int m = 1; m <= Narray; m++){	
    for (int n = 1; n <= Narray; n++){
      int l = (m-1)*Narray + n;
      /* x-coordinate of left edge of heating element l (meters). */
      xL[l] = xr + (m)*(xl-xr)/(Narray);
      /* x-coordinate of right edge of heating element l (meters). */
      xR[l] = xr + ( m - 1.0)*(xl-xr)/(Narray);
      /* y-coordinate of upper edge of heating element l (meters). */
      yU[l] = yd + (n)*(yu-yd)/(Narray);
      /* y-coordinate of lower edge of heating element l (meters). */
      yD[l] = yd + (n - 1.0)*(yu-yd)/(Narray); 
      cout<<"l "<<l<<" xl "<<xL[l]<<" xr "<<xR[l]<<" yu"<< yU[l]<<" yd "<<yD[l]<<endl;
    }
  }

  /* x-coordinate of left edge of total heating area (meters). */
  XL = xl;
  /* x-coordinate of right edge of total heating area (meters). */
  XR = xr; 
  /* y-coordinate of upper edge of total heating area (meters). */
  YU = yu; 
  /* y-coordinate lower edge of total heating area (meters). */
  YD = yd;


  if (time_d) {
    // Check parameters
    if (!Ntimesteps || !dt)
      throw(getInstanceName() + 
	    ": Both Ntimesteps and dt must be specified if time_d=1");

    // Change element flags
    setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

    
    /* Thermal impedance vector Rth_i(t) is precomputed at all
       required timesteps prior to the iterative coupled
       electro-thermal solution. */
    Rth = newmatrix(0,Ntimesteps,1,Nsquared*Nsquared); /* Dynamically allocate
                                              storage for thermal
                                              impedance vector
                                              Rth_i(t) evaluated at
                                              specified times t. */
    for (int i = 1; i <= Nsquared*Nsquared; i++){
      Rth[0][i] = 0.0;		/* At t = 0 temperature rise is always 0. */
    }
    Rthtemp = newvector(1,Nsquared*Nsquared);
    Pstore = newmatrix(1,Ntimesteps,1,Nsquared);
    /* Dynamically allocate storage for the self-consistently 
       calculated power dissipation Pj (Watts) at each timestep. */
    Pstoretotal = newvector(1,Ntimesteps);
    storedDeltaT = new double[Nsquared];

    //----------- Allow storage and re-use of precomputed Rth(t).   
    if (!(read_input)) 
      fp = fopen("Kovar_Rth_elements.dat","w");
    if (read_input) 
      fp = fopen("Kovar_Rth_elements.dat","r");
    
    for (int n = 1; n <= Ntimesteps; n++){	
      double t = dt * n;
      if (!(read_input) ){
	if (n==1) printf("\nPrecomputation for the Rth(t) ...\n\n");
	// printf("%g\n",t);
	TimeDomain(t,Rthtemp); // Calculate Rth_i(t) as a function of time t.
	for (int i = 1; i <= Nsquared*Nsquared; i++){
	  Rth[n][i] = Rthtemp[i];		
	  //	  printf("%g ",Rth[n][i]);
	  fprintf(fp,"%f \n",Rth[n][i]);
	}
	//	printf("\n");
//	fprintf(fp,"\n");
      }
      else if (read_input) {
	if (n==1) 
	  printf("\nRead file input for the Rth(t) ...\n\n");
	for (int i = 1; i <= Nsquared*Nsquared; i++){
	  fscanf(fp,"%lf ",&input);
	  Rth[n][i] = input;		
	}
      }
    }
    printf("Precomputation/read for the Rth_i(t) complete.\n\n");
    fclose(fp);
  }
}

  
unsigned ThermalHeatsink::getExtraRC(const unsigned& eqn_number, 
				  const MNAMType& type)
{
  // Keep the equation numbers assigned to this element
  for (int i=0; i < Nsquared; i++) 
    my_row[i] = eqn_number + i;

  // Add Nsquared extra RCs
  return Nsquared;
}

void ThermalHeatsink::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row[0]);
  first_eqn = my_row[0];
  n_rows = Nsquared;
}


void ThermalHeatsink::fillMNAM(FreqMNAM* mnam)
{
  /* Should be a safe `small' number, essentially zero. */
  const double SMALL=1.0e0;	
  double_complex rth;
  double_complex s(zero, twopi * mnam->getFreq());
  
  for (int i=0; i < Nsquared; i++) {
    // Uncomment this when you have TimeIndependent()
    if (abs(s) < SMALL) 
    	rth = - TimeIndependent(i+1);
    else if (abs(s) >= SMALL)
    rth = - LaplaceDomain(i+1, s);
    cout<<"Rth "<<rth<<endl;
    // Set structural ones
    mnam->setOnes(getTerminal(i)->getRC(), getTerminal(Nsquared)->getRC(),
		  my_row[i]);
    // Fill elements is row i of thermal block
    // All the elements in one row are the same, right?
    for (int j=0; j < Nsquared; j++)
      mnam->setElement(my_row[i], my_row[j], rth);
  }
}
  

// This routine is called at each time step
// (Actually, several times at each time step)
void ThermalHeatsink::svTran(TimeDomainSV* tdsv)
{
  if(tdsv->DC()) {
    for (int i=0; i < Nsquared; i++) {
      tdsv->i(i) = tdsv->getX(i);
      tdsv->u(i) = zero;
    }
  }
  else {
    int n = tdsv->getIndex();
    // Check consistency with svtr parameters
    assert(n <= Ntimesteps);
    assert(tdsv->getdt() == dt);

    // This is the whole thermal impedance calculation (after precomputation)
    Pstoretotal[n] = zero;
    for (int i = 1; i <= Nsquared; i++){ 
      tdsv->i(i-1) = tdsv->getX(i-1);
      // There seems to be no need for Pstore[][]
            Pstore[n][i] = tdsv->i(i-1);	
      /* The self-consistently obtained value of P_i at the given
         timestep would be stored in Transim. */
      Pstoretotal[n] += tdsv->getX(i-1);
      /* Only total power dissipation is required for temperature rise
         calculation. */
    }
    /* Calculate temperature rise of each grid array heating element. */
    for (int  i = 1; i <= Nsquared; i++){
      // Avoid doing the convolution at each nonlinear iteration.
      if (tdsv->firstTime()) {
	storedDeltaT[i-1] = zero; 
	for (int m = n; m > 1; m--){
          for (int  j = 1; j <= Nsquared; j++)
	  storedDeltaT[i-1] += (Rth[m][(i-1)*Nsquared+j] - Rth[m-1][(i-1)*Nsquared+j])*Pstore[n-m+1][j-1];
	}
      }
      double deltaTi = storedDeltaT[i-1];
       for (int  j = 1; j <= Nsquared; j++)
	deltaTi += (Rth[1][(i-1)*Nsquared+j] - Rth[0][(i-1)*Nsquared+j]) * Pstore[n][j-1];
      /* Given Ptotal at specified timestep, the values of the deltaTi
         are the final results. */
      /* deltaTi is given in Kelvin (and t in seconds). */

      //Implement inverse Kirchhoff transformation
      deltaTi = Tambient*(pow(one-(b-one)*deltaTi/Tambient,-one/(b-one))-one);
      
      tdsv->u(i-1) = deltaTi;
    }
  }
}

void ThermalHeatsink::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->cleanJac();
  // Ji is an identity matrix
  if(tdsv->DC()) {
    for (int  i = 0; i < Nsquared; i++)
      tdsv->getJi()(i,i) = one;
  }
  else {
    double pst = zero;
    for (int i = 1; i <= Nsquared; i++)
      pst += tdsv->getX(i-1);

    assert(!tdsv->firstTime());
    for(int i =1; i<=Nsquared; i++)
	for(int j =1; j<=Nsquared; j++){
	tdsv->getJu()(i-1,j-1) = Rth[1][(i-1)*Nsquared+j];
	}
    for (int  i = 0; i < Nsquared; i++) {
        tdsv->getJi()(i,i) = one;
        double deltaTi = storedDeltaT[i] ;
        for(int j =0 ; j < Nsquared; j++)
		deltaTi += (Rth[1][i*Nsquared+j+1] - Rth[0][i*Nsquared+j+1]) * tdsv->getX(j);
      
      	double dJi = pow(one-(b-one)/Tambient*deltaTi, b/(one-b)) *(Rth[1][i+1] - Rth[0][i+1]);
      
        for (int  j = 0; j < Nsquared; j++) {
	tdsv->getJu()(i,j) *= dJi;
      }
    }
  }
}


/* Function LaplaceDomain returns the Laplace s-space thermal
   impedance vector elements. */
/* This subroutine contains all the geometrical and material details
   required for analytical construction of the s-space thermal
   impedance. */
double ThermalHeatsink::LaplaceDomain(int i, int j, double s)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imni,I00i,Imnj,I00j;

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body radiation and natural convection to order (T - T0). */

  k = Ks/(rho*C);	/* Diffusivity. */
P=1.0;
  //P = 1.0/((XL-XR)*(YU-YD)); /* W.m-2 */ 
  /* Areal power dissipation over the whole grid array heating area (total power set at 1 Watt for normalisation). */

  mMAX = 50;	/* Number of terms in the m and n summations of the double series for the thermal impedance. */	
  nMAX = 50;	
  result = 0.0;
  for (m = 0; m <= mMAX; m++){
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++){
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n + s/k);

      /* Imni are area integrals over the individual heating elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xL[i]) - sin(lambda_m*xR[i]))*(sin(mu_n*yU[i]) - sin(mu_n*yD[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xL[i] - xR[i])*(sin(mu_n*yU[i]) - sin(mu_n*yD[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xL[i]) - sin(lambda_m*xR[i]))*(yU[i] - yD[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xL[i] - xR[i])*(yU[i] - yD[i]);
			
      /* Imnj are the area integrals over the whole grid array heating area. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xL[j]) - sin(lambda_m*xR[j]))*(sin(mu_n*yU[j]) - sin(mu_n*yD[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xL[j] - xR[j])*(sin(mu_n*yU[j]) - sin(mu_n*yD[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xL[j]) - sin(lambda_m*xR[j]))*(yU[j] - yD[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xL[j] - xR[j])*(yU[j] - yD[j]);

	//Imni = Imni/I00i;
	Imnj = Imnj/I00j;
      /* Insert into analytical expression for thermal impedance. */
      
result += (4.0/((1.0 + delta_m0)*(1.0 + delta_n0)))*(Imni*Imnj/I00i)*tanh(gamma_mn*D)/(gamma_mn*s*Ks);
//      result += (4.0/((1.0 + delta_m0)*(1.0 + delta_n0)))*(Imni*Imnj/I00i)*((H*tanh(gamma_mn*D)+Ks*gamma_mn)/s)/
//	((H*H + gamma_mn*gamma_mn*Ks*Ks)*(1.0/tanh(gamma_mn*D)) + 2.0*H*Ks*gamma_mn);
    }
  }

  result = result*P/(L*W);
  return result;
}




double ThermalHeatsink::TimeIndependent(int i)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imni,I00i,Imnj,I00j;

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body radiation and natural convection to order (T - T0). */

  k = Ks/(rho*C);	/* Diffusivity. */

  P = 1.0/((XL-XR)*(YU-YD)); /* W.m-2 */ 
  /* Areal power dissipation over the whole grid array heating area (total power set at 1 Watt for normalisation). */

  mMAX = nMAX = 35;	/* Number of terms in the m and n summations of the double series for the thermal impedance. */	
	
  result = 0.0;
  for (m = 0; m <= mMAX; m++){
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++){
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n);

      /* Imni are area integrals over the individual heating elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xL[i]) - sin(lambda_m*xR[i]))*(sin(mu_n*yU[i]) - sin(mu_n*yD[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xL[i] - xR[i])*(sin(mu_n*yU[i]) - sin(mu_n*yD[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xL[i]) - sin(lambda_m*xR[i]))*(yU[i] - yD[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xL[i] - xR[i])*(yU[i] - yD[i]);
			
      /* Imnj are the area integrals over the whole grid array heating area. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*XL) - sin(lambda_m*XR))*(sin(mu_n*YU) - sin(mu_n*YD))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (XL - XR)*(sin(mu_n*YU) - sin(mu_n*YD))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*XL) - sin(lambda_m*XR))*(YU - YD)/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (XL - XR)*(YU - YD);

      /* Insert into analytical expression for thermal impedance. */

      if (!(m==0 && n==0)) result += (4.0/((1.0 + delta_m0)*(1.0 + delta_n0)))*(Imni*Imnj/I00i)*(H*tanh(gamma_mn*D)+Ks*gamma_mn)/
	((H*H + gamma_mn*gamma_mn*Ks*Ks)*tanh(gamma_mn*D) + 2.0*H*Ks*gamma_mn);
    }
  }

  result += ((H*D/Ks)+1.0)*(I00j/H)/((H*D/Ks)+2.0);
  result = result*P/(L*W);

  return result;
}


double_complex ThermalHeatsink::LaplaceDomain(int i,  double_complex s)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double P, lambda_m, mu_n, H, T0, k, Imni,I00i,Imnj,I00j;

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body
                                                   radiation and
                                                   natural convection
                                                   to order (T -
                                                   T0). */
  k = Ks/(rho*C);	/* Diffusivity. */

  P = one/((XL-XR)*(YU-YD)); /* W.m-2 */ 
  /* Areal power dissipation over the whole grid array heating area
     (total power set at 1 Watt for normalisation). */

  mMAX = nMAX = 35; /* Number of terms in the m and n summations of
                       the double series for the thermal impedance. */
  double_complex result(zero);
  for (m = 0; m <= mMAX; m++){
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++){
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      double_complex gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n + s/k);

      /* Imni are area integrals over the individual heating elements. */

      if (m != 0 && n != 0) 
	Imni = (sin(lambda_m*xL[i]) - sin(lambda_m*xR[i]))*(sin(mu_n*yU[i]) - sin(mu_n*yD[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) 
	Imni = (xL[i] - xR[i])*(sin(mu_n*yU[i]) - sin(mu_n*yD[i]))/mu_n;
      if (m != 0 && n == 0) 
	Imni = (sin(lambda_m*xL[i]) - sin(lambda_m*xR[i]))*(yU[i] - yD[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xL[i] - xR[i])*(yU[i] - yD[i]);
			
      /* Imnj are the area integrals over the whole grid array heating area. */

      if (m != 0 && n != 0) 
	Imnj = (sin(lambda_m*XL) - sin(lambda_m*XR))*(sin(mu_n*YU) - sin(mu_n*YD))/(lambda_m*mu_n);
      if (m == 0 && n != 0) 
	Imnj = (XL - XR)*(sin(mu_n*YU) - sin(mu_n*YD))/mu_n;
      if (m != 0 && n == 0) 
	Imnj = (sin(lambda_m*XL) - sin(lambda_m*XR))*(YU - YD)/lambda_m;
      if (m == 0 && n == 0) 
	Imnj = I00j = (XL - XR)*(YU - YD);

	    // Calculate exp(-gamma_mn*D)
      double_complex expx = -gamma_mn*D, imag_i(zero,one);
      expx = exp(expx.real())*(cos(expx.imag())+imag_i*sin(expx.imag()));
	

      /* Insert into analytical expression for thermal impedance. */
	      result += (one/((one + delta_m0)*(one + delta_n0)))
	* (Imni*Imnj/I00i)*(H*(one-expx*expx/(one+expx*expx))+Ks*gamma_mn) /
	((H*H + gamma_mn*gamma_mn*Ks*Ks)*(one-expx*expx/(one+expx*expx)) 
	 + 2.0*H*Ks*gamma_mn) ; 
    }
  }
  result = result*4.0*P/(L*W);

  return result;
}



/* Function TimeDomain performs the inverse Laplace transform
   numerically based on Stehfest's algorithm. */
/* In fact the inversion is easy to do analytically, but the numerical
   approach is general and computationally cheap. */
void ThermalHeatsink::TimeDomain(double t, double *Rthtemp)
{
  int v,i,j;
  double result, increment;
	
  for (i = 1; i <= Nsquared; i++){
  for (j = 1; j <= Nsquared; j++){
    increment = result = 0;		
    for (v=1;v<=Np;v++){
      increment = CalculateWeight(v) * LaplaceDomain(i,j,v*LN2/t);
      result = result + increment;
    }			
    result = result*LN2/t;
    Rthtemp[(i-1)*Nsquared+j] = result;

  }
  }
}

/* Function CalculateWeight constructs the weights for the numerical
   Laplace inversion. */
double ThermalHeatsink::CalculateWeight(int v)
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
	
      increment = pow((double)k,Np/2.0) * fact(2*k);
      increment = increment / (fact(Np/2-k) * fact(k) * fact(k-1) * fact(v-k) * fact(2*k-v));
      result = result + increment;
		
    }
  result = result*pow((double)-1,(double)(Np/2+v));
  return result;
}


/* Function factorial returns x! */
double ThermalHeatsink::fact(int x)
{	
  int i;
  double cumul = 1;
  for(i=1;i<=x;i++) cumul = cumul *  i;
  if (x==0) cumul = one;
  if (x < 0) {
    printf("Argument less than zero in factorial\n");
    exit(1);
  }
  return cumul;
}	


/* nrerror prints an error message to stderr on failure. */

void ThermalHeatsink::nrerror(const char *error_text)
{
  fprintf(stderr,"Run time error....\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"... now exiting to system ..\n");
  exit(1);
}

/* Dynamic memory allocation and deallocation routines. */

/* free a double matrix allocated by matrix() */
void ThermalHeatsink::free_matrix(double** m, int nrl, int nrh, int ncl, int nch)
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}


/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double** ThermalHeatsink::newmatrix(int nrl, int nrh, int ncl, int nch)
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

double* ThermalHeatsink::newvector(int nl,int nh)	
{ 
  double *v;
  v=(double *) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("Allocation failure in newvector()"); 
  return v-nl+NR_END;
}

void ThermalHeatsink::free_vector(double *v,int nl,int nh)	
{
  free((char*) (v+nl-NR_END));
}
