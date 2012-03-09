#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "ThermalNportNL.h"
#include <cstdio>

#define Np 6 
// Number of terms used in Stehfest algorithm for
// numerical Laplace inversion. 
#define NR_END 0 /* For use in dynamic allocation of arrays. */

#define FREE_ARG char*

#define LN2 0.6931471806

// Static members
// Number of netlist paramenters (remember to update!)
const unsigned ThermalNportNL::n_par = 18;
// Stefan-Boltzmann constant W.m-2.K-4.
const double ThermalNportNL::sigma = 5.67051e-8;

// Element information
ItemInfo ThermalNportNL::einfo =
{
  "thermalnportnl",
  "N-port grid array substrate with NxN surface heating elements, MxM surface discretisation and non linear boundary conditions",
  "Bill Batty, Carlos Christoffersen mailto:w.batty@elec-eng.leeds.ac.uk",
  DEFAULT_ADDRESS"category:thermal",
  "2001_10_10"
};

// Parameter information
// true means the paramenter is required and false means it
// can be omitted in the netlist.
// See ../network/NetListItem.h for a list of types (TR_INT, TR_DOUBLE, etc.)
ParmInfo ThermalNportNL::pinfo[] =
{
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
  {"ks", "Thermal conductivity (W/m.K).", TR_DOUBLE, false},
  {"rho", "Density  (kg.m-3).", TR_DOUBLE, false},
  {"c", "Specific heat (J/kg.K).", TR_DOUBLE, false},
  {"xi","Adjustment for T**4 non linearity", TR_DOUBLE, false},
  {"eta", "Adjustment for natural convection", TR_DOUBLE, false},
  {"epsilon", "Emissivity", TR_DOUBLE, false},
  {"narray", "Order of NxN grid array", TR_INT, true},
  {"ndevices", "Number of heat dissipating devices", TR_INT, false},
  {"msubstrate", "Order of MxM substrate surface discretisation", TR_INT, false},
  {"b","Exponent in power law temperature dependence of thermal conductivity",TR_DOUBLE,false}
};


ThermalNportNL::ThermalNportNL(const string& iname) : Element(&einfo, pinfo, n_par, iname)
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
  paramvalue[8] = &(Ks = .294);
  paramvalue[9] = &(rho = 1900.);
  paramvalue[10] = &(C = 1150.);
  paramvalue[11] = &(xi = 1.3);
  paramvalue[12] = &(eta = 3.);
  paramvalue[13] = &(epsilon = .7);
  paramvalue[14] = &(Narray);
  paramvalue[15] = &(Ndevices = 1);
  paramvalue[16] = &(Msubstrate = 3);
  paramvalue[17] = &(b = 0.0);

  // Set flags
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);

  // Set some unallocated pointers to 0
  Rth00 = Rth0D = RthD0 = RthDD = 0;
  Rth00temp = Rth0Dtemp = RthD0temp = RthDDtemp = 0;
  Pstoretotal = storedDeltaT = 0;

  Nsquared = 0;
  Msquared = 0;
}


ThermalNportNL::~ThermalNportNL()
{
  if (Rth00)
    free_tensor(Rth00,0,Ntimesteps,0,Nsquared+Msquared,0,Nsquared+Msquared);
  if (Rth0D)
    free_tensor(Rth0D,0,Ntimesteps,0,Nsquared+Msquared,1,Msquared);
  if (RthD0)
    free_tensor(RthD0,0,Ntimesteps,1,Msquared,0,Nsquared+Msquared);
  if (RthDD)
    free_tensor(RthDD,0,Ntimesteps,1,Msquared,1,Msquared);
  if (Rth00temp)
    free_matrix(Rth00temp,0,Nsquared+Msquared,0,Nsquared+Msquared);
  if (Rth0Dtemp)
    free_matrix(Rth0Dtemp,0,Nsquared+Msquared,1,Msquared);
  if (RthD0temp)
    free_matrix(RthD0temp,1,Msquared,0,Nsquared+Msquared);
  if (RthDDtemp)
    free_matrix(RthDDtemp,1,Msquared,1,Msquared);
  if (Pstoretotal)
    free_vector(Pstoretotal,1,Ntimesteps);
  if (storedDeltaT)
    delete [] storedDeltaT;
  if (Nsquared) {
    delete [] xL0;
    delete [] xR0;
    delete [] yU0;
    delete [] yD0;
  }
  if (Msquared) {
    delete [] xLD;
    delete [] xRD;
    delete [] yUD;
    delete [] yDD;
  }

}

// This is the initialization routine.
void ThermalNportNL::init() throw(string&)
{
  FILE *fp;
  double input;

  if (Narray < 1)
    throw(getInstanceName() +
	": Narray must be greater than zero.");
  if (Msubstrate < 1)
    throw(getInstanceName() +
	": Msubstrate must be greater than zero.");
  Nsquared = Narray * Narray;
  Msquared = Msubstrate * Msubstrate;
  // Allocate memory for vectors

  xL0 = new double[Nsquared+Msquared+1];
  xR0 = new double[Nsquared+Msquared+1];
  yU0 = new double[Nsquared+Msquared+1];
  yD0 = new double[Nsquared+Msquared+1];

  xLD = new double[Msquared+1];
  xRD = new double[Msquared+1];
  yUD = new double[Msquared+1];
  yDD = new double[Msquared+1];

  // Set number of terminals and state variables
  if (getNumTerms() != unsigned(Nsquared + 2*Msquared + 1))
    throw(getInstanceName() + ": Incorrect number of terminals.");
  setNumTerms(Nsquared + 2*Msquared + 1);
  // resize vector with MNAM equation row numbers.
  my_row.resize(Nsquared + 2*Msquared);
  // Set number of state variables
  setNumberOfStates(Nsquared + 2*Msquared);

  /* Element coordinates could easily be read from a data file, to
	allow easy input of different grid array layouts. */

  /* Use a single averaged heating element, to represent the effects
	of heat spreading by surface metallisation, between power
	dissipating elements. */

  /* x-coordinate of left edge of total heating area (meters). */
  xL0[0] = L/2.0 + 0.5*L*(Narray - 1.0/3.0)/(Narray+1.0);
  /* x-coordinate of right edge of total heating area (meters). */
  xR0[0] = L/2.0 - 0.5*L*(Narray - 1.0/3.0)/(Narray+1.0);
  /* y-coordinate of upper edge of total heating area (meters). */
  yU0[0] = W/2.0 + 0.5*W*(Narray - 1.0/3.0)/(Narray+1.0);
  /* y-coordinate lower edge of total heating area (meters). */
  yD0[0] = W/2.0 - 0.5*W*(Narray - 1.0/3.0)/(Narray+1.0);

  /* Surface heating element base coordinates */

  for (int m = 1; m <= Narray; m++){
    for (int n = 1; n <= Narray; n++){
      int l = (m-1)*Narray + n;
      /* x-coordinate of left edge of heating element l (meters). */
      xL0[l] =  ( 1.0/3.0 + m)*L/(Narray+1.0);
      /* x-coordinate of right edge of heating element l (meters). */
      xR0[l] =  (-1.0/3.0 + m)*L/(Narray+1.0);
      /* y-coordinate of upper edge of heating element l (meters). */
      yU0[l] =  ( 1.0/3.0 + n)*W/(Narray+1.0);
      /* y-coordinate of lower edge of heating element l (meters). */
      yD0[l] =  (-1.0/3.0 + n)*W/(Narray+1.0);
    }
  }

  /* Discretised surface element coordinates */

  for (int m = 1; m <= Msubstrate; m++){
    for (int n = 1; n <= Msubstrate; n++){
      int l = Nsquared + (m-1)*Msubstrate + n;
      /* x-coordinate of left  edge of z=0 and z=D surface element l (meters). */
      xL0[l] = xLD[l-Nsquared] =  m*L/(Msubstrate*1.0);
      /* x-coordinate of right edge of z=0 and z=D surface element l (meters). */
      xR0[l] = xRD[l-Nsquared] = (m-1)*L/(Msubstrate*1.0);
      /* y-coordinate of upper edge of z=0 and z=D surface element l (meters). */
      yU0[l] = yUD[l-Nsquared] =  n*W/(Msubstrate*1.0);
      /* y-coordinate of lower edge of z=0 and z=D surface element l (meters). */
      yD0[l] = yDD[l-Nsquared] = (n-1)*W/(Msubstrate*1.0);

    }
  }

  if (time_d) {
    // Check parameters
    if (!Ntimesteps || !dt)
      throw(getInstanceName() +
		": Both Ntimesteps and dt must be specified if time_d=1");

    // Change element flags
    setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

    /* Thermal impedance matrices Rth00_ij(t), Rth0D_ij(t), RthD0_ij(t), RthDD_ij(t) are precomputed at all
		required timesteps prior to the iterative coupled electro-thermal solution. */

    Rth00 = newtensor(0,Ntimesteps,0,Nsquared+Msquared,0,Nsquared+Msquared);
    Rth0D = newtensor(0,Ntimesteps,0,Nsquared+Msquared,1,Msquared);
    RthD0 = newtensor(0,Ntimesteps,1,Msquared,0,Nsquared+Msquared);
    RthDD = newtensor(0,Ntimesteps,1,Msquared,1,Msquared);

    /* Dynamically allocate storage for thermal impedance matrices evaluated at specified times t. */

    /* At t = 0 temperature rise is always 0. */

    for (int i = 0; i <= Nsquared+Msquared; i++){
      for (int j = 0; j <= Nsquared+Msquared; j++){
				Rth00[0][i][j] = 0.0;
      }
    }

    for (int i = 0; i <= Nsquared+Msquared; i++){
      for (int j = 1; j <= Msquared; j++){
				Rth0D[0][i][j] = 0.0;
      }
    }

    for (int i = 1; i <= Msquared; i++){
      for (int j = 0; j <= Nsquared+Msquared; j++){
				RthD0[0][i][j] = 0.0;
      }
    }

    for (int i = 1; i <= Msquared; i++){
      for (int j = 1; j <= Msquared; j++){
				RthDD[0][i][j] = 0.0;
      }
    }


    Rth00temp = newmatrix(0,Nsquared+Msquared,0,Nsquared+Msquared);
    Rth0Dtemp = newmatrix(0,Nsquared+Msquared,1,Msquared);
    RthD0temp = newmatrix(1,Msquared,0,Nsquared+Msquared);
    RthDDtemp = newmatrix(1,Msquared,1,Msquared);

    Pstoretotal = newvector(1,Ntimesteps);
    storedDeltaT = new double[Nsquared+2*Msquared];

    //----------- Allow storage and re-use of precomputed Rth(t).
    if (!(read_input))
      fp = fopen("Grid_Rth_elements.dat","w");
    if (read_input)
      fp = fopen("Grid_Rth_elements.dat","r");

    for (int n = 1; n <= Ntimesteps; n++){
      double t = dt * n;
      if (!(read_input) ){
				if (n==1) printf("\nPrecomputation for the Rth(t) ...\n\n");
        printf ("Timestep %d of %d\n",n,Ntimesteps);
				// printf("%g\n",t);
				TimeDomain00(t,Rth00temp); // Calculate the Rth_ij(t) as a function of time t.
				TimeDomain0D(t,Rth0Dtemp);
				TimeDomainD0(t,RthD0temp);
				TimeDomainDD(t,RthDDtemp);

				for (int i = 0; i <= Nsquared+Msquared; i++){
					for (int j = 0; j <= Nsquared+Msquared; j++){
						Rth00[n][i][j] = Rth00temp[i][j];
						//	  printf("%g ",Rth00[n][i][j]);
						fprintf(fp,"%f ",Rth00[n][i][j]);
					}
					//	printf("\n");
					fprintf(fp,"\n");
				}
				//	printf("\n");
				fprintf(fp,"\n");

				for (int i = 0; i <= Nsquared+Msquared; i++){
					for (int j = 1; j <= Msquared; j++){
						Rth0D[n][i][j] = Rth0Dtemp[i][j];
						//	  printf("%g ",Rth0D[n][i][j]);
						fprintf(fp,"%f ",Rth0D[n][i][j]);
					}
					//	printf("\n");
					fprintf(fp,"\n");
				}
				//	printf("\n");
				fprintf(fp,"\n");

				for (int i = 1; i <= Msquared; i++){
					for (int j = 0; j <= Nsquared+Msquared; j++){
						RthD0[n][i][j] = RthD0temp[i][j];
						//	  printf("%g ",RthD0[n][i][j]);
						fprintf(fp,"%f ",RthD0[n][i][j]);
					}
					//	printf("\n");
					fprintf(fp,"\n");
				}
				//	printf("\n");
				fprintf(fp,"\n");


				for (int i = 1; i <= Msquared; i++){
					for (int j = 1; j <= Msquared; j++){
						RthDD[n][i][j] = RthDDtemp[i][j];
						//	  printf("%g ",RthDD[n][i][j]);
						fprintf(fp,"%f ",RthDD[n][i][j]);
					}
					//	printf("\n");
					fprintf(fp,"\n");
				}
				//	printf("\n");
				fprintf(fp,"\n");



      }
      else if (read_input) {
				if (n==1)
					printf("\nRead file input for the Rth(t) ...\n\n");
				for (int i = 0; i <= Nsquared+Msquared; i++){
					for (int j = 0; j <= Nsquared+Msquared; j++){
						fscanf(fp,"%lf ",&input);
						Rth00[n][i][j] = input;
					}
				}
				for (int i = 0; i <= Nsquared+Msquared; i++){
					for (int j = 1; j <= Msquared; j++){
						fscanf(fp,"%lf ",&input);
						Rth0D[n][i][j] = input;
					}
				}

				for (int i = 1; i <= Msquared; i++){
					for (int j = 0; j <= Nsquared+Msquared; j++){
						fscanf(fp,"%lf ",&input);
						RthD0[n][i][j] = input;
					}
				}

				for (int i = 1; i <= Msquared; i++){
					for (int j = 1; j <= Msquared; j++){
						fscanf(fp,"%lf ",&input);
						RthDD[n][i][j] = input;
					}
				}

      }
    }
    printf("Precomputation/read for the Rth_ij(t) complete.\n\n");
    fclose(fp);
  }
}


unsigned ThermalNportNL::getExtraRC(const unsigned& eqn_number,
const MNAMType& type)
{
  // Keep the equation numbers assigned to this element
  for (int i=0; i < Nsquared+2*Msquared; i++)
    my_row[i] = eqn_number + i;

  // Add Nsquared extra RCs
  return Nsquared+2*Msquared;
}

void ThermalNportNL::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  assert(my_row[0]);
  first_eqn = my_row[0];
  n_rows = Nsquared+2*Msquared;
}


void ThermalNportNL::fillMNAM(FreqMNAM* mnam)
{
  /* Should be a safe `small' number, essentially zero. */
  const double SMALL=1.0e0;
  double_complex rth;
  double_complex s(zero, twopi * mnam->getFreq());

  cout << "f = " << mnam->getFreq() << endl;

  if (abs(s) < SMALL) {
    for (int  i = 0; i < Nsquared + Msquared; i++){
      rth = - TimeIndependent00(i+1, 0);
      // Set structural ones
      mnam->setOnes(getTerminal(i)->getRC(),
			getTerminal(Nsquared+2*Msquared)->getRC(),
			my_row[i]);
      // Fill elements is row i of thermal block
      // All the elements in one row are the same
      for (int j=0; j < Nsquared; j++)
				mnam->setElement(my_row[i], my_row[j], rth);

      for (int j = Nsquared; j < Nsquared+Msquared; j++){
				rth = - TimeIndependent00(i+1, j+1);
				mnam->setElement(my_row[i], my_row[j], rth);
      }
      for (int j = 0; j < Msquared; j++){
				rth = - TimeIndependent0D(i+1, j+1);
				mnam->setElement(my_row[i], my_row[j+Nsquared+Msquared], rth);
      }
    }
    for (int  i = Nsquared+Msquared; i < Nsquared + 2*Msquared; i++){
      rth = - TimeIndependentD0(i+1-Nsquared-Msquared, 0);
      // Set structural ones
      mnam->setOnes(getTerminal(i)->getRC(),
			getTerminal(Nsquared+2*Msquared)->getRC(),
			my_row[i]);
      // Fill elements is row i of thermal block
      // All the elements in one row are the same
      for (int j=0; j < Nsquared; j++)
				mnam->setElement(my_row[i], my_row[j], rth);

      for (int j = Nsquared; j < Nsquared+Msquared; j++){
				rth = - TimeIndependentD0(i+1-Nsquared-Msquared, j+1);
				mnam->setElement(my_row[i], my_row[j], rth);
      }
      for (int j = 0; j < Msquared; j++){
				rth = - TimeIndependentDD(i+1-Nsquared-Msquared, j+1);
				mnam->setElement(my_row[i], my_row[j+Nsquared+Msquared], rth);
      }
    }
  }
  else {
    for (int  i = 0; i < Nsquared + Msquared; i++){
      rth = - LaplaceDomain00(i+1, 0, s);
      // Set structural ones
      mnam->setOnes(getTerminal(i)->getRC(),
			getTerminal(Nsquared+2*Msquared)->getRC(),
			my_row[i]);
      // Fill elements is row i of thermal block
      // All the elements in one row are the same
      for (int j=0; j < Nsquared; j++)
				mnam->setElement(my_row[i], my_row[j], rth);

      for (int j = Nsquared; j < Nsquared+Msquared; j++){
				rth = - LaplaceDomain00(i+1, j+1, s);
				mnam->setElement(my_row[i], my_row[j], rth);
      }
      for (int j = 0; j < Msquared; j++){
				rth = - LaplaceDomain0D(i+1, j+1, s);
				mnam->setElement(my_row[i], my_row[j+Nsquared+Msquared], rth);
      }
    }
    for (int  i = Nsquared+Msquared; i < Nsquared + 2*Msquared; i++){
      rth = - LaplaceDomainD0(i+1-Nsquared-Msquared, 0, s);
      // Set structural ones
      mnam->setOnes(getTerminal(i)->getRC(),
			getTerminal(Nsquared+2*Msquared)->getRC(),
			my_row[i]);
      // Fill elements is row i of thermal block
      // All the elements in one row are the same
      for (int j=0; j < Nsquared; j++)
				mnam->setElement(my_row[i], my_row[j], rth);

      for (int j = Nsquared; j < Nsquared+Msquared; j++){
				rth = - LaplaceDomainD0(i+1-Nsquared-Msquared, j+1, s);
				mnam->setElement(my_row[i], my_row[j], rth);
      }
      for (int j = 0; j < Msquared; j++){
				rth = - LaplaceDomainDD(i+1-Nsquared-Msquared, j+1, s);
				mnam->setElement(my_row[i], my_row[j+Nsquared+Msquared], rth);
      }
    }
  }


}

// This routine is called at each time step
// (Actually, several times at each time step)
void ThermalNportNL::svTran(TimeDomainSV* tdsv)
{
  if(tdsv->DC()) {
    for (unsigned i=0; i < getNumberOfStates(); i++) {
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
    for (int i = 0; i < Nsquared; i++){
      tdsv->i(i) = tdsv->getX(i);
      Pstoretotal[n] += tdsv->getX(i);
      /* Only total power dissipation is required for temperature rise
			calculation. */
    }

    // Avoid doing the convolution at each nonlinear iteration.
    if (tdsv->firstTime()) {
      /* Calculate temperature rise of each discretised heating
			element on top surface. */
      for (int  i = 0; i < Nsquared + Msquared; i++){
        storedDeltaT[i] = zero;
				for (int m = n; m > 1; m--){
					storedDeltaT[i] += (Rth00[m][i+1][0] - Rth00[m-1][i+1][0])
					*Pstoretotal[n-m+1];
					for (int j = Nsquared+1; j <= Nsquared+Msquared; j++){
						storedDeltaT[i] += (Rth00[m][i+1][j] - Rth00[m-1][i+1][j])
						*tdsv->getX(j-1);
					}
					for (int j = 1; j <= Msquared; j++){
						storedDeltaT[i] += (Rth0D[m][i+1][j] - Rth0D[m-1][i+1][j])
						*tdsv->getX(Nsquared+Msquared+j-1);
					}
        }
      }
      /* Calculate temperature rise of each discretised heating
			element on top surface. */
      for (int  i = Nsquared+Msquared; i < Nsquared + 2*Msquared; i++){
        storedDeltaT[i] = zero;
				for (int m = n; m > 1; m--){
					storedDeltaT[i] += (RthD0[m][i+1-Nsquared-Msquared][0]
					- RthD0[m-1][i+1-Nsquared-Msquared][0])
					*Pstoretotal[n-m+1];
					for (int j = Nsquared+1; j <= Nsquared+Msquared; j++){
						storedDeltaT[i] += (RthD0[m][i+1-Nsquared-Msquared][j]
						- RthD0[m-1][i+1-Nsquared-Msquared][j])
						*tdsv->getX(j-1);
					}
					for (int j = 1; j <= Msquared; j++){
						storedDeltaT[i] += (RthDD[m][i+1-Nsquared-Msquared][j]
						- RthDD[m-1][i+1-Nsquared-Msquared][j])
						*tdsv->getX(Nsquared+Msquared+j-1);
          }
        }
      }
    }

    double deltaTi;
    double& T0 = Tambient;		/* Ambient temperature. */
    /* Black body radiation and natural convection to order (T - T0). */
    double H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;

    /* Calculate temperature rise of each discretised heating element
		on top surface. */
    for (int  i = 0; i < Nsquared + Msquared; i++){
      deltaTi = storedDeltaT[i] + Rth00[1][i+1][0]*Pstoretotal[n];
      for (int j = Nsquared+1; j <= Nsquared+Msquared; j++){
				deltaTi += Rth00[1][i+1][j]*tdsv->getX(j-1);
      }
      for (int j = 1; j <= Msquared; j++){
				deltaTi += Rth0D[1][i+1][j]*tdsv->getX(Nsquared+Msquared+j-1);
      }
      tdsv->i(i) = tdsv->getX(i) - /**********/ 0.0002777777*H*deltaTi;
      //Implement inverse Kirchhoff transformation
      deltaTi = Tambient*(pow(one-(b-one)*deltaTi/Tambient,-one/(b-one))-one);
      tdsv->u(i) = deltaTi;
    }
    /* Calculate temperature rise of each discretised heating element
		on top surface. */
    for (int  i = Nsquared+Msquared; i < Nsquared + 2*Msquared; i++){
      deltaTi = storedDeltaT[i]
			+ RthD0[1][i+1-Nsquared-Msquared][0]*Pstoretotal[n];
      for (int j = Nsquared+1; j <= Nsquared+Msquared; j++){
				deltaTi += RthD0[1][i+1-Nsquared-Msquared][j]*tdsv->getX(j-1);
      }
      for (int j = 1; j <= Msquared; j++){
        deltaTi += RthDD[1][i+1-Nsquared-Msquared][j]
				*tdsv->getX(Nsquared+Msquared+j-1);
      }
      tdsv->i(i) = - tdsv->getX(i) - /***********/ 0.0002777777*H*deltaTi;
      //Implement inverse Kirchhoff transformation
      deltaTi = Tambient*(pow(one-(b-one)*deltaTi/Tambient,-one/(b-one))-one);
      tdsv->u(i) = deltaTi;
    }

    /* Given Ptotal at specified timestep, the values of the deltaTi
		are the final results. */
    /* deltaTi is given in Kelvin (and t in seconds). */
  }
}

void ThermalNportNL::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->cleanJac();
  // Ji is an identity matrix
  if(tdsv->DC()) {
    for (unsigned i = 0; i <  getNumberOfStates(); i++)
      tdsv->getJi()(i,i) = one;
  }
  else {
    double pst = zero;
    for (int i = 0; i < Nsquared; i++)
      pst += tdsv->getX(i);

    assert(!tdsv->firstTime());
    double& T0 = Tambient;		/* Ambient temperature. */
    /* Black body radiation and natural convection to order (T - T0). */
    double H = /***********/
		0.0002777777*xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;

    // Build J_theta and J_i
    // Elements from Rth00
    for (int i=0; i< Nsquared+Msquared; i++)
      for (int j=0; j< Nsquared; j++)
				tdsv->getJu()(i,j) = Rth00[1][i+1][0];
    for (int i=0; i< Nsquared+Msquared; i++)
		for (int j=Nsquared; j< Nsquared+Msquared; j++) {
			tdsv->getJu()(i,j) = Rth00[1][i+1][j+1];
			tdsv->getJi()(i,j) = - H * tdsv->getJu()(i,j);
		}
    // Elements from RthD0
    for (int i=Nsquared+Msquared; i< Nsquared+2*Msquared; i++)
      for (int j=0; j< Nsquared; j++)
				tdsv->getJu()(i,j) = RthD0[1][i+1-Nsquared-Msquared][0];
    for (int i=Nsquared+Msquared; i< Nsquared+2*Msquared; i++)
		for (int j=Nsquared; j< Nsquared+Msquared; j++) {
			tdsv->getJu()(i,j) = RthD0[1][i+1-Nsquared-Msquared][j+1];
			tdsv->getJi()(i,j) = - H * tdsv->getJu()(i,j);
		}
    for (int i=Nsquared; i< Nsquared+2*Msquared; i++)
      for (int j=0; j< Nsquared; j++)
				tdsv->getJi()(i,j) = - H * tdsv->getJu()(i,j);
    // Elements from Rth0D
    for (int i=0; i< Nsquared+Msquared; i++)
		for (int j=Nsquared+Msquared; j< Nsquared+2*Msquared; j++) {
			tdsv->getJu()(i,j) = Rth0D[1][i+1][j+1-Nsquared-Msquared];
			tdsv->getJi()(i,j) = - H * tdsv->getJu()(i,j);
		}
    // Elements from RthDD
    for (int i=Nsquared+Msquared; i< Nsquared+2*Msquared; i++)
		for (int j=Nsquared+Msquared; j< Nsquared+2*Msquared; j++) {
			tdsv->getJu()(i,j) =
			RthDD[1][i+1-Nsquared-Msquared][j+1-Nsquared-Msquared];
			tdsv->getJi()(i,j) = - H * tdsv->getJu()(i,j);
		}
    // Add identity matrix to J_i
    for (int i = 0; i < Nsquared+Msquared; i++)
      tdsv->getJi()(i,i) += one;
    for (int i = Nsquared+Msquared; i < Nsquared+2*Msquared; i++)
      tdsv->getJi()(i,i) -= one;


    /* Calculate temperature rise of each discretised heating element
		on top surface. */
    double deltaTi;
    for (int  i = 0; i < Nsquared + Msquared; i++){
      deltaTi = storedDeltaT[i] + Rth00[1][i+1][0]*pst;
      for (int j = Nsquared+1; j <= Nsquared+Msquared; j++){
				deltaTi += Rth00[1][i+1][j]*tdsv->getX(j-1);
      }
      for (int j = 1; j <= Msquared; j++){
				deltaTi += Rth0D[1][i+1][j]*tdsv->getX(Nsquared+Msquared+j-1);
      }
      double dJi = pow(one-(b-one)/Tambient*deltaTi, b/(one-b));
      for (int  j = 0; j < Nsquared+2*Msquared; j++) {
				tdsv->getJu()(i,j) *= dJi;
      }
    }
    /* Calculate temperature rise of each discretised heating element
		on top surface. */
    for (int  i = Nsquared+Msquared; i < Nsquared + 2*Msquared; i++){
      deltaTi = storedDeltaT[i] + RthD0[1][i+1-Nsquared-Msquared][0]*pst;
      for (int j = Nsquared+1; j <= Nsquared+Msquared; j++){
				deltaTi += RthD0[1][i+1-Nsquared-Msquared][j]*tdsv->getX(j-1);
      }
      for (int j = 1; j <= Msquared; j++){
        deltaTi += RthDD[1][i+1-Nsquared-Msquared][j]*tdsv->getX(Nsquared+Msquared+j-1);
      }
      double dJi = pow(one-(b-one)/Tambient*deltaTi, b/(one-b));
      for (int  j = 0; j < Nsquared+2*Msquared; j++) {
				tdsv->getJu()(i,j) *= dJi;
      }
    }
  }
}


/* Function LaplaceDomain returns the Laplace s-space thermal
impedance matrix elements. */
/* This subroutine contains all the geometrical and material details
required for analytical construction of the s-space thermal
impedance matrix. */
double ThermalNportNL::LaplaceDomain00(int i, int j, double s)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imni,I00i,Imnj,I00j;

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body radiation and natural convection to order (T - T0). */

  k = Ks/(rho*C);	/* Diffusivity. */

  P = 1.0/((xL0[0]-xR0[0])*(yU0[0]-yD0[0])); /* W.m-2 */
  /* Areal power dissipation over the whole grid array heating area (total power set at 1 Watt for normalisation). */

  mMAX = nMAX = 35;	/* Number of terms in the m and n summations of the double series for the thermal impedance. */

  result = 0.0;
  for (m = 0; m <= mMAX; m++){
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++){
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n + s/k);

      /* Imni are area integrals over z=0 temperature response elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xL0[i]) - sin(lambda_m*xR0[i]))*(sin(mu_n*yU0[i]) - sin(mu_n*yD0[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xL0[i] - xR0[i])*(sin(mu_n*yU0[i]) - sin(mu_n*yD0[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xL0[i]) - sin(lambda_m*xR0[i]))*(yU0[i] - yD0[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xL0[i] - xR0[i])*(yU0[i] - yD0[i]);

      /* Imnj are the area integrals z=0 heating elements. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xL0[j]) - sin(lambda_m*xR0[j]))*(sin(mu_n*yU0[j]) - sin(mu_n*yD0[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xL0[j] - xR0[j])*(sin(mu_n*yU0[j]) - sin(mu_n*yD0[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xL0[j]) - sin(lambda_m*xR0[j]))*(yU0[j] - yD0[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xL0[j] - xR0[j])*(yU0[j] - yD0[j]);

      /* Insert into analytical expression for thermal impedance. */

      result += (4.0/((1.0 + delta_m0)*(1.0 + delta_n0)))*(Imni*Imnj/I00i)*((H*tanh(gamma_mn*D)+Ks*gamma_mn)/s)/
			((H*H + gamma_mn*gamma_mn*Ks*Ks)*tanh(gamma_mn*D) + 2.0*H*Ks*gamma_mn);
    }
  }

  result = result*P/(L*W);

  return result;
}


double ThermalNportNL::TimeIndependent00(int i, int j)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imni,I00i,Imnj,I00j;

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body radiation and natural convection to order (T - T0). */

  k = Ks/(rho*C);	/* Diffusivity. */

  P = 1.0/((xL0[0]-xR0[0])*(yU0[0]-yD0[0])); /* W.m-2 */
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

      /* Imni are area integrals over z=0 temperature response elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xL0[i]) - sin(lambda_m*xR0[i]))*(sin(mu_n*yU0[i]) - sin(mu_n*yD0[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xL0[i] - xR0[i])*(sin(mu_n*yU0[i]) - sin(mu_n*yD0[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xL0[i]) - sin(lambda_m*xR0[i]))*(yU0[i] - yD0[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xL0[i] - xR0[i])*(yU0[i] - yD0[i]);

      /* Imnj are the area integrals z=0 heating elements. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xL0[j]) - sin(lambda_m*xR0[j]))*(sin(mu_n*yU0[j]) - sin(mu_n*yD0[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xL0[j] - xR0[j])*(sin(mu_n*yU0[j]) - sin(mu_n*yD0[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xL0[j]) - sin(lambda_m*xR0[j]))*(yU0[j] - yD0[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xL0[j] - xR0[j])*(yU0[j] - yD0[j]);

      /* Insert into analytical expression for thermal impedance. */

      if (!(m==0 && n==0)) result += (4.0/((1.0 + delta_m0)*(1.0 + delta_n0)))*(Imni*Imnj/I00i)*(H*tanh(gamma_mn*D)+Ks*gamma_mn)/
				((H*H + gamma_mn*gamma_mn*Ks*Ks)*tanh(gamma_mn*D) + 2.0*H*Ks*gamma_mn);
    }
  }

  result += ((H*D/Ks)+1.0)*(I00j/H)/((H*D/Ks)+2.0);
  result = result*P/(L*W);

  return result;
}

double_complex ThermalNportNL::LaplaceDomain00(int i, int j, double_complex s)
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

  P = one/((xL0[0]-xR0[0])*(yU0[0]-yD0[0])); /* W.m-2 */
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

      /* Imni are area integrals over z=0 temperature response elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xL0[i]) - sin(lambda_m*xR0[i]))*(sin(mu_n*yU0[i]) - sin(mu_n*yD0[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xL0[i] - xR0[i])*(sin(mu_n*yU0[i]) - sin(mu_n*yD0[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xL0[i]) - sin(lambda_m*xR0[i]))*(yU0[i] - yD0[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xL0[i] - xR0[i])*(yU0[i] - yD0[i]);

      /* Imnj are the area integrals z=0 heating elements. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xL0[j]) - sin(lambda_m*xR0[j]))*(sin(mu_n*yU0[j]) - sin(mu_n*yD0[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xL0[j] - xR0[j])*(sin(mu_n*yU0[j]) - sin(mu_n*yD0[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xL0[j]) - sin(lambda_m*xR0[j]))*(yU0[j] - yD0[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xL0[j] - xR0[j])*(yU0[j] - yD0[j]);

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

double ThermalNportNL::LaplaceDomain0D(int i, int j, double s)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imni,I00i,Imnj,I00j;

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body radiation and natural convection to order (T - T0). */

  k = Ks/(rho*C);	/* Diffusivity. */

  P = 1.0/((xL0[0]-xR0[0])*(yU0[0]-yD0[0])); /* W.m-2 */
  /* Areal power dissipation over the whole grid array heating area (total power set at 1 Watt for normalisation). */

  mMAX = nMAX = 35;	/* Number of terms in the m and n summations of the double series for the thermal impedance. */

  result = 0.0;
  for (m = 0; m <= mMAX; m++){
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++){
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n + s/k);

      /* Imni are area integrals over z=0 temperature response elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xL0[i]) - sin(lambda_m*xR0[i]))*(sin(mu_n*yU0[i]) - sin(mu_n*yD0[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xL0[i] - xR0[i])*(sin(mu_n*yU0[i]) - sin(mu_n*yD0[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xL0[i]) - sin(lambda_m*xR0[i]))*(yU0[i] - yD0[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xL0[i] - xR0[i])*(yU0[i] - yD0[i]);

      /* Imnj are the area integrals z=D heating elements. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xLD[j]) - sin(lambda_m*xRD[j]))*(sin(mu_n*yUD[j]) - sin(mu_n*yDD[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xLD[j] - xRD[j])*(sin(mu_n*yUD[j]) - sin(mu_n*yDD[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xLD[j]) - sin(lambda_m*xRD[j]))*(yUD[j] - yDD[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xLD[j] - xRD[j])*(yUD[j] - yDD[j]);

      /* Insert into analytical expression for thermal impedance. */

      result += (4.0/((1.0 + delta_m0)*(1.0 + delta_n0)))*(Imni*Imnj/I00i)*((Ks*gamma_mn)/s)/
			((H*H + gamma_mn*gamma_mn*Ks*Ks)*sinh(gamma_mn*D) + 2.0*H*Ks*gamma_mn*cosh(gamma_mn*D));
    }
  }

  result = -result*P/(L*W);
  return result;
}

double ThermalNportNL::TimeIndependent0D(int i, int j)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imni,I00i,Imnj,I00j;

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body radiation and natural convection to order (T - T0). */

  k = Ks/(rho*C);	/* Diffusivity. */

  P = 1.0/((xL0[0]-xR0[0])*(yU0[0]-yD0[0])); /* W.m-2 */
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

      /* Imni are area integrals over z=0 temperature response elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xL0[i]) - sin(lambda_m*xR0[i]))*(sin(mu_n*yU0[i]) - sin(mu_n*yD0[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xL0[i] - xR0[i])*(sin(mu_n*yU0[i]) - sin(mu_n*yD0[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xL0[i]) - sin(lambda_m*xR0[i]))*(yU0[i] - yD0[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xL0[i] - xR0[i])*(yU0[i] - yD0[i]);

      /* Imnj are the area integrals z=D heating elements. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xLD[j]) - sin(lambda_m*xRD[j]))*(sin(mu_n*yUD[j]) - sin(mu_n*yDD[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xLD[j] - xRD[j])*(sin(mu_n*yUD[j]) - sin(mu_n*yDD[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xLD[j]) - sin(lambda_m*xRD[j]))*(yUD[j] - yDD[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xLD[j] - xRD[j])*(yUD[j] - yDD[j]);

      /* Insert into analytical expression for thermal impedance. */

      if (!(m==0 && n==0)) result += (4.0/((1.0 + delta_m0)*(1.0 + delta_n0)))*(Imni*Imnj/I00i)*(Ks*gamma_mn)/
				((H*H + gamma_mn*gamma_mn*Ks*Ks)*cosh(gamma_mn*D) + 2.0*H*Ks*gamma_mn*sinh(gamma_mn*D));
    }
  }

  result += (I00j/H)/((H*D/Ks)+2.0);
  result = -result*P/(L*W);

  return result;
}


double_complex ThermalNportNL::LaplaceDomain0D(int i, int j, double_complex s)
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

  P = one/((xL0[0]-xR0[0])*(yU0[0]-yD0[0])); /* W.m-2 */
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

      /* Imni are area integrals over z=0 temperature response elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xL0[i]) - sin(lambda_m*xR0[i]))*(sin(mu_n*yU0[i]) - sin(mu_n*yD0[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xL0[i] - xR0[i])*(sin(mu_n*yU0[i]) - sin(mu_n*yD0[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xL0[i]) - sin(lambda_m*xR0[i]))*(yU0[i] - yD0[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xL0[i] - xR0[i])*(yU0[i] - yD0[i]);

      /* Imnj are the area integrals z=D heating elements. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xLD[j]) - sin(lambda_m*xRD[j]))*(sin(mu_n*yUD[j]) - sin(mu_n*yDD[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xLD[j] - xRD[j])*(sin(mu_n*yUD[j]) - sin(mu_n*yDD[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xLD[j]) - sin(lambda_m*xRD[j]))*(yUD[j] - yDD[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xLD[j] - xRD[j])*(yUD[j] - yDD[j]);

      // Calculate exp(-gamma_mn*D)
      double_complex expx = -gamma_mn*D, imag_i(zero,one);
      expx = exp(expx.real())*(cos(expx.imag())+imag_i*sin(expx.imag()));


      /* Insert into analytical expression for thermal impedance. */
      result += (one/((one + delta_m0)*(one + delta_n0)))
			* (Imni*Imnj/I00i)*(Ks*gamma_mn) /
			((H*H + gamma_mn*gamma_mn*Ks*Ks)*(expx-one/expx)/2.0
			+ 2.0*H*Ks*gamma_mn*((expx+one/expx)/2.0)) ;
    }
  }
  result = result*4.0*P/(L*W);

  return result;
}

double ThermalNportNL::LaplaceDomainD0(int i, int j, double s)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imni,I00i,Imnj,I00j;

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body radiation and natural convection to order (T - T0). */

  k = Ks/(rho*C);	/* Diffusivity. */

  P = 1.0/((xL0[0]-xR0[0])*(yU0[0]-yD0[0])); /* W.m-2 */
  /* Areal power dissipation over the whole grid array heating area (total power set at 1 Watt for normalisation). */

  mMAX = nMAX = 35;	/* Number of terms in the m and n summations of the double series for the thermal impedance. */

  result = 0.0;
  for (m = 0; m <= mMAX; m++){
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++){
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n + s/k);

      /* Imni are area integrals over z=D temperature response elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xLD[i]) - sin(lambda_m*xRD[i]))*(sin(mu_n*yUD[i]) - sin(mu_n*yDD[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xLD[i] - xRD[i])*(sin(mu_n*yUD[i]) - sin(mu_n*yDD[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xLD[i]) - sin(lambda_m*xRD[i]))*(yUD[i] - yDD[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xLD[i] - xRD[i])*(yUD[i] - yDD[i]);

      /* Imnj are the area integrals z=0 heating elements. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xL0[j]) - sin(lambda_m*xR0[j]))*(sin(mu_n*yU0[j]) - sin(mu_n*yD0[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xL0[j] - xR0[j])*(sin(mu_n*yU0[j]) - sin(mu_n*yD0[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xL0[j]) - sin(lambda_m*xR0[j]))*(yU0[j] - yD0[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xL0[j] - xR0[j])*(yU0[j] - yD0[j]);

      /* Insert into analytical expression for thermal impedance. */

      result += (4.0/((1.0 + delta_m0)*(1.0 + delta_n0)))*(Imni*Imnj/I00i)*((Ks*gamma_mn)/s)/
			((H*H + gamma_mn*gamma_mn*Ks*Ks)*sinh(gamma_mn*D) + 2.0*H*Ks*gamma_mn*cosh(gamma_mn*D));
    }
  }

  result = result*P/(L*W);

  return result;
}

double ThermalNportNL::TimeIndependentD0(int i, int j)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imni,I00i,Imnj,I00j;

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body radiation and natural convection to order (T - T0). */

  k = Ks/(rho*C);	/* Diffusivity. */

  P = 1.0/((xL0[0]-xR0[0])*(yU0[0]-yD0[0])); /* W.m-2 */
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

      /* Imni are area integrals over z=D temperature response elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xLD[i]) - sin(lambda_m*xRD[i]))*(sin(mu_n*yUD[i]) - sin(mu_n*yDD[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xLD[i] - xRD[i])*(sin(mu_n*yUD[i]) - sin(mu_n*yDD[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xLD[i]) - sin(lambda_m*xRD[i]))*(yUD[i] - yDD[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xLD[i] - xRD[i])*(yUD[i] - yDD[i]);

      /* Imnj are the area integrals z=0 heating elements. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xL0[j]) - sin(lambda_m*xR0[j]))*(sin(mu_n*yU0[j]) - sin(mu_n*yD0[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xL0[j] - xR0[j])*(sin(mu_n*yU0[j]) - sin(mu_n*yD0[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xL0[j]) - sin(lambda_m*xR0[j]))*(yU0[j] - yD0[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xL0[j] - xR0[j])*(yU0[j] - yD0[j]);

      /* Insert into analytical expression for thermal impedance. */

      if (!(m==0 && n==0)) result += (4.0/((1.0 + delta_m0)*(1.0 + delta_n0)))*(Imni*Imnj/I00i)*(Ks*gamma_mn)/
				((H*H + gamma_mn*gamma_mn*Ks*Ks)*sinh(gamma_mn*D) + 2.0*H*Ks*gamma_mn*cosh(gamma_mn*D));
    }
  }

  result += (I00j/H)/((H*D/Ks)+2.0);
  result = result*P/(L*W);

  return result;
}

double_complex ThermalNportNL::LaplaceDomainD0(int i, int j, double_complex s)
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

  P = one/((xL0[0]-xR0[0])*(yU0[0]-yD0[0])); /* W.m-2 */
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

      /* Imni are area integrals over z=D temperature response elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xLD[i]) - sin(lambda_m*xRD[i]))*(sin(mu_n*yUD[i]) - sin(mu_n*yDD[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xLD[i] - xRD[i])*(sin(mu_n*yUD[i]) - sin(mu_n*yDD[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xLD[i]) - sin(lambda_m*xRD[i]))*(yUD[i] - yDD[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xLD[i] - xRD[i])*(yUD[i] - yDD[i]);

      /* Imnj are the area integrals z=0 heating elements. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xL0[j]) - sin(lambda_m*xR0[j]))*(sin(mu_n*yU0[j]) - sin(mu_n*yD0[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xL0[j] - xR0[j])*(sin(mu_n*yU0[j]) - sin(mu_n*yD0[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xL0[j]) - sin(lambda_m*xR0[j]))*(yU0[j] - yD0[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xL0[j] - xR0[j])*(yU0[j] - yD0[j]);

      // Calculate exp(-gamma_mn*D)
      double_complex expx = -gamma_mn*D, imag_i(zero,one);
      expx = exp(expx.real())*(cos(expx.imag())+imag_i*sin(expx.imag()));

      /* Insert into analytical expression for thermal impedance. */
      result += (one/((one + delta_m0)*(one + delta_n0)))
			* (Imni*Imnj/I00i)*(Ks*gamma_mn) /
			((H*H + gamma_mn*gamma_mn*Ks*Ks)*((expx-one/expx)/2.0)
			+ 2.0*H*Ks*gamma_mn*((expx+one/expx)/2.0)) ;
    }
  }
  result = result*4.0*P/(L*W);

  return result;
}

double ThermalNportNL::LaplaceDomainDD(int i, int j, double s)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imni,I00i,Imnj,I00j;

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body radiation and natural convection to order (T - T0). */

  k = Ks/(rho*C);	/* Diffusivity. */

  P = 1.0/((xL0[0]-xR0[0])*(yU0[0]-yD0[0])); /* W.m-2 */
  /* Areal power dissipation over the whole grid array heating area (total power set at 1 Watt for normalisation). */

  mMAX = nMAX = 35;	/* Number of terms in the m and n summations of the double series for the thermal impedance. */

  result = 0.0;
  for (m = 0; m <= mMAX; m++){
    if (m == 0) delta_m0 = 1; else delta_m0 = 0;
    lambda_m = m*pi/L;
    for (n = 0; n <= nMAX; n++){
      if (n == 0) delta_n0 = 1; else delta_n0 = 0;
      mu_n = n*pi/W;
      gamma_mn = sqrt(lambda_m*lambda_m + mu_n*mu_n + s/k);

      /* Imni are area integrals over z=D temperature response elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xLD[i]) - sin(lambda_m*xRD[i]))*(sin(mu_n*yUD[i]) - sin(mu_n*yDD[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xLD[i] - xRD[i])*(sin(mu_n*yUD[i]) - sin(mu_n*yDD[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xLD[i]) - sin(lambda_m*xRD[i]))*(yUD[i] - yDD[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xLD[i] - xRD[i])*(yUD[i] - yDD[i]);

      /* Imnj are the area integrals z=D heating elements. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xLD[j]) - sin(lambda_m*xRD[j]))*(sin(mu_n*yUD[j]) - sin(mu_n*yDD[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xLD[j] - xRD[j])*(sin(mu_n*yUD[j]) - sin(mu_n*yDD[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xLD[j]) - sin(lambda_m*xRD[j]))*(yUD[j] - yDD[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xLD[j] - xRD[j])*(yUD[j] - yDD[j]);

      /* Insert into analytical expression for thermal impedance. */

      result += (4.0/((1.0 + delta_m0)*(1.0 + delta_n0)))*(Imni*Imnj/I00i)*((H*tanh(gamma_mn*D)+Ks*gamma_mn)/s)/
			((H*H + gamma_mn*gamma_mn*Ks*Ks)*tanh(gamma_mn*D) + 2.0*H*Ks*gamma_mn);
    }
  }

  result = -result*P/(L*W);

  return result;
}

double ThermalNportNL::TimeIndependentDD(int i, int j)
{
  int m,n,delta_m0,delta_n0,mMAX,nMAX;
  double result, P, lambda_m,mu_n,gamma_mn, H, T0, k, Imni,I00i,Imnj,I00j;

  T0 = Tambient;		/* Ambient temperature. */
  H = xi*(1.0+eta)*4.0*T0*T0*T0*sigma*epsilon;	/* Black body radiation and natural convection to order (T - T0). */

  k = Ks/(rho*C);	/* Diffusivity. */

  P = 1.0/((xL0[0]-xR0[0])*(yU0[0]-yD0[0])); /* W.m-2 */
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

      /* Imni are area integrals over z=D temperature response elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xLD[i]) - sin(lambda_m*xRD[i]))*(sin(mu_n*yUD[i]) - sin(mu_n*yDD[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xLD[i] - xRD[i])*(sin(mu_n*yUD[i]) - sin(mu_n*yDD[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xLD[i]) - sin(lambda_m*xRD[i]))*(yUD[i] - yDD[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xLD[i] - xRD[i])*(yUD[i] - yDD[i]);

      /* Imnj are the area integrals z=D heating elements. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xLD[j]) - sin(lambda_m*xRD[j]))*(sin(mu_n*yUD[j]) - sin(mu_n*yDD[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xLD[j] - xRD[j])*(sin(mu_n*yUD[j]) - sin(mu_n*yDD[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xLD[j]) - sin(lambda_m*xRD[j]))*(yUD[j] - yDD[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xLD[j] - xRD[j])*(yUD[j] - yDD[j]);

      /* Insert into analytical expression for thermal impedance. */

      if (!(m==0 && n==0)) result += (4.0/((1.0 + delta_m0)*(1.0 + delta_n0)))*(Imni*Imnj/I00i)*(H*tanh(gamma_mn*D)+Ks*gamma_mn)/
				((H*H + gamma_mn*gamma_mn*Ks*Ks)*tanh(gamma_mn*D) + 2.0*H*Ks*gamma_mn);
    }
  }

  result += ((H*D/Ks)+1.0)*(I00j/H)/((H*D/Ks)+2.0);
  result = -result*P/(L*W);

  return result;
}

double_complex ThermalNportNL::LaplaceDomainDD(int i, int j, double_complex s)
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

  P = one/((xL0[0]-xR0[0])*(yU0[0]-yD0[0])); /* W.m-2 */
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

      /* Imni are area integrals over z=D temperature response elements. */

      if (m != 0 && n != 0) Imni = (sin(lambda_m*xLD[i]) - sin(lambda_m*xRD[i]))*(sin(mu_n*yUD[i]) - sin(mu_n*yDD[i]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imni = (xLD[i] - xRD[i])*(sin(mu_n*yUD[i]) - sin(mu_n*yDD[i]))/mu_n;
      if (m != 0 && n == 0) Imni = (sin(lambda_m*xLD[i]) - sin(lambda_m*xRD[i]))*(yUD[i] - yDD[i])/lambda_m;
      if (m == 0 && n == 0) Imni = I00i = (xLD[i] - xRD[i])*(yUD[i] - yDD[i]);

      /* Imnj are the area integrals z=D heating elements. */

      if (m != 0 && n != 0) Imnj = (sin(lambda_m*xLD[j]) - sin(lambda_m*xRD[j]))*(sin(mu_n*yUD[j]) - sin(mu_n*yDD[j]))/(lambda_m*mu_n);
      if (m == 0 && n != 0) Imnj = (xLD[j] - xRD[j])*(sin(mu_n*yUD[j]) - sin(mu_n*yDD[j]))/mu_n;
      if (m != 0 && n == 0) Imnj = (sin(lambda_m*xLD[j]) - sin(lambda_m*xRD[j]))*(yUD[j] - yDD[j])/lambda_m;
      if (m == 0 && n == 0) Imnj = I00j = (xLD[j] - xRD[j])*(yUD[j] - yDD[j]);

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
  result = -result*4.0*P/(L*W);

  return result;
}

/* Function TimeDomain performs the inverse Laplace transform
numerically based on Stehfest's algorithm. */
/* In fact the inversion is easy to do analytically, but the numerical
approach is general and computationally cheap. */
void ThermalNportNL::TimeDomain00(double t, double **Rth00temp)
{
  int v,i,j;
  double result, increment;

  for (i = 0; i <= Nsquared+Msquared; i++){
    for (j = 0; j <= Nsquared+Msquared; j++){
      increment = result = 0;
      for (v=1;v<=Np;v++){
				increment = CalculateWeight(v) * LaplaceDomain00(i,j,v*LN2/t);
				result = result + increment;
      }
      result = result*LN2/t;
      Rth00temp[i][j] = result;
    }
  }
}

void ThermalNportNL::TimeDomain0D(double t, double **Rth0Dtemp)
{
  int v,i,j;
  double result, increment;
  for (i = 0; i <= Nsquared+Msquared; i++){
    for (j = 1; j <= Msquared; j++){
      increment = result = 0;
      for (v=1;v<=Np;v++){
				increment = CalculateWeight(v) * LaplaceDomain0D(i,j,v*LN2/t);
				result = result + increment;
      }
      result = result*LN2/t;
      Rth0Dtemp[i][j] = result;
    }
  }
}

void ThermalNportNL::TimeDomainD0(double t, double **RthD0temp)
{
  int v,i,j;
  double result, increment;

  for (i = 1; i <= Msquared; i++){
    for (j = 0; j <= Nsquared+Msquared; j++){
      increment = result = 0;
      for (v=1;v<=Np;v++){
				increment = CalculateWeight(v) * LaplaceDomainD0(i,j,v*LN2/t);
				result = result + increment;
      }
      result = result*LN2/t;
      RthD0temp[i][j] = result;
    }
  }
}

void ThermalNportNL::TimeDomainDD(double t, double **RthD0temp)
{
  int v,i,j;
  double result, increment;

  for (i = 1; i <= Msquared; i++){
    for (j = 1; j <= Msquared; j++){
      increment = result = 0;
      for (v=1;v<=Np;v++){
				increment = CalculateWeight(v) * LaplaceDomainDD(i,j,v*LN2/t);
				result = result + increment;
      }
      result = result*LN2/t;
      RthDDtemp[i][j] = result;
    }
  }
}

/* Function CalculateWeight constructs the weights for the numerical
Laplace inversion. */
double ThermalNportNL::CalculateWeight(int v)
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
double ThermalNportNL::fact(int x)
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
void ThermalNportNL::nrerror(const char *error_text)
{
  fprintf(stderr,"Run time error....\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"... now exiting to system ..\n");
  exit(1);
}

/* Dynamic memory allocation and deallocation routines. */
/* free a double matrix allocated by matrix() */
void ThermalNportNL::free_matrix(double** m, int nrl, int nrh, int ncl, int nch)
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double** ThermalNportNL::newmatrix(int nrl, int nrh, int ncl, int nch)
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

double* ThermalNportNL::newvector(int nl,int nh)
{
  double *v;
  v=(double *) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("Allocation failure in newvector()");
  return v-nl+NR_END;
}

void ThermalNportNL::free_vector(double *v,int nl,int nh)
{
  free((char*) (v+nl-NR_END));
}


double*** ThermalNportNL::newtensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;

  /* allocate pointers to pointers to rows */
  t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
  if (!t) nrerror("allocation failure 1 in newtensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
  if (!t[nrl]) nrerror("allocation failure 2 in newtensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in newtensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void ThermalNportNL::free_tensor(double ***t, int nrl, int nrh,
int ncl, int nch, int ndl, int ndh)
/* free a double d3tensor allocated by newtensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}

