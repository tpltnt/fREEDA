#include "SVTran.h"
using std::cout;
using std::endl;


// This is required for the allocTimeV_P function
extern "C"
{
#include "../inout/ftvec.h"
#include "../inout/report.h"
}

#define MAX(X,Y) ((X) > (Y) ?  (X) : (Y))
#define MIN(X,Y)  ((X) < (Y) ?  (X) : (Y))


// This function develops the Incidence Matrix T. It is defined in Tincidence.cc
// It also returns elem_vec, number of state variables and maximum number of state
// variables 

void buildTIncidence(ElemFlag mask, Circuit*& my_circuit,
		IntDenseMatrix& T, ElementVector& elem_vec,
		int& n_states, int& max_n_states);

// number of analysis parameters, a static variable
const unsigned SVTran::n_par = 15;


// Element information
ItemInfo SVTran::ainfo =
{
	"SVTran",
	"State-Variable-Based Time-Marching Variable Time Step Transient Analysis with Newton Iterations",
	"Shivam Priyadarshi, Chris Saunders, Eric Wyers, Sonali Luniya, Carlos E. Christoffersen",
	DEFAULT_ADDRESS
};



// Parameter description and types
ParmInfo SVTran::pinfo[] =
{
	{"tstop", "Stop time (s)", TR_DOUBLE, true},
	{"tstep", "Time step (s)", TR_DOUBLE, true},
	{"nst", "No save time (s)", TR_DOUBLE, false},
	{"deriv", "Approximate derivatives or use automatic diff.", TR_INT, false},
	{"im", "Integration method", TR_INT, false},  // im=1 Trapezoidal else Backward Euler
	{"out_steps", "Number of steps skipped for output simulation progress",
		TR_INT, false},
	{"gcomp", "Compensation network conductance (S)", TR_DOUBLE, false},
	{"opt", "run a DC analysis up front", TR_BOOLEAN, false},
	{"reltol", "Relative Error" , TR_DOUBLE, false},
	{"abstol", "Absolute Error" , TR_DOUBLE, false},
	{"tsteptol", "Time step tolerance parameter" , TR_DOUBLE, false}, // Used in calculating upper bound of time step
	{"tm", "Time step method", TR_INT, false},  // if tm=2 it takes precedence over im and follow Backward Euler->Trapezoidal else follow Integration method(im)
	{"max_iter" , "Maximum number of iterations ", TR_INT, false},
	{"tsmin" , "Minimum time step", TR_DOUBLE , false},
	{"tsmax", "Maximum time step", TR_DOUBLE, false}
};




// note some parameters are obtained from the netlist only
SVTran::SVTran() : Analysis(&ainfo, pinfo, n_par), ls_size(0)
{
	// Parameter stuff
	paramvalue[0] = &(tf);
	paramvalue[1] = &(h);
	paramvalue[2] = &(nst = zero);
	paramvalue[3] = &(deriv = 0);
	paramvalue[4] = &(int_method = 1);
	paramvalue[5] = &(out_steps = 200);
	paramvalue[6] = &(gcomp = 1e-5);
	paramvalue[7] = &(opt = false);
	paramvalue[8] = &(reltol=1e-2);
	paramvalue[9] = &(abstol=1e-4);
	paramvalue[10] = &(tsteptol=1e-2);
	paramvalue[11] = &(tstep_method=2);
	paramvalue[12] = &(max_iter=1000);
	paramvalue[13] = &(tsmin=0);
	paramvalue[14] = &(tsmax=0);
}


void SVTran::run(Circuit * cir)
{
	this->cir = cir;
	char msg[80];
	int Newton_Raphson_Converge;


	// Build time domain MNAM, adding all linear circuit elements
	ElemFlag mnam_mask(LINEAR);
	TimeMNAM mnam(cir, mnam_mask);

	// Build incidence matrix, T, and nonlinear element vector
	n_states = 0;
	max_n_states = 0;
	ElemFlag mask(NONLINEAR);
	// this routine builds the incidence matrix T based on the
	// number of nonlinear elements, and assigns the appropriate
	// values to n_states (the number of state variables
	// T is of dimension (2 x n_states) and is represents
	// how the linear and nonlinear elements are connected
	// at the "interface"
	buildTIncidence(mask, cir, T, elem_vec, n_states, max_n_states);

	// Add the compensation resistors, if desired, to the MNAM
	if (gcomp)
	{
		for (int i = 0; i < n_states; i++)
			mnam.setMAdmittance(T(0,i), T(1,i) , gcomp);
	}


	// now the MNAM is ready
	mnam.FillComplete();
	//mnam.print();    

	// handle all the convolution-based elements  
	mask = CONVOLUTION;
	// Loop through all selected elements and build convolution vector if
	// needed.
	ElementVector conv_elem_vec;
	// Reserve a fixed amount since we do not know in advance
	elem_vec.reserve(100);
	cir->setFirstElement(mask);
	Element* elem = cir->nextElement();
	while (elem)
	{
		conv_elem_vec.push_back(elem);
		// get next element pointer
		elem = cir->nextElement();
	}

	// Loop through all the nonlinear elements
	int n_elem = elem_vec.size();

	// cns contains the number of state variables per nonlinear element
	cns.resize(n_elem);

	// cns contains the number of secondary state variables per nonlinear element
	nss.resize(n_elem);
	n_sec_states = 0;
	for (int k=0; k < n_elem; k++)
	{
		cns[k] = elem_vec[k]->getNumberOfStates();
		nss[k] = elem_vec[k]->getNumberOfSecStates();
		n_sec_states += nss[k];
	}

	// Setup simulation variables (use circular vectors)
	ls_size = mnam.getDim();  // call mnam size as nm


	// list_size is the size of the circular vectors
	// each circular vector stores variable values for all time steps
	// so the size of the circ vectors is the the number of time steps
	if(tstep_method > 0)
		list_size = int(((tf - nst) / h + 3)*16);
	else
		list_size = int((tf - nst) / h + 3);


	// initialize all circular vectors
	cU = new CircVector(list_size, ls_size);
	cX = new CircVector(list_size, n_states);
	cVnl = new CircVector(list_size, n_states);
	cInl = new CircVector(list_size, n_states);
	cY = new CircVector(list_size, n_sec_states);
	cTime = new CircVector(list_size, 1);
	cTimestep = new CircVector(list_size,1);

	sf = DenseDoubleVector(ls_size);
	vL = DenseDoubleVector(n_states);

	// Choose an integration method
	// for now, these can be either Trapezoidal or Euler integration
	// Trapezoidal is the default
	// If variable time step with 2 integration techniques is selected

	if(tstep_method == 2)
	{  
		int_method = 0; // then Euler will be followed by Trapezoidal
		itypes = 2; // 2 integration technqiues
	}
	else
		itypes = 1; // only 1 integration technique

	// If only 1 integration type, select between Euler and Trapezoidal
	if (int_method == 1)
	{
		l_im = new LTrapezoidal(cU, &mnam);
		nl_im = new NLTrapezoidal(cX, cY, cU, cTime, h);
		//nl_im = new NLTrapezoidal(cX, cY, h);
	}
	else
	{
		l_im = new LEuler(cU, &mnam);
		nl_im = new NLEuler(cX, cY,cU, cTime,  h);
		//nl_im = new NLEuler(cX, cY, h);
		if(tstep_method == 2)
		{
			l_im2 = new LTrapezoidal(cU, &mnam);
			nl_im2 = new NLTrapezoidal(cX, cY, cU, cTime, h);
			//nl_im2 = new NLEuler(cX, cY, h);
		}
	}




	// Create a sparse matrix, M, which contains the contributions
	// of the static as well as dynamic portions of the MNAM
	// in other words, M = G + a*C
	// M changes as time step changes. M is Sparse Matrix
	SerialCommunicator Comm;

	// Epetra_Map is for partioning vectors and matrices
	ProcessorMap Map(ls_size, 0, Comm);
	ProcessorMap Map2(ls_size + n_states, 0, Comm);
	DoubleSparseColMatrix *M;

	// Build the MNAM for 4 timesteps: h, h/2, h/4 and h/8
	// Build MSVs T*(G + a*C)^-1 * T^T for each of these MNAMS
	// Msv stores the pointers to these 4 MSVs and Msv_ptr 
	// is the pointer to the selected MSV


	//Code to generate Tc - conventional T incidence matrix
	lhs_v = new DistributedDoubleVector(Map);

	(*lhs_v).PutScalar(zero);

	// now generate Tc, using lhs_v as a temporary vector
	// Tc_matrix is Incidence matrix in conventional format
	Tc_matrix = new DoubleDenseMatrix(n_states, ls_size); 
	for (int i = 0; i < n_states; i++)
	{
		if (T(0,i))
			(*lhs_v)[T(0,i) - 1] += one;
		if (T(1,i))
			(*lhs_v)[T(1,i) - 1] -= one;
		// copy contents of lhs_v into i_th row of Tc
		for (int j = 0; j < ls_size; j++)
			(*Tc_matrix)(i,j) = (*lhs_v)[j];
		// Clear vector for next iteration
		(*lhs_v).PutScalar(zero);
	}

	delete lhs_v;

	//create Tc transpose Tc^T
	Tct_matrix = new DoubleDenseMatrix(*Tc_matrix, Teuchos::TRANS);


	// M_map is vector containing a pair of  M matrix (for Euler Integration) and corresponding time_step_exponent 
	// M_map2 is vector containing a pair of M matrix (for Trapezoidal Integration) and corresponding time_step_exponent
	t_step_exp = 0;


	// Initial call the build the MNAM
	// It calculates G + a*C

	build_M(M,Map);

	// Create the interface class
	//tdsv = new TimeDomainSV(nl_im, max_n_states);

	tdsv = new TimeDomainSV(&(cX->getCurrent()[0]),&(cVnl->getCurrent()[0]),
			&(cInl->getCurrent()[0]),max_n_states,nl_im);

	if (tdsv->DC())
		cout <<"\n*** DC \n***\n";

	// --- Optional DC analysis setup
	// Initial condition vectors
	DenseDoubleVector X0(n_states), Vnl0(n_states);
	DenseDoubleVector Inl0(n_states), Vl0(n_states);
	if (opt && n_states)
	{
		// run a DC analysis
		DC * dca = new DC();
		double res = dca->getBias(cir, LINEAR, NONLINEAR,	X0, Vnl0, Inl0, Vl0);
		sprintf(msg, "DC analysis residual = %g", res);
		report(MESSAGE, msg);

		// Initialize the circular vectors
		for (int j = 0; j < list_size; j++)
		{
			cX->getCurrent() = X0;
			cX->advance();
			cVnl->getCurrent() = Vnl0;
			cVnl->advance();
			cInl->getCurrent() = Inl0;
			cInl->advance();
			for (int i = 0; i < n_states; i++)
				cU->getCurrent()[i] = Vl0[i];
			cU->advance();
		}
		delete dca;
	}



	// Begin main time stepping loop
	// the time steps are fixed
	residual = zero;
	report(MESSAGE, "--- Starting transient simulation ...\n");
	sprintf(msg, "Number of Nonlinear state variables: %d",n_states);
	report(MESSAGE,msg);
	sepLine();
	report(MESSAGE, "|   Step   |    Time (s)    |   Residual (V)   |");
	sepLine();

	// Set zero initial conditions
	// ctime is the current time
	double ctime = zero;  
	double prev_time = zero;
	int rc = 0;


	U_matrix = new DoubleDenseMatrix(ls_size, 1);
	X_matrix = new DoubleDenseMatrix(n_states, 1);

	X_vector = DenseDoubleVector(n_states);
	X_vector_last_iter = DenseDoubleVector(n_states);
	X_vector_last_tstep = DenseDoubleVector(n_states);
	X_vector_first_int = DenseDoubleVector(n_states); // Corresponding to first integration method if itypes=2

	U_vector = DenseDoubleVector(ls_size);
	U_vector_last_iter = DenseDoubleVector(ls_size);
	U_vector_last_tstep = DenseDoubleVector(ls_size); 
	U_vector_first_int = DenseDoubleVector(ls_size); // Corresponding to first integration method if itypes=2

	Integ_predX = new double[n_states];  // Stores X value predicted by forward Euler
	Integ_predU = new double[ls_size];   // Stored U value predicted by forward Euler

	// Number of convolution-based elements
	int n_c_elem = conv_elem_vec.size();


	//cTime->getCurrent()[0]=-2*h;
	//cTime->advance();
	cTime->getCurrent()[0]=-1*h;
	cTime->advance();
	cTime->getCurrent()[0]=zero;
	cTime->advance();
	cTimestep->getCurrent()[0]=h;
	cTimestep->advance();
	cTimestep->getCurrent()[0]=h;
	cTimestep->advance();
	//cTimestep->getCurrent()[0]=h;
	//cTimestep->advance();

	cU->advance();
	cU->advance();

	cX->advance();
	cX->advance();

	cVnl->advance();
	cVnl->advance();

	cInl->advance();
	cInl->advance();

	cY->advance();
	cY->advance();


	ctime = h;
	prev_time = zero;
	prev_timestep = h;
	orig_timestep  = h;

	// Here, we set the minimum and maximum allowable timesteps
	if( !tsmin)
		tsmin = orig_timestep / pow(2.0,4);

	if( !tsmax)
		tsmax = orig_timestep * pow(2.0,5);


	int tsmax_exp = int (floor( log(tsmax/orig_timestep) / log(2.0) ) );
	int tsmin_exp = int (ceil( log(tsmin/orig_timestep) / log(2.0) ) );

	// Time index. Printed in output as Step 
	nt = 1;

	int firsttime = 2;
	int skip_time_update = 0;
	cktorder = 1;

	// To keep the track of time step change
	changestep = 0;

	// Find 'a' corresponding to time step a = 1/h for Euler
	nl_im->changeStep(h);
	l_im->changeStep(h);

	if(tstep_method == 2)
	{
		// a = 2 / h for Trapezoidal
		nl_im2->changeStep(h);
		l_im2->changeStep(h);
	}


	// ------------------
	// the main time loop
	// ------------------
	while(ctime < tf)
	{
		// Update time in interface class
		//cTime->getCurrent()[0] = ctime;

		// update current time
		ctime = prev_time + h;
		// update current time and time step in Circular vectors storing them
		cTime->getCurrent()[0] = ctime;
		cTimestep->getCurrent()[0] = h;

		skip_time_update = 0;


		//Beginning of itype loop
		for(int c_itype = 0; c_itype < itypes; c_itype++)
		{
			tdsv->setTime(&(cX->getCurrent()[0]),
					&(cVnl->getCurrent()[0]),
					&(cInl->getCurrent()[0]),
					nt, ctime, h);


			//tdsv->clearflag();

			switch(c_itype)
			{
				case 0:    // Backward Euler case  

					// update 'a' and 'M' if h (time step) has changed
					if (changestep)
					{
						nl_im->changeStep(h);
						l_im->changeStep(h);
						build_M(M,Map);

					}

					if(tstep_method != 2)
						changestep = 0;

					// Select appropriate MNAM  
					// First search time_step_exponent in vector and assign the pointer of location to M_map_iter
					M_map_iter = M_map.find(t_step_exp);
					//  M_ptr is pointer to M matrix (G + aC) 
					M_ptr = &(*M_map_iter).second;


					// Update source vector (sf). It Changes on every time step so (sf) is updated on every time step.
					// sf = source_vector - M1p * a * b_(n-1) = source_vector - C * a * b_(n-1)

					l_im->buildSf(sf, ctime);


					//Predict the initial guess of X with FE
					nl_im->predictX(n_states,Integ_predX,cTimestep);
					//Predict the initial guess of U with FE
					nl_im->predictU(ls_size,Integ_predU,cTimestep);

					//Set intial value of state variables to a FE prediction.
					if(tstep_method >= 1)
					{
						for(int i = 0; i < n_states; i++)
						{
							X_vector[i] = Integ_predX[i];   
						}
					}
					break;

				case 1:	    

					// update 'a' and 'M' if h (time step) has changed
					if (changestep)
					{
						nl_im2->changeStep(h);
						l_im2->changeStep(h);

						build_M(M,Map);

						changestep = 0;
					}

					// First search time_step_exponent in vector and assign the pointer of location to M_map_iter  
					M_map_iter = M_map2.find(t_step_exp);
					//  M_ptr is pointer to M matrix (G + aC)
					M_ptr = &(*M_map_iter).second;

					// Update source vector (sf). It Changes on every time step so (sf) is updated on every time step
					//  sf = source_vector - M1p * a * b_(n-1) = source_vector - C * a * b_(n-1) 

					l_im2->buildSf(sf, ctime);

					break;
			}


			// Setting up left-hand and right-hand vectors for linear solve

			//ProcessorMap Map2(ls_size + n_states, 0, Comm);
			rhs_v = new DistributedDoubleVector(Map2);
			(*rhs_v).PutScalar(zero);
			lhs_v = new DistributedDoubleVector(Map2);
			(*lhs_v).PutScalar(zero);
			//const char * SolverType = "Klu";
			const char * SolverType = "Amesos_Klu";

			// Fills Larger Jacobian Matrix having contribution of all Non Linear elements. Ji and Jv are calculated

			eval_NL(&(X_vector[0]));


			// Call function to solve the linear problem Ax = b
			// It calls build_A_matrix() which build A matrix
			// It calles  build_b_vector(double* X) which build B matrix 
			// This function will do all necessary iterations until solution is found
			// The criteria of convergence is abs( U_vector[j] -  U_vector_last_iter[j] ) < reltol
			// and abs( X_vector[j] -  X_vector_last_iter[j] ) < reltol 
			Newton_Raphson_Converge = solve_L(Map2);

			(*lhs_v).PutScalar(zero);
			(*rhs_v).PutScalar(zero);

			delete lhs_v;
			delete rhs_v;



			// Update the linear terminal values after the solution is found
			eval_L();


			// Handle the situation of changing integration methods IF we use more than one
			switch(c_itype)
			{
				case 0:
					for(int i = 0 ; i < n_states; i++)
					{
						cX->getCurrent()[i]=  X_vector[i];

						if(itypes == 2)
							X_vector_first_int[i] = X_vector[i];
					}

					//nl_im->store();
					//l_im->store();

					if(itypes == 2)
					{

						for(int i = 0 ; i < ls_size; i++)
						{
							U_vector_first_int[i] = U_vector[i];
						}

						delete tdsv;
						tdsv = new TimeDomainSV(&(cX->getCurrent()[0]),&(cVnl->getCurrent()[0]),
								&(cInl->getCurrent()[0]),max_n_states,nl_im2);
						tdsv->clearflag();
					}
					break;

				case 1:
					for(int i =0 ; i < n_states; i++)
						cX->getCurrent()[i]= X_vector[i];

					//nl_im2->store();
					//l_im2->store();

					delete tdsv;
					tdsv = new TimeDomainSV(&(cX->getCurrent()[0]),&(cVnl->getCurrent()[0]),
							&(cInl->getCurrent()[0]),max_n_states,nl_im);
					tdsv->clearflag();
					break;
			} 
		} //END of itype loop





		//Predict new time step and bound the new time step to multiples of h

		if(Newton_Raphson_Converge)  // Proceed below only if Newton Raphson Converge : Find new step only if NR converge otherwise h = h/2
		{  
			if(!firsttime && (tstep_method > 0)) 
			{

				prev_timestep = h;
				new_timestep = h;
				PredictNewStep(X_vector,&new_timestep);


				if(new_timestep < 0.85 * prev_timestep)
				{        
					changestep = 1;
					//cout << "\nDecrementing time step exponent." << endl;
					t_step_exp--;

					skip_time_update = 1;
					throwaway++;



					if(new_timestep < tsmin )
					{
						t_step_exp = tsmin_exp;
						skip_time_update = 0;
					}

					h = orig_timestep * pow(2.0,t_step_exp);
					//cout << "\nNEW timestep: " << h << endl;


				}

				else if(new_timestep >= 1.5*prev_timestep)
				{
					changestep = 1;
					//cout << "\nInrementing time step exponent." << endl;
					t_step_exp++;


					if(new_timestep > tsmax )
						t_step_exp = tsmax_exp;    

					h = orig_timestep * pow(2.0,t_step_exp);
					//cout << "\nNEW timestep: " << h << endl;
				}

			}

		}

		else
		{
			changestep = 1;
			skip_time_update = 1;
			throwaway++;
			t_step_exp--;
			if(new_timestep < tsmin )
			{
				t_step_exp = tsmin_exp;
				skip_time_update = 0;
			}

			h = orig_timestep * pow(2.0,t_step_exp);
		}

		// Time update section
		// If we reduced the time step, throw away the solution and shrink the timestep
		// If timestep stayed the same or got larger, keep the solution and move forward

		if(skip_time_update == 0)
		{

			nl_im->store();
			l_im->store();

			if(itypes == 2)
			{
				nl_im2->store();
				l_im2->store();
			}

			if(firsttime)
			{
				firsttime--;
			}

			for(int i = 0 ; i < n_states; i++)
			{
				X_vector_last_tstep[i] = X_vector[i];
			}

			for(int i = 0 ; i < ls_size; i++)
			{
				U_vector_last_tstep[i] = U_vector[i];
			}  

			// Update convolution-based elements
			// Go through all convolution elements and send result vector
			for (int k = 0; k < n_c_elem; k++)
				conv_elem_vec[k]->setLastResult(cU->getCurrent(), ctime);


			// Advance the circular lists for the next step.
			cU->advance();
			cX->advance();
			cY->advance();
			cVnl->advance();
			cInl->advance();
			cTime->advance();
			cTimestep->advance();

			// Print table line
			if (!(nt % out_steps))
			{
				if (rc)
					sprintf(msg, "|  %6d  |  %e  |  %e *  |", nt, ctime, residual);
				else
					sprintf(msg, "|  %6d  |  %e  |  %e    |", nt, ctime, residual);

				report(MESSAGE, msg);
				rc = 0;
			}

			prev_time = ctime;
			nt++;

		} // End time step update if()

	} 


	// ------------------
	// end main time loop
	// ------------------


	//cout << "Throwaway solution is :" << throwaway << endl;
	//cout << "Maximum Newton Raphson Iteration :" << Max_NewtonRaphson_itr << endl;

	doOutput();

	// Erase allocated space

	delete tdsv;
	delete nl_im;
	delete nl_im2;
	delete l_im;
	delete l_im2;
	delete cY;
	delete cInl;
	delete cVnl;
	delete cX;
	delete cU;
	delete Tc_matrix;
	delete Tct_matrix;
	delete [] Integ_predX;
	delete [] Integ_predU;
	delete cTime;
	delete cTimestep;
	delete U_matrix;
	delete X_matrix;

	return;
}




// calculate the error function
void SVTran::func_ev(double * X, double * errFunc)
{

	// This function is required by inheritance, but not necessary to run, so only returns
	return;
}



// Used to update the linear node voltages

void SVTran::eval_L()
{
	for(int i = 0; i < ls_size; i++)
	{
		cU->getCurrent()[i] = U_vector[i];
	}


}







void SVTran::eval_NL(double * X)
{


	if(vnl_matrix!=NULL)
		delete vnl_matrix;
	if(inl_matrix!=NULL)
		delete inl_matrix;
	if(Jv_matrix!=NULL)
		delete Jv_matrix;
	if(Ji_matrix!=NULL)
		delete Ji_matrix ;

	// Evaluate first the time-domain elements.
	updateVInl(X);


	//for (int j=0; j < n_states; j++)
	//cX->getCurrent()[j] = X[j];


	// Tell tdsv that the first evaluation is already completed.
	//tdsv->clearflag();

	// get nonlinear current and voltages vector
	DenseDoubleVector& vnl = cVnl->getCurrent();
	DenseDoubleVector& inl = cInl->getCurrent();

	// Create "matrix" versions of vnl and inl so they are column vectors
	vnl_matrix = new DoubleDenseMatrix(vnl, Teuchos::NO_TRANS);
	inl_matrix = new DoubleDenseMatrix(inl, Teuchos::NO_TRANS);


	int n_elem = elem_vec.size();

	//tdsv->cleanJac();

	DoubleDenseMatrix& Ju = tdsv->getJu();
	DoubleDenseMatrix& Ji = tdsv->getJi();

	int ibase = 0 , jbase = 0;

	Jv_matrix = new DoubleDenseMatrix(n_states, n_states);
	Ji_matrix = new DoubleDenseMatrix(n_states, n_states);

	int rowcol_counter = 0;

	for (int j=0; j < n_elem; j++)
	{ 
		// Call element evaluation - removed and put in updateVInl()

		tdsv->setBase(ibase, cns[j], jbase, nss[j]);
		//elem_vec[j]->svTran(tdsv);
		elem_vec[j]->deriv_svTran(tdsv);

		DoubleDenseMatrix& Ju = tdsv->getJu();
		DoubleDenseMatrix& Ji = tdsv->getJi();


		/*///////
		  cout<<"\n\n\nPrint Temp Ju matrix...\n\n";
		  Ju.print(std::cout);

		  cout<<"\n\n\nPrint Temp Ji matrix...\n\n";
		  Ji.print(std::cout);
		/////// */

		// This fills the larger jacobian matrices Jv and Ji with the individual element Jacobians
		for(int k=0; k < cns[j]; k++)
			for(int l=0; l < cns[j]; l++)
			{
				(*Jv_matrix)(rowcol_counter+k,rowcol_counter+l) = 1*(Ju)(k,l);
				(*Ji_matrix)(rowcol_counter+k,rowcol_counter+l) = 1*(Ji)(k,l);
			}


		rowcol_counter += cns[j];

		//tdsv->cleanJac();


		ibase += cns[j];
		jbase += nss[j];


	} // End nonlinear element loop
	/*
	////////
	cout<<"\n\n\nPrint Jv matrix...\n\n";
	Jv_matrix->print(std::cout);

	cout<<"\n\n\nPrint Ji matrix...\n\n";
	Ji_matrix->print(std::cout);
	////////
	 */



}







// call the nonlinear element eval routines to get nonlinear I and V values
void SVTran::updateVInl(double* x_p)
{
	// Set current state variable values.
	for (int j=0; j < n_states; j++)
		cX->getCurrent()[j] = x_p[j];

	// Number of elements
	int n_elem = elem_vec.size();
	// Go through all the nonlinear elements
	for (int i=0, j=0, k=0; k < n_elem; k++)
	{
		// Set base index in interface object
		tdsv->setBase(i, cns[k], j, nss[k]);
		// Call element evaluation
		elem_vec[k]->svTran(tdsv);
		//elem_vec[k]->deriv_svTran(tdsv);
		i += cns[k];
		j += nss[k];



		DoubleDenseMatrix& Ju = tdsv->getJu();
		DoubleDenseMatrix& Ji = tdsv->getJi();

		/*
		///////
		cout<<"\n**************\n\nPrint Temp Ju matrix...\n\n";
		Ju.print(std::cout);

		cout<<"\n\n\nPrint Temp Ji matrix...\n\n";
		Ji.print(std::cout);
		cout << "********************\n";
		///////
		 */

	}
}



// This is the linear sparse matrix solver, it will also iterate until the solution is found within
// the specified error margin

int SVTran::solve_L(ProcessorMap& map)
{

	int iteration=0, eflag=0;
	eerror_X=zero;
	eerror_U=zero;

	do
	{
		A_matrix = new DoubleSparseColMatrix (Copy, map, int(ls_size/2));

		// Here, we build the A matrix of the Ax=b problem...
		(*A_matrix).PutScalar(zero);
		build_A_matrix();
		(*A_matrix).FillComplete();


		// and then build the b vector
		build_b_vector(&(X_vector[0]));

		LinProblem.SetOperator(A_matrix);
		LinProblem.SetLHS(lhs_v);
		LinProblem.SetRHS(rhs_v);
		const char * SolverType = "Klu";
		LinSolver = LinFactory.Create(SolverType, LinProblem);

		LinSolver->SymbolicFactorization();

		LinSolver->NumericFactorization();

		LinSolver->Solve();

		delete LinSolver;

		// Update U vector with most recent linear solution values
		for(int j= 0 ; j < ls_size ; j++)
		{      
			U_vector_last_iter[j] = U_vector[j];
			(*U_matrix)(j,0) = (*lhs_v)[j];
			U_vector[j] = (*U_matrix)(j,0);
		}

		// Update X vector with most recent nonlinear state variable values
		for(int j= 0 ; j < n_states ; j++)
		{
			X_vector_last_iter[j] = X_vector[j];
			(*X_matrix)(j,0) = (*lhs_v)[j+ls_size];
			X_vector[j] = (*X_matrix)(j,0);
		}

		// Use eval_NL here to update nonlinear voltages and currents based on state variables
		eval_NL(&(X_vector[0]));

		eerror_X = zero;
		eerror_U = zero;
		eflag = 0;

		for(int j= 0 ; j < ls_size && eflag==0; j++)
		{
			//eerror = abs( 1 - ( U_vector[j] /  U_vector_last_iter[j] ) );
			eerror_U = abs( U_vector[j] -  U_vector_last_iter[j] );
			if(eerror_U > reltol*0.01*MAX(abs(U_vector[j]),abs(U_vector_last_iter[j])+abstol))
				eflag = 1; 
		} // end linear error calcs

		if(n_states)
		{
			for(int j =0 ; j < n_states && eflag==0; j++)
			{
				eerror_X = abs( X_vector[j] -  X_vector_last_iter[j] );

				if(eerror_X > reltol*0.01*MAX(abs(X_vector[j]),abs(X_vector_last_iter[j]))+abstol)
					eflag = 1; 
			} // end nonlinear state variable error calcs


		} // end if(n_states)


		iteration++;


		delete  A_matrix;

	} while (iteration<max_iter && eflag==1);

	// END DO LOOP for single time step.

	Max_NewtonRaphson_itr = MAX(iteration,Max_NewtonRaphson_itr); 

	if (iteration>=max_iter)
	{
		char msg[80];
		sprintf(msg , "Newton Raphson  didn't converge: Reducing time step to half");
		report(MESSAGE, msg);
		return 0;
		//exit(1);
	}
	else 
		return 1; // NR converge

}






void SVTran::build_M(DoubleSparseColMatrix *M, ProcessorMap& map)
{

	M_map_iter = M_map.find(t_step_exp);

	if(M_map_iter != M_map.end() )
	{
		return;
	}


	M = new DoubleSparseColMatrix (Copy, map, int(ls_size/2));
	lhs_v = new DistributedDoubleVector(map);
	rhs_v = new DistributedDoubleVector(map);
	// Now compute the L and U factors of M  
	LinProblem.SetOperator(M);
	LinProblem.SetLHS(lhs_v);
	LinProblem.SetRHS(rhs_v);
	const char * SolverType = "Klu";
	LinSolver = LinFactory.Create(SolverType, LinProblem);

	// This function call calculates M = G + a C  which is M = M1 + a M1p
	// G or M1 represents conductances of MNAM
	// M1p represents capacitor/inductor part of MNAM a = 1/h
	l_im->buildMd(*M, h);


	// M is complete
	(*M).FillComplete();

	// factorize, to produce the L and U components
	LinSolver->NumericFactorization();

	delete LinSolver;

	// Insert this MNAM into the map, it will ignore if it is already present
	// M_map is storing pair of time_step_exponent and M matrix
	M_map.insert( pair<int, DoubleSparseColMatrix> (t_step_exp,*M) );

	(*lhs_v).PutScalar(zero);
	(*rhs_v).PutScalar(zero);
	(*M).PutScalar(zero);
	delete M;
	delete lhs_v;
	delete rhs_v;

	if(itypes == 2)
	{
		M = new DoubleSparseColMatrix (Copy, map, int(ls_size/2));
		lhs_v = new DistributedDoubleVector(map);
		rhs_v = new DistributedDoubleVector(map);
		// Now compute the L and U factors of M  
		LinProblem.SetOperator(M);
		LinProblem.SetLHS(lhs_v);
		LinProblem.SetRHS(rhs_v);

		LinSolver = LinFactory.Create(SolverType, LinProblem);

		// This function call calculates M = G + a C  which is M = M1 + a M1p
		// G or M1 represents conductances of MNAM
		// M1p represents capacitor/inductor part of MNAM.  a = 2/h
		l_im2->buildMd(*M, h);

		// M is complete
		(*M).FillComplete();

		// factorize, to produce the L and U components
		LinSolver->NumericFactorization();

		// M_map2 is vector for keeping MNAMS in case of itype=2
		// M_map2 is storing time_step_exponent and M matrix
		M_map2.insert( pair<int, DoubleSparseColMatrix> (t_step_exp,*M) );

		(*lhs_v).PutScalar(zero);
		(*rhs_v).PutScalar(zero);
		(*M).PutScalar(zero);
		delete M;
		delete lhs_v;
		delete rhs_v;
		delete LinSolver;
	}


}





void SVTran::build_A_matrix()
{

	//Store the size of the M-matrix
	int M_numrows = M_ptr->NumGlobalRows();
	int M_numcols = M_ptr->NumGlobalCols();

	int dimension = M_numrows;
	double * rowValExtract = new double[dimension];
	int * colIndExtract = new int[dimension];
	int num_non_zeros;

	for(int j = 0; j < M_numrows; j++)
	{
		M_ptr->ExtractGlobalRowCopy(j, dimension, num_non_zeros, rowValExtract, colIndExtract);
		A_matrix->InsertGlobalValues(j, num_non_zeros, rowValExtract, colIndExtract);
	}


	double* temp_double = new double(0.0);
	int* temp_int = new int(0);


	//Copy in Tc matrix
	for(int j = 0; j < Tc_matrix->numRows(); j++)
		for(int k = 0; k < Tc_matrix->numCols(); k++)
			if((*Tc_matrix)(j,k) == 1 || (*Tc_matrix)(j,k) == -1)
			{
				*temp_double = (*Tc_matrix)(j,k);
				*temp_int = k;
				A_matrix->InsertGlobalValues(M_numrows + j, 1, temp_double, temp_int);
			}



	//Calculate the upper-right corner of the A matrix: -1*Tc_transpose*Ji

	DoubleDenseMatrix* Upper_Right = new DoubleDenseMatrix(ls_size, n_states);

	if(n_states)
	{
		//// Changed this to +1....
		Upper_Right->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, (*Tct_matrix), (*Ji_matrix), zero);
	}

	//Copy in -1*Tc_transpose*Ji
	for(int j = 0; j < Upper_Right->numRows(); j++)
		for(int k = 0; k < Upper_Right->numCols(); k++)
			if((*Upper_Right)(j,k) != 0 )
			{
				*temp_double = (*Upper_Right)(j,k);
				*temp_int = k + M_numrows;
				A_matrix->InsertGlobalValues(j, 1, temp_double, temp_int);
			}


	//Insert -1*Jv into the lower right corner

	//Copy in -1*Jv
	for(int j = 0; j < Jv_matrix->numRows(); j++)
		for(int k = 0; k < Jv_matrix->numCols(); k++)
			if((*Jv_matrix)(j,k) != 0 )
			{
				*temp_double = -1*(*Jv_matrix)(j,k);
				*temp_int = k + M_numrows;
				A_matrix->InsertGlobalValues(j + M_numrows, 1, temp_double, temp_int);
			}

	delete [] rowValExtract;
	delete [] colIndExtract;
	delete  temp_double;
	delete  temp_int;
	delete   Upper_Right; 

}


void SVTran::build_b_vector(double* X)
{
	//convert the X vector to a single column matrix
	DenseDoubleVector temp_x;
	temp_x.size(n_states);

	for (int j=0; j < n_states; j++)
		temp_x[j] = X[j];


	DoubleDenseMatrix* x_matrix;
	x_matrix = new DoubleDenseMatrix(temp_x, Teuchos::NO_TRANS);

	//The upper portion of the b column vector is all related algebraically, so it
	//will be formed here as: sf + Tct*[inl(x) - Ji*x]

	DoubleDenseMatrix* upper_portion = new DoubleDenseMatrix(ls_size,1);
	DoubleDenseMatrix* temp_portion = new DoubleDenseMatrix(n_states,1);

	//Build -1*Ji*x
	if(n_states)
		temp_portion->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1, (*Ji_matrix), (*x_matrix), zero);


	//Build inl(x) - Ji*x
	(*temp_portion) += (*inl_matrix);


	//Build Tct*[inl(x) - Ji*x]
	//Changed this to -1
	if(n_states)
		upper_portion->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1, (*Tct_matrix), (*temp_portion), zero);

	DoubleDenseMatrix* sf_matrix;
	sf_matrix = new DoubleDenseMatrix(sf, Teuchos::NO_TRANS);

	(*upper_portion) += (*sf_matrix);


	//Now we will build the lower portion that is algebraically related
	//Defined as vnl - Jv*x

	DoubleDenseMatrix* lower_portion;
	lower_portion = new DoubleDenseMatrix(n_states,1);

	//DoubleDenseMatrix lower_portion(n_states,1);

	if(n_states)
		lower_portion->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1, (*Jv_matrix), (*x_matrix), zero);

	(*lower_portion) += (*vnl_matrix);

	//Fill in rhs_v with the upper and lower portions

	//Filling in upper portion

	double* temp_double = new double(0.0);
	int* temp_int = new int(0);

	for(int j = 0; j < ls_size; j++)
	{
		*temp_double = (*upper_portion)(j,0);
		*temp_int = j;
		rhs_v->ReplaceGlobalValues(1, temp_double, temp_int);
	}

	//Now we need to fill in the lower portion

	for(int j = 0; j < n_states; j++)
	{
		*temp_double = (*lower_portion)(j,0);
		*temp_int = j + ls_size;
		rhs_v->ReplaceGlobalValues(1, temp_double, temp_int);
	}


	delete x_matrix;
	delete temp_portion;
	delete upper_portion;
	delete lower_portion;
	delete sf_matrix;
	delete temp_double;
	delete temp_int;
}



// write output vectors
void SVTran::doOutput()
{
	// First check if the result matrices contain any data
	assert(cX);

	report(MESSAGE, "--- Writing output vectors ...");

	int out_size;
	if(nt>list_size)
		out_size=list_size;
	else
		out_size=nt;


	// Allocate and copy global time vector
	allocTimeV_P(out_size);
	for (int tindex = 0; tindex < out_size; tindex++)
		TimeV_P[tindex] = cTime->getPrevious(out_size-tindex)[0];


	// Temporary vectors
	DenseDoubleVector tmp_x(out_size);
	DenseDoubleVector tmp_i(out_size);
	DenseDoubleVector tmp_u(out_size);

	// Now fill currents and port voltages of time domain devices.
	int current_ns;
	int n_elem = elem_vec.size();   // Number of nonlinear elements
	int i=0;
	for (int k=0; k < n_elem; k++)
	{
		current_ns = elem_vec[k]->getNumberOfStates();
		for(int j=0; j< current_ns; j++)
		{
			// For the current, decompensate while copying.
			for (int tindex=0; tindex < out_size; tindex++)
			{
				tmp_x[tindex] = cX->getPrevious(out_size - tindex)[j+i];
				tmp_i[tindex] = cInl->getPrevious(out_size - tindex)[j+i];
				tmp_u[tindex] = cVnl->getPrevious(out_size - tindex)[j+i];
			}
			elem_vec[k]->getElemData()->setRealX(j, tmp_i);
			elem_vec[k]->getElemData()->setRealI(j, tmp_i);
			elem_vec[k]->getElemData()->setRealU(j, tmp_u);
		}
		i += current_ns;
	}

	// For each terminal, assign voltage vector
	Terminal* term = NULL;
	cir->setFirstTerminal();
	while((term = cir->nextTerminal()))
	{
		// Get MNAM index
		if (term->getRC())
		{
			// Remember that MNAM indices begin at 1
			i = term->getRC() - 1;
			for (int tindex=0; tindex < out_size; tindex++)
				tmp_u[tindex] = cU->getPrevious(out_size - tindex)[i];
		}
		else
		{
			// This is a reference terminal
			int tmp_u_length = tmp_u.length();
			for (int i = 0; i < tmp_u_length; i++)
				tmp_u[i] = zero;
		}

		// Set terminal vector
		term->getTermData()->setRealV(tmp_u);
	}

	ElemFlag mask = LINEAR;
	// Loop throw all selected elements in circuit and fill current
	// vector if needed.
	cir->setFirstElement(mask);
	Element* elem = cir->nextElement();
	unsigned first_eqn, no_eqn,num_terms;
	while(elem)
	{
		elem->getExtraRC(first_eqn, no_eqn);
		//num_terms = elem->getNumTerms();
		if (first_eqn)
		{
			// Fill current vector(s) in element
			//for (unsigned i = 0; i < num_terms; i++)
			for (unsigned i = 0; i < no_eqn; i++)
			{
				// Remember that MNAM indices begin at 1
				unsigned row = first_eqn-1 + i;
				for (int tindex=0; tindex < out_size; tindex++)
					tmp_u[tindex] = cU->getPrevious(out_size - tindex)[row];
				elem->getElemData()->setRealI(i, tmp_u);
			}
		}
		// get next element pointer
		elem = cir->nextElement();
	}
}





void SVTran::PredictNewStep(DenseDoubleVector& X, double *timestep)
{
	double tol=0, diff=0, tmp=0, timetemp = 1e18;

	int tstep_red_flag = 0;  //flag to indicate the time step reduction

	// Start Spice like time step control

	double LTEX,LTEX_max,X_max,timestep_Xtemp=0,timestep_Xtemp_min=0,temp=0;
	//double tol=0;

	for(int i=0 ; i < n_states ; i++)
	{
		if (tstep_method == 2)
		{
			if(abs(X_vector_first_int[i]) > abs(X_vector[i]))
				X_max = abs(X_vector_first_int[i]);
			else
				X_max = abs(X_vector[i]); 

			LTEX = abs((X_vector[i])-(X_vector_first_int[i]));

			if(i==0) LTEX_max = LTEX;
			else if(LTEX > LTEX_max)
				LTEX_max = LTEX;
		}   

		else
		{
			if(abs(Integ_predX[i]) > abs(X_vector[i]))
				X_max = abs(Integ_predX[i]);
			else
				X_max = abs(X_vector[i]);

			LTEX = abs((Integ_predX[i])-(X_vector[i]));
			if(i==0) LTEX_max = LTEX;
			else if(LTEX > LTEX_max)
				LTEX_max = LTEX;
		}

		temp = sqrt(1*(X_max * reltol + abstol)/(LTEX_max));
		timestep_Xtemp =  prev_timestep * sqrt(temp*10);

		if(i==0) timestep_Xtemp_min = timestep_Xtemp;
		else
			if(timestep_Xtemp < timestep_Xtemp_min)
				timestep_Xtemp_min = timestep_Xtemp; 
	}

	double LTEU,LTEU_max, U_max, timestep_Utemp, timestep_Utemp_min = 0;

	for(int i=0; i < ls_size; i++)
	{

		if (tstep_method == 2)
		{
			if(abs(U_vector_first_int[i]) > abs(U_vector[i]))
				U_max = abs(U_vector_first_int[i]);
			else
				U_max = abs(U_vector[i]); 

			LTEU = abs((U_vector[i])-(U_vector_first_int[i]));


			if(i==0) LTEU_max = LTEU;
			else if(LTEU > LTEU_max)
				LTEU_max = LTEU;
		}

		else
		{
			if(abs(Integ_predU[i]) > abs(U_vector[i]))
				U_max = abs(Integ_predU[i]);
			else
				U_max = abs(U_vector[i]);

			LTEU = abs((Integ_predU[i])-(U_vector[i]));
			if(i==0) LTEU_max = LTEU;
			else if(LTEU > LTEU_max)
				LTEU_max = LTEU;
		}

		temp = sqrt(1*(U_max * reltol + abstol)/(LTEU_max));
		timestep_Utemp =  prev_timestep * sqrt(temp*10);

		residual = LTEU_max;  // Residual is difference between actual solution and predicted solution  
		if(i==0) timestep_Utemp_min = timestep_Utemp;
		else
			if(timestep_Utemp < timestep_Utemp_min)
				timestep_Utemp_min = timestep_Utemp; 
	}


	temp = MIN(timestep_Utemp,timestep_Xtemp);

	*timestep = temp;

	// End Spice like time step control
	return;
} 
