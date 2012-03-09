#include "SVHB.h"
#include "FreqMNAM.h"

extern "C"
{
	#include "../inout/ftvec.h"
	#include "../inout/report.h"
}

#include "FreqDomainSV.h"

// NOX Objects
#include "NOX.H"
#include "NOX_Epetra.H"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"


#include "Teuchos_ParameterList.hpp"
// Static members
const unsigned SVHB::n_par = 12;

// Element information
ItemInfo SVHB::ainfo =
{
  "svhb",
  "State Variable-Based Harmonic Balance Analysis",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS
};

// Parameter information
ParmInfo SVHB::pinfo[] =
{
  {"n_freqs", "Maximum index for first tone (DC not included)",
	TR_INT, true},
  {"fundamental", "Fundamental frequency (Hz)", TR_DOUBLE, true},
  {"oversample", "Use oversample in the FFT", TR_INT, false},
  {"steps", "source stepping", TR_INT, false},
  {"deriv", "Approximate derivatives or use automatic diff.", TR_INT, false},
  {"verbosity", "Amount of output to print", TR_INT, false},
  {"fundamental2", "Frequency of the second tone.", TR_DOUBLE, false},
  {"n_freqs2", "Maximum index for second tone (DC not included)",
	TR_INT, false},
  {"regrowth", "Indicates whether a spectral regrowth analysis is performed",
	TR_BOOLEAN, false},
  {"n_fund", "Number of fundamental frequencies for spectral regrowth",
	TR_INT, false},
  {"f_step", "Frequency step for spectral regrowth (Hz)", TR_DOUBLE, false},
  {"step_init", "Initial step in source stepping", TR_INT, false},
  {"max_cfreq",
		"Maximum center frequency to consider in spectral regrowth analysis",
	TR_DOUBLE, false}
};


SVHB::SVHB() : Analysis(&ainfo, pinfo, n_par)
{
  paramvalue[0] = &n_freqs;
  paramvalue[1] = &fundamental;
  paramvalue[2] = &(oversample_par = 1);
  paramvalue[3] = &(steps_par = 0);
  //paramvalue[4] = &(diff_par = 0);
  paramvalue[4] = &(verbosity_par = 1);
  paramvalue[5] = &(fundamental2 = zero);
  paramvalue[6] = &(n_freqs2 = 0);
  paramvalue[7] = &(regrowth = false);
  paramvalue[8] = &(n_fund = 0);
  paramvalue[9] = &(f_step = zero);
  paramvalue[10] = &(step_init = 0);
  paramvalue[11] = &(max_cfreq = zero);

  my_cir = NULL;
  tape_flag = false;
}


SVHB::~SVHB()
{
}

void SVHB::run(Circuit* cir)
{
  my_cir = cir;

  // Char array to store output messages
  //  char msg[80];
  report(MESSAGE, "*** State Variable-Based Harmonic Balance Analysis ***");
  sepLine();
  newLine();

  if (!isSet(&n_freqs) || !isSet(&fundamental))
    report(FATAL, "Both n_freqs and fundamental must be specified");
  //Check that the number of tones and index for each are consistent
  if (isSet(&n_freqs2) && !isSet(&fundamental2) ||
		!isSet(&n_freqs2) && isSet(&fundamental2))
	report(FATAL,
	"Both n_freqs2 and fundamental2 must be set");
  // Check consistency of other paramenters
  if (n_freqs2 || regrowth)
	{
    if (oversample_par != 1)
      report(WARNING, "Setting oversample to 1 for multitone analysis");
    oversample_par = 1;
  }

  // First create the options structure so the code can run unchanged
  // Assign the parameter values to the elements to the structure svHB_options
  svHBopt = (svHB_options_t *) malloc(sizeof(svHB_options_t));

  svHBopt->numtones = (n_freqs2) ? 2 : 1;;
  svHBopt->oversample = oversample_par;
  svHBopt->steps = steps_par;
  svHBopt->diff = diff_par;
  svHBopt->verbosity = verbosity_par;

  if (regrowth)
	{
    // Handle spectral regrowth analysis
    if (n_fund < 2)
      report(FATAL, "n_fund must be greater than 2.");
    if (!f_step)
      report(FATAL, "f_step must be greater than zero.");
    report(WARNING,
		"n_fund must be odd by now. Not sure what happens otherwise.");
    // Increase numtones to avoid trouble
    svHBopt->numtones = 10;
  }

  // Allocate memory for the fundamental vectors
  svHBopt->tone = Mlib_DNewVec(svHBopt->numtones);
  svHBopt->h = Mlib_INewVec(svHBopt->numtones);
  // Set the tones by hand until we solve the "sweep problem"
  svHBopt->h[0] = n_freqs; //h is the number of frequencies
  svHBopt->tone[0] = fundamental;
  if (svHBopt->numtones == 2)
	{
    svHBopt->h[1] = n_freqs2;
    svHBopt->tone[1] = fundamental2;
  }
  /* Initial time */
  float t0;
  time_t tot_time;
  double residual;
  /* Vector of unknowns (state variables) */
  doublev_t X;

  /*
	* set the initial time:
	*                       t0: set-up time measurement (cpu seconds)
	*                       tot_time: total time (real seconds => depends on
	*                               processor load).
	*/
  t0 = (float) clock() / CLOCKS_PER_SEC;
  tot_time = time(NULL);

  /*
	* Create the compressed source vectors
	*
	* svSource = T M^(-1) S
	*
	* and the compressed MNAM matrices
	*
	* svMatrix = T M^(-1) T'
	*
	* M: MNAM
	* S: source vector
	*
	*/
  create_svHBdata();

  /* Create the state vector. We prefer to create it as a long real
	* vector to make it easier to interface with external equation
	* solvers. The dimension was calculated in Create_svHBdata().
	*/
  X = Mlib_DNewVec(svHBdata->sysdim);

  /* Generate the initial guess for X
	* Initial guess = 0.0
	*/
  for(int i=0; i < svHBdata->sysdim; i++)
    X[i] = zero;

  /* Create work space for calculations, in SVHB_phys.cc */
  createWs();

  /* output the set-up time */
  fprintf(output_F, "\n  ** Set-up CPU time: %6.2f s\n\n",
	(float)clock()/CLOCKS_PER_SEC - t0);
  fflush(output_F);

  /* After this point we are done with the initialization. Now we send
	* the rest of the problem to the nonlinear equation solver. It is
	* called from the interface routine svHBsolve().
	*/
  residual = solve(X);
  //Copy over the X vector from the interface class
  //As harmonice balance needs this vector to compute
  //the output.
  for(int i=0; i < svHBdata->sysdim; i++)
	  X[i] = interface->x_svtran[i];
  /* Now get the voltages at all the circuit terminals and write
	* results in the terminal and node tables.
	*/
  calc_output(X, residual);

  /* Free vectors, structures, etc. The frequency vector is not freed
	* because it is used by transim output routines.
	*/

  /* Free calculation work space */
  destroyWs();

  Mlib_DFreeVec(X);

  destroy_svHBdata();

  // Free memory allocated in this routine
  Mlib_DFreeVec(svHBopt->tone);
  Mlib_IFreeVec(svHBopt->h);
  free(svHBopt);

  /* output total time */
  tot_time = time(NULL) - tot_time;
  fprintf(output_F, " *** Total analysis time: %d min %d sec\n\n",
	int(tot_time/60), int(tot_time%60));
  fflush(output_F);
}

void SVHB::create_svHBdata()
{
  int i, findex;
  dcxv_t Source, Xtemp;

  /* Allocate memory for svHBdata. The fields of this structure should
	* be set in this routine only. Other routines can read the structure
	* but they should not change anything.
	*/
  svHBdata = (svHB_data_t*) malloc(sizeof(svHB_data_t));

  /* First, set up the frequency vector.
	* We need the frequency vector to create the MNAM.
	*/
  createFvec();

  // Create frequency vector as a reference to existing vector
  DenseDoubleVector freq_vector(Teuchos::View, svHBdata->FreqV_P, svHBdata->NoFreqPoints);
  // Build MNAM
  // Set flags to linear
  ElemFlag mask = LINEAR;
  FreqMNAM* mnam = new FreqMNAM(freq_vector, my_cir, mask);
  svHBdata->mnam = mnam;
  // this step is required
  svHBdata->mnam->FillComplete();

  /* Fill Omega vector. Also get the value of NoSamples.
	* Memory is allocated inside the routine.
	*/
  svHBdata->Omega = createOmega(&(svHBdata->NoSamples));

  /* Create T matrix. The routine returns the number of states of our
	* circuit.  Allocate space for matrix T and W. T is the transpose
	* of W (W=T'), so we store only one matrix. They are the same for
	* all frequencies.  Call the routine at the physical level to
	* create the matrix.
	*
	* Size of T is:
	*               n_states x mnam_st->dimension
	*/
  svHBdata->n_states = createT();

  fprintf(output_F,"\nTotal number of frequencies: %d\n",
	svHBdata->NoFreqPoints);

  fprintf(output_F,"\nNumber of state variables: %d\n", svHBdata->n_states);

  /* Allocate memory for vectors and matrices inside svHBdata.
	*/
  svHBdata->svSource = Mlib_CNewMat(svHBdata->NoFreqPoints,
	svHBdata->n_states);
  svHBdata->svMatrix =
	Mlib_CNew3DMat(svHBdata->NoFreqPoints, svHBdata->n_states,
	svHBdata->n_states);

  /* System Dimension:
	*           n_states * ( 2 * svHBdata->NoFreqPoints - 1 )
	*
	* We don't store the imaginary part of the DC component.
	*/
  svHBdata->sysdim = (2*svHBdata->NoFreqPoints - 1) * svHBdata->n_states;

  fprintf(output_F,"\nSystem dimension: %d\n", svHBdata->sysdim);

  /* Allocate memory for temporary vectors */
  Source = Mlib_CNewVec(svHBdata->mnam->getDim());
  Xtemp = Mlib_CNewVec(svHBdata->mnam->getDim());

  /* Temporary vectors */
  DenseComplexVector tempvec1(svHBdata->mnam->getDim());
  DenseComplexVector tempvec2(svHBdata->mnam->getDim());

  /* Create a compressed vector and matrix for each frequency.
	*/
  for(findex = 0; findex < svHBdata->NoFreqPoints; findex++)
	{
    /*
		* We don't invert the MNA explicitly. Instead:
		* - First, we re-order the MNA to prevent having zeros in the diagonal
		* - Second, perform LU decomposition
		* - Then, solve M x = s
		* - Then we solve M x = (T')i for each column of T' (= row of T).
		* - Finally, we invert the sign of the compressed matrix.
		*/

    // Solve for default source vector
    svHBdata->mnam->solve(findex, tempvec2);

    // Copy tempvec2 to Xtemp
    for (unsigned j=0; j<svHBdata->mnam->getDim(); j++)
		{
      Xtemp[j].re = tempvec2[j].real();
      Xtemp[j].im = tempvec2[j].imag();
    }

    /* Now we should multiply T and Xtemp to get vector svSource */
    /* arguments: result, matrix, vector, dimensions */
    Mlib_CMultVecRec(svHBdata->svSource[findex], svHBdata->T, Xtemp,
		svHBdata->n_states, svHBdata->mnam->getDim());

		/* We do the same for each column of W (= row of T)
		*/
    for (i=0; i<svHBdata->n_states; i++)
		{
      // Copy row of T into tempvec
      for (unsigned j=0; j<svHBdata->mnam->getDim(); j++)
			{
        tempvec1[j].real() = svHBdata->T[i][j].re;
        tempvec1[j].imag() = svHBdata->T[i][j].im;
			}

      svHBdata->mnam->solve(findex, tempvec1, tempvec2);

      // Copy tempvec2 to Xtemp
      for (unsigned j=0; j<svHBdata->mnam->getDim(); j++)
			{
				Xtemp[j].re = tempvec2[j].real();
				Xtemp[j].im = tempvec2[j].imag();
      }
      Mlib_CMultVecRec(svHBdata->svMatrix[findex][i], svHBdata->T, Xtemp,
			svHBdata->n_states, svHBdata->mnam->getDim());
    }
  }

  /* Free temporary vectors */
  Mlib_CFreeVec(Xtemp);
  Mlib_CFreeVec(Source);

  fprintf(output_F,
	"\n*****************************************************************************\n\n");
}

void SVHB::destroy_svHBdata()
{
  // free MNAM
  delete svHBdata->mnam;

  // Free matrices and vectors first.
  Mlib_DFreeVec(svHBdata->FreqV_P);

  Mlib_DFreeVec(svHBdata->Omega);
  Mlib_CFreeMat(svHBdata->svSource);
  Mlib_CFree3DMat(svHBdata->svMatrix);

  destroyT();

  free(svHBdata);
}

double SVHB::solve(doublev_t X)
{
  int retcod, max_s, i, j, k;
  double residual = 0;
  dcxm_t tempm;

  if (svHBopt->steps > 0)
	{
    max_s = svHBopt->steps;
    /* Allocate temporary space to save the original s.v. source
		* vector.
		*/
    tempm= Mlib_CNewMat(svHBdata->NoFreqPoints, svHBdata->n_states);

    /* Save the original vector */
    for(j=0; j < svHBdata->NoFreqPoints; j++)
		{
      for(i=0; i < svHBdata->n_states; i++)
			{
				tempm[j][i].re = svHBdata->svSource[j][i].re;
				tempm[j][i].im = svHBdata->svSource[j][i].im;
      }
    }

    if (svHBopt->numtones > 1)
		{
      // ********** Begin Setup Nonlinear Solver *******************
      SerialCommunicator Comm1;

      // Get the process ID and the total number of processors
      int MyPID = Comm1.MyPID();
      int NumProc = Comm1.NumProc();

      // Get the number of elements
      int NumGlobalElements = svHBdata->sysdim; 

      // The number of unknowns must be at least equal to the 
      // number of processors.
      if (NumGlobalElements < NumProc) {
        cout << "numGlobalBlocks = " << NumGlobalElements << " cannot be < number of processors = " << NumProc << endl;
        cout << "Test failed!" << endl;
        throw "NOX Error";
      }
      // Create the interface between NOX and the application
      // This object is derived from NOX::Epetra::Interface
      interface = Teuchos::rcp(new Nox_Interface(this, NumGlobalElements, Comm1));
      //interface = Teuchos::rcp(new Interface(this, NumGlobalElements, Comm));

      // Get the vector from the Problem
      Teuchos::RCP<Epetra_Vector> InitialGuess = interface->getSolution();

      // Set the initial guess 
      InitialGuess->PutScalar(0.0);

      // Create the top level parameter list
      Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
        Teuchos::rcp(new Teuchos::ParameterList);
      Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

      // Set the nonlinear solver method
      nlParams.set("Nonlinear Solver", "Line Search Based");

      // Set the printing parameters in the "Printing" sublist
      Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
      printParams.set("MyPID", MyPID); 
      printParams.set("Output Precision", 3);
      printParams.set("Output Processor", 0);
      printParams.set("Output Information", NOX::Utils::Error + NOX::Utils::TestDetails);

      // Create a print class for controlling output below
      NOX::Utils printing(printParams);

      // Sublist for line search 
      Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
      searchParams.set("Method", interface->method);
      Teuchos::ParameterList dirParams;
      Teuchos::ParameterList directionParams;

      switch (interface->direction)
      {
        case 1:
          // Sublist for direction
          dirParams = nlParams.sublist("Direction");
          dirParams.set("Method", "Newton");
          directionParams = dirParams.sublist("Newton");
          directionParams.set("Forcing Term Method", "Constant");
          break;
        default:
          dirParams = nlParams.sublist("Direction");
          dirParams.set("Method", "Newton");
          directionParams = dirParams.sublist("Newton");
          directionParams.set("Forcing Term Method", "Constant");
      }

      // Sublist for linear solver for the Newton method
      Teuchos::ParameterList& lsParams = directionParams.sublist("Linear Solver");
      lsParams.set("Aztec Solver", "GMRES");  
      lsParams.set("Max Iterations", 800);  
      lsParams.set("Tolerance", 1e-4);
      lsParams.set("Aztec Precondictioner", "ilu");

      Teuchos::RCP<Epetra_CrsMatrix> A = interface->getJacobian();

      NOX::Epetra::Vector noxInitGuess(InitialGuess, NOX::DeepCopy);

      // Create the linear system
      Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
      Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
        Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
              lsParams, iReq, iJac, A, noxInitGuess));

      // Create the Group
      Teuchos::RCP<NOX::Epetra::Group> grpPtr = 
        Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, noxInitGuess, linSys));
      NOX::Epetra::Group& grp = *grpPtr;

      // Set up the status tests
      Teuchos::RCP<NOX::StatusTest::NormF> testNormF = 
        Teuchos::rcp(new NOX::StatusTest::NormF(interface->ftol));
      Teuchos::RCP<NOX::StatusTest::MaxIters> testMaxIters = 
        Teuchos::rcp(new NOX::StatusTest::MaxIters(interface->maxit));
      // this will be the convergence test to be used
      Teuchos::RCP<NOX::StatusTest::Combo> combo = 
        Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
              testNormF, testMaxIters));

      // Create the convergence tests
      Teuchos::RCP<NOX::StatusTest::NormF> absresid = 
        Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
      Teuchos::RCP<NOX::StatusTest::NormF> relresid = 
        Teuchos::rcp(new NOX::StatusTest::NormF(grp, 1.0e-2));
      Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
        Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
      Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
        Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
      Teuchos::RCP<NOX::StatusTest::Combo> converged =
        Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
      converged->addStatusTest(absresid);
      converged->addStatusTest(relresid);
      converged->addStatusTest(wrms);
      converged->addStatusTest(update);
      Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
        Teuchos::rcp(new NOX::StatusTest::MaxIters(interface->maxit));
      Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
        Teuchos::rcp(new NOX::StatusTest::FiniteValue);
      combo->addStatusTest(fv);
      combo->addStatusTest(converged);
      combo->addStatusTest(maxiters);

      // Create the solver
      solver = NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);

      NOX::StatusTest::StatusType solvStatus;
      NOX::Epetra::Vector* guess;
      const NOX::Epetra::Group* soln_group = 0;
      const Epetra_Vector* next_guess = 0;			
      for (k=step_init; k <= max_s; k++) 
      {
        for(j=0; j < svHBdata->NoFreqPoints; j++)
        {
          for(i=0; i < svHBdata->n_states; i++)
          {
            svHBdata->svSource[j][i].re = tempm[j][i].re * (float)k / max_s;
            svHBdata->svSource[j][i].im = tempm[j][i].im * (float)k / max_s;
          }
        }

        printf("\n ------------------------------ Step = %d \n", k);

        Teuchos::RCP<Epetra_Vector> InitialGuess = interface->getSolution();
        // Set the initial guess 
        InitialGuess->PutScalar(0.0);
        NOX::Epetra::Vector noxInitGuess(InitialGuess, NOX::DeepCopy);
        //Reset the solver.
        solver->reset(noxInitGuess); 
        // Solve the nonlinear equation (call the nonlinear solver).
        int solvStatus = interface->evaluate(solver);				 
        //Get the residual
        const NOX::Epetra::Group& soln_group = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
        residual = soln_group.getNormF();

        printf("\n\n *** Residual = %g\n", residual);
        if (svHBopt->verbosity > 0)
        {
          fprintf(output_F, "\n --- Step = %d, Residual = %g\n\n", k, residual);
          fflush(output_F);
        }
      }
    }
    else
    {
      // Single tone source stepping => use also tone stepping
      // Save the original number of tones
      int nACtones = svHBopt->h[0];

      /* Solve the system at each step. The result from one step is the
       * initial guess for the next.
       */
      for (k=step_init; k <= max_s; k++)
      {
        // Use tone-stepping with source stepping
        // This is a quick hack to reduce the computation
        // time for the first steps. We solve only for 2 tones
        // at DC (the minimun would be one, but that number may cause
        // problems in some parts of the program). Then we increment the
        // number of tones at each step. Note that we can use always
        // the same vectors and matrices, since they are dimensioned for
        // the maximum number of tones. For the multi-tone case, things
        // would be a little more difficult.
        svHBopt->h[0] = 1 + k * (nACtones - 1) / max_s;
        svHBdata->NoFreqPoints = svHBopt->h[0]+1;
        // Set number of used frequencies (save time building the Jacobian).
        fdsv->setNFreqs(svHBdata->NoFreqPoints);
        svHBdata->sysdim = (2*svHBdata->NoFreqPoints - 1) *
          svHBdata->n_states;
        // ********** Begin Setup Nonlinear Solver *******************
        SerialCommunicator Comm1;

        // Get the process ID and the total number of processors
        int MyPID = Comm1.MyPID();
        int NumProc = Comm1.NumProc();

        // Get the number of elements
        int NumGlobalElements = svHBdata->sysdim; 

        // The number of unknowns must be at least equal to the 
        // number of processors.
        if (NumGlobalElements < NumProc) {
          cout << "numGlobalBlocks = " << NumGlobalElements << " cannot be < number of processors = " << NumProc << endl;
          cout << "Test failed!" << endl;
          throw "NOX Error";
        }
        // Create the interface between NOX and the application
        // This object is derived from NOX::Epetra::Interface
        interface = Teuchos::rcp(new Nox_Interface(this, NumGlobalElements, Comm1));

        // Get the vector from the Problem
        Teuchos::RCP<Epetra_Vector> InitialGuess = interface->getSolution();

        // Set the initial guess 
        InitialGuess->PutScalar(0.0);

        // Create the top level parameter list
        Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
          Teuchos::rcp(new Teuchos::ParameterList);
        Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

        // Set the nonlinear solver method
        nlParams.set("Nonlinear Solver", "Line Search Based");

        // Set the printing parameters in the "Printing" sublist
        Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
        printParams.set("MyPID", MyPID); 
        printParams.set("Output Precision", 3);
        printParams.set("Output Processor", 0);
        printParams.set("Output Information", NOX::Utils::Error +
            NOX::Utils::TestDetails);

        // Create a print class for controlling output below
        NOX::Utils printing(printParams);

        // Sublist for line search 
        Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
        searchParams.set("Method", interface->method);
        Teuchos::ParameterList dirParams;
        Teuchos::ParameterList directionParams;

        switch (interface->direction)
        {
          case 1:
            // Sublist for direction
            dirParams = nlParams.sublist("Direction");
            dirParams.set("Method", "Newton");
            directionParams = dirParams.sublist("Newton");
            directionParams.set("Forcing Term Method", "Constant");
            break;
          default:
            dirParams = nlParams.sublist("Direction");
            dirParams.set("Method", "Newton");
            directionParams = dirParams.sublist("Newton");
            directionParams.set("Forcing Term Method", "Constant");
        }  	
        // Sublist for linear solver for the Newton method
        Teuchos::ParameterList& lsParams = directionParams.sublist("Linear Solver");
        lsParams.set("Aztec Solver", "GMRES");  
        lsParams.set("Max Iterations", 800);  
        lsParams.set("Tolerance", 1e-4);
        lsParams.set("Output Frequency", 50);
        lsParams.set("Aztec Precondictioner", "ilu");

        Teuchos::RCP<Epetra_CrsMatrix> A = interface->getJacobian();

        NOX::Epetra::Vector noxInitGuess(InitialGuess, NOX::DeepCopy);

        // Create the linear system
        Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
        Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
        Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
          Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                lsParams, iReq, iJac, A, noxInitGuess));

        // Create the Group
        Teuchos::RCP<NOX::Epetra::Group> grpPtr = 
          Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, noxInitGuess, linSys));
        NOX::Epetra::Group& grp = *grpPtr;

        // Set up the status tests
        Teuchos::RCP<NOX::StatusTest::NormF> testNormF = 
          Teuchos::rcp(new NOX::StatusTest::NormF(interface->ftol));
        Teuchos::RCP<NOX::StatusTest::MaxIters> testMaxIters = 
          Teuchos::rcp(new NOX::StatusTest::MaxIters(interface->maxit));
        // this will be the convergence test to be used
        Teuchos::RCP<NOX::StatusTest::Combo> combo = 
          Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
                testNormF, testMaxIters));

        // Create the convergence tests
        Teuchos::RCP<NOX::StatusTest::NormF> absresid = 
          Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
        Teuchos::RCP<NOX::StatusTest::NormF> relresid = 
          Teuchos::rcp(new NOX::StatusTest::NormF(grp, 1.0e-2));
        Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
          Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
        Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
          Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
        Teuchos::RCP<NOX::StatusTest::Combo> converged =
          Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
        converged->addStatusTest(absresid);
        converged->addStatusTest(relresid);
        converged->addStatusTest(wrms);
        converged->addStatusTest(update);
        Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
          Teuchos::rcp(new NOX::StatusTest::MaxIters(interface->maxit));
        Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
          Teuchos::rcp(new NOX::StatusTest::FiniteValue);
        combo->addStatusTest(fv);
        combo->addStatusTest(converged);
        combo->addStatusTest(maxiters);

        // Create the solver
        solver = NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);
        NOX::StatusTest::StatusType solvStatus;
        NOX::Epetra::Vector* guess;
        //const NOX::Epetra::Group* soln_group = 0;
        const Epetra_Vector* next_guess = 0;

        //NLSInterface* nlsi =
        //new NLSInterface(this, svHBdata->sysdim, svHBopt->diff);

        // Do source stepping
        for(j=1; j < svHBdata->NoFreqPoints; j++)
        {
          for(i=0; i < svHBdata->n_states; i++)
          {
            svHBdata->svSource[j][i].re = tempm[j][i].re * (float)k / max_s;
            svHBdata->svSource[j][i].im = tempm[j][i].im * (float)k / max_s;
          }
        }

        printf("\n ------------------------------ Step = %d \n", k);
        printf("    Number of tones: %d\n   System dimension: %d\n",
            svHBdata->NoFreqPoints, svHBdata->sysdim);
        //Reset the solver.
        solver->reset(noxInitGuess); 
        // Solve the nonlinear equation (call the nonlinear solver).
        solvStatus = interface->evaluate(solver);

        //Get the residual
        const NOX::Epetra::Group& soln_group = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
        residual = soln_group.getNormF();

        printf("\n\n *** Residual = %g\n", residual);
        if (svHBopt->verbosity > 0)
        {
          fprintf(output_F, "\n --- Step = %d, Residual = %g\n\n", k, residual);
          fflush(output_F);
        }
      }
    }
    /* Free temporary vector */
    Mlib_CFreeMat(tempm);
  }
  else
  {
    // ********** Begin Setup Nonlinear Solver *******************		
    SerialCommunicator Comm1;

    // Get the process ID and the total number of processors
    int MyPID = Comm1.MyPID();
    int NumProc = Comm1.NumProc();

    // Get the number of elements
    int NumGlobalElements = svHBdata->sysdim; 

    // The number of unknowns must be at least equal to the 
    // number of processors.
    if (NumGlobalElements < NumProc) 
    {
      cout << "numGlobalBlocks = " << NumGlobalElements << " cannot be < number of processors = " << NumProc << endl;
      cout << "Test failed!" << endl;
      throw "NOX Error";
    }
    // Create the interface between NOX and the application
    // This object is derived from NOX::Epetra::Interface
    interface = Teuchos::rcp(new Nox_Interface(this, NumGlobalElements, Comm1));

    // Get the vector from the Problem
    Teuchos::RCP<Epetra_Vector> InitialGuess = interface->getSolution();

    // Set the initial guess 
    InitialGuess->PutScalar(0.0);

    // Create the top level parameter list
    Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

    // Set the nonlinear solver method
    nlParams.set("Nonlinear Solver", "Line Search Based");

    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
    printParams.set("MyPID", MyPID); 
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    printParams.set("Output Information", NOX::Utils::Error + NOX::Utils::TestDetails);

    // Create a print class for controlling output below
    NOX::Utils printing(printParams);

    // Sublist for line search 
    Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
    searchParams.set("Method", interface->method);
    Teuchos::ParameterList dirParams;
    Teuchos::ParameterList directionParams;

    switch (interface->direction)
    {
      case 1:
        // Sublist for direction
        dirParams = nlParams.sublist("Direction");
        dirParams.set("Method", "Newton");
        directionParams = dirParams.sublist("Newton");
        directionParams.set("Forcing Term Method", "Constant");
        break;
      default:
        dirParams = nlParams.sublist("Direction");
        dirParams.set("Method", "Newton");
        directionParams = dirParams.sublist("Newton");
        directionParams.set("Forcing Term Method", "Constant");
    }  

    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = directionParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");  
    lsParams.set("Max Iterations", 200);  
    lsParams.set("Tolerance", 1e-4);
    lsParams.set("Output Frequency", 50);
    lsParams.set("Aztec Precondictioner", "ilu");

    Teuchos::RCP<Epetra_CrsMatrix> A = interface->getJacobian();

    NOX::Epetra::Vector noxInitGuess(InitialGuess, NOX::DeepCopy);

    // Create the linear system
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
            lsParams, iReq, iJac, A, noxInitGuess));

    // Create the Group
    Teuchos::RCP<NOX::Epetra::Group> grpPtr = 
      Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, noxInitGuess, linSys));
    NOX::Epetra::Group& grp = *grpPtr;

    // Set up the status tests
    Teuchos::RCP<NOX::StatusTest::NormF> testNormF = 
      Teuchos::rcp(new NOX::StatusTest::NormF(interface->ftol));
    Teuchos::RCP<NOX::StatusTest::MaxIters> testMaxIters = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(interface->maxit));
    // this will be the convergence test to be used
    Teuchos::RCP<NOX::StatusTest::Combo> combo = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
            testNormF, testMaxIters));

    // Create the convergence tests
    Teuchos::RCP<NOX::StatusTest::NormF> absresid = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    Teuchos::RCP<NOX::StatusTest::NormF> relresid = 
      Teuchos::rcp(new NOX::StatusTest::NormF(grp, 1.0e-2));
    Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
      Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
    Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
      Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
    Teuchos::RCP<NOX::StatusTest::Combo> converged =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    converged->addStatusTest(absresid);
    converged->addStatusTest(relresid);
    converged->addStatusTest(wrms);
    converged->addStatusTest(update);
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(interface->maxit));
    Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
      Teuchos::rcp(new NOX::StatusTest::FiniteValue);
    combo->addStatusTest(fv);
    combo->addStatusTest(converged);
    combo->addStatusTest(maxiters);

    // Create the solver
    solver = NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);


    NOX::StatusTest::StatusType solvStatus;
    NOX::Epetra::Vector* guess;
    const Epetra_Vector* next_guess = 0;

    // Set the initial guess 
    InitialGuess->PutScalar(0.0);
    //Reset the solver.
    solver->reset(noxInitGuess); 
    // Solve the nonlinear equation (call the nonlinear solver).
    solvStatus = interface->evaluate(solver);

    //Get the residual
    const NOX::Epetra::Group& soln_group = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    residual = soln_group.getNormF();
  }

  return(residual);
}


void SVHB::func_ev(double * X, double * errFunc)
{
  int i,findex;
  float t0;
  /* Voltages and currents of the nonlinear devices */
  dcxm_t U_NL, I_NL, U_L;

  /* Print a dot to indicate normal function evaluation. */
  printf("."); fflush(stdout);

  /* Set zero time */
  t0 = (float)clock()/CLOCKS_PER_SEC;

  /* First, call the physical interface to calculate U_NL and I_NL.
	* Also, the function sets the pointer U_L to a previously allocated
	* space.
	*/
  get_U_and_I(X, &U_NL, &I_NL, &U_L);

  /* Calculate the linear voltages */
  for (findex=0; findex < svHBdata->NoFreqPoints; findex++)
  {
    Mlib_CMultVec(U_L[findex],
		svHBdata->svMatrix[findex],
		I_NL[findex],
		svHBdata->n_states);
    /* Now add svSource to U_L.
		*/
    for (i=0; i< svHBdata->n_states; i++)
    {
      U_L[findex][i].re += svHBdata->svSource[findex][i].re;
      U_L[findex][i].im += svHBdata->svSource[findex][i].im;
    }
  }

  /* Calculate the DC error components.
	*/
  for (i=0; i< svHBdata->n_states; i++)
    errFunc[i] = U_L[0][i].re - U_NL[0][i].re;

  /* Calculate the rest of the error function. The F index calculation
	* is somewhat complicated because we have to convert a complex
	* matrix to a single real vector.
	*/
  for (findex=1; findex < svHBdata->NoFreqPoints; findex++)
  {
    for (i=0; i< svHBdata->n_states; i++)
    {
      /* real part */
      errFunc[(2*findex-1)*(svHBdata->n_states) + 2*i] =
        U_L[findex][i].re - U_NL[findex][i].re;
      /* imaginary part */
      errFunc[(2*findex-1)*(svHBdata->n_states) + 2*i+1] =
        U_L[findex][i].im - U_NL[findex][i].im;
    }
  }

  /* output evaluation time */
  if (svHBopt->verbosity > 2)
	{
    fprintf(output_F, "        Function CPU time: %6.2f s \n",
		(float)clock()/CLOCKS_PER_SEC - t0);
    fflush(output_F);
  }
}

void SVHB::calc_output(doublev_t X, double residual)
{
  int findex;
  dcxv_t Source;
  /* Voltages and currents of the nonlinear devices */
  dcxm_t U_NL, I_NL, U_L, VI_P;

  /* Calculate the nonlinear currents from the state
	* variables.
	*/
  get_U_and_I(X, &U_NL, &I_NL, &U_L);

  /* Allocate memory for circuit voltage vector */
  VI_P = Mlib_CNewMat(svHBdata->NoFreqPoints, svHBdata->mnam->getDim());

  /* Allocate memory for temporary source vector */
  Source = Mlib_CNewVec(svHBdata->mnam->getDim());
  /* Temporary vectors */
	DenseComplexVector tempvec1(svHBdata->mnam->getDim());
  DenseComplexVector tempvec2(svHBdata->mnam->getDim());

  /* Solve voltages (and some currents) for all frequencies.
	*/
  for(findex=0; findex<svHBdata->NoFreqPoints; findex++)
	{
    /*
		* The MNAMs are already factorized.
		*/

    /* First calculate the contribution of the nonlinear devices to the
		* source vector.
		*/
    Mlib_CMultVecRecTrans(Source, svHBdata->T, I_NL[findex],
		svHBdata->n_states, svHBdata->mnam->getDim());

    /* Now add the contribution of the fixed sources.
		*/
    svHBdata->mnam->getSource(findex, tempvec1);
    // Add Source to temp1
    for (unsigned j=0; j<svHBdata->mnam->getDim(); j++)
    {
      tempvec1[j].real() += Source[j].re;
      tempvec1[j].imag() += Source[j].im;
    }

    svHBdata->mnam->solve(findex, tempvec1, tempvec2);

    // Copy tempvec2 to Xtemp
    for (unsigned j=0; j<svHBdata->mnam->getDim(); j++)
		{
      VI_P[findex][j].re = tempvec2[j].real();
      VI_P[findex][j].im = tempvec2[j].imag();
    }
  }

  /* Output results */
  doOutput(VI_P, I_NL, X, residual);

  /* Free temporary vectors */
  Mlib_CFreeVec(Source);
  Mlib_CFreeMat(VI_P);
}

void SVHB::createFvec()
{
  int i, j, jini;
  int findex;

  if (!regrowth)
	{
    if (svHBopt->numtones == 1)
		{
      svHBdata->NoFreqPoints = svHBopt->h[0]+1;
      // allocate memory for the frequency vector
      svHBdata->FreqV_P = Mlib_DNewVec(svHBdata->NoFreqPoints);
      // fill vector
      for (findex=0; findex < svHBdata->NoFreqPoints; findex++)
				svHBdata->FreqV_P[findex] = svHBopt->tone[0] * findex;
    }
    else
		{
      // numtones=2, map artificial frequencies.
      svHBdata->NoFreqPoints = (svHBopt->h[0] + 1)*(2*svHBopt->h[1] + 1)
			- svHBopt->h[1];
      // allocate memory for the frequency vector
      svHBdata->FreqV_P = Mlib_DNewVec(svHBdata->NoFreqPoints);

      for (i=0, findex=0; i<= svHBopt->h[0]; i++)
			{
				jini = (i==0) ? 0 : -(svHBopt->h[1]);
				for (j=jini; j<= svHBopt->h[1]; findex++, j++)
					svHBdata->FreqV_P[findex] =
				fabs(svHBopt->tone[0] * i + svHBopt->tone[1] * j);
      }
    }
  }
  else
	{
    // Create vector using method from Borges de Carvalho and Pedro.
    // First calculate bandwidths
    int bw1 = n_freqs * (n_fund - 1) + 1;
    int bw2 = bw1 - n_fund + 1;
    int ndc = (n_freqs % 2) ? (bw2 + 1) / 2 : (bw1 + 1) / 2;
    double fund_center = fundamental + f_step * (n_fund>>1);
    int max_nf = n_freqs;
    if (isSet(&max_cfreq))
		{
      max_nf = int(max_cfreq / fund_center);
      if (!max_nf || max_nf > n_freqs)
			{
				max_nf = n_freqs;
				report(WARNING, "max_cfreq parameter ignored");
      }
    }
    int nbw1 = (max_nf + 1) / 2;
    int nbw2 = max_nf / 2;

    svHBdata->NoFreqPoints = ndc + bw1 * nbw1 + bw2 * nbw2;
    // allocate memory for the frequency vector
    svHBdata->FreqV_P = Mlib_DNewVec(svHBdata->NoFreqPoints);

    // Fill vector
    for (findex = 0; findex < ndc; findex++)
      svHBdata->FreqV_P[findex] = f_step * findex;
    int nf=1;
    int cw = (n_freqs % 2) ? bw1 : bw2;
    while (nf <= max_nf) {
      double fcenter = fund_center * nf;
      for (int nsf = -cw/2; nsf <= cw/2; nsf++)
				svHBdata->FreqV_P[findex++] = fcenter + f_step * nsf;
      cw = (cw == bw1) ? bw2 : bw1;
      nf++;
    }
    assert(findex == svHBdata->NoFreqPoints);
  }
}

doublev_t SVHB::createOmega(int *n)
{
  int i, j, jini, nsamples;
  int findex, NoFreq, order, n_min_freq;
  doublev_t omega;

  NoFreq = 0;
  if (svHBopt->numtones == 1)
	{
		/* Search the minimum power of 2 greater than svHBdata->NoFreqPoints */
    for (order=1; NoFreq < svHBopt->h[0]+1; order++)
      NoFreq = (1 << order);

    NoFreq *= svHBopt->oversample;
    nsamples = 2 * NoFreq;

    /* allocate memory for the vector */
    omega = Mlib_DNewVec(NoFreq);

    /* fill vector */
    for (findex=0; findex < NoFreq; findex++)
		{
      omega[findex] = svHBopt->tone[0] * twopi * findex;
		}
  }
  else if (svHBopt->numtones == 2)
	{
    /* Search the minimum power of 2 greater than the number of
		* freq. required for rectangular truncation.
		*/
    n_min_freq = (svHBopt->h[0] + 1)*(2*svHBopt->h[1] + 1)
		- svHBopt->h[1];
    for (order=1; NoFreq < n_min_freq; order++)
      NoFreq = (1 << order);

    nsamples = 2 * NoFreq;

    /* allocate memory for the vector */
    omega = Mlib_DNewVec(NoFreq);

    for (i=0, findex=0; findex < NoFreq; i++)
		{
      jini = (i==0) ? 0 : -(svHBopt->h[1]);
      for (j=jini; j<= svHBopt->h[1]; findex++, j++)
			{
				if (findex == NoFreq)
					break;
				omega[findex] = fabs(svHBopt->tone[0] * i + svHBopt->tone[1] * j)
				* twopi;
      }
    }
  }
  else
	{
    // Assume regrowth
    assert(regrowth);
    n_min_freq = svHBdata->NoFreqPoints;
    for (order=1; NoFreq < n_min_freq; order++)
      NoFreq = (1 << order);
    nsamples = 2 * NoFreq;
    /* allocate memory for the vector */
    omega = Mlib_DNewVec(NoFreq);
    for (findex = 0; findex < n_min_freq; findex++)
      omega[findex] = svHBdata->FreqV_P[findex] * twopi;
    // put dummy frequencies in the rest of the vector
    for (findex = n_min_freq; findex < NoFreq; findex++)
      omega[findex] = (fundamental * (n_freqs+1) + f_step * findex) * twopi;
  }

  *n = nsamples;
  return(omega);
}

