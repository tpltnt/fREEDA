#include "SVTr.h"
using std::cout;
using std::endl;

#include "DC.h"
#include "FreqMNAM.h"
#include "Euler.h"
#include "TimeDomainSV.h"
#include "rfftw.h"

extern "C"
{
#include "../inout/ftvec.h"
#include "../inout/report.h"
}

// NOX objects
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

#include "Nox_Interface.h"
#include "Teuchos_ParameterList.hpp"

void buildTIncidence(ElemFlag mask, Circuit*& my_circuit,
    IntDenseMatrix& T, ElementVector& elem_vec,
    int& n_states, int& max_n_states);

// Static members
const unsigned SVTr::n_par = 15;

// Element information
ItemInfo SVTr::ainfo =
{
  "svtr",
  "State-Variable-Based Convolution Transient Analysis",
  "Carlos E. Christoffersen, Mete Ozkar",
  DEFAULT_ADDRESS
};

// Parameter information
ParmInfo SVTr::pinfo[] =
{
  {"verbosity", "Output verbosity", TR_INT, false},
  {"tstop", "Stop time (s)", TR_DOUBLE, true},
  {"tstep", "Time step size (s)", TR_DOUBLE, true},
  {"n_freqs", "Number of frequencies considered", TR_INT, false},
  {"gcomp", "Compensation conductance (mhos)", TR_DOUBLE, false},
  {"tolerance", "Tolerance used to evaluate errors", TR_DOUBLE, false},
  {"filter_freq", "Corner frequency for artificial filtering (Hz)",
	TR_DOUBLE, false},
  {"n_samples", "Number of time samples for impulse response", TR_INT, false},
  {"ntest", "Number of points used to test the impulse response",
	TR_INT, false},
  {"imp_tol", "Contribution percentage of last ntest samples of the impulse response",
	TR_DOUBLE, false},
  {"check_imp", "Flag to signal if impulse checking is desired", TR_BOOLEAN,
	false},
  {"out_steps", "Number of steps skipped for output simulation progress",
	TR_INT, false},
  {"opt", "Start transient assuming a biased circuit", TR_BOOLEAN, false},
  {"adjust", "Adjust the impulse response to produce a correct DC level",
	TR_BOOLEAN, false},
  {"deriv", "Approximate derivatives or use automatic diff.", TR_INT, false}
};


SVTr::SVTr() : Analysis(&ainfo, pinfo, n_par), mnam_mask(TR_FREQ_DOMAIN),
nle_mask(TR_TIME_DOMAIN)
{
  // Parameter stuff
  paramvalue[0] = &(verbosity = 1);
  paramvalue[1] = &tstop;
  paramvalue[2] = &tstep;
  paramvalue[3] = &(n_freqs_net = 1024);
  paramvalue[4] = &(gcomp = 0.005);
  paramvalue[5] = &(tolerance = 1e-8);
  paramvalue[6] = &(filter_freq = 0);
  paramvalue[7] = &(n_samples = 0);
  paramvalue[8] = &(ntest = 0);
  paramvalue[9] = &(imp_tol = one);
  paramvalue[10] = &(check_imp = true);
  paramvalue[11] = &(out_steps = 100);
  paramvalue[12] = &(opt = false);
  paramvalue[13] = &(adjust = false);
  paramvalue[14] = &(deriv = 0);
}

SVTr::~SVTr()
{
}

void SVTr::run(Circuit* cir)
{
  this->cir = cir;

  // Char array to store output messages
  char msg[80];
  report(MESSAGE, "*** State Variable-Based Transient Analysis ***");
  sepLine();
  newLine();

  // ----------------------------------------------------------------------
  // -------------- Set up and check parameters
  // ----------------------------------------------------------------------
  checkParams();

  // Add one to n_freqs so the value in the netlist is a more intuitive
  // power of two.
  n_freqs = n_freqs_net + 1;

  // Number of time samples
  if (!isSet(&n_samples))
    n_samples = 7 * (n_freqs - 1) / 4;
  else if (n_samples > 2 * (n_freqs - 1))
    n_samples = 2 * (n_freqs - 1);

  sprintf(msg, "n_samples = %d", n_samples);
  report(MESSAGE, msg);

  // Check the last 50th of the impulse response by default
  if (!isSet(&ntest))
    ntest = n_samples / 50;

  // ----------------------------------------------------------------------
  // -------------- Set up time and frequency vectors.
  // ----------------------------------------------------------------------

  // Create frequency vector
  fstep = .5 / tstep / n_freqs;
  f_vec.resize(n_freqs);
  for (int i = 0; i < n_freqs; i++)
    f_vec[i] = fstep * i;

  sprintf(msg, "Frequency step = %g Hz", fstep);
  report(MESSAGE, msg);
  sprintf(msg, "Maximum frequency = %g Hz", f_vec[n_freqs-1]);
  report(MESSAGE, msg);

  // Number of time steps
  n_tsteps = int(tstop / tstep) + 1;
  // Fill the time vector
  t_vec.resize(n_tsteps);
  for (int i = 0; i < n_tsteps; i++)
    t_vec[i] = tstep * i;

  sprintf(msg, "Number of time steps = %d", n_tsteps);
  report(MESSAGE, msg);

  // ----------------------------------------------------------------------
  // -------------- Build one Msv matrix at each frequency.
  // ----------------------------------------------------------------------
  buildMsv();
  sprintf(msg, "Number of state variables = %d", n_states);
  report(MESSAGE, msg);

  // ----------------------------------------------------------------------
  // -------------- Calculate DC bias if requested
  // ----------------------------------------------------------------------
  // Initial condition vectors.
  DenseDoubleVector X0(n_states), U0(n_states), I0(n_states);
  if (opt || adjust)
	{
    // Run a DC analysis
    DC* dca = new DC();
    double res = dca->getBias(cir, mnam_mask, nle_mask, X0, U0, I0);
    sprintf(msg, "DC analysis residual = %g", res);
    report(MESSAGE, msg);
    delete dca;
  }

	// ----------------------------------------------------------------------
  // -------------- Filter Msv in the frequency domain.
  // ----------------------------------------------------------------------
  filterMsv();

  // ----------------------------------------------------------------------
  // -------------- Convert Msv to the time domain.
  // ----------------------------------------------------------------------
  iDFTMsv();

  // ----------------------------------------------------------------------
  // ---------- Allocate memory for the simulation.
  // ----------------------------------------------------------------------

  // resize the convolution accumulator
  p_conv.resize(n_states);
  // Vector to hold DC contribution
  DenseDoubleVector DC_conv(n_states);

  // Allocate memory for the state variable vector
  double* X_P = new double[n_states];
  DenseDoubleVector X(Teuchos::View, X_P, n_states);

  // These matrices hold the result of the simulation.
  // MI_NL is also used in the convolution.
  MI_NL.reshape(n_tsteps, n_states);
  MU_NL.reshape(n_tsteps, n_states);

  // Circular lists
  const int list_size = 30;
  cX = new CircVector(list_size, n_states);
  cY = new CircVector(list_size, n_sec_states);
  cU = new CircVector(list_size, n_states);
  cI = new CircVector(list_size, n_states);
  // Create nonlinear integration method class
  NLIntegMethod *nl_im = new NLEuler(cX, cY, tstep);
  // Create the interface class
  tdsv = new TimeDomainSV(nl_im, max_n_states);

  // The residual is the 2-norm of the resiual vector.
  double residual = zero;

  // ----------------------------------------------------------------------
  // ----------------- Run the transient
  // ----------------------------------------------------------------------

  // scale the impulse response to produce the correct DC level if
  // requested.
  if (adjust)
	{
    double DC_acc;
    // Initialize DC contribution to convolution
    // Also rescale impulse response to get the correct DC value
    double imp_error = zero;
    for (int i = 0; i < n_states; i++)
		{
      DC_conv[i] = U0[i];
      DC_acc = zero;
      // First calculate contribution from this row
      for (int j = 0; j < n_states; j++)
			{
				double acc = zero;
				for (int ntau = 0; ntau < n_samples; ntau++)
					acc += Msv[i * n_states + j][ntau];
				acc *= (- I0[j] + U0[j] * gcomp);
				DC_acc += acc;
      }
      // Now rescale the row of Msv
      double c_factor = U0[i] / DC_acc;
      if (c_factor > .5)
			{
				double error = abs(one - c_factor);
				if (error > imp_error)
					imp_error = error;
				for (int j = 0; j < n_states; j++)
					for (int ntau = 0; ntau < n_samples; ntau++)
						Msv[i * n_states + j][ntau] *= c_factor;
      }
      else
			{
				// Since the final value should be close to zero, compensate
				// by adding a correcting term
				int valid_states = 0;
				for (int j = 0; j < n_states; j++)
				{
					double div = - I0[j] + U0[j] * gcomp;
					// Perform correction only if the divider is nonzero
					if (div > tolerance / 10.)
						valid_states++;
				}
				double c_term = (U0[i] - DC_acc) / valid_states / (n_samples - 1);
				for (int j = 0; j < n_states; j++)
				{
					double div = - I0[j] + U0[j] * gcomp;
					// Perform correction only if the divider is nonzero
					if (div > tolerance / 10.0)
						for (int ntau = 1; ntau < n_samples; ntau++)
							Msv[i * n_states + j][ntau] += c_term / div;
				}
      }
    }
    sprintf(msg, "Impulse response DC error: %4.2f %%", 100.0*imp_error);
    report(MESSAGE, msg);
  }
  else
	{
    for (int i = 0; i < n_states; i++)
		{
      // Calculate contribution from this row
      for (int j = 0; j < n_states; j++)
			{
				double acc = zero;
				for (int ntau = 0; ntau < n_samples; ntau++)
					acc += Msv[i * n_states + j][ntau];
				acc *= (- I0[j] + U0[j] * gcomp);
				DC_conv[i] += acc;
      }
    }
  }

  // Initialize variables
  if (opt)
	{
    // Just be sure that the DC analysis is doing what we want.
    assert(int(X0.length()) == n_states);
    X = X0;
    // Also initialize the circular vectors!
    for (int j = 0; j < list_size; j++)
		{
      for (int i = 0; i < X0.length(); i++)
        cX->getCurrent()[i] = X0[i];
      cX->advance();
      for (int i = 0; i < DC_conv.length(); i++)
        cU->getCurrent()[i] = DC_conv[i];
      cU->advance();
      for (int i = 0; i < I0.length(); i++)
        cI->getCurrent()[i] = I0[i];
      cI->advance();
    }
    // Initialize MI_NL(0,:)
    for (int j = 0; j < n_states; j++)
		{
      MI_NL(0,j) = - I0[j] + DC_conv[j] * gcomp;
      MU_NL(0,j) = DC_conv[j];
    }
  }
  else
	{
    for (int i = 0; i < X.length(); i++)
      X[i] = zero;
		for (int k = 0; k <= 0; k++)
		{
			for (int l = 0; l <= n_states-1; l++)
			{
				MI_NL(k,l) = zero;
				MU_NL(k,l) = zero;
			}
		}
  }

  SerialCommunicator Comm;

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // Get the number of elements from the command line
  int NumGlobalElements = n_states;

  // The number of unknowns must be at least equal
  // to the number of processors
  if (NumGlobalElements < NumProc)
  {
    cout << "numGlobalBlocks = " << NumGlobalElements << 
      " cannot be < number of processors = " << NumProc << endl;
    cout << "Test failed!" << endl;
    throw "NOX Error";
  }

  // ********** Begin Setup Nonlinear Solver *******************
  // Create the interface between NOX and the application
  // This object is derived from NOX::Epetra::Interface
  Teuchos::RCP<Nox_Interface> interface = 
    Teuchos::rcp(new Nox_Interface(this, NumGlobalElements, Comm));
  
  // Get the vector from the Problem
  Teuchos::RCP<DistributedDoubleVector> InitialGuess = interface->getSolution();
 
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

  // Update time in interface class
  tdsv->setTime(&(cX->getCurrent()[0]), &(cU->getCurrent()[0]),
    &(cI->getCurrent()[0]), 0, 0.0, tstep);
    
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
  Teuchos::RCP<NOX::Solver::Generic> solver = 
    NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);
  
  NOX::StatusTest::StatusType solvStatus;
  // this is the guess vector which will be updated every time step
  NOX::Epetra::Vector* guess = &noxInitGuess;
  const NOX::Epetra::Group* soln_group = 0;
  const Epetra_Vector* next_guess = 0;
  // *********** End Setup Nonlinear Solver *********************
  
  report(MESSAGE, "--- Starting transient simulation ...\n");
  sepLine();
  report(MESSAGE, "|   Step   |    Time (s)    |   Residual (V)   |");
  sepLine();
  
  for (nt = 1; nt < n_tsteps; nt++)
  {
    // Perform main part of the convolution sum
    // Clear convolution accumulator
    for (int i = 0; i < p_conv.length(); i++)
      p_conv[i] = zero;

    // nmax is the upper limit of the summation
    int nmax;
    if (nt < n_samples)
		{
      nmax = nt;
      // Add DC contribution if requested
      if (opt)
			{
				for (int i = 0; i < n_states; i++)
				{
					for (int j = 0; j < n_states; j++)
					{
						DC_conv[i] -= Msv[i * n_states + j][nt - 1] * MI_NL(0, j);
					}
					p_conv[i] = DC_conv[i];
				}
      }
    }
    else
      nmax = n_samples;

    for (int i = 0; i < n_states; i++)
		{
      for (int j = 0; j < n_states; j++)
			{
				for (int ntau = 1; ntau < nmax; ntau++)
				{
					p_conv[i] += Msv[i * n_states + j][ntau] * MI_NL(nt - ntau, j);
				}
      }
    }

    // Update time in interface class
    tdsv->setTime(&(cX->getCurrent()[0]),	&(cU->getCurrent()[0]),	&(cI->getCurrent()[0]),
        nt, nt * tstep, tstep);
    
    //*** NOX Solver Start****
    // Reset the solver so it isn't converged anymore and pass the 
    // result from the previous time step
    solver->reset(*guess);
    // Solve the nonlinear equation (call the nonlinear solver).
    solvStatus = interface->evaluate(solver);
    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& soln_group = 
      dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    const NOX::Epetra::Vector& next_guess = 
      (dynamic_cast<const NOX::Epetra::Vector&>(soln_group.getX()));
    residual = soln_group.getNormF();
    *guess = next_guess;
    //*******NOX Solver End******
    
    // Save variables.
    // Recalculate U and I since there is no warranty from
    // the nonlinear solver that the values are updated.
    guess->getEpetraVector().ExtractCopy(&X_P);
    updateUI(X);
    const DenseDoubleVector& I_nl = cI->getCurrent();
    const DenseDoubleVector& U_nl = cU->getCurrent();
    for (int i=0; i < n_states; i++)
		{
      MU_NL(nt, i) = U_nl[i];
      // The values of the current need to be compensated since
      // we added the compensation resistor in the freq. domain.
      MI_NL(nt, i) = - I_nl[i] + U_nl[i] * gcomp;
    }
    // Advance the circular lists for the next step.
    cX->advance();
    cU->advance();
    cI->advance();

    if (!(nt % out_steps))
		{
      sprintf(msg, "|  %6d  |  %e  |  %e *  |", nt, t_vec[nt], sqrt(residual));
      report(MESSAGE, msg);
    }
  }
  residual = sqrt(residual);

  sprintf(msg, "--- Residual: %g", residual);
  sepLine();
  newLine();
  report(MESSAGE, msg);
  
  // Write output
  doOutput();

  // ----------------------------------------------------------------------
  // ------------- Free unused memory.
  // ----------------------------------------------------------------------
  delete tdsv;
  delete nl_im;
  delete cI;
  delete cU;
  delete cX;
  delete [] X_P;
  delete [] Msv;
}


void SVTr::buildMsv()
{
  newLine();
  report(MESSAGE, "--- Building Msv(f) ...");

  FreqMNAM* mnam = new FreqMNAM(f_vec, cir, mnam_mask);

  // Get T matrix
  buildTIncidence(nle_mask, cir, T, elem_vec, n_states, max_n_states);

  // Go through all the nonlinear elements
  int n_elem = elem_vec.size();
  cns.resize(n_elem);
  nss.resize(n_elem);
  n_sec_states = 0;
  for (int k=0; k < n_elem; k++)
	{
    cns[k] = elem_vec[k]->getNumberOfStates();
    nss[k] = elem_vec[k]->getNumberOfSecStates();
    n_sec_states += nss[k];
  }

  // Allocate memory for Msv. It is implemented as an array of
  // vectors, so we can quickly perform the IFFT of each vector.
  Msv = new DenseDoubleVector[n_states * n_states];
	for (int i = 0; i < n_states*n_states; i++)
		Msv[i] = DenseDoubleVector(n_freqs << 1);
  
  // Temporary vectors (do we need them?)
  DenseComplexVector tempvec1(mnam->getDim());
  DenseComplexVector tempvec2(mnam->getDim());

  // Add compensation admittance to the MNAMs at all frequencies
  for (int findex = 0; findex < n_freqs; findex++)
  {
    mnam->setFreqIndex(findex);
    for (int i = 0; i < n_states; i++)
      mnam->setAdmittance(T(0,i), T(1,i), double_complex(gcomp));
  }
  
  mnam->FillComplete();

  // Loop at all the frequencies
  for (int findex=0; findex < n_freqs; findex++)
	{
    // Solve for each row of T and
    for (int i=0; i < n_states; i++)
		{
      // Set the 1 and -1 in tempvec tempvec1
      if (T(0,i))
        tempvec1[T(0,i)-1].real() += one;
      if (T(1,i))
        tempvec1[T(1,i)-1].real() -= one;
      
      mnam->solve(findex, tempvec1, tempvec2);

      // Restore tempvec1 to zero
      if (T(0,i))
        tempvec1[T(0,i)-1].real() = zero;
      if (T(1,i))
        tempvec1[T(1,i)-1].real() = zero;

      // Multiply T times tempvec2
      for (int j=0; j < n_states; j++)
			{
				double_complex c1(zero);
				if (T(0,j))
					c1 = tempvec2[T(0,j) - 1];
				if (T(1,j))
					c1 -= tempvec2[T(1,j) - 1];

				Msv[n_states * i + j][2 * findex] = c1.real();
				Msv[n_states * i + j][2 * findex + 1] = c1.imag();
      }
    }
  }

  // Since the MNAM is no longer used, free the memory associated
  delete mnam;
}


void SVTr::filterMsv()
{
  if (filter_freq > one)
	{
    int si = int(filter_freq / fstep);

    if (si < n_freqs)
		{
      report(MESSAGE, "--- Filtering ...");

      for(int i = si + 1; i < n_freqs; i++)
			{
				double factor = double(n_freqs - i) / (n_freqs - si);
				factor *= factor;
				for(int j = 0; j < n_states; j++)
				{
					for(int k = 0; k < n_states; k++)
					{
						Msv[n_states * j + k][2 * i] *= factor;
						Msv[n_states * j + k][2 * i] +=
						Msv[n_states * j + k][0] * (one - factor);
						Msv[n_states * j + k][2 * i + 1] *= factor;
					}
				}
      }
    }
  }
}


void SVTr::iDFTMsv()
{
  char msg[80];

  report(MESSAGE, "--- Converting Msv to the time domain ...");

  int nfft = 2 * (n_freqs - 1);

  // Allocate space for a temporary vector to do the transformation
  double* fd_p = new double[nfft];
  double* td_p = new double[nfft];

  rfftw_plan p = rfftw_create_plan(nfft,
	FFTW_COMPLEX_TO_REAL,
	FFTW_MEASURE);

  // Transform each matrix element
  for (int j = 0; j < n_states; j++)
	{
    for (int i = 0; i < n_states; i++)
		{
      // Copy the frequency elements to a double vector with the
      // order used by fftw (halfcomplex array).
      fd_p[0] = Msv[n_states * i + j][0];
      fd_p[nfft>>1] = Msv[n_states * i + j][2 * n_freqs - 2];
      for (int findex=1; findex < ((nfft+1) / 2); findex++)
			{
				fd_p[findex] = Msv[n_states * i + j][2 * findex];
				fd_p[nfft - findex] = Msv[n_states * i + j][2 * findex + 1];
      }

      // Perform the FFT
      rfftw_one(p, fd_p, td_p);

      // Copy the time domain vector in the convetional way
      double mDC = zero;
      for (int tindex = 0; tindex < n_samples; tindex++)
			{
				Msv[n_states * i + j][tindex] = td_p[tindex] / nfft;
				mDC += abs(td_p[tindex]);
      }

      if (check_imp)
			{
				// Check that the last part of the impulse response is small enough.
				double res = zero;
				const double itol = imp_tol / 100.;
				for (int tindex = n_samples - ntest; tindex < n_samples; tindex++)
				{
					res += abs(td_p[tindex]);
				}
				res *= (n_samples / mDC / ntest);
				if (abs(res) > itol)
				{
					sprintf(msg, "Last %d samples impulse contribution: %4.2f %%",
					ntest, abs(res) * 100);
					report(WARNING, msg);
				}
      }
    }
  }
  // Destroy the plan to save memory.
  rfftw_destroy_plan(p);

  delete [] fd_p;
  delete [] td_p;
}


void SVTr::func_ev(double * X, double * errFunc)
{
  // Evaluate first the time-domain elements.
  DenseDoubleVector X_vec(Teuchos::View, X, n_states);
  updateUI(X_vec);
  // Tell tdsv that the first evaluation is already completed.
  tdsv->clearflag();

  const DenseDoubleVector& I_nl = cI->getCurrent();
  const DenseDoubleVector& U_nl = cU->getCurrent();
  // The values of the current need to be compensated since
  // we added the compensation resistor in the freq. domain.
  for (int i = 0; i < n_states; i++)
    MI_NL(nt, i) = -I_nl[i] + U_nl[i] * gcomp;

  // Complete the convolution sum with the last term
  for (int i = 0; i < n_states; i++)
	{
    double alpha = p_conv[i];
    for (int j = 0; j < n_states; j++)
		{
      // Since the result from the FFT is already divided by two,
      // we do not need the factor for the trapezoidal rule.
      alpha += Msv[i * n_states + j][0] *  MI_NL(nt, j);
    }
    errFunc[i] = alpha - U_nl[i];
  }
}

void SVTr::updateUI(const DenseDoubleVector& X)
{
  // Set current state variable values.
	assert(X.length() == cX->getCurrent().length());
  int len = X.length();
  for (int i = 0; i < len; i++)
    cX->getCurrent()[i] = X[i];

  // Number of elements
  int n_elem = elem_vec.size();
  // Go through all the nonlinear elements
  for (int i=0, j1=0, k=0; k < n_elem; k++)
	{
    // Set base index in interface object
    tdsv->setBase(i, cns[k], j1, nss[k]);
    // Call element evaluation
    elem_vec[k]->svTran(tdsv);

    i += cns[k];
    j1 += nss[k];
  }
}

void SVTr::jacobian(double * X, DoubleSparseColMatrix& Jacobian)
{
  DoubleDenseMatrix J(n_states, n_states);

  // Evaluate first the time-domain elements.
  DenseDoubleVector X_vec(Teuchos::View, X, n_states);
  // Set current state variable values.
	assert(X_vec.length() == cX->getCurrent().length());
  for (int i = 0; i < X_vec.length(); i++)
    cX->getCurrent()[i] = X_vec[i];

  // Number of elements
  int n_elem = elem_vec.size();
  DoubleDenseMatrix& Ju = tdsv->getJu();
  DoubleDenseMatrix& Ji = tdsv->getJi();
  J.putScalar(zero);
  // Go through all the nonlinear elements
  for (int i=0, j1=0, k=0; k < n_elem; k++)
	{
    // Set base index in interface object
    tdsv->setBase(i, cns[k], j1, nss[k]);
    // Call element evaluation
    elem_vec[k]->deriv_svTran(tdsv);
    // correct Ji with gcomp
    for (int j=0; j < cns[k]; j++)
			for (int n=0; n < cns[k]; n++)
		{
			Ji(j, n) = Ju(j,n) * gcomp - Ji(j,n);
			J(i+j, i+n) = - Ju(j, n);
		}
    for (int l=0; l < n_states; l++)
      for (int j=0; j < cns[k]; j++)
				for (int n=0; n < cns[k]; n++)
					J(l, i+j) += Msv[l * n_states + i+n][0] * Ji(n, j);

    i += cns[k];
    j1 += nss[k];
  }
  
  // copy contents of dense jacobian J into a sparse equivalent Jacobian
  // required by NOX
  std::vector<int> indices(n_states, 0);
  std::vector<double> values(n_states, 0.0);

  int k;
  // row loop
  for (int i = 0; i < n_states; i++)
  {
    k = 0;
    // column loop
    for (int j = 0; j < n_states; j++)
    {
      if (J(i,j) != 0.0)
      {
        // we have a non-zero value, so save it
        values[k] = J(i,j);
        indices[k] = j; // save the column the non zero value was found in
        k++;
      }
    }
    if (k != 0)
    {
      // non zero values were found in row i
      Jacobian.ReplaceMyValues(i, k, &values[0], &indices[0]);
    }
  }
}

void SVTr::doOutput()
{
  // First check if the result matrices contain any data
  assert(MI_NL.numCols() == n_states);
  assert(MU_NL.numRows() == n_tsteps);

  report(MESSAGE, "--- Writing output vectors ...");

  // Allocate and copy global time vector
  allocTimeV_P(n_tsteps);
  for (int tindex = 0; tindex < n_tsteps; tindex++)
    TimeV_P[tindex] = t_vec[tindex];

  // Temporary vectors
  DenseDoubleVector tmp_i(n_tsteps);
  DenseDoubleVector tmp_u(n_tsteps);

  // Now fill currents and port voltages of time domain devices.
  int current_ns;
  int n_elem = elem_vec.size();   // Number of elements
  int i=0;
  for (int k=0; k < n_elem; k++)
	{
    current_ns = elem_vec[k]->getNumberOfStates();
    for(int j=0; j< current_ns; j++)
		{
      // For the current, decompensate while copying.
      for (int tindex=0; tindex < n_tsteps; tindex++)
			{
				tmp_i[tindex] = - MI_NL(tindex, j+i) + MU_NL(tindex, j+i) * gcomp;
				tmp_u[tindex] = MU_NL(tindex, j+i);
      }
      elem_vec[k]->getElemData()->setRealI(j, tmp_i);
      elem_vec[k]->getElemData()->setRealU(j, tmp_u);
    }
    i += current_ns;
  }

  // Once the output is made, discard the memory used by the matrices
  MI_NL.reshape(0,0);
  MU_NL.reshape(0,0);
}

