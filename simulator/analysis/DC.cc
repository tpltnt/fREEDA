#include "DC.h"
#include "FreqMNAM.h"
#include "TimeDomainSV.h"

extern "C"
{
	#include "../inout/ftvec.h"
	#include "../inout/report.h"
}

#include <cstdio>

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
 
#include "Nox_Interface.h"
#include "Teuchos_ParameterList.hpp"

void buildTIncidence(ElemFlag mask, Circuit*& my_circuit,
IntDenseMatrix& T, ElementVector& elem_vec,
int& n_states, int& max_n_states);

// Static members
const unsigned DC::n_par = 6;

// Element information
ItemInfo DC::ainfo =
{
  "dc",
  "State Variable DC Analysis",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS
};

// Parameter information
ParmInfo DC::pinfo[] =
{
  {"sweep", "Vsource to be swept", TR_STRING, false},
  {"start", "Start sweep value (V)", TR_DOUBLE, false},
  {"stop", "Stop sweep value (V)", TR_DOUBLE, false},
  {"step", "Sweep step (V)", TR_DOUBLE, false},
  {"gcomp", "Compensation network conductance (S)", TR_DOUBLE, false},
  {"deriv", "Use numerical derivatives", TR_INT, false},
};


DC::DC() : Analysis(&ainfo, pinfo, n_par)
{
  paramvalue[0] = &(sweep = "none");
  paramvalue[1] = &(vstart = zero);
  paramvalue[2] = &(vstop = zero);
  paramvalue[3] = &(vstep = zero);
  paramvalue[4] = &(gcomp = 1e-4);
  paramvalue[5] = &(deriv = 0);

  cir = NULL;
  mnam = NULL;
  tdsv = NULL;
  allocflag = false;
}


DC::~DC()
{
  if (allocflag)
    clean();
}

void DC::run(Circuit* cir)
{
  // Char array to store output messages
  char msg[80];
  report(MESSAGE, "*** State Variable-Based DC Analysis ***");
  sepLine();
  newLine();

  DenseComplexVector sourcev;
  DenseComplexVector result;

  if (isSet(&sweep))
  {
    report(MESSAGE, "I hope your DC parameters make sense (no checking made)");
    nsteps = int((vstop-vstart)/vstep) + 1;
    Usweep = new DenseDoubleVector[nsteps];
    Xsweep = new DenseDoubleVector[nsteps];
    Vsweep = new DenseDoubleVector[nsteps];
    Isweep = new DenseDoubleVector[nsteps];
    // Perform voltage source sweep
    double residual =
      getBias(cir, LINEAR, NONLINEAR, Xsweep[0], Vsweep[0], Isweep[0]);
    sprintf(msg, "Number of state variables = %d", n_states);
    report(MESSAGE, msg);
    sprintf(msg, "Residual = %g", residual);
    report(MESSAGE, msg);
    // Build full source vector (sourcev1)
    mnam->getSource(0, sourcev);
    DenseComplexVector sourcev1(sourcev.length());
    sourcev1 = sourcev;

    // Add contributions from nonlinear devices
    // Remember current has sign of passive devices
    for (int j=0; j < n_states; j++)
    {
      double itmp= Iw[j] - gcomp * Vw[j];
     	if (T(0,j))
        sourcev1[T(0,j)-1].real() -= itmp;
     	if (T(1,j))
        sourcev1[T(1,j)-1].real() += itmp;
    }
    // resize nodal vector.
    Usweep[0].resize(mnam->getDim());
    result.resize(mnam->getDim());
    // Solve for sourcev1 vector
    mnam->solve(0, sourcev1, result);
    for (unsigned j=0; j < mnam->getDim(); j++)
      Usweep[0][j] = result[j].real();

    report(MESSAGE, "--- Starting DC simulation ...\n");
    sepLine();
    report(MESSAGE, "|   Step  (V)    |   Residual (V) |");
    sepLine();
    for (int i=1; i < nsteps; i++)
    {
      double cvolt = vstart + vstep * i;
      // vsource->setParam("vdc", &cvolt, TR_DOUBLE);
      // Remember that MNAM indices begin at 0
      sourcev[nrow-1].real() = cvolt; 
      // Initial guess is last result
      double residual = solveDC(sourcev);
			Xsweep[i].resize(Xw.length());
      Xsweep[i] = Xw;
			Vsweep[i].resize(Vw.length());
      Vsweep[i] = Vw;
			Isweep[i].resize(Iw.length());
      Isweep[i] = Iw;
      sprintf(msg, "|  %e  |  %e  |", cvolt, residual);
      report(MESSAGE, msg);
      // Save results and use as a guess for next iteration
      // Build full source vector (sourcev1)
      sourcev1 = sourcev;
      // Add contributions from nonlinear devices
      // Remember current has sign of passive devices
      for (int j=0; j < n_states; j++)
      {
       	double itmp = Iw[j] - gcomp * Vw[j];
       	if (T(0,j))
          sourcev1[T(0,j)-1].real() -= itmp;
       	if (T(1,j))
          sourcev1[T(1,j)-1].real() += itmp;
      }
      // resize nodal vector.
      Usweep[i].resize(mnam->getDim());
      // Solve for sourcev1 vector
      mnam->solve(0, sourcev1, result);
      for (unsigned j=0; j < mnam->getDim(); j++)
				Usweep[i][j] = result[j].real();
    }
    sepLine();
    newLine();
  }
  else
  {
    nsteps = 1;
    Usweep = new DenseDoubleVector[nsteps];
    Xsweep = new DenseDoubleVector[nsteps];
    Vsweep = new DenseDoubleVector[nsteps];
    Isweep = new DenseDoubleVector[nsteps];
    double residual =
      getBias(cir, LINEAR, NONLINEAR, Xsweep[0], Vsweep[0], Isweep[0]);
    sprintf(msg, "Number of state variables = %d", n_states);
    report(MESSAGE, msg);
    sprintf(msg, "Residual = %g", residual);
    report(MESSAGE, msg);
    // Build full source vector (sourcev1)
    mnam->getSource(0, sourcev);
    DenseComplexVector sourcev1(sourcev.length());
    sourcev1 = sourcev;

    // Add contributions from nonlinear devices
    // Remember current has sign of passive devices
    for (int j=0; j < n_states; j++)
    {
      double itmp= Iw[j] - gcomp * Vw[j];
     	if (T(0,j))
        sourcev1[T(0,j)-1].real() -= itmp;
     	if (T(1,j))
        sourcev1[T(1,j)-1].real() += itmp;
    }
    // resize result vector.
    Usweep[0].resize(mnam->getDim());
    result.resize(mnam->getDim());
    // Solve for sourcev1 vector
    mnam->solve(0, sourcev1, result);
    for (unsigned j=0; j < mnam->getDim(); j++)
      Usweep[0][j] = result[j].real();
  }

  // Write output
  doOutput();
  // Clean space
  clean();
  if (nsteps)
  {
    delete [] Isweep;
    delete [] Vsweep;
    delete [] Xsweep;
    delete [] Usweep;
  }
}

double DC::getBias(Circuit* cir,
		   const ElemFlag& mnam_mask, const ElemFlag& nonlinear_mask,
		   DenseDoubleVector& X, DenseDoubleVector& U, DenseDoubleVector& I)
{
  char msg[80];
  assert(!allocflag);
  this->cir = cir;

  // ----------------------------------------------------------------------
  // -------------- Set up and check parameters
  // ----------------------------------------------------------------------

  // checkParams();
  // Create frequency vector
  DenseDoubleVector f_vec(1);

  // ----------------------------------------------------------------------
  // -------------- Build one Msv matrix
  // ----------------------------------------------------------------------
  // Create the MNAM
  mnam = new FreqMNAM(f_vec, cir, mnam_mask);
  mnam->setFreqIndex(0);

  // Get T matrix
  buildTIncidence(nonlinear_mask, cir, T, elem_vec, n_states, max_n_states);
  if (gcomp) // Add the compensation resistors to the MNAM
    for (int i = 0; i < n_states; i++)
      mnam->setAdmittance(T(0,i), T(1,i), gcomp);
  
  // this step is required
  mnam->FillComplete();

  // Go through all the nonlinear elements
  int n_elem = elem_vec.size();
  cns.resize(n_elem);
  nss.resize(n_elem);
  n_sec_states = 0;
  for (int k = 0; k < n_elem; k++)
  {
    cns[k] = elem_vec[k]->getNumberOfStates();
    nss[k] = elem_vec[k]->getNumberOfSecStates();
    n_sec_states += nss[k];
  }

  DenseComplexVector sourcev;
  mnam->getSource(0, sourcev);
  
  // This is not ideal, but it is the easiest way for now
  if (isSet(&sweep))
  {
    Element* vsource = cir->getElement(sweep);
    if (vsource->getName() != "vsource")
		{
			sprintf(msg, "%s: not a vsource element", sweep.c_str());
			report(FATAL, msg);
    }
    // First time
    unsigned neq;
    // Get the position in the source vector
    vsource->getExtraRC(nrow, neq);
    // Remember that MNAM indices begin at 0
    sourcev[nrow-1] = vstart;
  }
  
  // Allocate memory for Msv and Ssv
  Msv.reshape(n_states, n_states);
  Ssv.resize(n_states);

  // Temporary vectors
  DenseComplexVector tempvec1(mnam->getDim());
  DenseComplexVector tempvec2(mnam->getDim());

  // Solve for each row of T and
  for (int i=0; i < n_states; i++)
  {
    // Set the 1 and -1 in tempvec tempvec1
    if (T(0,i))
      tempvec1[T(0,i)-1].real() += one;
    if (T(1,i))
      tempvec1[T(1,i)-1].real() -= one;

    mnam->solve(0, tempvec1, tempvec2);

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
      Msv(i, j) = c1.real();
    }
  }
  
  // ----------------------------------------------------------------------
  // ---------- Allocate memory for the simulation.
  // ----------------------------------------------------------------------
  Xw.resize(n_states);
  Xw1.resize(n_states);
  Vw.resize(n_states);
  Iw.resize(n_states);

  // Create the interface class using DC constructor
  tdsv = new TimeDomainSV(&(Xw1[0]), &(Vw[0]), &(Iw[0]), max_n_states);

  // Set initial guess to zero
  for (int j = 0; j < Xw.length(); j++)
    Xw[j] = zero;

  // ----------------------------------------------------------------------
  // ----------------- Run the nonlinear solver
  // ----------------------------------------------------------------------
  SerialCommunicator Comm1;

  // Get the process ID and the total number of processors
  int MyPID = Comm1.MyPID();
  int NumProc = Comm1.NumProc();

  // Get the number of elements from the command line
  int NumGlobalElements = n_states;

  // The number of unknowns must be at least equal to the 
  // number of processors.
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
  interface = Teuchos::rcp(new Nox_Interface(this, NumGlobalElements, Comm1));
  
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
  // *********** End Setup Nonlinear Solver *********************

  // Call solving routine
  allocflag = true;
  double residual = solveDC(sourcev);

  X.resize(Xw.length());
  for (int j = 0; j < Xw.length(); j++)
    X[j] = Xw[j];
  U.resize(Vw.length());
  for (int j = 0; j < Vw.length(); j++)
    U[j] = Vw[j];
  I.resize(Iw.length());
  for (int j = 0; j < Iw.length(); j++)
    I[j] = Iw[j];

  return residual;
}

// overloaded version of getBias. This function also sets the
// nodal voltages at the interface.
// no DC sweep is possible here. This function is used to perform
// a DC analysis before starting tran2.
double DC::getBias(Circuit* cir, const ElemFlag& mnam_mask,
           const ElemFlag& nonlinear_mask, DenseDoubleVector& X,
           DenseDoubleVector& U, DenseDoubleVector& I, DenseDoubleVector & Ul)
{
  // char msg[80];
  assert(!allocflag);
  this->cir = cir;

  // Create frequency vector
  //double f_vec_temp[] = {0.0};
  //DenseDoubleVector f_vec(Teuchos::View, f_vec_temp, 1);
  DenseDoubleVector f_vec(1);

  // ----------------------------------------------------------------------
  // -------------- Build one Msv matrix
  // ----------------------------------------------------------------------
  // Create the MNAM
  mnam = new FreqMNAM(f_vec, cir, mnam_mask);
  mnam->setFreqIndex(0);

  // Get T matrix
  buildTIncidence(nonlinear_mask, cir, T, elem_vec, n_states, max_n_states);
  if (gcomp) // Add the compensation resistors to the MNAM
    for (int i = 0; i < n_states; i++)
      mnam->setAdmittance(T(0,i), T(1,i) , gcomp);

  // this step is required
  mnam->FillComplete();

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

  DenseComplexVector sourcev;
  mnam->getSource(0, sourcev);

  // Allocate memory for Msv and Ssv
  Msv.reshape(n_states, n_states);
  Ssv.resize(n_states);

  // Temporary vectors
  DenseComplexVector tempvec1(mnam->getDim());
  DenseComplexVector tempvec2(mnam->getDim());

  // Solve for each row of T and
  for (int i=0; i < n_states; i++)
  {
    // Set the 1 and -1 in tempvec tempvec1
    if (T(0,i))
      tempvec1[T(0,i) - 1] += one;
    if (T(1,i))
      tempvec1[T(1,i) - 1] -= one;

    mnam->solve(0, tempvec1, tempvec2);

    // Restore tempvec1 to zero
    if (T(0,i))
      tempvec1[T(0,i) - 1] = zero;
    if (T(1,i))
      tempvec1[T(1,i) - 1] = zero;

    // Multiply T times tempvec2
    for (int j=0; j < n_states; j++)
    {
      double_complex c1(zero);
      if (T(0,j))
	      c1 = tempvec2[T(0,j) - 1];
      if (T(1,j))
	      c1 -= tempvec2[T(1,j) - 1];
      Msv(i, j) = c1.real();
    }
  }

  // ----------------------------------------------------------------------
  // ---------- Allocate memory for the simulation.
  // ----------------------------------------------------------------------
  Xw.resize(n_states);
  Xw1.resize(n_states);
  Vw.resize(n_states);
  Iw.resize(n_states);

  // Create the interface class using DC constructor
  tdsv = new TimeDomainSV(&(Xw1[0]), &(Vw[0]), &(Iw[0]), max_n_states);

  // Set initial guess to zero
  Xw.putScalar(zero);

  // ----------------------------------------------------------------------
  // ----------------- Run the nonlinear solver
  // ----------------------------------------------------------------------
  SerialCommunicator Comm1;

  // Get the process ID and the total number of processors
  int MyPID = Comm1.MyPID();
  int NumProc = Comm1.NumProc();

  // Get the number of elements from the command line
  int NumGlobalElements = n_states;

  // The number of unknowns must be at least equal to the 
  // number of processors.
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
  const NOX::Epetra::Group* soln_group = 0;
  const Epetra_Vector* next_guess = 0;

  // Call solving routine
  allocflag = true;
  double residual = solveDC(sourcev);

  // now we have Ssv, Msv and Iw.
  // Ul = Ssv + Msv*Iw;

  DenseDoubleVector tempProd;
  tempProd.resize(n_states);
  // Multiply Msv with Iw (matrix vector prod)
  for (int i = 0; i < n_states; i++)
    for (int j = 0; j < n_states; j++)
      tempProd[i] += Msv(i,j)*Iw[j];

  // add the elements of Ssv and tempProd and put it into
  // Ul
  for (int j = 0; j < Ssv.length(); j++)
    Ul[j] = Ssv[j] + tempProd[j];

  X.resize(Xw.length());
  X = Xw;
  U.resize(Vw.length());
  U = Vw;
  I.resize(Iw.length());
  I = Iw;

  return residual;
}

double DC::solveDC(DenseComplexVector& sourcev)
{
  assert(allocflag);
  // Solve for default source vector.
  DenseComplexVector tempvec2(mnam->getDim());
  mnam->solve(0, sourcev, tempvec2);
  // Calculate Ssv
  for (int j=0; j < n_states; j++)
  {
    double_complex c1(zero);
    if (T(0,j))
      c1 = tempvec2[T(0,j) - 1];
    if (T(1,j))
      c1 -= tempvec2[T(1,j) - 1];
    Ssv[j] = c1.real();
  }

  // ----------------------------------------------------------------------
  // ----------------- Run the nonlinear solver
  // ----------------------------------------------------------------------
  double residual = 0.0;
  Teuchos::RCP<DistributedDoubleVector> InitialGuess = interface->getSolution();
  // Set the initial guess 
  InitialGuess->PutScalar(0.0);
  NOX::Epetra::Vector noxInitGuess(InitialGuess, NOX::DeepCopy);
  //Reset the solver.
  solver->reset(noxInitGuess); 
  // Solve the nonlinear equation (call the nonlinear solver).
  int solvStatus = interface->evaluate(solver);
  //Get the residual
  const NOX::Epetra::Group& soln_group = 
    dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  residual = soln_group.getNormF();
  
  if (residual > 1e-3)
  {
    // Implement rudimentary source stepping
    report(MESSAGE, "Performing source stepping.");
    DenseDoubleVector Ssv1;
    Ssv1.resize(Ssv.length());
    Ssv1 = Ssv;
    for (int j = 0; j < Xw.length(); j++)
      Xw[j] = zero;
    int nsteps1 = 10;
    for (int k=1; k <= nsteps1; k++)
    {
      for (int i=0; i < n_states; i++)
	       Ssv[i] = Ssv1[i] * double(k) / double(nsteps1);
      int solvStatus = interface->evaluate(solver);
      const NOX::Epetra::Group& soln_group = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
      residual = soln_group.getNormF();
      char msg[200];
      sprintf(msg,"*** Step: %d    Residual: %g\n", k, residual);
      report(MESSAGE, msg);
    }
  }
  // The nonlinear solver interface must guarantee the correct
  // evaluation of V and I
	return residual;
}

void DC::clean()
{
  assert(allocflag);
  delete tdsv;
  delete mnam;
  allocflag = false;
}

void DC::func_ev(double * X, double * errFunc)
{
  // Evaluate first the time-domain elements.
  DenseDoubleVector X_vec(Teuchos::View, X, n_states);
  updateUI(X_vec);
  // Tell tdsv that the first evaluation is already completed.
  tdsv->clearflag();
  
  // Calculate error function
  for (int i = 0; i < n_states; i++)
  {
    double alpha = zero;
    for (int j = 0; j < n_states; j++)
      alpha -= Msv(i, j) * (Iw[j] - gcomp * Vw[j]);
    errFunc[i] = Ssv[i] + alpha - Vw[i];
  }
}

void DC::jacobian(double* X, DoubleSparseColMatrix& Jacobian)
{
  DoubleDenseMatrix J(n_states, n_states);
  
  DenseDoubleVector X_vec(Teuchos::View, X, n_states);
  assert(Xw1.length() == (unsigned int)n_states);
  int len = X_vec.length();
  for (int j = 0; j < len; j++)
    Xw1[j] = X_vec[j];

  // Number of elements
  int n_elem = elem_vec.size();
  DoubleDenseMatrix& Ju = tdsv->getJu();
  DoubleDenseMatrix& Ji = tdsv->getJi();

  // Go through all the nonlinear elements
  for (int i=0, j1=0, k=0; k < n_elem; k++)
  {
    // Set base index in interface object
    tdsv->setBase(i, cns[k], j1, nss[k]);
    // Call element evaluation
    elem_vec[k]->deriv_svTran(tdsv);

    for (int l=0; l < n_states; l++)
    {
      for (int j=0; j < cns[k]; j++)
      {
        J(l, i+j) = zero;
        for (int n=0; n < cns[k]; n++)
          J(l, i+j) -= Msv(l, i+n) * (Ji(n, j) - gcomp * Ju(n,j));
      }
    }

    for (int j=0; j < cns[k]; j++)
      for (int n=0; n < cns[k]; n++)
        J(i+j, i+n) -= Ju(j, n);

    i += cns[k];
    j1 += nss[k];
  }

  // copy contents of dense jacobian J into a sparse equivalent Jacobian
  // required by NOX
  std::vector<int> indices(n_states, 0);
  std::vector<double> values(n_states, 0.0);

  int k;
  // row loop
  for(int i = 0; i < n_states; i++)
  {
    k = 0;
    // column loop
    for(int j = 0; j < n_states; j++)
    {
      if (J(i,j) != 0)
      { 
        //We have a non-zero value save it
        values[k] = J(i,j);
        indices[k] = j; //Save the column the value was found in.
        k++;
      }
    }
    if (k != 0)
    {
      //Non zero values were found for row i
      Jacobian.ReplaceMyValues(i, k, &values[0], &indices[0]);
    }
  }  
}

void DC::updateUI(const DenseDoubleVector& X)
{
  // Set current state variable values.
  assert(Xw1.length() == X.length());
  int len = X.length();
  for (int j = 0; j < len; j++)
    Xw1[j] = X[j];

  // Number of elements
  int n_elem = elem_vec.size();
  // Go through all the nonlinear elements
  for (int i=0, j1=0, k=0; k < n_elem; k++)
  {
    // Set base index in interface object
    tdsv->setBase(i, cns[k], j1, nss[k]);

    // Call element evaluation
    // Always use time domain function evaluation.
    elem_vec[k]->svTran(tdsv);

    i += cns[k];
    j1 += nss[k];
  }
}

void DC::doOutput()
{
  // First check that the run() routine has been performed
  assert(nsteps);

  report(MESSAGE, "--- Writing output voltages and currents ...");

  // Allocate and copy global time vector
  allocTimeV_P(nsteps);
  for (int i=0; i<nsteps; i++)
    TimeV_P[i] = vstart + vstep * i;

  // Use time domain vectors for output.  Some day we may define a
  // custom output variable for DC analysis.
  int current_ns;
  int n_elem = elem_vec.size();   // Number of elements
  int i=0;
  for (int k=0; k < n_elem; k++)
  {
    current_ns = elem_vec[k]->getNumberOfStates();
    for(int j=0; j< current_ns; j++)
    {
      DenseDoubleVector tmp_x(nsteps);
      DenseDoubleVector tmp_i(nsteps);
      DenseDoubleVector tmp_u(nsteps);

      for (int l=0; l<nsteps; l++)
      {
	       tmp_x[l] = Xsweep[l][j+i];
	       tmp_i[l] = Isweep[l][j+i];
	       tmp_u[l] = Vsweep[l][j+i];
      }
      elem_vec[k]->getElemData()->setRealX(j, tmp_x);
      elem_vec[k]->getElemData()->setRealI(j, tmp_i);
      elem_vec[k]->getElemData()->setRealU(j, tmp_u);
    }
    i += current_ns;
  }

  // Set nodal voltages
  // For each terminal, assign voltage
  Terminal* term = NULL;
  cir->setFirstTerminal();
  while((term = cir->nextTerminal()))
  {
    DenseDoubleVector tmp_v(nsteps);
    // Get MNAM index
    if (term->getRC())
    {
      // Remember that MNAM indices begin at 0
      for (int l = 0; l < nsteps; l++)
	      tmp_v[l] = Usweep[l][term->getRC()-1];
    }
    else
    {
      // This is a reference terminal
      for (int j = 0; j < tmp_v.length(); j++)
        tmp_v[j] = zero;
    }

    // Set terminal vector
    term->getTermData()->setRealV(tmp_v);
  }

  // Temporary vectors (initialized to zero)
  DenseDoubleVector tmp_i(nsteps);

  ElemFlag mask = LINEAR;
  // Loop throw all selected elements in circuit and fill current
  // vector if needed.
  cir->setFirstElement(mask);
  Element* elem = cir->nextElement();
  unsigned first_eqn, no_eqn;
  while(elem)
  {
    elem->getExtraRC(first_eqn, no_eqn);
    if (first_eqn)
    {
      // Fill current vector(s) in element
      for (unsigned i = 0; i < no_eqn; i++)
      {
	      // Remember that MNAM indices begin at 0
	      unsigned row = first_eqn - 1 + i;
	      for (int l = 0; l < nsteps; l++)
	        tmp_i[l] = Usweep[l][row];
	      elem->getElemData()->setRealI(i, tmp_i);
      }
    }
    // get next element pointer
    elem = cir->nextElement();
  }
}

