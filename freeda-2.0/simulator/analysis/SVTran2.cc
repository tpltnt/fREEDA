#include "SVTran2.h"
using std::cout;
using std::endl;

#include <stdio.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <fstream>
#include <iomanip>


// This is required for the allocTimeV_P function
extern "C"
{
#include "../inout/ftvec.h"
#include "../inout/report.h"
}

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


// number of analysis parameters, a static variable
const unsigned SVTran2::n_par = 10;


// Element information
ItemInfo SVTran2::ainfo =
{
  "SVTran2",
  "State-Variable-Based Time-Marching Transient Analysis",
  "Carlos E. Christoffersen, Nikhil M. Kriplani",
  DEFAULT_ADDRESS
};


// Parameter description and types
ParmInfo SVTran2::pinfo[] =
{
  {"tstop", "Stop time (s)", TR_DOUBLE, true},
  {"tstep", "Time step (s)", TR_DOUBLE, true},
  {"nst", "No save time (s)", TR_DOUBLE, false},
  {"deriv", "Approximate derivatives or use automatic diff.", TR_INT, false},
  {"im", "Integration method (1=trapezoidal) (2=Backward Euler)", TR_INT, false},
  {"out_steps", "Number of steps skipped for output simulation progress", TR_INT, false},
  {"gcomp", "Compensation network conductance (S)", TR_DOUBLE, false},
  {"opt", "Run a DC analysis up front", TR_BOOLEAN, false},
  {"verbose", "Print verbose nonlinear solver information", TR_BOOLEAN, false},
  {"memscheme", "Allocate all required memory up front or do it incrementally", TR_INT, false}
};


// default values for all parameters
// note some parameters are obtained from the netlist only
SVTran2::SVTran2() : Analysis(&ainfo, pinfo, n_par), ls_size(0)
{
  // Parameter stuff
  paramvalue[0] = &(tf);
  paramvalue[1] = &(h);
  paramvalue[2] = &(nst = zero);
  paramvalue[3] = &(deriv = 0);
  paramvalue[4] = &(int_method = 1);
  paramvalue[5] = &(out_steps = 200);
  paramvalue[6] = &(gcomp = 1e-4);
  paramvalue[7] = &(opt = false);
  paramvalue[8] = &(verbose = false);
  paramvalue[9] = &(memscheme = 1);
}


void SVTran2::run(Circuit * cir)
{
  this->cir = cir;
  char msg[80];

  // Build time domain MNAM, adding all linear circuit elements
  ElemFlag mnam_mask(LINEAR);
  TimeMNAM mnam(cir, mnam_mask);

  // Build incidence matrix, T, and nonlinear element vector
  n_states = 0;
  max_n_states = 0;
  ElemFlag mask(NONLINEAR);
  // this routine builds the incidence matrix T based on the
  // number of nonlinear elements, and assigns the appropriate
  // values to n_states (the number of state variables).
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
    if (elem->getName() == "yee") // look for yee element
      conv_elem_vec[conv_elem_vec.size() - 1]->getTimeInfo(h, tf);
    // get next element pointer
    elem = cir->nextElement();
  }

  // Loop through all the nonlinear elements
  int n_elem = elem_vec.size();

  // cns contains the number of state variables per nonlinear element
  cns = DenseIntVector(n_elem);

  // nss contains the number of secondary state variables per nonlinear element
  nss = DenseIntVector(n_elem);
  n_sec_states = 0;
  for (int k=0; k < n_elem; k++)
  {
    cns[k] = elem_vec[k]->getNumberOfStates();
    nss[k] = elem_vec[k]->getNumberOfSecStates();
    n_sec_states += nss[k];
  }

  if (memscheme == 2)
  {
    // this vector contains all the time delays specified in the netlist
    std::vector<double> tempTimeDelays;
    tempTimeDelays.push_back(0.0); // set at least one element

    // Mask to search through nonlinear elements for delays
    mask = NONLINEAR;
    cir->setFirstElement(mask);
    Element* nl_elem = cir->nextElement();

    // Declare a temporary delay vector storage container
    DenseDoubleVector tempDelayVec;

    // Parse the NL elements and push the delay vector contents into vector<double> tempTimeDelays
    while (nl_elem)
    {
      tempDelayVec = nl_elem->getDelayVec();
      for (int k = 0; k < tempDelayVec.length(); k++)
        tempTimeDelays.push_back(tempDelayVec[k]);
      nl_elem = cir->nextElement();
    }

    // this is the max of all the time delays
    double maxTimeDelay = *max_element(tempTimeDelays.begin(), tempTimeDelays.end());

    // Setup simulation variables (use circular vectors)
    ls_size = mnam.getDim();
    int out_size = int((tf - nst) / h + 1.49999);

    // list_size is the size of the circular vectors
    // each circular vector stores variable values for all time steps
    // steps so the size of the circ vectors is the the number of time steps
    // list_size is reduced to the bare minimum number of time steps
    // in order to conserve memory in case of long transient simulations.
    // list_size's value is used in the allocation of all the circular
    // vectors below.
    list_size = int(2 + maxTimeDelay/h);
  }
  else
  {
    ls_size = mnam.getDim();

    // list_size is the size of the circular vectors
    // each circular vector stores variable values for all time steps
    // so the size of the circ vectors is the the number of time steps
    list_size = int((tf - nst) / h + 2);
  }

  // initialize all circular vectors
  cU = new CircVector(list_size, ls_size);
  cX = new CircVector(list_size, n_states);
  cVnl = new CircVector(list_size, n_states);
  cInl = new CircVector(list_size, n_states);
  cY = new CircVector(list_size, n_sec_states);

  Ssv = DenseDoubleVector(n_states);
  sf = DenseDoubleVector(ls_size);
  vL = DenseDoubleVector(n_states);

  // Choose an integration method
  // for now, these can be either Trapezoidal or Euler integration
  // Trapezoidal is the default
  LIntegMethod * l_im;
  NLIntegMethod * nl_im;
  if (int_method == 1)
  {
    l_im = new LTrapezoidal(cU, &mnam);
    nl_im = new NLTrapezoidal(cX, cY, h);
  }
  else
  {
    l_im = new LEuler(cU, &mnam);
    nl_im = new NLEuler(cX, cY, h);
  }

  // Create a sparse matrix, M, which contains the contributions
  // of the static as well as dynamic portions of the MNAM
  // in other words, M = G + a*C
  SerialCommunicator Comm;
  ProcessorMap Map(ls_size, 0, Comm);
  DoubleSparseColMatrix M(Copy, Map, int(ls_size/2));
  l_im->buildMd(M, h);
  // M is complete
  M.FillComplete();

  // Now compute the L and U factors of M
  lhs_v = new DistributedDoubleVector(Map);
  rhs_v = new DistributedDoubleVector(Map);
  LinProblem.SetOperator(&M);
  LinProblem.SetLHS(lhs_v);
  LinProblem.SetRHS(rhs_v);
  const char * SolverType = "Klu";
  LinSolver = LinFactory.Create(SolverType, LinProblem);

  // factorize, to produce the L and U components
  LinSolver->NumericFactorization();

  // Create Msv matrix
  buildMsv();
  // now Tc and Minv are generated

  // Create the interface class
  tdsv = new TimeDomainSV(nl_im, max_n_states);

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
  double residual = zero;
  double residual_max = zero;
  report(MESSAGE, "--- Starting transient simulation ...\n");
  sepLine();
  report(MESSAGE, "|   Step   |    Time (s)    |   Residual    | Max Residual |");
  sepLine();
  double ctime = zero;

  // Number of convolution-based elements
  int n_c_elem = conv_elem_vec.size();

  // set the number of time steps
  n_tsteps = int(tf / h + 0.499999);

  // Ssv_temp = Minv * Tc
  DoubleDenseMatrix Ssv_temp(n_states, ls_size);
  if (n_states)
    Ssv_temp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, one, Tc, Minv, zero);

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

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
  if (verbose)
    printParams.set("Output Information",
			     NOX::Utils::OuterIteration +
			     NOX::Utils::OuterIterationStatusTest +
			     NOX::Utils::InnerIteration +
			     NOX::Utils::LinearSolverDetails +
			     NOX::Utils::Parameters +
			     NOX::Utils::Details +
			     NOX::Utils::Warning +
           NOX::Utils::Debug +
			     NOX::Utils::TestDetails +
			     NOX::Utils::Error);
  else
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
  Teuchos::RCP<NOX::Epetra::Group> grpPtr = Teuchos::rcp(new NOX::Epetra::Group(printParams,
    iReq,	noxInitGuess, linSys));
  NOX::Epetra::Group& grp = *grpPtr;

  // Update time in interface class
  tdsv->setTime(&(cX->getCurrent()[0]), &(cVnl->getCurrent()[0]),
     &(cInl->getCurrent()[0]), 0, ctime, h);

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
  NOX::Epetra::Vector* guess = &noxInitGuess;
  const NOX::Epetra::Group* soln_group = 0;
  const Epetra_Vector* next_guess = 0;
  // *********** End Setup Nonlinear Solver *********************

  // ------------------
  // the main time loop
  // ------------------
  for (int nt = 0; nt <= n_tsteps; nt++)
  {
    // update current time variable
    ctime = h * nt;

    // Update time in interface class
    tdsv->setTime(&(cX->getCurrent()[0]),
                  &(cVnl->getCurrent()[0]),
                  &(cInl->getCurrent()[0]),
                  nt, ctime, h);

    // Update source vector sf
    // sf = source_vector - Mp * b_(n-1)
    l_im->buildSf(sf, ctime);

    // build the Ssv
    // Ssv = Tc * Minv * sf
    if (n_states)
      Ssv.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -one, Ssv_temp, sf, zero);

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
    if (residual > residual_max)
      residual_max = residual;
    *guess = next_guess;
    //*******NOX Solver End******

    updateU(sf);

    // Update convolution-based elements
    // Go through all convolution elements and send result vector
    for (int k = 0; k < n_c_elem; k++)
      conv_elem_vec[k]->setLastResult(cU->getCurrent(), ctime);

    // Update previous derivatives in integration method (if required)
    nl_im->store();
    l_im->store();

    // Print table line
    if (!(nt % out_steps))
    {
      sprintf(msg, "|  %6d  |  %e  |  %e  |  %e  |", nt, ctime, residual, residual_max);
      report(MESSAGE, msg);
    }

    if (memscheme == 2)
    {
      // PROCESS OUTPUT at each time step instead of doing it at the
      // end of the time loop. This is done to prevent the analysis
      // routine from allocating a large chunk memory of memory before
      // the start of the loop.

      // Fill currents and port voltages of time domain devices.
      int current_ns;
      int n_elem = elem_vec.size();   // Number of nonlinear elements
      int i = 0;
      double tmp_u, tmp_x, tmp_i;

      for (int k = 0; k < n_elem; k++)
      {
        current_ns = elem_vec[k]->getNumberOfStates();
        for (int j = 0; j < current_ns; j++)
        {
          // For the current, decompensate while copying.
          tmp_x = cX->getCurrent()[j+i];
          tmp_i = cInl->getCurrent()[j+i];
          tmp_u = cVnl->getCurrent()[j+i];

          elem_vec[k]->getElemData()->setRealX(j, tmp_x);
          elem_vec[k]->getElemData()->setRealI(j, tmp_i);
          elem_vec[k]->getElemData()->setRealU(j, tmp_u);
        }
        i += current_ns;
      }

      Terminal * term = NULL;
      cir->setFirstTerminal();
      while ((term = cir->nextTerminal()))
      {
        // Get MNAM index
        if (term->getRC())
        {
          // Remember that MNAM indices begin at 1
          i = term->getRC() - 1;
          tmp_u = cU->getCurrent()[i];
        }
        else // This is a reference terminal
          tmp_u = 0.0;

        term->getTermData()->setRealV(tmp_u);
      }

      // now for linear element currents
      ElemFlag mask = LINEAR;
      // Loop through all selected elements in circuit and fill current
      // vector if needed.
      cir->setFirstElement(mask);
      Element* elem = cir->nextElement();
      unsigned first_eqn, no_eqn;
      while (elem)
      {
        elem->getExtraRC(first_eqn, no_eqn);
        if (first_eqn)
        {
          // Fill current vector(s) in element
          for (unsigned i = 0; i < no_eqn; i++)
          {
            // Remember that MNAM indices begin at 1
            unsigned row = first_eqn-1 + i;
            tmp_u = cU->getCurrent()[row];
            elem->getElemData()->setRealI(i, tmp_u);
          }
        }
        // get next element pointer
        elem = cir->nextElement();
      }
    }

    // Advance the circular lists for the next step.
    cU->advance();
    cX->advance();
    cY->advance();
    cVnl->advance();
    cInl->advance();
  }
  // ------------------
  // end main time loop
  // ------------------

  residual = sqrt(residual);

  sprintf(msg, "--- Residual: %g", residual);
  sepLine();
  newLine();
  report(MESSAGE, msg);

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup =
    dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector& finalSolution =
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  // Output the parameter list
  if (verbose)
  {
    if (printing.isPrintType(NOX::Utils::Parameters))
    {
      printing.out() << endl << "Final Parameters" << endl
        << "****************" << endl;
      solver->getList().print(printing.out());
      printing.out() << endl;
    }
  }

  // Tests
  int status = 0; // Converged
  // 1. Convergence
  if (solvStatus != NOX::StatusTest::Converged)
  {
    status = 1;
    if (printing.isPrintType(NOX::Utils::Error))
      printing.out() << "Nonlinear solver failed to converge!" << endl;
  }

  if (memscheme == 2)
  {
    int out_size = int((tf - nst) / h + 1.49999);
    // Allocate and copy global time vector
    allocTimeV_P(out_size);
    for (int tindex = 0; tindex < out_size; tindex++)
      TimeV_P[tindex] = nst + h*tindex;
  }

  if (memscheme != 2) // this is the default case
    doOutput();

  // Erase allocated space
  delete tdsv;
  delete rhs_v;
  delete lhs_v;
  delete nl_im;
  delete l_im;
  delete cY;
  delete cInl;
  delete cVnl;
  delete cX;
  delete cU;

  return;
}


void SVTran2::buildMsv()
{
  // Two matrices are needed here:
  // 1: The inverse of the mnam M, say Minv. M is [G + aC] in the theory.
  // and
  // 2: A conventional incidence matrix, Tc, which can be extracted
  // from the the T matrix above. Tc will contain 1, -1, and 0s only.
  // dimension of Tc is [n_states x mnam_dim]
  // Finally, the Msv = Tc * Minv * Tc_transpose
  Minv = DoubleDenseMatrix(ls_size, ls_size);

  // perform a Solve() successively setting the RHS to each column of the
  // identity matrix and copy each LHS vector into matrix Minv,
  // column by column. Minv ends up being dense.
  for (int i = 0; i < ls_size; i++)
  {
    (*rhs_v)[i] = one;
    // the solution is now in lhs_v
    LinSolver->Solve();
    // copy lhs_v to the i_th column of Minv
    for (int j = 0; j < ls_size; j++)
      Minv(j,i) = (*lhs_v)[j];
    // now clear the vector rhs_v for the next iteration
    (*rhs_v)[i] = zero;
  }
  // clear the LHS
  (*lhs_v).PutScalar(zero);

  // now generate Tc, using lhs_v as a temporary vector
  Tc = DoubleDenseMatrix(n_states, ls_size);
  for (int i = 0; i < n_states; i++)
  {
    // Build by columns
    if (T(0,i))
      (*lhs_v)[T(0,i) - 1] -= one;
    if (T(1,i))
      (*lhs_v)[T(1,i) - 1] += one;
    // copy contents of lhs_v into i_th row of Tc
    for (int j = 0; j < ls_size; j++)
      Tc(i,j) = (*lhs_v)[j];
    // Clear vector for next iteration
    (*lhs_v).PutScalar(zero);
  }

  // Finally, create Msv
  Msv = DoubleDenseMatrix(n_states, n_states);
  if (n_states)
  {
    DoubleDenseMatrix Tc_temp(n_states, ls_size);
    // Tc_temp is a temporary matrix product
    Tc_temp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, one, Tc, Minv, zero);
    Msv.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, -one, Tc_temp, Tc, zero);
  }
  return;
}


// update the U vector
// from theory:
// U = [G + aC]^-1 * [sf + T^T * i_NL]
// or
// U = Minv * [sf + T^T * i_NL]
void SVTran2::updateU(const DenseDoubleVector & sourceVec)
{
  // init a temporary vector
  DenseDoubleVector tempVec(n_states);

  // make shortcuts to the circular vectors
  // for nonlinear current and voltage
  DenseDoubleVector & inl = cInl->getCurrent();
  DenseDoubleVector & vnl = cVnl->getCurrent();

  // obtain the nonlinear element contributions
  for (int i = 0; i < n_states; i++)
    tempVec[i] = inl[i] - gcomp * vnl[i];

  DenseDoubleVector tempVec2(ls_size);

  if (n_states)
    tempVec2.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, one, Tc, tempVec, zero);

  // add the source vector contribution
  for (int i = 0; i < ls_size; i++)
    tempVec2[i] += sourceVec[i];

  DenseDoubleVector u_v(Teuchos::View, &(cU->getCurrent()[0]), ls_size);
  u_v.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, one, Minv, tempVec2, zero);

  return;
}


// calculate the error function
// this is called from the nonlinear interface
void SVTran2::func_ev(double * X, double * errFunc)
{
  // Evaluate first the time-domain elements.
  updateVInl(X);

  // Tell tdsv that the first evaluation is already completed.
  tdsv->clearflag();

  // get nonlinear current and voltages vector
  DenseDoubleVector & vnl = cVnl->getCurrent();
  DenseDoubleVector & inl = cInl->getCurrent();

  DenseDoubleVector tempVec(n_states);

  for (int i = 0; i < n_states; i++)
    tempVec[i] = inl[i] - vnl[i]*gcomp;

  // perform vL = Msv * (inl - vnl*gcomp)
  if (n_states)
    vL.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, one, Msv, tempVec, zero);

  for (int i = 0; i < n_states; i++)
    errFunc[i] = Ssv[i] + vL[i] - vnl[i];
}


// call the nonlinear element eval routines to get nonlinear I and V values
void SVTran2::updateVInl(double* x_p)
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
    i += cns[k];
    j += nss[k];
  }
}


// calculate the jacobian
void SVTran2::jacobian(double* X, DoubleSparseColMatrix& Jacobian)
{
  DoubleDenseMatrix J(n_states, n_states);

  // Set current state variable values.
  for (int j=0; j < n_states; j++)
    cX->getCurrent()[j] = X[j];

  // Number of elements
  int n_elem = elem_vec.size();
  DoubleDenseMatrix & Ju = tdsv->getJu();
  DoubleDenseMatrix & Ji = tdsv->getJi();

  // Go through all the nonlinear elements
  for (int i=0, j1=0, k=0; k < n_elem; k++)
  {
    // Set base index in interface object
    tdsv->setBase(i, cns[k], j1, nss[k]);

    // Call element evaluation
    elem_vec[k]->deriv_svTran(tdsv);
    for (int l=0; l < n_states; l++)
      for (int j=0; j < cns[k]; j++)
      {
        J(l, i+j) = zero;
        for (int n=0; n < cns[k]; n++)
          J(l, i+j) += Msv(l, i+n) * (Ji(n, j) - gcomp * Ju(n,j));
      }

    for (int j=0; j < cns[k]; j++)
      for (int n=0; n < cns[k]; n++)
      {
        J(i+j, i+n) -= Ju(j, n);
      }

    i += cns[k];
    j1 += nss[k];
  }

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
      if (J(i,j) != 0.0)
      {
        //We have a non-zero value save it
        values[k] = J(i,j);
        indices[k] = j; // Save the column the value was found in.
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


// write output vectors
void SVTran2::doOutput()
{
  // First check if the result matrices contain any data
  assert(cX);

  report(MESSAGE, "--- Writing output vectors ...");

  int out_size = list_size - 1;
  // Allocate and copy global time vector
  allocTimeV_P(out_size);
  for (int tindex = 0; tindex < out_size; tindex++)
    TimeV_P[tindex] = nst + h * tindex;

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
        unsigned row = first_eqn - 1 + i;
        for (int tindex = 0; tindex < out_size; tindex++)
          tmp_u[tindex] = cU->getPrevious(out_size - tindex)[row];
        elem->getElemData()->setRealI(i, tmp_u);
      }
    }
    // get next element pointer
    elem = cir->nextElement();
  }
  return;
}

