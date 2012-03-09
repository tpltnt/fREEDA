// NOX include (for iostream, cmath, etc...)
#include "NOX_Common.H"
#include "Teuchos_ParameterList.hpp"

// Class Definition
#include "Nox_Interface.h"

// Epetra includes
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"


extern "C"
{
// include this to use the symbol table
#include <cstdio>
#include "../compat/stuff.h"
//#include "../compat/st.h"
// use the report() function for messages
#include "../inout/report.h"
}


// Constructor - creates the Epetra objects (maps and vectors) 
Nox_Interface::Nox_Interface(OFunction_Nox* of, int numGlobalElements, Epetra_Comm& comm, double xmin_,
                     double xmax_) :
  NumGlobalElements(numGlobalElements),
  NumMyElements(0),  // gets set after map creation
  MyPID(comm.MyPID()),
  NumProc(comm.NumProc()),
  xmin(xmin_),
  xmax(xmax_),
  factor(1.0),
  Comm(&comm),
  StandardMap(0),
  OverlapMap(0),
  Importer(0),
  rhs(0),
  Graph(0)
{
  this->of = of;
  
  n_states = numGlobalElements;
  
  // Construct a Source Map that puts approximately the same 
  // Number of equations on each processor in uniform global ordering
  StandardMap = new Epetra_Map(NumGlobalElements, 0, *Comm);

  // Get the number of elements owned by this processor
  NumMyElements = StandardMap->NumMyElements();
  
  // Construct an overlaped map for the finite element fill *************
  // For single processor jobs, the overlap and standard map are the same
  if (NumProc == 1) 
  {
    OverlapMap = new Epetra_Map(*StandardMap);
  } 
  else 
  {
    int OverlapNumMyElements;
    int OverlapMinMyGID;
    OverlapNumMyElements = NumMyElements + 2;
    if ((MyPID == 0) || (MyPID == NumProc - 1)) 
      OverlapNumMyElements --;
    
    if (MyPID==0) 
      OverlapMinMyGID = StandardMap->MinMyGID();
    else 
      OverlapMinMyGID = StandardMap->MinMyGID() - 1;
    
    int* OverlapMyGlobalElements = new int[OverlapNumMyElements];
    
    for (int i = 0; i < OverlapNumMyElements; i ++) 
      OverlapMyGlobalElements[i] = OverlapMinMyGID + i;
    
    OverlapMap = new Epetra_Map(-1, OverlapNumMyElements, 
			    OverlapMyGlobalElements, 0, *Comm);

    delete [] OverlapMyGlobalElements;

  } // End Overlap map construction *************************************

  // Construct Linear Objects  
  Importer = new Epetra_Import(*OverlapMap, *StandardMap);
  initialSolution = Teuchos::rcp(new DistributedDoubleVector(*StandardMap));
  
  // Construct a matrix
  jacobian = Teuchos::rcp(CreateJacobian()); 

  // Clean-up 
  jacobian->FillComplete();

  initializeSoln();
  
  // Variables to read symbol table
  generic_t gval;
  int type;
  // Set method
  if (St_GetSym(capOptionsT_P, "method", &type, &gval) == ST_SYM_FOUND)
  {
    if (type != GEN_TYPE_INT)
    {
      gval.i = 2;//Use default
      report(WARNING, "Option \"method\" should be an integer.");
    }
  }
  else
    gval.i = 2;//Use default

  switch (gval.i) 
  {
    case 1:
      // Set line search to full step
      report(MESSAGE, "Using Line Search Method - (Full Step).");
      method = "Full Step";
      break;
    case 2:
      // Set line search to back track
      report(MESSAGE, "Using Line Search Method - (Back Track).");
      method = "Backtrack";
      break;
    default:
      // Set line search (default)
      report(MESSAGE, "Using Line Search Method - (Back Track).");
      method = "Backtrack";
  }

  ftol = 1e-3; //default
  // set the stop tolerances
  if (St_GetSym(capOptionsT_P, "ftol", &type, &gval) == ST_SYM_FOUND)
    if (type == GEN_TYPE_DOUBLE)
      ftol = gval.d;
    else if (type == GEN_TYPE_FLOAT)
      ftol = gval.f;
    else
      report(WARNING, "Option \"ftol\" should be float.");

  maxit = 200; //default
  // set the maximun number of iterations
  if (St_GetSym(capOptionsT_P, "maxit", &type, &gval) == ST_SYM_FOUND)
    if (type == GEN_TYPE_INT)
      maxit = gval.i;
    else
      report(WARNING, "Option \"maxit\" should be an integer.");

  // Set direction
  direction = 0; //Default
  if (St_GetSym(capOptionsT_P, "direction", &type, &gval) == ST_SYM_FOUND)
    if (type != GEN_TYPE_INT)
      report(WARNING, "Option \"direction\" should be an integer.");
    else
      direction = gval.i;

  // set verbose level
  int output = 0;
  if (St_GetSym(capOptionsT_P, "output", &type, &gval) == ST_SYM_FOUND)
    if (type == GEN_TYPE_INT)
      output = (gval.i > -1) ? gval.i : 0;
    else
      report(WARNING, "Option \"output\" should be an integer.");
}

// Destructor
Nox_Interface::~Nox_Interface()
{
  delete Graph;
  delete Importer;
  delete OverlapMap;
  delete StandardMap;
}

bool Nox_Interface::computeF(const DistributedDoubleVector& x, 
		      DistributedDoubleVector& FVec, 
		      NOX::Epetra::Interface::Required::FillType fillType)
{
  x_svtran = x.Values();
  FVec_svtran = FVec.Values();
  this->of->func_ev(x_svtran, FVec_svtran);
  return true;
}

bool Nox_Interface::computeJacobian(const DistributedDoubleVector& x,
				Epetra_Operator& Jac)
{
  DoubleSparseColMatrix * J;
  J = dynamic_cast<DoubleSparseColMatrix*>(&Jac);
  if (J == NULL) 
  {
    cout << "*ERR* Problem_Interface::computeJacobian() - The supplied" << endl;
    cout << "*ERR* Epetra_Operator is NOT an Epetra_CrsMatrix!" << endl;
    throw;
  }
 
  x_svtran = x.Values();
  this->of->jacobian(x_svtran, *J);
  
  return true;
}

bool Nox_Interface::computePreconditioner(const DistributedDoubleVector & x,
				      Epetra_Operator& Prec,
				      Teuchos::ParameterList* precParams)
{
  cout << "ERROR: Interface::preconditionVector() - "
       << "Use Explicit Jacobian only for this test problem!" << endl;
  throw "Interface Error";
}

// Matrix and Residual Fills
NOX::StatusTest::StatusType Nox_Interface::evaluate(Teuchos::RCP<NOX::Solver::Generic>  &nox_solver)
{
  NOX::StatusTest::StatusType rc;
  rc = nox_solver->solve();
 
  return rc;
}

Teuchos::RCP<DistributedDoubleVector> Nox_Interface::getSolution()
{
  return initialSolution;
}
  
Teuchos::RCP<DoubleSparseColMatrix> Nox_Interface::getJacobian()
{
  return jacobian;
}

bool Nox_Interface::createGraph()
{}

// Set initialSolution to desired initial condition
bool Nox_Interface::initializeSoln()
{
  initialSolution->PutScalar(0); // Default initialization
  return true;
}

DoubleSparseColMatrix * Nox_Interface::CreateJacobian()
{
  int n_states = NumGlobalElements;
  int NumGlobalElements_1 = n_states; //n_states * n_states;
    
  // create a map
  Epetra_Map * Map = new Epetra_Map(NumGlobalElements_1,0,*Comm);
  // local number of rows
  int NumMyElements = Map->NumMyElements();
  // get update list
  int * MyGlobalElements = Map->MyGlobalElements();
  
  DoubleSparseColMatrix * A = new DoubleSparseColMatrix(Copy,*Map,n_states);
    
  // Add  rows one-at-a-time
  double *Values = new double[n_states];
  int *Indices = new int[n_states];
  
  for (int i=0; i<n_states; i++)
  {
      MyGlobalElements[i]=i;
      Indices[i] = i;
      Values[i]  = 0;
  }   
  int NumEntries = n_states;  
  for( int i=0 ; i<n_states; ++i ) 
  {
    A->InsertGlobalValues(i, NumEntries, Values, Indices);
  }
  // put matrix in local ordering
  A->FillComplete();

  delete [] Indices;
  delete [] Values;
  delete    Map;

  return A;
} 

