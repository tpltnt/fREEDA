#ifndef SVTran2_h
#define SVTran2_h 1

#include <fstream>
#include <iomanip>
#include <string>
#include <iostream>

#include "Analysis.h"
#include "OFunction_Nox.h"
#include "CircVector.h"
#include "TimeMNAM.h"
#include "TimeDomainSV.h"
#include "Euler.h"
#include "Trapezoidal.h"
#include "DC.h"

// The system of equations for the analysis routine is
// G*u + C*du/dt = sf + T^T*iNL
// vL = T*u = ssv + Msv*iNL
// where
// Ssv = T[G + a*C]^-1 [sf - C*b] and
// Msv = T[G + a*C]^-1 T^T
// --- finally, the error function
// f = vL - vNL = 0
class Nox_Interface;

// Main class definition follows
class SVTran2 : public Analysis, public OFunction_Nox
{
public:
  SVTran2();
  ~SVTran2() { }

  // the workhorse routine: here is where all the vectors and
  // matrices are set up and the time stepping takes place
  void run(Circuit * cir);

  //Return the name of the element
  static const char* getNetlistName()
  {
     return ainfo.name;
  }

  // evaluate the error function
  virtual void func_ev(double * X, double * errFunc);
  
  // Evaluate Jacobian
  virtual void jacobian(double * X, DoubleSparseColMatrix& J);

  // the output routine, which fills all the vectors
  // required for printing the simulation output variables
  void doOutput();

private:

  // build the Msv matrix
  void buildMsv();
  
  // update vector U
  void updateU(const DenseDoubleVector & sourceVec);
  
  // call the nonlinear eval routines to get nonlinear I and V values
  void updateVInl(double* x_p);

  // a pointer to the circuit being simulated
  Circuit * cir;

  // ls_size is the dimension of the MNAM (the linear system)
  int ls_size;

  // list_size is the number of time steps, and sets the size
  // of the cuircular vectors
  int list_size;

  // number of state variables
  int n_states, n_sec_states, max_n_states;

  // The incidence matrix
  IntDenseMatrix T;

  // The element vector, containing the nonlinear elements
  ElementVector elem_vec;

  // Vectors for number of primary and secondary SV.
  DenseIntVector cns, nss;
  
  // time domain interface pointer
  TimeDomainSV * tdsv;
  
  // the number of time steps
  int n_tsteps;

  // declare all circular vectors
  CircVector *cU; // linear elements contribution
  CircVector *cX; // state variables
  CircVector *cY; // secondary state variables
  CircVector *cVnl; // non-linear voltages
  CircVector *cInl; // non-linear currents

  DenseDoubleVector Ssv, sf, vL;
  DoubleDenseMatrix Msv;
  DoubleDenseMatrix Minv, Tc, Ssv_temp;

  // lhs_v and rhs_v are used as vectors for the LU factorization
  // so, for a sparse matrix M (the MNAM), the linear system
  // of equations is M * lhs_v = rhs_v
  DistributedDoubleVector * lhs_v;
  DistributedDoubleVector * rhs_v;

  // LU factorization variables
  Epetra_LinearProblem LinProblem;
  Amesos_BaseSolver * LinSolver;
  Amesos LinFactory;

  // Parameter-related variables
  // Analysis information
  static ItemInfo ainfo;

  // Number of parameters of this analysis
  static const unsigned n_par;

  // Analysis parameters
  double h, tf, nst, gcomp;
  int deriv, int_method, out_steps;
  bool opt, verbose;
  // memscheme = 1: allocate memory for all time steps up front
  // memscheme = 2: incrementally allocate memory as simulation progresses
  // The first scheme requires more memory up front, but computation at
  // each step executes more quickly than the second scheme.
  int memscheme;
  string lsm; //linear solver string

  // Parameter information
  static ParmInfo pinfo[];

  //Teuchos::RCP<Interface> interface;
};

#endif

