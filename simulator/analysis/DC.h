// Class definition for the DC analysis
// 
// Author:
//          Carlos E. Christoffersen

#ifndef DC_h
#define DC_h 1

#include "Analysis.h"
#include "OFunction_Nox.h"
#include "NOX.H"

class TimeDomainSV;
class Nox_Interface;

// Main class definition follows

class DC : public Analysis, public OFunction_Nox
{
public:

  DC();
  ~DC();

  // The main analysis routine. 
  virtual void run(Circuit* cir);

  //Return the name of the element
  static const char* getNetlistName()
  {
     return ainfo.name;
  }

  // Evaluate error function vector F given X
  virtual void func_ev(double * X, double * errFunc);
  virtual void jacobian(double * X, DoubleSparseColMatrix& J);

  // Run the DC analysis and get solution vectors and residual
  double getBias(Circuit* cir,
		 const ElemFlag& mnam_mask, const ElemFlag& nonlinear_mask,
		 DenseDoubleVector& X, DenseDoubleVector& U, DenseDoubleVector& I);
  double getBias(Circuit* cir,
		 const ElemFlag& mnam_mask, const ElemFlag& nonlinear_mask,
		 DenseDoubleVector& X, DenseDoubleVector& U, DenseDoubleVector& I,
		 DenseDoubleVector& Ul);

  // Solve the DC problem for a given source vector
  double solveDC(DenseComplexVector& sourcev);

  // Clean the space allocated with getBias()
  void clean();

private:

  // Write results to network.
  void doOutput();

  // Update cX, cU and cI by evaluating the time-domain elements
  // setting the derivatives to zero.
  void updateUI(const DenseDoubleVector& X);

  // ------------------------------  Variables
  // Circuit pointer
  Circuit* cir;

  // Total number of state variables
  int n_states, n_sec_states;
  DenseIntVector cns, nss;    // Vectors for number of primary and secondary SV.
  // Maximum number of state variables aported by an element
  int max_n_states;
  // Vector to hold nonlinear elements
  ElementVector elem_vec;
  // Pointer to MNAM
  FreqMNAM* mnam;
  // Incidence matrix
  IntDenseMatrix T;
  // UL(X) = Ssv + Msv INL(X)
  DoubleDenseMatrix Msv;
  DenseDoubleVector Ssv;
  // Result vectors
  DenseDoubleVector* Usweep;
  DenseDoubleVector* Xsweep;
  DenseDoubleVector* Vsweep;
  DenseDoubleVector* Isweep;
  int nsteps;

  // Work vectors
  DenseDoubleVector Xw, Xw1;
  DenseDoubleVector Vw;
  DenseDoubleVector Iw;
  // Element interface object
  TimeDomainSV* tdsv;
 
  bool allocflag;
  unsigned nrow;
  // ------------------- Parameter-related variables
  string sweep;
  double vstart, vstop, vstep, gcomp;
  int deriv;

  // Analysis information
  static ItemInfo ainfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter information
  static ParmInfo pinfo[];
  
  Teuchos::RCP<Nox_Interface> interface;
  Teuchos::RCP<NOX::Solver::Generic> solver;
};

#endif
