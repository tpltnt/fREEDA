// Class definition for the new state variable transient analysis
//
// Author:
//          Carlos E. Christoffersen
//          Based on transient code by Mete Ozkar.
//

#ifndef SVTr_h
#define SVTr_h 1

#include "Analysis.h"
#include "OFunction_Nox.h"
#include "NOX.H"

class TimeDomainSV;
class CircVector;

// Main class definition follows

class SVTr : public Analysis, public OFunction_Nox
{
public:

  SVTr();

  ~SVTr();

  // The main analysis routine.
  virtual void run(Circuit* cir);

  //Return the name of the element
  static const char* getNetlistName()
  {
     return ainfo.name;
  }

  // Evaluate error function vector F given X
  virtual void func_ev(double * X, double * errFunc);

  // Evaluate Jacobian of F
  virtual void jacobian(double * X, DoubleSparseColMatrix& Jacobian);
  
private:

  // Write results to network.
  void doOutput();

  // Build one Msv matrix at each frequency.
  void buildMsv();

  // Filter Msv in the frequency domain.
  void filterMsv();

  // Convert Msv to the time domain.
  void iDFTMsv();

  // Update cX, cU and cI by evaluating the time-domain elements
  void updateUI(const DenseDoubleVector& X);

  // ------------------------------  Variables
  const ElemFlag mnam_mask;
  const ElemFlag nle_mask;

  // Circuit pointer
  Circuit* cir;

  // Real number of frequencies (the parameter + 1)
  int n_freqs;
  // Frequency vector
  DenseDoubleVector f_vec;
  // Frequency step
  double fstep;
  // Time vector
  DenseDoubleVector t_vec;
  // Time step is a parameter of the analysis.
  // Number of time steps
  int n_tsteps;
  // Current time step (index).
  int nt;

  // Total number of state variables
  int n_states, n_sec_states;
  DenseIntVector cns, nss;    // Vectors for number of primary and secondary SV.
  // Maximum number of state variables aported by an element
  int max_n_states;
  // Vector to hold time-domain elements
  ElementVector elem_vec;
  // Impulse matrix
  DenseDoubleVector* Msv;
  // Convolution accumulator
  DenseDoubleVector p_conv;
  // Incidence Matrix T in compressed format
  IntDenseMatrix T;
  // Circular vectors
  CircVector* cX;
  CircVector* cY;
  CircVector* cU;
  CircVector* cI;
  // These matrices hold the result of the simulation.
  // MI_NL is also used in the convolution.
  DoubleDenseMatrix MI_NL;
  DoubleDenseMatrix MU_NL;

  // Element interface object
  TimeDomainSV* tdsv;

  // ------------------- Parameter-related variables

  // Analysis information
  static ItemInfo ainfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Analysis parameters
  int verbosity;
  double tstop;
  double tstep;
  int n_freqs_net;
  double gcomp;
  double tolerance;
  double filter_freq;
  int n_samples;
  int ntest;
  double imp_tol;
  bool check_imp;
  int out_steps;
  bool opt;
  bool adjust;
  int deriv;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
