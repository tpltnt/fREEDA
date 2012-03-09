// Class definition for the state-variable oriented Harmonic Balance
// analysis.
//
// Please note that this is only an adaptation taken from the old transim.
//
// Author:
//         Carlos E. Christoffersen
//

#ifndef SvHB_h
#define SvHB_h 1

// C headers for compatibility
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

extern "C"
{
	#include "../compat/mltypes.h"
	#include "../compat/dcx.h"
	#include "../compat/ml.h"
}

#include "../network/CircuitManager.h"
#include "Analysis.h"
#include "FreqMNAM.h"
#include "OFunction_Nox.h"
#include "Nox_Interface.h"
#include "NOX.H"

//class TimeDomainSV;
class Nox_Interface;
class FreqDomainSV;

extern FILE *output_F; // Defined in init_clean.cc

// This definitions are for compatibility
#define TRUE 1
#define FALSE 0

// Main class definition follows

class SVHB : public Analysis, public OFunction_Nox
{
	public:

  SVHB();
  ~SVHB();

  // The main analysis routine. It also performs output.
  virtual void run(Circuit* main_cir);

  //Return the name of the element
  static const char* getNetlistName()
  {
     return ainfo.name;
  }
  
  virtual void func_ev(double * X, double * errFunc);
  virtual void jacobian(double * X, DoubleSparseColMatrix& J);
  /* svHB_Func_ev:
	*
	*              Calculates  F(X) = U_L(X) - U_NL(X)
	*
	*   U_L(X) = svSource + svMatrix * I_NL(X)            (for each frequency)
	*
	*   U_NL(X) and I_NL(X) are the voltages and currents calculated
	* from the nonlinear element models.
	*/

  /* svHB_Jacobian
	*
	* Calculate the full Jacobian:
	*
	*      J = expanded(svMatrix) * Ji - Jv
	*/
  
  Teuchos::RCP<Nox_Interface> interface;
  Teuchos::RCP<NOX::Solver::Generic> solver;
	private:

  /*
	* Allocate and fill the T matrix
	*/
  int createT();

  // free T matrix
  void destroyT();

  /* svHB_Create_ws: Initialize the work structures for the evaluation
	* routine.
	*/
  void createWs();

  /* svHB_Destroy_ws: free memory.
	*/
  void destroyWs();

  /* This routine creates the main data structure containing the compressed
	* matrices and the fundamental parameters for the evaluation of the
	* error function.
	*/
  void create_svHBdata();

  // Create frequency vector
  void createFvec();

  /* Free memory of the svHBdata structure.
	*/
  void destroy_svHBdata();

  /* Calculate the voltages at all the circuit nodes. Also write the
	* terminal and node tables.
	*/
  void calc_output(doublev_t X, double residual);

  /* Solve the system equations
	*/
  double solve(doublev_t X);

  /*
	* Output results of the analysis to the graph structures.
	*/
  void doOutput(dcxm_t VI, dcxm_t I_NL,
	doublev_t X, double residual);

  /* svHB_Get_U_and_I: Call the nonlinear element routines and get the
	* nonlinear voltages and currents given the states.
	*/
  void get_U_and_I(doublev_t X, dcxm_t *U_NL, dcxm_t *I_NL, dcxm_t *U_L);

	//    /* Get one column of the expanded Jacobian (with complex DC rows
	//     * and columns).
	//     *
	//     * Arguments:
	//     *              X     : state variable vector
	//     *              reset : if FALSE continue from the previus call,
	//     *                      otherwise reset everything and exit.
	//     *              JvecU : Voltage Jacobian column.
	//     *              JvecI : Current Jacobian column.
	//     *              Jtmp  : Temporary Work Jacobian column.
	//     *              id_beg: begin index for vector mult.
	//     *              id_end: end index for vector mult.
	//     *
	//     * Return code:    2 - reset performed.
	//     *                 1 - normal calculation.
	//     *                 0 - end of calculation (last column reached).
	//     *                -1 - attempt to calculate past the last column.
	//     */
	//    int getJacCol(doublev_t X, int reset,
	//  		     dcxm_t *JvecU, dcxm_t *JvecI, dcxm_t *Jtmp,
	//  		     int *id_beg, int *id_end);

  /* This routine creates the omega vector. This vector is used by the
	* time-derivative routine. Also, it returns the number of time samples
	* required for the fft.
	*/
  doublev_t createOmega(int *n);

  //-------------------------------------------------------------------

  // Main circuit pointer
  Circuit* my_cir;

  // Vector to hold nonlinear elements
  ElementVector elem_vec;

  // Maximum number of state variables aported by an element
  int max_n_states;

  // Flag to signal if the tapes have been generated
  bool tape_flag;

  // Interface object with elements
  FreqDomainSV* fdsv;

  // This structure holds all the option information from the input
  // deck.  We keep this structure definition for compatibility with
  // existing code.
  typedef struct svHB_options
	{
    int numtones;                   /* Number of exciting tones */
    doublev_t tone;                 /* Vector of exiting tones */
    integerv_t h;                   /* Maximum index for each tone */
    int oversample;                 /* Use oversample in the FFT */
    int steps;                      /* source stepping */
    int diff;                       /* Approximate derivatives or use a.d. */
    int verbosity;                  /* Amount of output to print */
  } svHB_options_t;

  /* This structure holds all the information needed by the HB error
	* function evaluation routine.
	*/
  typedef struct svHB_data
	{
    int n_states;                   /* Number of states */
    int NoFreqPoints;               /* Number of frequencies */
    int NoSamples;                  /* Number of time samples considered */
    doublev_t Omega;                /* Angular frequency vector */
    doublev_t FreqV_P;              /* Frequency vector */
    FreqMNAM* mnam;                 /* Modified nodal admittance matrices */
    dcxm_t T;                       /* Transformation matrix */
    int sysdim;                     /* Nonlinear system dimension */
    dcxm_t svSource;                /* Vector of compressed source vectors */
    dcxmv_t svMatrix;               /* Vector of compressed MNAMs */
  } svHB_data_t;

  typedef struct svHB_work_space
	{
    dcxm_t JvecU,    /* Jacobian voltage column */
		JvecI,           /* Jacobian current column */
		Jtmp,            /* temporary space */
		U_NL,            /* Nonlinear voltage vectors */
		I_NL,            /* Nonlinear current vectors */
		U_L;             /* Linear voltage vectors */
  } svHB_work_space_t;

  svHB_data_t *svHBdata;

  svHB_options_t *svHBopt;

  svHB_work_space_t *svHBws;

  // ------------------- Parameter-related variables

  // Analysis information
  static ItemInfo ainfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  int n_freqs, n_freqs2;
  double fundamental, fundamental2;
  int oversample_par, steps_par, diff_par, verbosity_par, step_init;
  bool regrowth;
  int n_fund;
  double f_step, max_cfreq;

  // Parameter information
  static ParmInfo pinfo[];
  
};

#endif

