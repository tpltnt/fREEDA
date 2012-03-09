#ifndef SVTran_h
#define SVTran_h 1

#include <fstream>
#include <iomanip>
#include <string>
#include <iostream>

#include "Analysis.h"
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
// --- Finally, the error function
// f = vL - vNL = 0
// At each time step, the Newton Raphson iteration results in following set of linear equations
//|G+Ca   -T^TJi | |un^j+1| = | sf,n-Cabn-1 + T^T(inl(xn^j)-Ji xn^j )|
//|T      -Jv    | |xn^j+1|   | vnl(xn^j) - Jv xn^j                  |
//           Ax             =    b

// Main class definition follows
class SVTran : public Analysis
{
	public:
		SVTran();
		~SVTran() { }

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
		//virtual void jacobian(double * X, DoubleDenseMatrix & J);

		// the output routine, which fills all the vectors
		// required for printing the simulation output variables
		void doOutput();


	private:

		// build the Msv matrix
		//void buildMsv();

		// update vector U
		//void updateU(DenseDoubleVector & sourceVec);

		// call the nonlinear eval routines to get nonlinear I and V values
		void updateVInl(double* x_p);

		//call the function to predict new time step
		void PredictNewStep(DenseDoubleVector& X, double *timestep);

		//call this to generate the vnl and inl column vectors, as well as Jv and Ji
		void eval_NL(double * X);

		//call this to update u vector for doOutput
		void eval_L();

		//call this to build the M = G +aC matrix for a given timestep size
		void build_M(DoubleSparseColMatrix *M, ProcessorMap& map);

		//call this to build the A matrix for solving Ax=b
		void build_A_matrix();

		//call this to build the b vector
		void build_b_vector(double* X);

		//call this to solve linear system of equations with newton iteration
		int solve_L(ProcessorMap& map);

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
		int nt;



		LIntegMethod * l_im;
		NLIntegMethod * nl_im;
		LIntegMethod * l_im2;
		NLIntegMethod * nl_im2;

		// declare all circular vectors
		CircVector *cU; // linear elements contribution
		CircVector *cX; // state variables
		CircVector *cY; // secondary state variables
		CircVector *cVnl; // non-linear voltages
		CircVector *cInl; // non-linear currents
		CircVector *cTimestep;
		CircVector *cTime;

		DenseDoubleVector sf, vL;

		DoubleDenseMatrix * Tc_matrix;
		DoubleDenseMatrix * Tct_matrix;   
		DoubleSparseColMatrix * M_array[4];
		DoubleSparseColMatrix * M_array2[4];
		map<int, DoubleSparseColMatrix> M_map;
		map<int, DoubleSparseColMatrix> M_map2;

		map<int,DoubleSparseColMatrix>::iterator M_map_iter;

		DoubleSparseColMatrix * M_ptr;

		// lhs_v and rhs_v are used as vectors for the LU factorization
		// so, for a sparse matrix M (the MNAM), the linear system
		// of equations is M * lhs_v = rhs_v
		DoubleSparseColMatrix * A_matrix;
		DistributedDoubleVector * lhs_v;
		DistributedDoubleVector * rhs_v;

		// LU factorization variables
		// Epetra provides wrappers for selected BLAS and Lapack routines
		// Amesos is direct solver class supports solvers like SuperLU
		Epetra_LinearProblem LinProblem;
		Amesos_BaseSolver * LinSolver;
		Amesos LinFactory;

		// Parameter-related variables
		// Analysis information
		static ItemInfo ainfo;

		// Number of parameters of this analysis
		static const unsigned n_par;

		//Added in these variables!!!
		DoubleDenseMatrix* vnl_matrix;
		DoubleDenseMatrix* inl_matrix;

		DoubleDenseMatrix* Jv_matrix;
		DoubleDenseMatrix* Ji_matrix;

		//adding u,x matrices... it was easier to use 'multiply' with matrix type
		//in addtion to X_vector
		DoubleDenseMatrix* U_matrix;
		DoubleDenseMatrix* X_matrix;


		// typedef Teuchos::SerialDenseVector<int,double> DenseDoubleVector; int represents data type of length of vector
		// double represents data type of stored value 
		DenseDoubleVector X_vector;
		DenseDoubleVector X_vector_last_iter;
		DenseDoubleVector X_vector_last_tstep;
		DenseDoubleVector X_vector_first_int;

		DenseDoubleVector U_vector;
		DenseDoubleVector U_vector_last_iter;
		DenseDoubleVector U_vector_last_tstep;
		DenseDoubleVector U_vector_first_int;


		double* Integ_predX;
		double* Integ_predU;

		// Error/residual
		double eerror_X,eerror_U,residual;
		int throwaway;
		int Max_NewtonRaphson_itr;


		// Analysis parameters
		double h, tf, nst, gcomp,abstol, reltol, tsteptol, tsmin, tsmax;
		int deriv, int_method, out_steps, tstep_method;
		bool opt;
		int changestep, cktorder, itypes, t_step_exp;
		//int M_index;
		double prev_timestep, new_timestep, orig_timestep;

		//added an analysis parameter ==> max_iter
		int max_iter;
		// Parameter information
		static ParmInfo pinfo[];

};

#endif

