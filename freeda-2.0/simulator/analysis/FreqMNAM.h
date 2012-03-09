// Sparse implementation of FreqMNAM using Amesos

#ifndef FreqMNAM_h
#define FreqMNAM_h 1

#include "../network/Circuit.h"
#include <cassert>
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Amesos.h"
#include "Teuchos_RCP.hpp"

class FreqMNAM
{
	public:

  FreqMNAM(const DenseDoubleVector& fvec, Circuit* cir, const ElemFlag& mask);
  ~FreqMNAM();
  
  // all analysis routines using FreqMNAMs must call this
  // after the MNAM is created and before it is solved
  void FillComplete();

  //--------------------------------------------------------------  
  //Methods related to the frequency vector
  //--------------------------------------------------------------  
  // Get the dimension of frequency vector
  inline unsigned getFVecDim()
	{
		return unsigned(freq_vec.length());
	}

  // Get one frequency from freq. vector.
  inline double getFreq(const unsigned& findex)
	{
		return freq_vec[findex];
	}

  // Get one copy of the frequency vector.
  inline DenseDoubleVector getFreqVec()
	{
		return freq_vec;
	}

  // Method to set the current frequency
  void setFreqIndex(const unsigned& findex);

  // Get the current frequency
  inline const double& getFreq()
	{
		return freq_vec[current_findex];
	}

  //--------------------------------------------------------------
  // Methods to operate the matrices
  //--------------------------------------------------------------
  // Get matrix dimension
  inline unsigned getDim()
  {
    return dimension;
  }

  // Get the source vector at a frequency
  void getSource(const unsigned& findex, DenseComplexVector& sourceV);

  // Solve one matrix using the fixed source vector
  void solve(const unsigned& findex, DenseComplexVector& result);

  // Solve one matrix for a given source vector
  void solve(const unsigned& findex, const DenseComplexVector& sourceV,
    DenseComplexVector& result);

  //--------------------------------------------------------------
  // Methods to fill the matrices used by the elements
  //--------------------------------------------------------------
  void setElement(const unsigned& row, const unsigned& col,	
    const double_complex& val);

  void setAdmittance(const unsigned& term1, const unsigned& term2,
    const double_complex& val);

  void setQuad(const unsigned& row1, const unsigned& row2, const unsigned& col1,
    const unsigned& col2,	const double_complex& val);

  void setOnes(const unsigned& pos, const unsigned& neg, const unsigned& eqn);

  //--------------------------------------------------------------
  // Methods to fill the source vectors
  //--------------------------------------------------------------
  // Set one element of the current source vector
  inline void setSource(const unsigned& row, const double_complex& val)
  {
    assert(row <= dimension);
    if (row)
      (*sourceV_vec[current_findex])[row-1] = val;
  }

  // Add value to a pair of elements of the current source vector
  void addToSource(const unsigned& pos, const unsigned& neg,
    const double_complex& val);

	private:
  
  // create a vector of smart pointers to Sparse Matrix objects
  // this allows the creation of an array of matrices. For freq domain,
  // there is a MNAM matrix required for each frequency
  // For now, since epetra cannot hold complex data, and the
  // amesos solver can solve real matrices only, the input data
  // is split into a real matrix and a complex matrix and an
  // augmented real system of twice the size is constructed and solved.
  // In the future, with the advent of tpetra, only a single complex
  // sparse matrix of regular dimension will be necessary
  std::vector< Teuchos::RCP< Epetra_CrsMatrix > > Areal;
  std::vector< Teuchos::RCP< Epetra_CrsMatrix > > Aimag;
  std::vector< Teuchos::RCP< Epetra_CrsMatrix > > A;
  
  // build the auxillary Matrix A
  void buildA();

  // utility function used to fill the matrices
  inline void addElem(const unsigned& i, const unsigned& j,
    DoubleSparseColMatrix*& M1, double val)
	{
		// This is needed because we do not know if the element exists
		// in the matrix.
		double currval((*M1)[i-1][j-1]);
		(*M1)[i-1][j-1] = currval + val;
	}

  // The frequency vector
  DenseDoubleVector freq_vec;

  // The current frequency index
  unsigned current_findex;

  // Dimension of the system
  unsigned dimension;

  // Mask to select elements to be included in matrix
  ElemFlag mymask;

  // The circuit from which to extract the matrices.
  Circuit* my_circuit;

  // Use an element vector to hold the selected elements and
  // save time each time an update is requested.
  ElementVector elem_vec;

  // The source vectors
  DenseComplexVector** sourceV_vec;
  
  // Amesos info
  Epetra_LinearProblem Problem;
  Amesos_BaseSolver * Solver;
  Amesos Factory;
};

#endif

