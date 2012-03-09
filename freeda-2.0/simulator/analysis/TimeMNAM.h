// Sparse implementation of TimeMNAM using the MTL
// by Carlos E. Christoffersen

#ifndef TimeMNAM_h
#define TimeMNAM_h 1

#include "../network/Circuit.h"
#include "../containers.h"
#include "slist.h"

class TimeMNAM 
{
	public:
  TimeMNAM(Circuit* cir, const ElemFlag& mask);
	
  ~TimeMNAM();
  Sllist Lhead, clist, nlist;
  //--------------------------------------------------------------
  // Methods to operate the matrices
  //--------------------------------------------------------------
	void print() const;
  
  // Get matrix dimension
  inline const unsigned& getDim()
	{
		return dimension;
	}
	
  // Get matrices
  inline void getMatrices(DoubleSparseColMatrix*& M, DoubleSparseColMatrix*& Mp)
	{
		(*M) = *(this->M);
		(*Mp) = *(this->Mp);
	}
	
  // Get the source vector at a time point
  void getSource(const double &time, DenseDoubleVector& ansV);
	
  //--------------------------------------------------------------
  // Methods to fill the matrices used by the elements
  //--------------------------------------------------------------
  
  // Methods for filling M
  inline void setMElement(const unsigned& row, const unsigned& col, 
	const double& val)
	{
		if (row && col) 
			addElem(row, col, M, val);
	}
  
  inline void setMAdmittance(const unsigned& term1, const unsigned& term2,
	const double& val)
	{
		if (term1)
			addElem(term1, term1, M, val);
		if (term2)
			addElem(term2, term2, M, val);
		if (term1 && term2) {
			addElem(term1, term2, M, -val);
			addElem(term2, term1, M, -val);
		}
	}
  
  inline void setMQuad(const unsigned& row1, const unsigned& row2,
	const unsigned& col1, const unsigned& col2,
	const double& val)
	{
		if (row1 && col1)
			addElem(row1, col1, M, val);
		if (row2 && col2)
			addElem(row2, col2, M, val);
		if (row1 && col2)
			addElem(row1, col2, M, val);
		if (row2 && col1)
			addElem(row2, col1, M, val);
	}
  
  inline void setMOnes(const unsigned& pos, const unsigned& neg,
	const unsigned& eqn)
	{
		if (pos) 
		{
      rowVal[0] = one;
      colInd[0] = eqn-1;
      M->InsertGlobalValues(pos-1, 1, rowVal, colInd);
      rowVal[0] = one;
      colInd[0] = pos-1;
      M->InsertGlobalValues(eqn-1, 1, rowVal, colInd);
		}
		if (neg) 
		{
      rowVal[0] = -one;
      colInd[0] = eqn-1;
      M->InsertGlobalValues(neg-1, 1, rowVal, colInd);
      rowVal[0] = -one;
      colInd[0] = neg-1;
      M->InsertGlobalValues(eqn-1, 1, rowVal, colInd);
		}
	}
  
  // Methods for filling Mp
  inline void setMpElement(const unsigned& row, const unsigned& col, 
	const double& val)
	{
		if (row && col) 
		{
			addElem(row, col, Mp, val);
			Lhead=insert_sllist(row,col,val,Lhead);
		}
	}
  
  inline void setMpAdmittance(const unsigned& term1, const unsigned& term2,
	const double& val)
	{
		if (term1)
		{
			addElem(term1, term1, Mp, val);
			Lhead=insert_sllist(term1,term1,val,Lhead);
		}
		if (term2)
		{
			addElem(term2, term2, Mp, val);
			Lhead=insert_sllist(term2,term2,val,Lhead);
		}
		if (term1 && term2) 
		{
			addElem(term1, term2, Mp, -val);
			Lhead=insert_sllist(term1,term2,-val,Lhead);
			addElem(term2, term1, Mp, -val);
			Lhead=insert_sllist(term2,term1,-val,Lhead);
      
		}
	}
  
  inline void setMpQuad(const unsigned& row1, const unsigned& row2,
	const unsigned& col1, const unsigned& col2,
	const double& val)
	{
		if (row1 && col1)
			addElem(row1, col1, Mp, val);
		if (row2 && col2)
			addElem(row2, col2, Mp, val);
		if (row1 && col2)
			addElem(row1, col2, Mp, val);
		if (row2 && col1)
			addElem(row2, col1, Mp, val);
	}
  
  inline void setMpOnes(const unsigned& pos, const unsigned& neg,
	const unsigned& eqn)
	{
		if (pos) 
		{
      rowVal[0] = one;
      colInd[0] = eqn-1;
      Mp->InsertGlobalValues(pos-1, 1, rowVal, colInd);
      rowVal[0] = one;
      colInd[0] = pos-1;
      Mp->InsertGlobalValues(eqn-1, 1, rowVal, colInd);
		}
		if (neg) 
		{
      rowVal[0] = -one;
      colInd[0] = eqn-1;
      Mp->InsertGlobalValues(neg-1, 1, rowVal, colInd);
      rowVal[0] = -one;
      colInd[0] = neg-1;
      Mp->InsertGlobalValues(eqn-1, 1, rowVal, colInd);
		}
	}
	
  //--------------------------------------------------------------
  // Methods to fill the source vector
  //--------------------------------------------------------------
  // Reminder: reference source element is not stored
	
  // Set one element of the source vector
  inline void setSource(const unsigned& row, const double& val)
	{
		assert(row <= dimension);
		if (row)
			(*sourceV)[row-1] = val;
	}
  
  // Add value to a pair of elements of the source vector
  inline void addToSource(const unsigned& pos, const unsigned& neg, 
	const double& val)
	{
		assert(pos <= dimension);
		assert(neg <= dimension);
		if (pos)
			(*sourceV)[pos-1] += val;
		if (neg)
			(*sourceV)[neg-1] -= val;
	}
	
  // Get the current time
  inline const double& getTime() 
	{
		return curr_time;
	}
  
  // Indicate that the matrices have been filled
  inline void FillComplete()
  {
    M->FillComplete();
    Mp->FillComplete();

  }
  
	private:
	
  // utility function used to fill the matrices
  inline void addElem(const unsigned& i, const unsigned& j, 
	DoubleSparseColMatrix*& M1, double val)
	{
    rowVal[0] = val;
    colInd[0] = j-1;
    M1->InsertGlobalValues(i-1, 1, rowVal, colInd);
	}
  
  // Dimension of the system
  unsigned dimension;
	
  // Current time
  double curr_time;
	
  // Mask to select elements to be included in matrix
  ElemFlag mymask;
	
  // The circuit from which to extract the matrices.
  Circuit* my_circuit;
	
  // Use element vectors to hold the selected elements and
  // save time each time an update is requested.
  ElementVector elem_vec;
  ElementVector source_elem_vec;
	
  // The source vector
  DenseDoubleVector *sourceV;
	
  // M and M' matrices
  DoubleSparseColMatrix *M, *Mp;
  // arrays to help fill the sparse matrices M and Mp
  double * rowVal;
  int * colInd;
};

#endif

