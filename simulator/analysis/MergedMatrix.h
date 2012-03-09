// This class is to handle merged sparse matrices efficiently
//
// Author:
//          Carlos E. Christoffersen
//

#ifndef MergedMatrix_h
#define MergedMatrix_h 1
#include "../freeda_mtl.h"

class MergedMatrix
{
	public:

  MergedMatrix(DoubleSparseColMatrix& M1, DoubleSparseColMatrix& M2);

  ~MergedMatrix();

  // Get total number of nonzeros of combination matrix
  inline const int& getNNZ() const
  {
    return nnz;
  }

  // Get nnzs due only to M1 or M2.
  inline void getSizes(int& nnz1, int& nnz2)
  {
    nnz1 = this->nnz1;
    nnz2 = this->nnz2;
  }

  // Selects one column and set index to zero.
  // Returns nnz for current column.
  inline const int setColumn(const int& col)
  {
    assert(col < cols);
    current_idx = 0;
    current_col = col;
    int numNonZeroElems = 0;
    for (int i = 0; i < iv[current_col]->Length(); i++)
    {
      if (iv[current_col][i] == 0)
        numNonZeroElems++;
    }
    return numNonZeroElems;
  }

  // Reset row index to zero (for multiple sweeps in one column).
  inline void resetRow()
  {
    current_idx = 0;
  }

  // Returns the type of element:
  //                             0 - M1 and M2
  //                             1 - M1
  //                             1 - M2
  inline int getNext(double& val1, double& val2, int& index)
  {
    assert(unsigned(current_idx) < v1[current_col]->nnz());
    val1 = (*v1[current_col])[current_idx];
    val2 = (*v2[current_col])[current_idx];
    index = (*iv[current_col])[current_idx];
    return (*typev[current_col])[current_idx++];
  }

	private:
  DenseDoubleVector **v1;
  DenseDoubleVector **v2;
  DenseIntVector **iv;
  DenseIntVector **typev;
  int nnz;
  int nnz1;
  int nnz2;
  int cols;
  int current_col;
  int current_idx;
};

#endif

