#include "MergedMatrix.h"

MergedMatrix::MergedMatrix(DoubleSparseColMatrix& M1,
DoubleSparseColMatrix& M2)
: nnz(0), nnz1(0), nnz2(0), cols(0), current_col(0), current_idx(0)
{
  assert(M1.ncols() == M2.ncols());
  assert(M1.nrows() == M2.nrows());

  cols = M1.ncols();
  v1 = new DenseDoubleVector*[cols];
  v2 = new DenseDoubleVector*[cols];
  iv = new DenseIntVector*[cols];
  typev = new DenseIntVector*[cols];

  for (int k=0; k < cols; k++)
	{
    // First pass: Count the number of elements in this column
    DoubleSparseColMatrix::iterator i;
    DoubleSparseColMatrix::iterator i1;
    DoubleSparseColMatrix::OneD::iterator j, jend;
    DoubleSparseColMatrix::OneD::iterator j1, j1end;
    i = M1.begin() + k;
    i1 = M2.begin() + k;
    j = (*i).begin(); jend = (*i).end();
    j1 = (*i1).begin(); j1end = (*i1).end();
    int counter(0);
    while (j != jend && j1 != j1end)
		{
      if (j.row() < j1.row())
				j++;
      else if (j.row() > j1.row())
				j1++;
      else
			{
				j++;
				j1++;
      }
      counter++;
    }
    // add remaining elements
    if (j == jend)
		while (j1 != j1end)
		{
			counter++;
			j1++;
		}
    else
		while (j != jend)
		{
			counter++;
			j++;
		}

    // Now we know how many elements this column has
    v1[k] = new DenseDoubleVector(counter);
    v2[k] = new DenseDoubleVector(counter);
    iv[k] = new DenseIntVector(counter);
    typev[k] = new DenseIntVector(counter);

    nnz += counter;

    // Second pass: Fill k column
    counter = 0;
    i = M1.begin() + k;
    i1 = M2.begin() + k;
    j = (*i).begin(); jend = (*i).end();
    j1 = (*i1).begin(); j1end = (*i1).end();
    while (j != jend && j1 != j1end)
		{
      if (j.row() < j1.row())
			{
				(*v1[k])[counter] = (*j);
				(*v2[k])[counter] = 0.;
				(*iv[k])[counter] = j.row();
				(*typev[k])[counter] = 1;
				nnz1++;
				j++;
      }
      else if (j.row() > j1.row())
			{
				(*v1[k])[counter] = 0.;
				(*v2[k])[counter] = (*j1);
				(*iv[k])[counter] = j1.row();
				(*typev[k])[counter] = 2;
				nnz2++;
				j1++;
      }
      else
			{
				(*v1[k])[counter] = (*j);
				(*v2[k])[counter] = (*j1);
				(*iv[k])[counter] = j.row();
				(*typev[k])[counter] = 0;
				j++;
				j1++;
      }
      counter++;
    }
    // add remaining elements
    if (j == jend)
		while (j1 != j1end)
		{
			(*v1[k])[counter] = 0.;
			(*v2[k])[counter] = (*j1);
			(*iv[k])[counter] = j1.row();
			(*typev[k])[counter] = 2;
			counter++;
			nnz2++;
			j1++;
		}
    else
		while (j != jend)
		{
			(*v1[k])[counter] = (*j);
			(*v2[k])[counter] = 0.;
			(*iv[k])[counter] = j.row();
			(*typev[k])[counter] = 1;
			counter++;
			nnz1++;
			j++;
		}
  }
}

MergedMatrix::~MergedMatrix()
{
  for (int k=0; k < cols; k++)
	{
    delete v1[k];
    delete v2[k];
    delete iv[k];
    delete typev[k];
  }
  delete [] v1;
  delete [] v2;
  delete [] iv;
  delete [] typev;
}

