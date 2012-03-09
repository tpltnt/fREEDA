#include "CircVector.h"

CircVector::CircVector(const unsigned& list_size, const unsigned& vec_size)
{
  this->list_size = list_size;
  curr_pos = 0;
  circ_vec = new DenseDoubleVector[list_size];
  // Initialize vectors
  for (unsigned i = 0; i < list_size; i++)
  {
    circ_vec[i].resize(vec_size);
    for (int j = 0; j < vec_size; j++)
      circ_vec[i][j] = zero;
  }
}

CircVector::~CircVector()
{
  delete [] circ_vec;
}

DenseDoubleVector& CircVector::getCurrent()
{
	return circ_vec[curr_pos];
}

DenseDoubleVector& CircVector::getPrevious(const unsigned& index)
{
  assert(index < list_size);
  if (curr_pos < index)
    return circ_vec[list_size + curr_pos - index];
  else
    return circ_vec[curr_pos - index];
}

void CircVector::advance()
{
  curr_pos++;
  if (curr_pos == list_size)
    curr_pos = 0;
}
