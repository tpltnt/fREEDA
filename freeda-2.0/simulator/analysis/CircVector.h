// This class contains a circular set of vectors. It is used to hold data
// from the current time step up to a fixed number of previous time steps.
// tn, tn-1, tn-2, ... , tn-size+1, tn-size
// by Carlos E. Christoffersen

#ifndef CircVector_h
#define CircVector_h 1

#include "../containers.h"
#include <cassert>

class CircVector
{
	public:
  CircVector(const unsigned& list_size, const unsigned& vec_size);
  ~CircVector();

	// Elements in the current vectors can be modified
  DenseDoubleVector& getCurrent();
  // Previous vectors
  DenseDoubleVector& getPrevious(const unsigned& index);

  // Advance one position in the list
  void advance();

	// The current position
  unsigned curr_pos;

	private:
  // We use a vector of vectors to hold the circular list
  DenseDoubleVector* circ_vec;

  // The size of the list
  unsigned list_size;
};

#endif

