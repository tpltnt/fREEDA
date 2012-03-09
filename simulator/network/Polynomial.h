#ifndef Polynomial_h
#define Polynomial_h 1

#include "ADInterface.h"

class Polynomial
{
  public:

  Polynomial(const int& dimension, const int& no_of_coeff, double* coefficient) ;

  ~Polynomial() {}

  inline int getDimension() ;

  // computes the value of the Polynomial.
  AD findvalue(AD *) ;

  private:

  AD subop(AD *, int) ;

  // stores the size of the controlling voltage vector.
  int contvoltvect_size ;

  // stores the size of the coefficient vector.
  int coeffvect_size ;

  // stores the coefficients.
  double* coefficient ;

  int f_internal,k_count;
  // stores the value of the Polynomial.
  AD value ;
};

#endif
