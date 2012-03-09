// This are the interfaces for the different integration methods
// by Carlos E. Christoffersen

#ifndef integMethod_h
#define integMethod_h 1

#include "../containers.h"
#include "CircVector.h"

class NLIntegMethod
{
	public:
  // virtual destructor
	virtual ~NLIntegMethod() {}
  // Returns approximate time derivative of x[i]
  virtual double derivX(const int& index) = 0;
  // Returns approximate time second derivative of x[i]
  virtual double deriv2X(const int& index) = 0;
  // Returns d(dx[i]/dt)/dx[i]
  virtual double getdx_dtFactor() = 0;
  // Returns d(d2x[i]/dt2)/dx[i]
  virtual double getd2x_dt2Factor() = 0;
  // Returns x[i](t-tau)
  virtual double delayX(const int& index, const double& t) = 0;
  // Returns d(x[i](t-tau))/dx[i]
  virtual double getDelayXFactor(const double& t) = 0;
  // Returns number of samples of x and dx/dt that must be stored to
  // approximate derivatives (e.g. Trapezoidal: 1,1).
  virtual void getNSamples(int &nx, int &ndx) = 0;
  virtual double derivY(const int& index, const double& currval) = 0;
  // Returns the derivative of the variable passed in x:
  // x[0]: current sample
  // x[1]: previous sample, etc.
  // dx[0]: not used
  // dx[1]: previous time derivative approximation, etc.
  virtual double deriv(double* x, double* dx) = 0;
  // Change time step value.
  // WARNING: use with caution with multi-step methods.
  virtual void changeStep(const double& h) = 0;
  // Store anything required before going to the next time step.
  virtual void store() {}
  virtual void predictX(int n_states, double *Integ_predX, CircVector*& cTimestep) {}
  virtual void predictU(int ls_size, double *Integ_predU, CircVector*& cTimestep) {}
};

class LIntegMethod
{
	public:
	// virtual destructor
	virtual ~LIntegMethod() {}
  // Build Md matrix ( = G + a C)
  virtual void buildMd(DoubleSparseColMatrix& M, const double& h) = 0;
  // Build right-hand-side vector
  virtual void buildSf(DenseDoubleVector& s1, const double& ctime) = 0;
  virtual void buildSf(double* s1, const double& ctime) {}
  // WARNING: use with caution with multi-step methods.
  virtual void changeStep(const double& h) = 0;
  // Store anything required before going to the next time step.
  virtual void store() {}
};

#endif

