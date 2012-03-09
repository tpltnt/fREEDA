// This is the standard interface of time domain state variable
// analysis with the elements.
// by Carlos E. Christoffersen

#ifndef TimeDomainSV_h
#define TimeDomainSV_h 1

#include "integMethod.h"
#include "../containers.h"

class TimeDomainSV
{
	public:
  // DC constructor
  TimeDomainSV(double *x, double *vnl, double *inl, 
	const int& max_n_states);
  
	// Transient constructor
  TimeDomainSV(double *x, double *vnl, double *inl, 
	const int& max_n_states,NLIntegMethod*& nl_im);
	
  TimeDomainSV(NLIntegMethod*& nl_im, const int& max_n_states);
	
  ~TimeDomainSV() {}
	
  // Methods used by the elements
	
  inline const bool& DC() const
  {
    return dc;
  }
	
  inline const bool& firstTime() const
  {
    return first_time;
  }
	
  inline const double& getdt() const
  {
    return tstep;
  }
	
  inline const double& getCurrentTime() const
  {
    return currentTime;
  }
	
  inline const int& getIndex() const
  {
    return nt;
  }
	
  inline const double& getX(const int& index) const
  {
    return x_b[index];
  }
	
  inline const double getdX_dt(const int& index) const
  {
    if (dc)
      return zero;
    else
      return nl_im->derivX(index + ibase);
  }
	
  inline const double getd2X_dt2(const int& index) const
  {
    if (dc)
      return zero;
    else
      return nl_im->deriv2X(index + ibase);
  }
	
  inline const double getDelayedX(const int& index, const double& t) const
  {
    if (dc)
      return x_b[index];
    else
      return nl_im->delayX(index + ibase, t);
  }
	
  inline double getdY_dt(const int& index, const double& currval)
  {
    if (dc)
      return zero;
    else
      return nl_im->derivY(index + jbase, currval);
  }
	
  inline double& u(const int& index)
  {
    return vnl_b[index];
  }
  
  inline double& i(const int& index)
  {
    return inl_b[index];
  }
  
  // Clean the Jacobian block associated with the current element.
  void cleanJac();
	
  // get factor to multiply derivatives with respect to (dx/dt)
  // (a in the time marching method, e.g. (1/h for B. Euler).
  inline const double& getdx_dtFactor()
  {
    return a;
  }
  
	// Same for second derivative
  inline const double& getd2x_dt2Factor()
  {
    return a2;
  }
  
	// Same for delays
  inline double getDelayFactor(double& t)
  {
    if (dc)
      return one;
    else
      return nl_im->getDelayXFactor(t);
  }
  
	// Get time-domain Jacobian
  // The complete Jacobian is not stored here because some analysis
  // may not need the Jacobian in that way.
  inline DoubleDenseMatrix& getJu()
  {
    return Ju;
  }
  
	inline DoubleDenseMatrix& getJi()
  {
    return Ji;
  }
	
  // Get integration method pointer
  inline NLIntegMethod* getIM()
  {
    return nl_im;
  }
	
  // Sets current time and sets first_time flag
  void setTime(double *x, double *vnl, double *inl,
	const int& nt, const double& time, const double& tstep);
	
  // Resets first_time flag
  inline void clearflag()
  {
    first_time = false;
  }
	
  inline void setIBase(const int& ibase, const int& cns)
  {
    this->ibase = ibase;
    this->cns = cns;
    x_b = x + ibase;
    vnl_b = vnl + ibase;
    inl_b = inl + ibase;
  }
	
  void setBase(const int& ibase, const int& cns, const int& jbase, const int& nss);
	
	private:
  // Factor from nl_im repeated here for the sake of efficiency
  double a;
  double a2;
	
  // Flag to signal that DC behavior is expacted
  bool dc;
	
  // Flag to tell elements whether the evaluation is for the first time
  // or not
  bool first_time;
	
  // Time index
  int nt;
	
  // Time step
  double tstep;
	
  // Current time
  double currentTime;
	
  // Current element index base
  int ibase;
  // Current element secondary state variable index base
  int jbase;
  // Current number of state variables
  int cns;
  // Current number of secondary state variables
  int nss;
	
  // Pointers to the current vectors
  double *x;
  double *vnl;
  double *inl;
  // Pointers to current base state variable
  double *x_b, *vnl_b, *inl_b;
	
  // Jacobian matrices
  DoubleDenseMatrix Ju, Ji;
	
  // Integration method class to calculate derivatives and delays.
  NLIntegMethod* nl_im;
};

#endif

