// AutiDiffElement is the base class for nonlinear elements which
// use Sacado. The derived class supplies a suitable eval() routine.

#ifndef ADInterface_h
#define ADInterface_h 1

#include "Element.h"
//#include "fadiff.h"
//using namespace fadbad;
#include "Sacado.hpp"

#include <cassert>

class TimeDomainSV;
class FreqDomainSV;

// synonym for automatic differentiation variable
//typedef F<double> AD;
typedef Sacado::Fad::DFad<double> AD;

class ADInterface : public Element
{
	public:

  ADInterface(ItemInfo* einfo, ParmInfo* param_desc, const int& numparms,
	const string& iname);

  virtual ~ADInterface();

  // ---------------------------------------
  // Analysis Methods
  // ---------------------------------------

  // State variable HB.
  virtual void svHB(FreqDomainSV* fdsv);
  virtual void deriv_svHB(FreqDomainSV* fdsv);

  // State variable transient.
  virtual void svTran(TimeDomainSV* tdsv);
  virtual void deriv_svTran(TimeDomainSV* tdsv);

  // Return the element's delay vector
  virtual DenseDoubleVector getDelayVec();

  // return the current time of simulation
  inline double getCurrentTime()
  {
    return currentTime;
  }

  // Generic state variable evaluation routine
  // The meaning of each x component is assigned in
  void initializeAD(const DenseIntVector& var);
  void initializeAD(const DenseIntVector& var, const DenseIntVector& dvar);
  void initializeAD(const DenseIntVector& var, const DenseIntVector& dvar,
                    const DenseIntVector& d2var, const DenseIntVector& t_var,
                    const DenseDoubleVector& delay);

  virtual void eval(AD * x, AD * effort, AD * flow) = 0;

	private:

  // Number of states
  unsigned my_nstates;
  DenseIntVector var, dvar, d2var, t_var;
  DenseDoubleVector delay;

  // Aux variables to save additions
  int d2const, delconst;

  // Auxiliary arrays for built-in analysis evaluation routines
  static double *x, *f;

  // Default size of allocated memory
  static int aux_size;

  // Flag to signal if the function is taped (ADOL-C tape)
  bool taped;
  int num_indep_ADvars, num_dep_ADvars;

  // This flag signal if memory is already allocated.
  static bool SVWSstatus;
  void allocSVWS();

  // Auxiliary functions for HB methods
  void svHB_timeX(FreqDomainSV*& fdsv);
  void svHB_getx(FreqDomainSV*& fdsv, const int& idx);

  // current time in transient simulation
  double currentTime;
};

#endif
