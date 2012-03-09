// This are the classes for backwards Trapezoidal integration
// by Carlos E. Christoffersen

#ifndef Trapezoidal_h
#define Trapezoidal_h 1

#include "integMethod.h"
#include "CircVector.h"
#include "TimeMNAM.h"

class NLTrapezoidal : public NLIntegMethod
{
	public:
  NLTrapezoidal(CircVector*& cX, CircVector*& cY, CircVector*& cTime, 
	double& h);
  NLTrapezoidal(CircVector*& cX,CircVector*& cY, CircVector*& cU, 
	CircVector*& cTime, double& h);
  NLTrapezoidal(CircVector*& cX,CircVector*& cY, double& h);
  virtual ~NLTrapezoidal() {}
	
  virtual double derivX(const int& index);
  virtual double deriv2X(const int& index);
  virtual double getdx_dtFactor();
  virtual double getd2x_dt2Factor();
  virtual double delayX(const int& index, const double& t);
  virtual double getDelayXFactor(const double& t);
  virtual double derivY(const int& index, const double& currval);
  virtual void getNSamples(int &nx, int &ndx);
  virtual double deriv(double* x, double * dx);
  virtual void changeStep(const double& h);
  virtual void store();
  virtual void predictX(int n_states, double *Integ_predX,
	CircVector*& cTimestep);
  virtual void predictU(int ls_size, double *Integ_predU,
	CircVector*& cTimestep);
	
	private:
  CircVector *cX, *cY;
  CircVector *cTime;
  CircVector *cU;
  DenseDoubleVector xpn, ypn;
  DenseDoubleVector xsn;
  // time step
  double a, a2, da;
};

class LTrapezoidal : public LIntegMethod
{
	public:
	
  LTrapezoidal(CircVector*& cU, TimeMNAM* mnam);
	
  virtual ~LTrapezoidal();
	
  virtual void buildMd(DoubleSparseColMatrix& M, const double& h);
  virtual void buildSf(DenseDoubleVector& s1, const double& ctime);
  virtual void buildSf(double* s1, const double& ctime);
  virtual void changeStep(const double& h);
  virtual void store();
  void multiply(DenseDoubleVector& u_n, double* s1) ;
	private:
	
  TimeMNAM *mnam;
  DoubleSparseColMatrix * M1;
  DoubleSparseColMatrix * M1p;
  double * rowValExtract;
  int * colIndExtract;
  CircVector *cU;
  DenseDoubleVector upn;
  DenseDoubleVector s2;
  int size;
  // Time step
  double a,temp; 
  int col_count,row_index;
};

#endif

