// This are the classes for backwards Euler integration
// by Carlos E. Christoffersen

#ifndef Euler_h
#define Euler_h 1

#include "integMethod.h"
#include "CircVector.h"
#include "TimeMNAM.h"

class NLEuler : public NLIntegMethod
{
	public:
  NLEuler(CircVector*& cX, CircVector*& cY, CircVector*& cU, 
	CircVector*& cTime, double& h);
  NLEuler(CircVector*& cX, CircVector*& cY, CircVector*& cTime, double& h);
  NLEuler(CircVector*& cX, CircVector*& cY, double& h);
  virtual ~NLEuler() {}
	
  virtual double derivX(const int& index);
  virtual double deriv2X(const int& index);
  virtual double getdx_dtFactor();
  virtual double getd2x_dt2Factor();
  virtual double delayX(const int& index, const double& t);
  virtual double getDelayXFactor(const double& t);
  virtual void getNSamples(int &nx, int &ndx);
  virtual double deriv(double* x, double * dx);
  virtual double derivY(const int& index, const double& currval);
  virtual void changeStep(const double& h);
  virtual void predictX(int n_states, double *Integ_predX, 
	CircVector*& cTimestep);
  virtual void predictU(int n_states, double *Integ_predU, 
	CircVector*& cTimestep);
	
	private:
  CircVector *cX, *cU, *cY;
  CircVector *cTime;
  // 1 / time step
  double a, a2;
};

class LEuler : public LIntegMethod
{
	public:
  LEuler(CircVector*& cU, TimeMNAM* mnam);
  virtual ~LEuler();
	
  virtual void buildMd(DoubleSparseColMatrix& M, const double& h);
  virtual void buildSf(DenseDoubleVector& s1, const double& ctime);
  virtual void buildSf(double* s1, const double& ctime);
  virtual void changeStep(const double& h);
  void multiply(DenseDoubleVector& u_n, double* s1) ;
	
	private:
  TimeMNAM *mnam;
  CircVector *cU;
  DoubleSparseColMatrix * M1;
  DoubleSparseColMatrix * M1p;
  double * rowValExtract;
  int * colIndExtract;
  DenseDoubleVector s2; 
  int size;
  // Time step
  double a,temp; 
  int col_count,row_index;
};

#endif

