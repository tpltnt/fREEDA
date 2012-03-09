// This is the Objective function interface definition for the Nox solver

#ifndef OFunction_Nox_h
#define OFunction_Nox_h 1

#include "../containers.h"

class OFunction_Nox
{
  public:
  OFunction_Nox() {}
  virtual ~OFunction_Nox() {}
  virtual void func_ev(double * X, double * errFunc) = 0;
  virtual void jacobian(double * X, DoubleSparseColMatrix& J){};
};

#endif

