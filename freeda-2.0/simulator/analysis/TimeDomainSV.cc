#include "TimeDomainSV.h"

TimeDomainSV::TimeDomainSV(NLIntegMethod*& nl_im, const int& max_n_states)
: Ju(max_n_states, max_n_states), Ji(max_n_states, max_n_states)
{
  this->nl_im = nl_im;
  a = nl_im->getdx_dtFactor();
  a2 = nl_im->getd2x_dt2Factor();
  dc = false;
  this->x = this->vnl = this->inl = NULL;
  this->x_b = this->vnl_b = this->inl_b = NULL;
  first_time = true;
  ibase = 0;
}

TimeDomainSV::TimeDomainSV(double* x, double *vnl, double *inl,
const int& max_n_states,NLIntegMethod*& nl_im)
: Ju(max_n_states, max_n_states), Ji(max_n_states, max_n_states)
{
  //used in tran and tran3
  this->nl_im = nl_im;
  a = nl_im->getdx_dtFactor();
  a2 = nl_im->getd2x_dt2Factor();
  dc = false;
  x_b = this->x = x;
  vnl_b = this->vnl = vnl;
  inl_b = this->inl = inl;
  //  a = zero;
  //a2 = zero;
  first_time = true;
  ibase = 0;
}

TimeDomainSV::TimeDomainSV(double* x, double *vnl, double *inl,
const int& max_n_states)
: Ju(max_n_states, max_n_states), Ji(max_n_states, max_n_states)
{
  //used in tran2
  dc = true;
  x_b = this->x = x;
  vnl_b = this->vnl = vnl;
  inl_b = this->inl = inl;

  // Set default values for unused variables
  nt = 0;
  tstep = currentTime = zero;
  first_time = true;
  nl_im = NULL;
  ibase = 0;
  a = a2 = zero;
}

void TimeDomainSV::setBase(const int& ibase, const int& cns,
const int& jbase, const int& nss)
{
  this->ibase = ibase;
  this->cns = cns;
  this->jbase = jbase;
  this->nss = nss;
  x_b = x + ibase;
  vnl_b = vnl + ibase;
  inl_b = inl + ibase;
}

void TimeDomainSV::setTime(double *x, double *vnl, double *inl,
const int& nt, const double& currentTime,
const double& tstep)
{
  this->x = x;
  this->vnl = vnl;
  this->inl = inl;
  this->nt = nt;
  this->currentTime = currentTime;
  this->tstep = tstep;
  // Set flag
  first_time = true;
}

//Chris: 04/14/09 Changed cleanJac()

void TimeDomainSV::cleanJac()
{
  for (int i = 0; i < Ju.numRows(); i++)
  {
    for (int j = 0; j < Ju.numCols(); j++)
    {
      Ju(i,j) = zero;
      Ji(i,j) = zero;
    }
  }
}



/*
void TimeDomainSV::cleanJac()
{
  for (int i = 0; i <= cns-1; i++)
  {
    for (int j = 0; j <= cns-1; j++)
    {
      Ju(i,j) = zero;
      Ji(i,j) = zero;
    }
  }
}
*/
