#include "Open.h"

// Static members
const unsigned Open::n_par = 0;

// Element information
ItemInfo Open::einfo =
{
  "open",
  "Open circuit to be used as a voltage meter",
  "Carlos Christoffersen",
  DEFAULT_ADDRESS"category:lumped",
  "2000_07_20"
};

// Parameter information
// ParmInfo Open::pinfo[] = {
//   {"rcomp", "Svtr compensation resistance (Ohms)", DOUBLE, false}
// };

Open::Open(const string& iname) : ADInterface(&einfo, NULL, n_par, iname)
{

  // Set the number of terminals
  setNumTerms(2);

  // Set flags
  setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

  // Set number of states
  setNumberOfStates(1);
  DenseIntVector var(1);
  DenseIntVector novar;
  DenseDoubleVector nodelay;
  initializeAD(var, novar, novar, novar, nodelay);
}

void Open::eval(AD * x, AD * effort, AD * flow)
{
  effort[0] = x[0];
  flow[0] = zero;
}

void Open::svHB(FreqDomainSV* fdsv)
{
  // The voltage is the state variable
  fdsv->getVp(0).resize(fdsv->getX(0).length());
  fdsv->getVp(0) = fdsv->getX(0);
  // The current is always zero
  for (int i = 0; i < fdsv->getIp(0).length(); i++)
    fdsv->getIp(0)[i] = zero;
}

void Open::deriv_svHB(FreqDomainSV* fdsv)
{
  fdsv->getJacVp(0,0).putScalar(zero);
  const int& nosamples = fdsv->noSamples();
  for (int i=0; i<nosamples; i++)
    fdsv->getJacVp(0,0)(i,i) = one;
  fdsv->getJacIp(0,0).putScalar(zero);
}

void Open::svTran(TimeDomainSV* tdsv)
{
  tdsv->u(0) = tdsv->getX(0);
  tdsv->i(0) = zero;
}

void Open::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = one;
  tdsv->getJi()(0,0) = zero;
}
