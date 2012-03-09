#include "ADInterface.h"
#include "../analysis/TimeDomainSV.h"
#include "../analysis/FreqDomainSV.h"

bool ADInterface::SVWSstatus = false;
int ADInterface::aux_size = 20;
double * ADInterface::x = NULL;
double * ADInterface::f = NULL;

ADInterface::ADInterface(ItemInfo* einfo, ParmInfo* param_desc,
                        const int& numparms, const string& iname) :
Element(einfo, param_desc, numparms, iname), my_nstates(0),
        num_indep_ADvars(0), num_dep_ADvars(0)
{
  taped = false;
  currentTime=0.0;
}

ADInterface::~ADInterface()
{
  if (SVWSstatus)
  {
    delete[] f;
    delete[] x;
  }
  SVWSstatus = false;
}

// State-variable related methods
void ADInterface::allocSVWS()
{
  if (SVWSstatus)
  {
    delete[] f;
    delete[] x;
  }
  x = new double[aux_size];
  f = new double[aux_size];
  // Set flag
  SVWSstatus=true;
}

// initialize automatic differentiation
// Each nonlinear element calls this function before control is
// transferred to the analysis routines
void ADInterface::initializeAD(const DenseIntVector& var)
{
  DenseIntVector novar;
  DenseDoubleVector nodelay;
  initializeAD(var, novar, novar, novar, nodelay);
}

void ADInterface::initializeAD(const DenseIntVector& var,
                               const DenseIntVector& dvar)
{
  DenseIntVector novar;
  DenseDoubleVector nodelay;
  initializeAD(var, dvar, novar, novar, nodelay);
}

void ADInterface::initializeAD(const DenseIntVector& var,
                               const DenseIntVector& dvar,
                               const DenseIntVector& d2var,
                               const DenseIntVector& t_var,
                               const DenseDoubleVector& delay)
{
  my_nstates=getNumberOfStates();
  assert(my_nstates);

  // Check sizes of variable arrays
  assert(var.length() <= my_nstates);
  assert(dvar.length() <= my_nstates);
  assert(d2var.length() <= my_nstates);
  assert(t_var.length() <= my_nstates);
  assert(t_var.length() == delay.length());

  this->var.resize(var.length());
  this->var = var;
  this->dvar.resize(dvar.length());
  this->dvar = dvar;
  this->d2var.resize(d2var.length());
  this->d2var = d2var;
  this->t_var.resize(t_var.length());
  this->t_var = t_var;
  this->delay.resize(delay.length());
  this->delay = delay;

  // set number of variables used by automatic differentiation package
  num_indep_ADvars = var.length() + dvar.length() + d2var.length() + t_var.length();
  num_dep_ADvars = 2 * my_nstates;
  d2const = var.length() + dvar.length();
  delconst = var.length() + dvar.length() + d2var.length();

  // Find the greater
  int max_dim = num_dep_ADvars;
  if (num_indep_ADvars > num_dep_ADvars)
    max_dim = num_indep_ADvars;
  // Reserve memory if needed
  if (!SVWSstatus || (max_dim > aux_size))
  {
    aux_size = max_dim;
    // realloc space
    allocSVWS();
  }

  // Initialize the inputs to zero
  for (int i=0; i < num_indep_ADvars; i++)
    x[i] = zero;

  // convert inputs to have form suitable for automatic differentiation
  AD ax[num_indep_ADvars];
  for (int i = 0; i < num_indep_ADvars; i++)
  {
    ax[i] = x[i];
    ax[i].diff(i, num_indep_ADvars);
  }

  // convert outputs to have form suitable for automatic differentiation
  AD effort[my_nstates];
  AD flow[my_nstates];

  // call the respective element's eval function
  eval(ax, effort, flow);

  // save evaluated values into conventional double form
  for (int i = 0; i < my_nstates; i++)
  {
    f[i] = effort[i].val();
    f[i + my_nstates] = flow[i].val();
  }
}

void ADInterface::svTran(TimeDomainSV* tdsv)
{
  // Prepare x
  for (unsigned i = 0; i < var.length(); i++)
    x[i] = tdsv->getX(var[i]);

  for (unsigned i = 0; i < dvar.length(); i++)
    x[var.length() + i] = tdsv->getdX_dt(dvar[i]);

  for (unsigned i = 0; i < d2var.length(); i++)
    x[d2const + i] = tdsv->getd2X_dt2(d2var[i]);

  for (unsigned i = 0; i < t_var.length(); i++)
    x[delconst + i] = tdsv->getDelayedX(t_var[i], delay[i]);

  // convert inputs to have form suitable for automatic differentiation
  AD ax[num_indep_ADvars];
  for (int i = 0; i < num_indep_ADvars; i++)
  {
    ax[i] = x[i];
    ax[i].diff(i, num_indep_ADvars);
  }

  // convert outputs to have form suitable for automatic differentiation
  AD effort[my_nstates];
  AD flow[my_nstates];

  // call the respective element's eval function
  eval(ax, effort, flow);

  // save evaluated values into conventional double form
  for (int i=0; i < my_nstates; i++)
  {
    f[i] = effort[i].val();
    f[i + my_nstates] = flow[i].val();
  }

  // set the current time
  currentTime = tdsv->getCurrentTime();

  for (unsigned i = 0; i < my_nstates; i++)
  {
    tdsv->u(i) = f[i];
    tdsv->i(i) = f[my_nstates + i];
  }

  return;
}

void ADInterface::deriv_svTran(TimeDomainSV* tdsv)
{
  // Prepare x
  for (unsigned i = 0; i < var.length(); i++)
    x[i] = tdsv->getX(var[i]);

  for (unsigned i = 0; i < dvar.length(); i++)
    x[var.length() + i] = tdsv->getdX_dt(dvar[i]);

  for (unsigned i = 0; i < d2var.length(); i++)
    x[d2const + i] = tdsv->getd2X_dt2(d2var[i]);

  for (unsigned i = 0; i < t_var.length(); i++)
    x[delconst + i] = tdsv->getDelayedX(t_var[i], delay[i]);

  // Evaluate derivatives
  tdsv->cleanJac();

  DoubleDenseMatrix& Ju = tdsv->getJu();
  DoubleDenseMatrix& Ji = tdsv->getJi();

  AD ax[num_indep_ADvars];
  AD effort[my_nstates];
  AD flow[my_nstates];

  for (int i = 0; i < num_indep_ADvars; i++)
  {
    ax[i] = x[i];
    ax[i].diff(i, num_indep_ADvars);
  }

  eval(ax, effort, flow);

  for (int j = 0; j < my_nstates; j++)
  {
    // df/dx
    for (int i = 0; i < var.length(); i++)
    {
      Ju(j, var[i]) = effort[j].fastAccessDx(i);
      Ji(j, var[i]) = flow[j].fastAccessDx(i);
    }

    // df/d(dx/dt)
    for (int i = 0; i < dvar.length(); i++)
    {
      Ju(j, dvar[i]) += tdsv->getdx_dtFactor() * effort[j].fastAccessDx(var.length() + i);
      Ji(j, dvar[i]) += tdsv->getdx_dtFactor() * flow[j].fastAccessDx(var.length() + i);
    }

    // df/d(d2x/dt2)
    for (int i = 0; i < d2var.length(); i++)
    {
      Ju(j, d2var[i]) += tdsv->getd2x_dt2Factor()
                         * effort[j].fastAccessDx(dvar.length() + var.length() + i);
      Ji(j, d2var[i]) += tdsv->getd2x_dt2Factor()
                         * flow[j].fastAccessDx(dvar.length() + var.length() + i);
    }

    // df/d(x-tau)
    for (int i = 0; i < t_var.length(); i++)
    {
      Ju(j, t_var[i]) += tdsv->getDelayFactor(delay[i])
                         * effort[j].fastAccessDx(d2var.length()
                         + dvar.length() + var.length() + i);
      Ji(j, t_var[i]) += tdsv->getDelayFactor(delay[i])
                         * flow[j].fastAccessDx(d2var.length()
                         + dvar.length() + var.length() + i);
    }
  }

  return;
}

void ADInterface::svHB_timeX(FreqDomainSV*& fdsv)
{
  // Tell fdsv to prepare the time-domain vectors
  for (unsigned i = 0; i < var.length(); i++)
  {
    fdsv->transformX(var[i]);
  }
  for (unsigned i = 0; i < dvar.length(); i++)
  {
    fdsv->evaldX_dt(dvar[i]);
    fdsv->transformdX_dt(dvar[i]);
  }
  for (unsigned i = 0; i < d2var.length(); i++)
  {
    fdsv->evald2X_dt2(d2var[i]);
    fdsv->transformd2X_dt2(d2var[i]);
  }
  for (unsigned i = 0; i < t_var.length(); i++)
  {
    fdsv->evalDelX(t_var[i], delay[i]);
    fdsv->transformDelX(t_var[i]);
  }
}

void ADInterface::svHB_getx(FreqDomainSV*& fdsv, const int& idx)
{
  // Prepare x
  for (unsigned i = 0; i < var.length(); i++)
  {
    x[i] = fdsv->getx(var[i])[idx];
  }
  for (unsigned i = 0; i < dvar.length(); i++)
  {
    x[var.length() + i] = fdsv->getdx_dt(dvar[i])[idx];
  }
  for (unsigned i = 0; i < d2var.length(); i++)
  {
    x[d2const + i] = fdsv->getd2x_dt2(d2var[i])[idx];
  }
  for (unsigned i = 0; i < t_var.length(); i++)
  {
    x[delconst + i] = fdsv->getDelx(t_var[i])[idx];
  }
}

void ADInterface::svHB(FreqDomainSV* fdsv)
{
  const int& noSamples = fdsv->noSamples();

  // Prepare time-domain inputs
  svHB_timeX(fdsv);

  // input vector in AD format
  AD ax[num_indep_ADvars];

  // output vectors in AD format
  AD effort[my_nstates];
  AD flow[my_nstates];

  // Loop through all time samples
  for (int idx = 0; idx < noSamples; idx++)
  {
    // prepare x
    svHB_getx(fdsv, idx);
    for (int i = 0; i < num_indep_ADvars; i++)
    {
      ax[i] = x[i];
      ax[i].diff(i, num_indep_ADvars);
    }

    // call the respective element's eval function
    eval(ax, effort, flow);

    // save evaluated values into conventional double form
    for (int i=0; i < my_nstates; i++)
    {
      f[i] = effort[i].val();
      f[i + my_nstates] = flow[i].val();
    }

    // store output
    for (unsigned i = 0; i < my_nstates; i++)
    {
      fdsv->vp(i)[idx] = f[i];
      fdsv->ip(i)[idx] = f[my_nstates + i];
    }
  }

  // Tell fdsv to convert the output back to the freq. domain
  for (unsigned i = 0; i < my_nstates; i++)
  {
    fdsv->transformvp(i);
    fdsv->transformip(i);
  }
}

void ADInterface::deriv_svHB(FreqDomainSV* fdsv)
{
  const int& noSamples = fdsv->noSamples();

  // Prepare time-domain inputs
  svHB_timeX(fdsv);

  AD ax[num_indep_ADvars];
  AD effort[my_nstates];
  AD flow[my_nstates];

  // Loop through all time samples
  for (int idx = 0; idx < noSamples; idx++)
  {
    // prepare x
    svHB_getx(fdsv, idx);

    for (int i = 0; i < num_indep_ADvars; i++)
    {
      ax[i] = x[i];
      ax[i].diff(i, num_indep_ADvars);
    }

    eval(ax, effort, flow);

    for (int j = 0; j < my_nstates; j++)
    {
      // df/dx
      for (int i = 0; i < var.length(); i++)
      {
        fdsv->Jvp(j, var[i])[idx] = effort[j].fastAccessDx(i);
        fdsv->Jip(j, var[i])[idx] = flow[j].fastAccessDx(i);
      }

      // df/d(dx/dt)
      for (int i = 0; i < dvar.length(); i++)
      {
        fdsv->Jvpdx_dt(j, dvar[i])[idx] = effort[j].fastAccessDx(var.length() + i);
        fdsv->Jipdx_dt(j, dvar[i])[idx] = flow[j].fastAccessDx(var.length() + i);
      }
      for (unsigned i = 0; i < d2var.length(); i++)
      {
        // df/d(d2x/dt2)
        fdsv->Jvpd2x_dt2(j, d2var[i])[idx] = effort[j].fastAccessDx(dvar.length()
                                             + var.length() + i);
        fdsv->Jipd2x_dt2(j, d2var[i])[idx] = flow[j].fastAccessDx(dvar.length()
                                             + var.length() + i);
      }
      for (unsigned i = 0; i < t_var.length(); i++)
      {
        fdsv->JvpDelx(j, t_var[i])[idx] = effort[j].fastAccessDx(d2var.length()
                                          + dvar.length() + var.length() + i);
        fdsv->JipDelx(j, t_var[i])[idx] = flow[j].fastAccessDx(d2var.length()
                                          + dvar.length() + var.length() + i);
      }
    }
  }

  // Clear Jacobians
  for (unsigned i = 0; i < my_nstates; i++)
  {
    for (unsigned j = 0; j < my_nstates; j++)
    {
      fdsv->getJacVp(i, j).putScalar(zero);
      fdsv->getJacIp(i, j).putScalar(zero);
    }
  }

  // Tell fdsv to buid the jacobian in the freq. domain.
  for (unsigned i = 0; i < var.length(); i++)
    fdsv->addJ(my_nstates, var[i]);
  for (unsigned i = 0; i < dvar.length(); i++)
    fdsv->addJdx_dt(my_nstates, dvar[i]);
  for (unsigned i = 0; i < d2var.length(); i++)
    fdsv->addJd2x_dt2(my_nstates, d2var[i]);
  for (unsigned i = 0; i < t_var.length(); i++)
    fdsv->addJDelx(my_nstates, t_var[i], delay[i]);
}

DenseDoubleVector ADInterface::getDelayVec()
{
  return delay;
}
