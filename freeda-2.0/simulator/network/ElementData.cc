#include "ElementData.h"

ElementData::ElementData(const unsigned& numterms)
{
  assert(numterms);
  this->numterms = numterms;

  // Allocate number of blank output vectors
  it_vv = new DenseDoubleVector[numterms];
  ut_vv = new DenseDoubleVector[numterms];
  xt_vv = new DenseDoubleVector[numterms];

  if_vv = new DenseComplexVector[numterms];
  xf_vv = new DenseComplexVector[numterms];
}

ElementData::~ElementData()
{
  // This memory is allocated when numterms is set.
  delete [] it_vv;
  delete [] ut_vv;
  delete [] xt_vv;

  delete [] if_vv;
  delete [] xf_vv;
}

void ElementData::setRealI(const unsigned& index, const DenseDoubleVector& i_v)
{
  assert(index < numterms);
	this->it_vv[index].resize(i_v.length());
  this->it_vv[index] = i_v;
}

void ElementData::setRealI(const unsigned& index, const double i_v)
{
  const int sz = this->it_vv[index].length();
  this->it_vv[index].resize(sz + 1);
  this->it_vv[index][sz] = i_v;
}

const DenseDoubleVector& ElementData::getRealI(const unsigned& index) const
{
  assert(index < numterms);
  return it_vv[index];
}

void ElementData::setRealU(const unsigned& index, const DenseDoubleVector& u_v)
{
  assert(index < numterms);
	this->ut_vv[index].resize(u_v.length());
	this->ut_vv[index] = u_v;
}

void ElementData::setRealU(const unsigned& index, const double u_v)
{
  const int sz = this->ut_vv[index].length();
  this->ut_vv[index].resize(sz + 1);
  this->ut_vv[index][sz] = u_v;
}

const DenseDoubleVector& ElementData::getRealU(const unsigned& index) const
{
  assert(index < numterms);
  return ut_vv[index];
}

void ElementData::setRealX(const unsigned& index, const DenseDoubleVector& x_v)
{
  assert(index < numterms);
  this->xt_vv[index].resize(x_v.length());
	this->xt_vv[index] = x_v;
}

void ElementData::setRealX(const unsigned& index, const double x_v)
{
  const int sz = this->xt_vv[index].length();
  this->xt_vv[index].resize(sz + 1);
  this->xt_vv[index][sz] = x_v;
}

const DenseDoubleVector& ElementData::getRealX(const unsigned& index) const
{
  assert(index < numterms);
  return xt_vv[index];
}

void ElementData::setPhasorI(const unsigned& index, const DenseComplexVector& i_v)
{
  assert(index < numterms);
	this->if_vv[index].resize(i_v.length());
  this->if_vv[index] = i_v;
}

const DenseComplexVector& ElementData::getPhasorI(const unsigned& index) const
{
  assert(index < numterms);
  return if_vv[index];
}

void ElementData::setPhasorX(const unsigned& index, const DenseComplexVector& x_v)
{
  assert(index < numterms);
	this->xf_vv[index].resize(x_v.length());
  this->xf_vv[index] = x_v;
}

const DenseComplexVector& ElementData::getPhasorX(const unsigned& index) const
{
  assert(index < numterms);
  return xf_vv[index];
}

