// ElementData contains the output information needed by one element
// by Carlos E. Christoffersen

#ifndef ElementData_h
#define ElementData_h 1

#include "../containers.h"
#include <cassert>

class ElementData
{
	public:
  
  ElementData(const unsigned& numterms);
  
  ~ElementData();
	
  // Set real current output vector
  void setRealI(const unsigned& index, const DenseDoubleVector& i_v);
  void setRealI(const unsigned& index, const double i_v);
	
  // Get real current output vector
  const DenseDoubleVector& getRealI(const unsigned& index) const;
	
  // Set real port voltage output vector
  void setRealU(const unsigned& index, const DenseDoubleVector& u_v);
  void setRealU(const unsigned& index, const double u_v);
	
  // Get real port voltage output vector
  const DenseDoubleVector& getRealU(const unsigned& index) const;
	
  // Set real port voltage output vector
  void setRealX(const unsigned& index, const DenseDoubleVector& x_v);
  void setRealX(const unsigned& index, const double x_v);
	
  // Get real port voltage output vector
  const DenseDoubleVector& getRealX(const unsigned& index) const;
	
  // Set phasor current output vector
  void setPhasorI(const unsigned& index, const DenseComplexVector& i_v);
	
  // Get Phasor current output vector
  const DenseComplexVector& getPhasorI(const unsigned& index) const;
	
  // Set Phasor state variable output vector
  void setPhasorX(const unsigned& index, const DenseComplexVector& x_v);
	
  // Get Phasor state variable output vector
  const DenseComplexVector& getPhasorX(const unsigned& index) const;
	
	private:
  
  // Number of terminals.
  unsigned numterms;
	
  // Currents in time domain
  DenseDoubleVector * it_vv;
  // Port voltages in time domain
  DenseDoubleVector * ut_vv;
  // State variable in time domain (or wavelet coefficients)
  DenseDoubleVector * xt_vv;
	
  // Currents in frequency domain
  DenseComplexVector * if_vv;
  // State variables in frequency domain
  DenseComplexVector * xf_vv;
};

#endif

