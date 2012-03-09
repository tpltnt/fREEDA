#ifndef TerminalData_h
#define TerminalData_h 1

#include "../containers.h"

class TerminalData
{
	public:
  TerminalData() {}
  ~TerminalData() {}

  // Set real voltage output vector
  inline void setRealV(const DenseDoubleVector& v_v)
	{
		this->vt_v.resize(v_v.length());
    this->vt_v = v_v;
	}

  inline void setRealV(const double v_v)
  {
    const int sz = this->vt_v.length();
    this->vt_v.resize(sz + 1);
    this->vt_v[sz] = v_v;
  }

  // Get real voltage output vector
  inline const DenseDoubleVector& getRealV() const
	{
		return vt_v;
	}

  // Set Phasor voltage output vector
  inline void setPhasorV(const DenseComplexVector& v_v)
	{
		this->vf_v.resize(v_v.length());
		this->vf_v = v_v;
	}

  // Get Phasor voltage output vector
  inline const DenseComplexVector& getPhasorV() const
	{
		return vf_v;
	}

	private:
  // Voltage in the time domain
  DenseDoubleVector vt_v;

  // Voltage in the freq domain
  DenseComplexVector vf_v;
};

#endif

