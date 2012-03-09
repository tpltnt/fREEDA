// This is the standard interface of frequency domain state variable
// analysis with the elements.
// by Carlos E. Christoffersen

#ifndef FreqDomainSV_h
#define FreqDomainSV_h 1

#include "../containers.h"
#include <rfftw.h>

class FreqDomainSV
{
	public:

  FreqDomainSV(DenseDoubleVector& omega,	const int& usedfreqs, const int& max_nstates);
  ~FreqDomainSV();

  inline void setNFreqs(const int& nuf)
	{
		assert(nuf <= max_calc);
		nuf2 = 2 * nuf;
	}

  inline const int noSamples() const
	{
		return nosamples;
	}

  // Methods to operate the frequency domain vectors.
  // Convert Xi to time domain
  inline void transformX(const int& i)
	{
		freq2Time(Xf[i], xt[i]);
	}

  // Calculate derivative of Xi
  inline void evaldX_dt(const int& i)
	{
    int len = mv_Xf[i]->length();
    // pointers temp1 and temp2 point to the
    // start of vectors mv_Xf[i] and mv_dXf[i]
    // this seems to be the only way copy values
    // from one vector to another, based on how they
    // are originally defined
    double * temp1 = mv_Xf[i]->values();
    double * temp2 = mv_dXf[i]->values();
    // both vectors are already of length "len"
    // so no check for that is necessary
    for (int k = 0; k < len; k++)
    {
      *temp2 = *temp1;
      temp1++;
      temp2++;
    }
		deriv(dXf[i]);
	}

  // Convert dXi/dt to time domain
  inline void transformdX_dt(const int& i)
	{
		freq2Time(dXf[i], dxt[i]);
	}

  // Calculate second derivative of Xi
  inline void evald2X_dt2(const int& i)
	{
    int len = mv_Xf[i]->length();
    // pointers temp1 and temp2 point to the
    // start of vectors mv_Xf[i] and mv_dXf[i]
    // this seems to be the only way copy values
    // from one vector to another, based on how they
    // are originally defined
    double * temp1 = mv_Xf[i]->values();
    double * temp2 = mv_d2Xf[i]->values();
    // both vectors are already of length "len"
    // so no check for that is necessary
    for (int k = 0; k < len; k++)
    {
      *temp2 = *temp1;
      temp1++;
      temp2++;
    }
		deriv(d2Xf[i]);
		deriv(d2Xf[i]);
	}

  // Convert d2Xi/dt2 to time domain
  inline void transformd2X_dt2(const int& i)
	{
		freq2Time(d2Xf[i], d2xt[i]);
	}

  // Calculate delayed X
  inline void evalDelX(const int& i, const double& t)
	{
    int len = mv_Xf[i]->length();
    // pointers temp1 and temp2 point to the
    // start of vectors mv_Xf[i] and mv_dXf[i]
    // this seems to be the only way copy values
    // from one vector to another, based on how they
    // are originally defined
    double * temp1 = mv_Xf[i]->values();
    double * temp2 = mv_delXf[i]->values();
    // both vectors are already of length "len"
    // so no check for that is necessary
    for (int k = 0; k < len; k++)
    {
      *temp2 = *temp1;
      temp1++;
      temp2++;
    }
		delay(delXf[i], t);
	}

  // Convert delayed X to time
  inline void transformDelX(const int& i)
	{
		freq2Time(delXf[i], delxt[i]);
	}

  // Get x_i (time domain vector)
  inline DenseDoubleVector& getx(const int& i)
	{
		return *(mv_xt[i]);
	}

  // Get time domain derivative
  inline DenseDoubleVector& getdx_dt(const int& i)
	{
		return *(mv_dxt[i]);
	}

  // Get time domain second derivative
  inline DenseDoubleVector& getd2x_dt2(const int& i)
	{
		return *(mv_d2xt[i]);
	}

  // Get time domain delayed x
  inline DenseDoubleVector& getDelx(const int& i)
	{
		return *(mv_delxt[i]);
	}

  // Get time domain port voltage
  inline DenseDoubleVector& vp(const int& i)
	{
		return *(mv_vpt[i]);
	}

  // Get time domain port current
  inline DenseDoubleVector& ip(const int& i)
	{
		return *(mv_ipt[i]);
	}

  // Convert vpt_i to frequency domain
  void transformvp(const int& i)
	{
		time2Freq(vpt[i], Vpf[i]);
	}

  // Convert ipt_i to frequency domain
  void transformip(const int& i)
	{
		time2Freq(ipt[i], Ipf[i]);
	}

  // Get time-domain partial Jacobian diagonals
  inline DenseDoubleVector& Jvp(const int& i,const int& j)
	{
		return mv_jvp[max_nstates*i+j];
	}

  inline DenseDoubleVector& Jip(const int& i,const int& j)
	{
		return mv_jip[max_nstates*i+j];
	}

  inline DenseDoubleVector& Jvpdx_dt(const int& i,const int& j)
	{
		return mv_jdvp[max_nstates*i+j];
	}

  inline DenseDoubleVector& Jipdx_dt(const int& i,const int& j)
	{
		return mv_jdip[max_nstates*i+j];
	}

  inline DenseDoubleVector& Jvpd2x_dt2(const int& i,const int& j)
	{
		return mv_jd2vp[max_nstates*i+j];
	}

  inline DenseDoubleVector& Jipd2x_dt2(const int& i,const int& j)
	{
		return mv_jd2ip[max_nstates*i+j];
	}

  inline DenseDoubleVector& JvpDelx(const int& i,const int& j)
	{
		return mv_jdelvp[max_nstates*i+j];
	}

  inline DenseDoubleVector& JipDelx(const int& i,const int& j)
	{
		return mv_jdelip[max_nstates*i+j];
	}

  // Build frequency-domain Jacobians
  inline void buildJ(const int& n_states)
	{
		buildJac(n_states, JacVp, mv_jvp, Gamma1);
		buildJac(n_states, JacIp, mv_jip, Gamma1);
	}

  // The following methods require the Jacobians already initialized
  inline void addJ(const int& n_states, const int& js)
	{
		addJac(n_states, js, JacVp, mv_jvp, Gamma1);
		addJac(n_states, js, JacIp, mv_jip, Gamma1);
	}

  inline void addJdx_dt(const int& n_states, const int& js)
	{
		addJac(n_states, js, JacVp, mv_jdvp, dGamma1);
		addJac(n_states, js, JacIp, mv_jdip, dGamma1);
	}

  inline void addJd2x_dt2(const int& n_states, const int& js)
	{
		addJac(n_states, js, JacVp, mv_jd2vp, d2Gamma1);
		addJac(n_states, js, JacIp, mv_jd2ip, d2Gamma1);
	}

  inline void addJDelx(const int& n_states, const int& js, const double& t)
	{
		// delGamma1 depends on t
		build_delGamma1(t);
		addJac(n_states, js, JacVp, mv_jdelvp, delGamma1);
		addJac(n_states, js, JacIp, mv_jdelip, delGamma1);
	}

  // Methods used internally or for custom element models
  // Input vector is in nr complex format.
  void freq2Time(double* ivf, double* ovt);
  // Output vector is in nr complex format.
  void time2Freq(double* ivt, double* ovf);

  // Calculate the derivative of a frequency-domain vector
  // with respect to time.
  const void deriv(double* x);
  // Delay the signal given in the frequency domain.
  const void delay(double* x, const double& Time);

  // Methods to be used by the HB analysis or for custom element models
  // State variable vectors
  inline DenseDoubleVector& getX(const int& i)
	{
		return *(mv_Xf[i]);
	}

  // This is intended just for custom elements
  inline const DenseDoubleVector& getdX_dt(const int& i)
	{
		return *(mv_dXf[i]);
	}

  inline const DenseDoubleVector& getd2X_dt2(const int& i)
	{
		return *(mv_d2Xf[i]);
	}

  inline const DenseDoubleVector& getDelX(const int& i)
	{
		return *(mv_delXf[i]);
	}

  // Output vectors
  inline DenseDoubleVector& getVp(const int& i)
	{
		return *(mv_Vpf[i]);
	}

  inline DenseDoubleVector& getIp(const int& i)
	{
		return *(mv_Ipf[i]);
	}

  // Frequency domain Jacobian matrices
  inline DoubleDenseMatrix& getJacVp(const int& i, const int& j)
	{
		return *(mv_JacVp[i * max_nstates + j]);
	}

  inline DoubleDenseMatrix& getJacIp(const int& i, const int& j)
	{
		return *(mv_JacIp[i * max_nstates + j]);
	}

	private:

  // Output vector is in nr complex format.
  // Accumulate the result in ovf intead of overwrite.
  // Handy to calculate the Jacobian.
  void time2FreqAdd(double* ivt, double* ovf);

  void build_delGamma1(const double& t);

  void buildJac(const int& n_states,
	double**& Jac, DenseDoubleVector*& mv_Jac, double*& Gamma);

  void addJac(const int& n_states, const int& js,
	double**& Jac, DenseDoubleVector*& mv_Jac, double*& Gamma);

  // Vector of angular frequencies
  DenseDoubleVector omega;
  // 2 * Number of actually used frequencies
  int max_calc;
  int nuf2;
  // Number of time samples
  int nosamples;
  // Normalization factor for the FFT
  double norm_f;
  // Maximun possible number of state variables for an element
  int max_nstates;

  // FFTW vars
  rfftw_plan pf2t, pt2f;
  double *fd_p;

  // Workspace
  double **Xf, **xt, **dXf, **dxt;
  double **d2Xf, **d2xt;
  double **delXf, **delxt;
  double **Vpf, **vpt, **Ipf, **ipt;

  // DenseDoubleVector versions
  // Use pointer because we need references to the workspace vectors
  DenseDoubleVector **mv_Xf, **mv_xt, **mv_dXf, **mv_dxt;
  DenseDoubleVector **mv_d2Xf, **mv_d2xt;
  DenseDoubleVector **mv_delXf, **mv_delxt;
  DenseDoubleVector **mv_Vpf, **mv_vpt, **mv_Ipf, **mv_ipt;

  DenseDoubleVector *mv_jvp, *mv_jip, *mv_jdvp, *mv_jdip;
  DenseDoubleVector *mv_jd2vp, *mv_jd2ip, *mv_jdelvp, *mv_jdelip;

  double **JacVp, **JacIp;
  DoubleDenseMatrix **mv_JacVp, **mv_JacIp;

  double *Gamma1, *dGamma1, *d2Gamma1, *delGamma1, *TempMat;

};

#endif

