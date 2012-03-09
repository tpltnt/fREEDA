#include "FreqDomainSV.h"

void allocArray(double**& x, DenseDoubleVector**& mv_x,
const int& nosamples, const int& max_nstates)
{
  int size = nosamples * max_nstates;
  x = new double*[max_nstates];
  mv_x = new DenseDoubleVector*[max_nstates];
  x[0] = new double[size];
  mv_x[0] = new DenseDoubleVector(Teuchos::View, x[0], nosamples);
  for (int i=1; i < max_nstates; i++)
	{
    x[i] = x[0] + i * nosamples;
    mv_x[i] = new DenseDoubleVector(Teuchos::View, x[i], nosamples);
  }
}


void deleteArray(double**& x, DenseDoubleVector**& mv_x,
const int& max_nstates)
{
  for (int i=1; i < max_nstates; i++)
	{
    delete mv_x[i];
	}
  delete [] x[0];
  delete [] mv_x;
  delete [] x;
}


FreqDomainSV::FreqDomainSV(DenseDoubleVector& omega,
const int& usedfreqs, const int& max_nstates)
: omega(omega), max_calc(2 * usedfreqs),
nuf2(max_calc), nosamples(2 * omega.length()),
norm_f(double(2)/nosamples), max_nstates(max_nstates)
{
  // Create FFTW plan
  // FFTW_ESTIMATE is good for small circuits.
  // FFTW_MEASURE produces an optimal plan, but it introduces
  // a noticeable delay at setup time.
  // This should probably be an user's option.
  pf2t = rfftw_create_plan(nosamples,
	FFTW_COMPLEX_TO_REAL,
	FFTW_ESTIMATE);
  pt2f = rfftw_create_plan(nosamples,
	FFTW_REAL_TO_COMPLEX,
	FFTW_ESTIMATE);
  // Allocate space for a temporary vector to do the transformation
  fd_p = new double[nosamples];

  // Transform routines are ready now
  // Allocate workspace
  // State variable vectors
  allocArray(Xf, mv_Xf, nosamples, max_nstates);
  allocArray(xt, mv_xt, nosamples, max_nstates);
  allocArray(dXf, mv_dXf, nosamples, max_nstates);
  allocArray(dxt, mv_dxt, nosamples, max_nstates);
  allocArray(d2Xf, mv_d2Xf, nosamples, max_nstates);
  allocArray(d2xt, mv_d2xt, nosamples, max_nstates);
  allocArray(delXf, mv_delXf, nosamples, max_nstates);
  allocArray(delxt, mv_delxt, nosamples, max_nstates);
  // Output port voltages and currents
  allocArray(Vpf, mv_Vpf, nosamples, max_nstates);
  allocArray(vpt, mv_vpt, nosamples, max_nstates);
  allocArray(Ipf, mv_Ipf, nosamples, max_nstates);
  allocArray(ipt, mv_ipt, nosamples, max_nstates);

  // Create workspace for jacobian calculation
  mv_jvp = new DenseDoubleVector[max_nstates * max_nstates];
	for (int i = 0; i < max_nstates*max_nstates; i++)
	{
		mv_jvp[i] = DenseDoubleVector(nosamples);
	}

  mv_jip = new DenseDoubleVector[max_nstates * max_nstates];
	for (int i = 0; i < max_nstates*max_nstates; i++)
	{
		mv_jip[i] = DenseDoubleVector(nosamples);
	}

  mv_jdvp = new DenseDoubleVector[max_nstates * max_nstates];
	for (int i = 0; i < max_nstates*max_nstates; i++)
	{
		mv_jdvp[i] = DenseDoubleVector(nosamples);
	}

  mv_jdip = new DenseDoubleVector[max_nstates * max_nstates];
	for (int i = 0; i < max_nstates*max_nstates; i++)
	{
		mv_jdip[i] = DenseDoubleVector(nosamples);
	}

  mv_jd2vp = new DenseDoubleVector[max_nstates * max_nstates];
	for (int i = 0; i < max_nstates*max_nstates; i++)
	{
		mv_jd2vp[i] = DenseDoubleVector(nosamples);
	}

  mv_jd2ip = new DenseDoubleVector[max_nstates * max_nstates];
	for (int i = 0; i < max_nstates*max_nstates; i++)
	{
		mv_jd2ip[i] = DenseDoubleVector(nosamples);
	}

  mv_jdelvp = new DenseDoubleVector[max_nstates * max_nstates];
	for (int i = 0; i < max_nstates*max_nstates; i++)
	{
		mv_jdelvp[i] = DenseDoubleVector(nosamples);
	}

  mv_jdelip = new DenseDoubleVector[max_nstates * max_nstates];
	for (int i = 0; i < max_nstates*max_nstates; i++)
	{
		mv_jdelip[i] = DenseDoubleVector(nosamples);
	}

  JacVp = new double*[max_nstates * max_nstates];
  JacIp = new double*[max_nstates * max_nstates];
  mv_JacVp = new DoubleDenseMatrix*[max_nstates * max_nstates];
  mv_JacIp = new DoubleDenseMatrix*[max_nstates * max_nstates];
  for (int i=0; i < max_nstates * max_nstates; i++)
	{
    JacVp[i] = new double[nosamples*max_calc];
    JacIp[i] = new double[nosamples*max_calc];
    mv_JacVp[i] = new DoubleDenseMatrix(Teuchos::View, JacVp[i], nosamples, nosamples, max_calc);
    mv_JacIp[i] = new DoubleDenseMatrix(Teuchos::View, JacIp[i], nosamples, nosamples, max_calc);
  }

  // Build Gamma1
  Gamma1 = new double[max_calc * nosamples];
  dGamma1 = new double[max_calc * nosamples];
  d2Gamma1 = new double[max_calc * nosamples];
  delGamma1 = new double[max_calc * nosamples];
  TempMat = new double[max_calc * nosamples];
  // build each column of gamma
  // Gamma1 is column-oriented
  // Use TempMat as scratch space
  for (int i=0; i < nosamples; i++)
	{
    TempMat[i] = zero;
	}
  for (int i=0; i < max_calc; i++)
	{
    TempMat[i] = one;
    freq2Time(TempMat, Gamma1 + i * nosamples);
    TempMat[i] = zero;
  }
  // dGamma1 is the same but with the derivative operator applied first.
  // The fist 2 columns are zeros because the derivatives of DC and
  // the last component must be zero to satisfy symmetry conditions.
  // The fist 2 colums are zero
  // d2Gamma1 is very similar too.
  for (int i=0; i < nosamples; i++)
	{
    dGamma1[i] = (dGamma1 + nosamples)[i] =
    d2Gamma1[i] = (d2Gamma1 + nosamples)[i] = zero;
  }

	// TempMat is zeroed from previous calculation.
  for (int i=2; i < max_calc; i+=2)
	{
    // Real part
    TempMat[i+1] = omega[i>>1];
    freq2Time(TempMat, dGamma1 + i * nosamples);
    TempMat[i+1] = zero;
    // Imaginary part
    TempMat[i] = -omega[i>>1];
    freq2Time(TempMat, dGamma1 + (i+1) * nosamples);
    TempMat[i] = zero;
  }

	// TempMat is zeroed from previous calculation.
  for (int i=2; i < max_calc; i++)
	{
    TempMat[i] = -pow(omega[i>>1], 2);
    freq2Time(TempMat, d2Gamma1 + i * nosamples);
    TempMat[i] = zero;
  }
}

FreqDomainSV::~FreqDomainSV()
{
  delete [] TempMat;
  delete [] delGamma1;
  delete [] d2Gamma1;
  delete [] dGamma1;
  delete [] Gamma1;

  for (int i=0; i < max_nstates * max_nstates; i++)
	{
    delete mv_JacIp[i];
    delete mv_JacVp[i];
    delete [] JacIp[i];
    delete [] JacVp[i];
  }

  delete [] mv_JacIp;
  delete [] mv_JacVp;
  delete [] JacIp;
  delete [] JacVp;

  delete [] mv_jdelip;
  delete [] mv_jdelvp;
  delete [] mv_jd2ip;
  delete [] mv_jd2vp;
  delete [] mv_jdip;
  delete [] mv_jdvp;
  delete [] mv_jip;
  delete [] mv_jvp;

  // delete vector workspace
  // Output port voltages and currents
  deleteArray(ipt, mv_ipt, max_nstates);
  deleteArray(Ipf, mv_Ipf, max_nstates);
  deleteArray(vpt, mv_vpt, max_nstates);
  deleteArray(Vpf, mv_Vpf, max_nstates);
  // State variable vectors
  deleteArray(delxt, mv_delxt, max_nstates);
  deleteArray(delXf, mv_delXf, max_nstates);
  deleteArray(d2xt, mv_d2xt, max_nstates);
  deleteArray(d2Xf, mv_d2Xf, max_nstates);
  deleteArray(dxt, mv_dxt, max_nstates);
  deleteArray(dXf, mv_dXf, max_nstates);
  deleteArray(xt, mv_xt, max_nstates);
  deleteArray(Xf, mv_Xf, max_nstates);

  // delete fftw workspace
  delete [] fd_p;
  rfftw_destroy_plan(pf2t);
  rfftw_destroy_plan(pt2f);
}


void FreqDomainSV::build_delGamma1(const double& t)
{
  // Build delGamma1
  // Not the most efficient way, but safe. Ok by now.
  // DC remains constant
  // Clear TempMat first
  TempMat[0] = one;
  for (int i=1; i < nosamples; i++)
	{
    TempMat[i] = zero;
    delGamma1[nosamples + i] = zero;
  }
  freq2Time(TempMat, delGamma1);
  TempMat[0] = zero;
  for (int i=2; i < max_calc; i+=2)
	{
    // Real part
    TempMat[i+1] = cos(omega[i>>1] * t);
    freq2Time(TempMat, delGamma1 + i * nosamples);
    TempMat[i+1] = zero;
    // Imaginary part
    TempMat[i] = -sin(omega[i>>1] * t);
    freq2Time(TempMat, delGamma1 + (i+1) * nosamples);
    TempMat[i] = zero;
  }
}


void FreqDomainSV::buildJac(const int& n_states,
double**& Jac, DenseDoubleVector*& mv_jac,
double*& Gamma)
{
  // For each matrix, multiply each row of Gamma (store in TempMat).
  // Then, apply FFT to each column
  for (int i=0; i < n_states; i++)
	{
    for (int j=0; j< n_states; j++)
	  {
		  // Row selector is k
	  	for (int k=0; k < nosamples; k++)
			{
			  for (int l=0; l < nuf2; l++)
				{
				  TempMat[l*nosamples + k] = Gamma[l*nosamples + k] *
					mv_jac[max_nstates*i+j][k];
				}
			}
		  // Now FFT TempMat and store in JacVp
		  for (int l=0; l < nuf2; l++)
			{
			  time2Freq(TempMat + l * nosamples,
				Jac[i * max_nstates + j] + l * nosamples);
			}
	  }
	}
}

void FreqDomainSV::addJac(const int& n_states, const int& js,
double**& Jac, DenseDoubleVector*& mv_jac,
double*& Gamma)
{
  // The same, but add result to Jacobian
  for (int i=0; i < n_states; i++)
	{
    // Row selector is k
    for (int k=0; k < nosamples; k++)
		{
      for (int l=0; l < nuf2; l++)
			{
				TempMat[l*nosamples + k] = Gamma[l*nosamples + k] *
				mv_jac[max_nstates*i+js][k];
			}
		}
    // Now FFT TempMat and store in JacVp
    for (int l=0; l < nuf2; l++)
		{
      time2FreqAdd(TempMat + l * nosamples,
			Jac[i * max_nstates + js] + l * nosamples);
		}
  }
}


void FreqDomainSV::freq2Time(double* ivf, double* ovt)
{
  // Input vector is in nr complex format.
  // Assume nosamples is a power of 2
  // Copy the frequency elements to a double vector with the
  // order used by fftw (halfcomplex array).
  // Normalize now if required.
  fd_p[0] = ivf[0];
  fd_p[nosamples>>1] = ivf[1] * 0.5;
  for (int idx = 1; idx < (nosamples>>1); idx++)
	{
    fd_p[idx] = ivf[idx<<1] * 0.5;
    fd_p[nosamples - idx] = ivf[(idx<<1) + 1] * 0.5;
  }
  // Perform the IFFT
  rfftw_one(pf2t, fd_p, ovt);
}

void FreqDomainSV::time2Freq(double* ivt, double* ovf)
{
  // Perform the FFT
  rfftw_one(pt2f, ivt, fd_p);

  // Assume nosamples is a power of 2
  // Copy the frequency elements to a double vector with the
  // order used by fftw (halfcomplex array).
  // Normalize now if required.
  ovf[0] = fd_p[0] / nosamples;
  ovf[1] = fd_p[nosamples>>1] * norm_f;
  for (int idx=1; idx < (nosamples>>1); idx++)
	{
    ovf[idx<<1] = fd_p[idx] * norm_f;
    ovf[(idx<<1) + 1] = fd_p[nosamples - idx] * norm_f;
  }
  // Output vector is in nr complex format.
}


void FreqDomainSV::time2FreqAdd(double* ivt, double* ovf)
{
  // Perform the FFT
  rfftw_one(pt2f, ivt, fd_p);

  // Assume nosamples is a power of 2
  // Copy the frequency elements to a double vector with the
  // order used by fftw (halfcomplex array).
  // Normalize now if required.
  ovf[0] += fd_p[0] / nosamples;
  ovf[1] += fd_p[nosamples>>1] * norm_f;
  for (int idx=1; idx < (nosamples>>1); idx++)
	{
    ovf[idx<<1] += fd_p[idx] * norm_f;
    ovf[(idx<<1) + 1] += fd_p[nosamples - idx] * norm_f;
  }
  // Output vector is in nr complex format.
}


const void FreqDomainSV::deriv(double *x)
{
  double xtmp;
  // Now multiply each phasor by jw. x[0] is the DC component.
  x[0] = 0.0;

  // x[1] is the real part of the highest frequency phasor. As the
	// imaginary part is zero, the real part of the derivative with
	// respect to time is also zero.
	x[1] = 0.0;
  for (int i = 1; i < nosamples/2 ; i++)
	{
    // x[(i<<1)] is the real part and x[(i<<1)+1] is the imaginary part
		// of the current phasor.
    xtmp = -omega[i] * x[(i<<1) + 1];
    x[(i<<1) + 1] = omega[i] * x[(i<<1)];
    x[(i<<1)] = xtmp;
  }
  return;
}


const void FreqDomainSV::delay(double* x, const double& Time)
{
  // We have to multiply each phasor by
	// exp(-jw Time) = cos(w Time) - j sin(w Time)
	// The DC component doesn't change.

  // The real part of the highest frequency x[1] should always be zero.
	// We have no way to delay it because we lack the imaginary component.
  x[1] = zero;

  // Now, calculate the rest of the frequency components.
  for (int i=1; i < nosamples/2 ; i++)
	{
    // x[(i<<1)] is the real part and x[(i<<1)+1] is the imaginary part
		// of the current phasor.
    double retmp = cos(omega[i] * Time);
    double imtmp = -sin(omega[i] * Time);
    double xtmp = x[(i<<1)];

    x[(i<<1)] = xtmp * retmp - x[(i<<1)+1] * imtmp;
    x[(i<<1)+1] = x[(i<<1)+1] * retmp + xtmp * imtmp;
  }
  return;
}

