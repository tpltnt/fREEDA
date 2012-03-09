
/*
 *  dsp.c--
 *
 */

#include <stdio.h>
#include <math.h>

#ifndef PI
#define PI 3.14159265358979323
#endif

#ifndef RAD2DEG
#define RAD2DEG 57.2957795131
#endif

#include "../compat/mltypes.h"
#include "dsp.h"
#include <rfftw.h>

/* hypot doesn't seems to work well */
#define hypot(x,y) sqrt((x)*(x) + (y)*(y))


/*
 * Copying and arithmetic for vectors of doubles
 */
void Dsp_DoublevCopy(int n, doublev_t r, doublev_t x)
{ while (n--) *r++ = *x++; }

/*
 * Copies the last specified number of elements of a double vector
 * into a new result vector 
 * int n The number of entries
 * int l The last 2^k number of entries 
 * doublev_t r The result vector
 * doublev_t x The vector whose last 2^k elements are copied into r
 */
void Dsp_DoublevLastCopy(int n, int l, doublev_t r, doublev_t x)
{
  x+=(n-l);
  while (l--) *r++ = *x++; 
}
/* Dsp_DoublevRepeatCopy
 *	This funtion concatenates a data into a file a specified number of times.
 *
 * Parameters:
 * 	n is the number of data lines
 *	l is the number of times the data will be concatenate in to the file
 *	doublev_t r The result vector
 *	c resets the number of entries after the while loop
 *	doublev_t x The vector whose elements are concatenated into r
 */ 
void Dsp_DoublevRepeatCopy(int n, int l, doublev_t r, doublev_t x)
{
  int i;
  int j;
  int t;
  t = 0;
  for(i=0;l>i;i++)
    {
      for(j=0;n>j;j++){
        r[j+t] = x[j];
      }
      t+=n;      
    }
}
/* Dsp_TimeRepeatCopy
 *	This funtion concatenates time data into a file a specified number 
 *	of times but increaments the value of each run by the number of
 *	lines aready written.
 *
 * Parameters:
 * 	n is the number of data lines
 *	l is the number of times the data will be concatenate in to the file
 *	doublev_t r The result vector
 *	c resets the number of entries after the while loop
 *	doublev_t x The vector whose elements are concatenated into r
 */ 
void Dsp_TimeRepeatCopy(int n, int l, doublev_t r, doublev_t x)
{
  int i;
  int j;
  
  
  
  for(i=0;l>i;i++)
    {
      for(j=0;n>j;j++)
	{
	  r[j+i*n] = x[j] + i*((x[1] - x[0]) *n);
	}
          
    }
}
void Dsp_Double2DoublevCopy(int n, doublev_t r, double x )
{   

  r[n+1] = x;
       
}

void Dsp_Doublev2Dcxv(int n, dcxv_t r, doublev_t x)
{ while (n--) r->re = *x++, r->im = 0.0, r++; }

void Dsp_DoublevAdd(int n, doublev_t r, doublev_t x, doublev_t y)
{ while (n--) *r++ = *x++ + *y++; }

void Dsp_DoublevSub(int n, doublev_t r, doublev_t x, doublev_t y)
{ while (n--) *r++ = *x++ - *y++; }

void Dsp_DoublevMult(int n, doublev_t r, doublev_t x, doublev_t y)
{ while (n--) *r++ = *x++ * *y++; }

void Dsp_DoublevDiv(int n, doublev_t r, doublev_t x, doublev_t y)
{ while (n--) *r++ = *x++ / *y++; }

/*
 * ...with scalars...
 */
void Dsp_DoublevSclAdd(int n, doublev_t r, doublev_t x, double s)
{ while (n--) *r++ = *x++ + s; }

void Dsp_DoublevSclSub(int n, doublev_t r, doublev_t x, double s)
{ while (n--) *r++ = *x++ - s; }

void Dsp_DoublevSclMult(int n, doublev_t r, doublev_t x, double s)
{ while (n--) *r++ = *x++ * s; }

void Dsp_DoublevSclDiv(int n, doublev_t r, doublev_t x, double s)
{ while (n--) *r++ = *x++ / s; }

/*
 * ...ordinary functions...
 */
void Dsp_DoublevRecip(int n, doublev_t r, doublev_t x)
{ while (n--) *r++ = 1.0 / *x++; }

void Dsp_DoublevNeg(int n, doublev_t r, doublev_t x)
{ while (n--) *r++ = -*x++; }

void Dsp_DoublevDB(int n, doublev_t r, doublev_t x)
{ while (n--) *r++ = log10(*x++) * 20.0; }

void Dsp_DoublevDB10(int n, doublev_t r, doublev_t x)
{ while (n--) *r++ = log10(*x++) * 10.0; }

void Dsp_DoublevRad2Deg(int n, doublev_t r, doublev_t x)
{ while (n--) *r++ = *x++ * RAD2DEG; }

void Dsp_DoublevDeg2Rad(int n, doublev_t r, doublev_t x)
{ while (n--) *r++ = *x++ / RAD2DEG; }

void Dsp_DoublevAbs(int n, doublev_t r, doublev_t x)
{ while (n--) *r++ = fabs(*x++); }

void Dsp_DoublevMinLmt(int n, doublev_t r, doublev_t x, double minlmt)
{ while (n--) *r++ = (*x > minlmt ? *x++ : (x++, minlmt)); }

void Dsp_DoublevMaxLmt(int n, doublev_t r, doublev_t x, double maxlmt)
{ while (n--) *r++ = (*x < maxlmt ? *x++ : (x++, maxlmt)); }

void Dsp_DoublevLog(int n, doublev_t r, doublev_t x)
{ while (n--) *r++ = log(*x++); }

void Dsp_DoublevLog10(int n, doublev_t r, doublev_t x)
{ while (n--) *r++ = log10(*x++); }

void Dsp_DoublevExp(int n, doublev_t r, doublev_t x)
{ while (n--) *r++ = exp(*x++); }

void Dsp_DoublevExp10(int n, doublev_t r, doublev_t x)
{ while (n--) *r++ = pow(*x++, 10.0); }


void Dsp_DoublevSum(int n, doublev_t r, doublev_t x)
{
  int i;

  r[0] = 0;
  for (i = 1; i < n; i++)
    r[i] = r[i - 1] + x[i - 1];
}


/*
 * Functions on complex vectors returning vectors of doubles
 *
 * Dsp_DcxvReal		return the real part of a complex vector
 * Dsp_DcxvImag		return the imaginary part of a complex vector
 * Dsp_DcxvMag		return the magnitude of a complex vector
 * Dsp_DcxvDB		return 20 * log10(magnitude) of a complex vector
 * Dsp_DcxvPrinPhase	return the principal value phase of a complex vector
 * Dsp_DcxvContPhase	return the continuous phase of a complex vector
 */

void Dsp_DcxvReal(int n, doublev_t r, dcxv_t x)
{ for (; n--; r++, x++) *r = x->re; }

void Dsp_DcxvImag(int n, doublev_t r, dcxv_t x)
{ for (; n--; r++, x++) *r = x->im; }

void Dsp_DcxvMag(int n, doublev_t r, dcxv_t x)
{ for (; n--; r++, x++) *r = hypot(x->re, x->im); }

void Dsp_DcxvDB(int n, doublev_t r, dcxv_t x)
{ for (; n--; r++, x++) *r = 20.0 * log10(hypot(x->re, x->im)); }

void Dsp_DcxvDB10(int n, doublev_t r, dcxv_t x)
{ for (; n--; r++, x++) *r = 10.0 * log10(hypot(x->re, x->im)); }

void Dsp_DcxvPrinPhase(int n, doublev_t r, dcxv_t x)
{ 
  for (; n--; r++, x++)
    if (x->re == 0.0)
      if (x->im == 0.0)
	*r = PI / 2.0;
      else
	*r = -PI / 2.0;
    else
      *r = atan2(x->im, x->re);
}

void Dsp_DcxvContPhase(int n, doublev_t r, dcxv_t x)
{
  int i;
  double diff, offset;

  offset = 0;
  for (i = 0; i < n; i++, x++, r++) {
    if (x->re == 0.0)
      if (x->im == 0.0)
	*r = PI / 2.0 + offset;
      else
	*r = -PI / 2.0 + offset;
    else
      *r = atan2(x->im, x->re) + offset;

    if (i) {
      if ((diff = (*r - *(r - 1))) > PI) {
	offset -= 2.0 * PI;
	*r -= 2.0 * PI;
      } else if (diff < -PI) {
	offset += 2.0 * PI;
	*r += 2.0 * PI;
      }
    }
			
  }
}


/*
 * Copying and arithmetic for complex vectors
 */
void Dsp_DcxvCopy(int n, dcxv_t r, dcxv_t x)
{ while (n--) *r++ = *x++; }

/*
 * Copies the last specified number of elements of a double 
 * complex vector into a new result vector
 */
void Dsp_DcxvLastCopy(int n, int l, dcxv_t r, dcxv_t x)
{
  x+= (n-l);
  while (n--) *r++ = *x++;
}
/* Dsp_DcxvRepeatCopy
 *	This funtion concatenates a data into a file a specified number of times.
 *
 * Parameters:
 * 	n is the number of data lines
 *	l is the number of times the data will be concatenate in to the file
 *	doublev_t r The result vector
 *	c resets the number of entries after the while loop
 *	doublev_t x The vector whose elements are concatenated into r
 */ 
void Dsp_DcxvRepeatCopy(int n, int l, dcxv_t r, dcxv_t x)
{   
  int i;
  int j;
  int t;
  t = 0;
  for(i=0;l>i;i++)
    {
      for(j=0;n>j;j++){
        r[j+t] = x[j];
      }
      t+=n;      
    }
}
void Dsp_Dcx2DcxvCopy(int n, dcxv_t r, dcx_t x)
{   

  r[n+1] = x;
       
}
void Dsp_DcxvNeg(int n, dcxv_t r, dcxv_t x)
{ for (; n--; r++, x++) { r->re = -x->re; r->im = -x->im; }}

void Dsp_DcxvRecip(int n, dcxv_t r, dcxv_t x)
{ double q; for (; n--; r++, x++) {
  q = x->re * x->re + x->im * x->im;
  r->re = x->re / q; 
  r->im = -x->im / q;
} }

void Dsp_DcxvConj(int n, dcxv_t r, dcxv_t x)
{ for (; n--; r++, x++) { r->re = x->re; r->im = -x->im; }}

void Dsp_DcxvAdd(int n, dcxv_t r, dcxv_t x, dcxv_t y)
{ for (; n--; r++, x++, y++) { r->re = x->re + y->re; r->im = x->im + y->im; }}

void Dsp_DcxvSub(int n, dcxv_t r, dcxv_t x, dcxv_t y)
{ for (; n--; r++, x++, y++) { r->re = x->re - y->re; r->im = x->im - y->im; }}

void Dsp_DcxvMult(int n, dcxv_t r, dcxv_t x, dcxv_t y)
{ for (; n--; r++, x++, y++) {
  r->re = x->re * y->re - x->im * y->im; 
  r->im = x->re * y->im + x->im * y->re;
} }

void Dsp_DcxvDiv(int n, dcxv_t r, dcxv_t x, dcxv_t y)
{ double q; for (; n--; r++, x++, y++) {
  q = y->re * y->re + y->im * y->im;
  r->re = (x->re * y->re + x->im * y->im) / q; 
  r->im = (x->im * y->re - x->re * y->im) / q;
} }

/*
 * ...with complex scalars...
 */
void Dsp_DcxvSclAdd(int n, dcxv_t r, dcxv_t x, dcx_t s)
{ for (; n--; r++, x++) { r->re = x->re + s.re; r->im = x->im + s.im; }}

void Dsp_DcxvSclSub(int n, dcxv_t r, dcxv_t x, dcx_t s)
{ for (; n--; r++, x++) { r->re = x->re - s.re; r->im = x->im - s.im; }}

void Dsp_DcxvSclMult(int n, dcxv_t r, dcxv_t x, dcx_t s)
{ for (; n--; r++, x++) {
  r->re = x->re * s.re - x->im * s.im; 
  r->im = x->re * s.im + x->im * s.re;
} }

void Dsp_DcxvSclDiv(int n, dcxv_t r, dcxv_t x, dcx_t s)
{ double q; for (; n--; r++, x++) {
  q = s.re * s.re + s.im * s.im;
  r->re = (x->re * s.re + x->im * s.im) / q; 
  r->im = (x->im * s.re - x->re * s.im) / q;
} }



/*
 * Fourier transforms
 */

/*
 * Real-data forward transform.  There should be 2^order real points in
 * x.  There will be 2^order complex points returned in r.
 */
void Dsp_DoublevFFT(int order, dcxv_t r, doublev_t x)
{
  int i, samples;
  double *tmpv1;
  rfftw_plan p;
  
  samples = 1 << order;
  
  tmpv1 = (double *) malloc(sizeof(double) * samples);
  
  p = rfftw_create_plan(samples,
			FFTW_REAL_TO_COMPLEX,
			FFTW_ESTIMATE);
  rfftw_one(p, x, tmpv1);
  rfftw_destroy_plan(p);
  
  /*
   * Copy the result from the temporary vector into the result vector.
   * Unpack the DC and Nyquist components. 
   * Multiply by 2 / samples (single-sided spectrum).
   */
  r[0].re = tmpv1[0] / samples;
  r[0].im = 0.;
  r[samples >> 1].re = tmpv1[samples >> 1] * 2. / samples;
  r[samples >> 1].im = 0.0;
  for (i = 1; i < (samples >> 1); i++) {
    r[i].re = tmpv1[i] * 2. / samples;
    r[i].im = tmpv1[samples-i] * 2. / samples;
  }
}

/*
 * Real-data inverse transform.  There should be 2^(order - 1) + 1 complex 
 * points in x.  There will be 2^order real points returned in r.
 */
void Dsp_DcxvInvFFT(int order, doublev_t r, dcxv_t x)
{
  int i, samples, nfreq;
  double *tmpv1;
  rfftw_plan p;
  
  samples = 1 << order;
  nfreq = samples >> 1;
  tmpv1 = (double *) malloc(sizeof(double) * samples);
  
  /*
   * Copy the response into the result vector in preparation for the
   * real-data fft ...
   * Divide by two since the fft assume double-sided spectrum.
   */
  tmpv1[0] = x[0].re;
  tmpv1[nfreq] = x[nfreq].re / 2.;
  for (i = 1; i < nfreq; i++) {
    tmpv1[i] = x[i].re / 2.;
    tmpv1[samples - i] = x[i].im / 2.;
  }
  
  p = rfftw_create_plan(samples,
			FFTW_COMPLEX_TO_REAL,
			FFTW_ESTIMATE);
  rfftw_one(p, tmpv1, r);
  rfftw_destroy_plan(p);

}




/*
 * Dsp_LPBwthRsp returns the butterworth LPF filter response
 */

void Dsp_LPBwthRsp(int n, dcxv_t response, doublev_t frq, 
		   double corner, int order)
{

  /*
   * n		- number of frequency points
   * response	- the response
   * frq		- frequency points (supplied)
   * corner	- corner frequency of filter
   * order	- order of filter (only order 1 and 2 supported)
   */

  /*
   *    Transfer function of a first order Butterworth Low Pass Filter
   *      H(f) = 1./( 1 + (jf/fc))      j = sqrt(-1)
   *      H(f) = c - j*b*c where b = f/fc
   *                             c = 1./(1. + b*b)
   *    Transfer function of a second order Butterworth Low Pass Filter
   *      H(f) = 1./( 1. + 1.414(jf/fc) + (jf/fc)**2 )      j = sqrt(-1)
   *      i.e. H(f) = 1./( 1 - (f/fc)**2 + 1.414(jf/fc) )      j = sqrt(-1)
   *      i.e. H(f) = a*c - j*b*c  where d = f/fc
   *                                     a = 1. - d*d
   *                                     b = 1.414*d
   *                                     c = 1./(a*a + b*b)
   */

  double a, b, c, d;	/* temporary variables */
  int i;


  if (corner == 0.0) {
    printf("Dsp_ButterworthLPRsp: zero corner frequency \
illegal.\n");
    exit(1);
  }

  switch (order) {
  case 1:
    for (i = 0; i < n; i++) {
      b = frq[i] / corner;
      c = 1.0 / (1.0 + b * b);
      response[i].re = c;
      response[i].im = -b * c;
    }
    break;
  case 2:
    for (i = 0; i < n; i++) {
      d = frq[i] / corner;
      a = 1.0 - d * d;
      b = 1.414 * d;
      c = 1.0 / (a * a + b * b);
      response[i].re = a * c;
      response[i].im = -b * c;
    }
    break;
  default:
    fprintf(stderr, "Only orders 1 and 2 supported for LP \
Butterworth filter\n");
    exit(1);
  }

  return;
}
