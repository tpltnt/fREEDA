/*
 * an_dsp_op.c
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "../compat/stuff.h"
#include "oper.h"
#include "../compat/ml.h"
#include "../compat/dsp.h"
#include "../compat/dcx.h"
#include "ftvec.h"


void Dsp_DcxvInvFFT(int order, doublev_t r, dcxv_t x);

/*
 * An_Op_SMPLTIME--
 *	return a vector of time points (in x and y) representing the
 * current analysis settings.
 */
int An_Op_SMPLTIME(capOutReq_t *result)
{
  int i, order, pts;
  double x, totLen;

  /*
   * Calculate the appropriate number of points (what INVFFT
   * would return).
   */
  for (x = NoFreqPoints - 1, order = -1; x > 0.5; x /= 2.0, order++)
    ;
  pts = 2 << order;

  /*
   * Allocate space and create two vectors of timebase points.
   */
  result->obj.dv = Mlib_DNewVec(2 << order);
  result->x = Mlib_DNewVec(2 << order);
  totLen = (double)(1 << order) / FreqV_P[1 << order];
  for (i = 0; i < pts; i++)
    result->obj.dv[i] = result->x[i] = totLen * i / pts;

  result->type = CAP_OBJ_DOUBLEV;
  result->size = pts;

  return 0;
}


/*
 * An_Op_SWEEPFREQ--
 *	return a vector of frequency points (in x and y) representing the
 * current analysis settings.
 */
int An_Op_SWEEPFREQ(capOutReq_t *result)
{
  int i;

  /*
   * Allocate space and create two vectors of timebase points.
   */
  result->obj.dv = Mlib_DNewVec(NoFreqPoints);
  result->x = Mlib_DNewVec(NoFreqPoints);
  for (i = 0; i < NoFreqPoints; i++)
    result->obj.dv[i] = result->x[i] = FreqV_P[i];

  result->type = CAP_OBJ_DOUBLEV;
  result->size = NoFreqPoints;

  return 0;
}


/*
 * An_Op_LPBWFRQ--
 *	given a vector w/ x data (assumed to be frequencies), returns
 * a complex vector w/ frequency response corresponding to the supplied
 * corner frequency and filter order.
 */
int An_Op_LPBWFRQ(capOutReq_t *sfarg, capOutReq_t *corner,
		  capOutReq_t *order, capOutReq_t *result)
{

  if (!sfarg->x)
    An_CalcError("Operator 'lpbwfrq' requires a vector argument \
with sweep frequencies");
  if (corner->obj.d <= 0.0)
    An_CalcError("Operator 'lpbwfrq' requires corner frequency \
> 0.0");
  if ((order->obj.i < 1) || (order->obj.i > 2))
    An_CalcError("Operator 'lpbwfrq' requires order = 1 or 2");

  result->obj.dcxv = Mlib_CNewVec(sfarg->size);
  Dsp_LPBwthRsp(sfarg->size, result->obj.dcxv, sfarg->x,
		corner->obj.d, order->obj.i);

  result->type = CAP_OBJ_DCXV;
  result->size = sfarg->size;
  An_CopyTimebase(result, sfarg);

  return 0;
}


/*
 * An_Op_MAKESWEEP--
 *	return a vector of frequency points (in both x and y).  The lowest
 * frequency point returned will be 0.0; the highest will be the value of
 * frqarg.  ptsarg points will be created.
 */
int An_Op_MAKESWEEP(capOutReq_t *frqarg, 
		    capOutReq_t *ptsarg, capOutReq_t *result)
{
  int i;

  if (frqarg->obj.d <= 0.0)
    An_CalcError("Operator 'MAKESWEEP' requires positive \
frequency argument");
  if (ptsarg->obj.i < 2)
    An_CalcError("Operator 'MAKESWEEP' require points \
argument > 1");

  result->obj.dv = Mlib_DNewVec(ptsarg->obj.i);
  result->x = Mlib_DNewVec(ptsarg->obj.i);

  for (i = 0; i < ptsarg->obj.i; i++)
    result->obj.dv[i] = result->x[i] =
      frqarg->obj.d * 
      (double) i / ((double) ptsarg->obj.i - 1);

  result->type = CAP_OBJ_DOUBLEV;
  result->size = ptsarg->obj.i;

  return 0;
}

/*
 * An_Op_MAKETIME--
 *	return a vector of time points (in both x and y).  The lowest
 * time point returned will be 0.0; the highest will be the value of
 * lenarg.  ptsarg points will be created.
 */
int An_Op_MAKETIME(capOutReq_t *lenarg, 
		   capOutReq_t *ptsarg, capOutReq_t *result)
{
  int i;

  if (lenarg->obj.d <= 0.0)
    An_CalcError("Operator 'makesweep' requires positive \
time argument");
  if (ptsarg->obj.i < 2)
    An_CalcError("Operator 'makesweep' require points \
argument > 1");

  result->obj.dv = Mlib_DNewVec(ptsarg->obj.i);
  result->x = Mlib_DNewVec(ptsarg->obj.i);

  for (i = 0; i < ptsarg->obj.i; i++)
    result->obj.dv[i] = result->x[i] =
      lenarg->obj.d * 
      (double) i / ((double) ptsarg->obj.i - 1);

  result->type = CAP_OBJ_DOUBLEV;
  result->size = ptsarg->obj.i;

  return 0;
}


/*
 * An_Op_SMPLCVT
 *	Sampling rate conversion.  
 */
int An_Op_SMPLCVT(capOutReq_t *srcarg, 
		  capOutReq_t *tbarg, capOutReq_t *result)
{
  if (!srcarg->x || !tbarg->x)
    An_CalcError("Both arguments to 'smplcvt' must have \
frequency data");

  if (srcarg->type == CAP_OBJ_DCXV) {	
    result->obj.dcxv = Mlib_CNewVec(tbarg->size);
    An_LinInterpCVec(srcarg->size, srcarg->obj.dcxv, srcarg->x,
		     tbarg->size, result->obj.dcxv, tbarg->x);
    result->type = CAP_OBJ_DCXV;
  } else {
    result->obj.dv = Mlib_DNewVec(tbarg->size);
    An_LinInterpVec(srcarg->size, srcarg->obj.dv, srcarg->x,
		    tbarg->size, result->obj.dv, tbarg->x);
    result->type = CAP_OBJ_DOUBLEV;
  }

  result->size = tbarg->size;
  An_CopyTimebase(result, tbarg);

  return 0;
}


/*
 * An_Op_SWEEPCVT--
 *	Sweep frequency conversion.  
 */
int An_Op_SWEEPCVT(capOutReq_t *srcarg, 
		   capOutReq_t *sfarg, capOutReq_t *result)
{
  if (!srcarg->x || !sfarg->x)
    An_CalcError("Both arguments to 'SWEEPCVT' must have \
frequency data");

  if (srcarg->type == CAP_OBJ_DCXV) {	
    result->obj.dcxv = Mlib_CNewVec(sfarg->size);
    An_LinInterpCVec(srcarg->size, srcarg->obj.dcxv, srcarg->x,
		     sfarg->size, result->obj.dcxv, sfarg->x);
    result->type = CAP_OBJ_DCXV;
  } else {
    result->obj.dv = Mlib_DNewVec(sfarg->size);
    An_LinInterpVec(srcarg->size, srcarg->obj.dv, srcarg->x,
		    sfarg->size, result->obj.dv, sfarg->x);
    result->type = CAP_OBJ_DOUBLEV;
  }

  result->size = sfarg->size;
  An_CopyTimebase(result, sfarg);

  return 0;
}

/*
 * An_Op_ZEROPAD--
 * Zero pad vector to to take size of vector to a power 2 in length.
 */
int An_Op_ZEROPAD(capOutReq_t *arg, capOutReq_t *result)
{
  int order, i;
  double x, deltax;

  if (arg->type != CAP_OBJ_DOUBLEV)
    An_CalcError("Operator 'ZEROPAD' requires a real vector argument");

  // Find size = a power of 2 points. Allocate result vector.
  for (x = arg->size, order = -1; x > 0.5; x /= 2.0, order++);

  result->size = 1 << order;
  result->obj.dv = Mlib_DNewVec(result->size);

  // Zero-pad if necessary.
  for (i = 0; i < arg->size; i++)
    result->obj.dv[i] = arg->obj.dv[i];
  for (i = arg->size; i < result->size; i++)
    result->obj.dv[i] = 0.0;

  // Construct x vector if an argument x vector was present.
  if (arg->x)
  {
    deltax = (arg->x[(arg->size)-1] - arg->x[0]) / ((arg->size)-1);
    result->x = Mlib_DNewVec(result->size);
    for (i = 0; i < arg->size; i++)
      result->x[i] = arg->x[i];
    for (i = arg->size; i <= result->size; i++)
      result->x[i] = result->x[i-1]+deltax;
  }

  result->type = CAP_OBJ_DOUBLEV;

  return 0;
}


/*
 * An_Op_GETLAST2N--
 * Keep the last 2^N points
 */
int An_Op_GETLAST2N(capOutReq_t *arg, capOutReq_t *result)
{
  int order, i, j;
  double x;

  if (arg->type != CAP_OBJ_DOUBLEV)
    An_CalcError("Operator 'GETLAST2N' requires a real vector argument");

  // Find size = a power of 2 points. Allocate result vector.
  for (x = arg->size, order = -1; x > 0.5; x /= 2.0, order++);

  if(arg->size >= 1<< (order))
    result->size = 1<<(order);
  else
    result->size = 1<<(order-1);

  result->obj.dv = Mlib_DNewVec(result->size);

  for(i = result->size , j = arg->size ; i ; i--, j--)
    result->obj.dv[i-1] = arg->obj.dv[j-1];

  if (arg->x)
  {
    result->x = Mlib_DNewVec(result->size);
    for(i = result->size , j = arg->size ; i ; i--, j--)
      result->x[i-1] = arg->x[j-1];
  }
  else
    result->x = 0;

  result->type = CAP_OBJ_DOUBLEV;

  return 0;
}



/*
 * An_Op_FFT--
 *	Calculate the forward fast Fourier transform of a vector of
 * real data.
 *	There is no check to insure that the timebase supplied is
 * evenly spaced.  You can expect some pretty weird results if it isn't!
 */
int An_Op_FFT(capOutReq_t *arg, capOutReq_t *result)
{
  doublev_t tmpv;
  dcxv_t dcxv;
  int order, i;
  double x, maxFrq;

  if (arg->type != CAP_OBJ_DOUBLEV)
    An_CalcError("Operator 'FFT' requires a real vector argument");

  /*
   * Find size = a power of 2 points.  Allocate temporary vector.
   */
  for (x = arg->size, order = -1; x > 0.5; x /= 2.0, order++)
    ;
  tmpv = Mlib_DNewVec(1 << order);
  Dsp_DoublevCopy(arg->size, tmpv, arg->obj.dv);

  /*
   * Zero-pad if necessary, printing a warning message.  We could
   * try interpolating...
   */	
  if (arg->size != 1 << order) {
    fprintf(output_F, "Note: zero-padding FFT data--wasn't power of \
2 in length--expect trouble!\n");
    for (i = arg->size; i < (1 << order); i++)
      tmpv[i] = 0.0;
  }

  /*
   * Allocate space for result and do the FFT...
   */
  dcxv = Mlib_CNewVec((1 << (order - 1)) + 1);
  Dsp_DoublevFFT(order, dcxv, tmpv);

  /*
   * Construct a vector of sweep points, if a timebase was present.
   */
  if (arg->x) {
    result->x = Mlib_DNewVec((1 << (order - 1)) + 1);
    maxFrq = 1.0 / (2.0 * (arg->x[1] - arg->x[0]));
    for (i = 0; i <= (1 << (order - 1)); i++)
      result->x[i] = maxFrq * i / (1 << (order - 1));
  } else
    result->x = NULL;

  /*
   * Save the result.
   */
  result->type = CAP_OBJ_DCXV;
  result->size = (1 << (order - 1)) + 1;
  result->obj.dcxv = dcxv;

  Mlib_DFreeVec(tmpv);

  return 0;
}


/*
 * An_Op_INVFFT
 */
int An_Op_INVFFT(capOutReq_t *arg, capOutReq_t *result)
{
  dcxv_t tmpv;
  int order, i;
  double x, totLen;
  doublev_t dv;

  if (arg->type != CAP_OBJ_DCXV)
    An_CalcError("Operator 'INVFFT' requires a \
complex vector argument");

  /*
   * Find size = a power of 2 points.  Allocate temporary vector.
   */
  for (x = arg->size - 1, order = -1; x > 0.5; x /= 2.0, order++)
    ;
  tmpv = Mlib_CNewVec((1 << order) + 1);
  Dsp_DcxvCopy(arg->size, tmpv, arg->obj.dcxv);

  /*
   * Zero-pad if necessary, printing a warning message.  We could
   * try interpolating...
   */	
  if (arg->size != ((1 << order) + 1)) {
    fprintf(output_F, "Note: zero-padding INVFFT data--wasn't \
2^n + 1 in length--expect trouble!\n");
    for (i = arg->size; i <= (1 << (order)); i++)
      tmpv[i].re = tmpv[i].im = 0.0;
  }

  /*
   * Allocate space for result and do the inverse FFT...
   */
  dv = Mlib_DNewVec(2 << order);
  Dsp_DcxvInvFFT(order + 1, dv, tmpv);

  /*
   * Construct a vector of timebase points, if a sweep was present.
   */
  if (arg->x) {
    result->x = Mlib_DNewVec(2 << order);
    totLen = (1 << order) / arg->x[1 << order];
    for (i = 0; i < (2 << order); i++)
      result->x[i] = totLen * i / (2 << order);
  } else
    result->x = NULL;

  /*
   * Save the result.
   */
  result->type = CAP_OBJ_DOUBLEV;
  result->size = 2 << order;
  result->obj.dv = dv;

  Mlib_CFreeVec(tmpv);

  return 0;
}


/*
 * An_Op_CCONV--
 *	Use FFT techniques to convolve two signals ("circular convolution").
 * The signal and response are zero-padded to insure that the convolution
 * does not wrap around.
 *	There is no check to insure that the timebase supplied is
 * evenly spaced.  You can expect some pretty weird results if it isn't!
 */
int An_Op_CCONV(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  doublev_t signal1v, signal2v;
  int order, i;
  double x; /*  , maxFrq; */
  dcxv_t frq1v, frq2v;

  /*
   * Both arguments should be the same size.  If the user wants to
   * convolve signals of different sizes, sampling rates, etc.,
   * he'll have to do the conversion himself.
   */
  An_CompareTimebase(arg1, arg2);

  /*
   * Find size = a power of 2 points.  Allocate temporary vector of
   * double length (to allow for padding).
   */
  for (x = arg1->size, order = -1; x > 0.5; x /= 2.0, order++)
    ;
  signal1v = Mlib_DNewVec(2 << order);
  signal2v = Mlib_DNewVec(2 << order);
  Dsp_DoublevCopy(arg1->size, signal1v, arg1->obj.dv);
  Dsp_DoublevCopy(arg2->size, signal2v, arg2->obj.dv);

  /*
   * Zero-pad to a power of two, then add that power of two more
   * zeros.  There is no warning for odd lengths, since it shouldn't
   * make any difference.
   */	
  if (arg1->size != 1 << order) {
    for (i = arg1->size; i < (1 << order); i++)
      signal1v[i] = 0.0;
    signal2v[i] = 0.0;
  }
  for (i = (1 << order); i < (2 << order); i++) {
    signal1v[i] = 0.0;
    signal2v[i] = 0.0;
  }

  /*
   * Allocate space for result and do two forward FFTs...
   */
  frq1v = Mlib_CNewVec((1 << order) + 1);
  frq2v = Mlib_CNewVec((1 << order) + 1);
  Dsp_DoublevFFT(order + 1, frq1v, signal1v);
  Dsp_DoublevFFT(order + 1, frq2v, signal2v);

  /*
   * Multiply in the frequency domain, putting result in frq1v...
   */
  Dsp_DcxvMult((1 << order) + 1, frq1v, frq1v, frq2v);

  /*
   * Now, compute the reverse transform of frq1v...
   */
  Dsp_DcxvInvFFT(order + 1, signal1v, frq1v);

  /*
   * Truncate and save the result.
   */
  result->obj.dv = Mlib_DNewVec(1 << order);
  Dsp_DoublevCopy(1 << order, result->obj.dv, signal1v);
  result->type = CAP_OBJ_DOUBLEV;
  result->size = 1 << order;
  return 0;
}


int An_Op_RECT2POLAR(capOutReq_t *arg1, capOutReq_t *result)
{
  int i;

  switch (arg1->type) {
  case CAP_OBJ_DCX:
    result->obj.dcx.re = Cx_DAbs(arg1->obj.dcx);
    result->obj.dcx.im = Cx_DAng(arg1->obj.dcx);
    result->type = CAP_OBJ_DCX;
    break;
  case CAP_OBJ_DCXV:
    An_CopyTimebase(result, arg1);
    result->obj.dcxv = Mlib_CNewVec(arg1->size);
    for (i=0; i<arg1->size;i++)  
      {
	result->obj.dcxv[i].re = Cx_DAbs(arg1->obj.dcxv[i]);
	result->obj.dcxv[i].im = Cx_DAng(arg1->obj.dcxv[i]);
      }
    result->type = CAP_OBJ_DCXV;
    result->size = arg1->size;
    break;
  }
  return(0);
}

int An_Op_POLAR2RECT(capOutReq_t *arg1, capOutReq_t *result)
{

  double theta;

  theta = arg1->obj.dcx.im;
  theta = 180.0/(3.1415926535897932384626433) * theta;

  switch (arg1->type) {
  case CAP_OBJ_DCX:
    result->obj.dcx.re = (arg1->obj.dcx.re)*cos(theta); 
    result->obj.dcx.im = (arg1->obj.dcx.re)*sin(theta); 
    break;
  }  
  return(0);
}


