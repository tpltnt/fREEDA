#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/*
 *  an_op_util.c--
 *
 */

#include "../compat/stuff.h"
#include "oper.h"
#include "../compat/ml.h"
#include "../compat/dsp.h"
#include "../compat/dcx.h"



/*
 * An_LinInterpPoint--
 *	Given a vector with a timebase (or other x data), and a single x
 * point which may or may not be in the timebase, find a corresponding
 * linearly-interpolated y value.  If the supplied value of x is outside
 * the limits of the known range, return the y value of the closer
 * endpoint.  n is the size of the vectors.
 *
 *	It is assumed that values in the timebase are in ascending
 * order, though they may not be evenly spaced.
 */
double An_LinInterpPoint(int n, doublev_t yv, doublev_t xv, double x)
{
  int i, h, l;

  /*
   * See if we're in range...
   */
  if (x <= xv[0])
    return 0.0;
  else if (x >= xv[n - 1])
    return 0.0;

  /*
   * Find the interpolation points in the vector using a binary search.
   */
  l = 0;
  h = n;
  i = (l + h) / 2;
  while ((h - l) > 1) {
    if (xv[i] == x)
      return yv[i];
    else if (x < xv[i])
      h = i;
    else
      l = i;
    i = (l + h) / 2;
  }

  /*
   * Now, interpolate.
   */
  return yv[l] + (yv[l + 1] - yv[l]) * (x - xv[l]) / (xv[l + 1] - xv[l]);
}

/*
 * An_LinInterpVec--
 *	Given a vector w/ timebase and a "new" timebase, create a new
 * interpolated vector fitted to the new timebase.  Linear interpolation
 * is used.
 */
void An_LinInterpVec(int n, doublev_t yv, doublev_t xv,
		     int newn, doublev_t newyv, doublev_t newxv)
{
  int i;

  for (i = 0; i < newn; i++)
    newyv[i] = An_LinInterpPoint(n, yv, xv, newxv[i]);
}


/*
 * An_LinInterpCPoint--
 *	Same as An_LinInterpPoint, but for complex y values.
 */
dcx_t An_LinInterpCPoint(int n, dcxv_t yv, doublev_t xv, double x)
{
  int i, h, l;
  dcx_t y;

  /*
   * See if we're in range...
   */
  if (x <= xv[0]) {
    y.re = y.im = 0.0;
    return y;
  }
  else if (x >= xv[n - 1]) {
    y.re = y.im = 0.0;
    return y;
  }

  /*
   * Find the interpolation points in the vector using a binary search.
   */
  l = 0;
  h = n;
  i = (l + h) / 2;
  while ((h - l) > 1) {
    if (xv[i] == x)
      return yv[i];
    else if (x < xv[i])
      h = i;
    else
      l = i;
    i = (l + h) / 2;
  }

  /*
   * Now, interpolate.
   */
  y.re = yv[l].re + (yv[l + 1].re - yv[l].re) * 
    (x - xv[l]) / (xv[l + 1] - xv[l]);
  y.im = yv[l].im + (yv[l + 1].im - yv[l].im) * 
    (x - xv[l]) / (xv[l + 1] - xv[l]);
  return y;
}

/*
 * An_LinInterpCVec--
 *	Same as An_LinInterpVec, except for complex y values.
 */
void An_LinInterpCVec(int n, dcxv_t yv, doublev_t xv,
		      int newn, dcxv_t newyv, doublev_t newxv)
{
  int i;

  for (i = 0; i < newn; i++)
    newyv[i] = An_LinInterpCPoint(n, yv, xv, newxv[i]);
}



/*
 * An_CopyTimebase--
 *	Copy a timebase from one record to another if one exists.
 */
void An_CopyTimebase(capOutReq_t *dest, capOutReq_t *src)
{
  if (src->x) {
    dest->x = Mlib_DNewVec(src->size);
    Dsp_DoublevCopy(src->size, dest->x, src->x);
  } else
    dest->x = NULL;
}


/*
 * An_CompareTimebase--
 *	Compare two timebases.
 */
int An_CompareTimebase(capOutReq_t *arg1, capOutReq_t *arg2)
{
  int i;
  double diff, sum;

  /*
    if ((arg1->x && !arg2->x) || (!arg1->x && arg2->x))
    An_CalcError("One argument has timebase, other doesn't");
  */
  if (arg1->size != arg2->size)
    An_CalcError("Timebases are of different length");
  if (!arg1->x)
    return 0;
  for (i = 0; i < arg1->size; i++) {
    if (arg1->x[i] != arg2->x[i])  {
      sum = arg1->x[i] + arg2->x[i];
      diff = fabs(arg1->x[i] - arg2->x[i]);
      if ((sum != 0.0) && ((diff / sum) > 1e-6))
	An_CalcError("Timebases are different \
in value");
    }
  }
  return 0;
}

/*
 * An_Promote--
 *	An_Promote one argument so that both are of the same type.
 */
void An_Promote(capOutReq_t *arg1, capOutReq_t *arg2)
{
  /*
   * Not very sophisticated but who cares...
	 */
  while (arg1->type != arg2->type) {

    if (arg1->type == CAP_OBJ_INT) {
      double x;
      x = arg1->obj.i;
      arg1->obj.d = x;
      arg1->type = CAP_OBJ_DOUBLE;
      continue;
    }

    if (arg2->type == CAP_OBJ_INT) {
      double x;
      x = arg2->obj.i;
      arg2->obj.d = x;
      arg2->type = CAP_OBJ_DOUBLE;
      continue;
    }

    if ((arg1->type == CAP_OBJ_DOUBLE) &&
	(arg2->type == CAP_OBJ_DOUBLEV)) {
      doublev_t dv;
      int i;

      dv = Mlib_DNewVec(arg2->size);
      for (i = 0; i < arg2->size; i++)
	dv[i] = arg1->obj.d;

      arg1->obj.dv = dv;
      arg1->type = CAP_OBJ_DOUBLEV;
      arg1->size = arg2->size;
      An_CopyTimebase(arg1, arg2);
      continue;
    }

    if ((arg2->type == CAP_OBJ_DOUBLE) &&
	(arg1->type == CAP_OBJ_DOUBLEV)) {
      doublev_t dv;
      int i;

      dv = Mlib_DNewVec(arg1->size);
      for (i = 0; i < arg1->size; i++)
	dv[i] = arg2->obj.d;

      arg2->obj.dv = dv;
      arg2->type = CAP_OBJ_DOUBLEV;
      arg2->size = arg1->size;
      An_CopyTimebase(arg2, arg1);
      continue;
    }

    if ((arg1->type == CAP_OBJ_DOUBLE) &&
	(arg2->type == CAP_OBJ_DCXV)) {
      dcxv_t dcxv;
      int i;

      dcxv = Mlib_CNewVec(arg2->size);
      for (i = 0; i < arg2->size; i++) {
	dcxv[i].re = arg1->obj.d;
	dcxv[i].im = 0.0;
      }

      arg1->obj.dcxv = dcxv;
      arg1->type = CAP_OBJ_DCXV;
      arg1->size = arg2->size;
      An_CopyTimebase(arg1, arg2);
      continue;
    }

    if ((arg2->type == CAP_OBJ_DOUBLE) &&
	(arg1->type == CAP_OBJ_DCXV)) {
      dcxv_t dcxv;
      int i;

      dcxv = Mlib_CNewVec(arg1->size);
      for (i = 0; i < arg1->size; i++) {
	dcxv[i].re = arg2->obj.d;
	dcxv[i].im = 0.0;
      }

      arg2->obj.dcxv = dcxv;
      arg2->type = CAP_OBJ_DCXV;
      arg2->size = arg1->size;
      An_CopyTimebase(arg2, arg1);
      continue;
    }

    if ((arg1->type == CAP_OBJ_DOUBLE) &&
	(arg2->type == CAP_OBJ_DCX)) {
      dcx_t dcx;

      dcx.re = arg1->obj.d;
      dcx.im = 0.0;
      arg1->obj.dcx = dcx;
      arg1->type = CAP_OBJ_DCX;
      continue;
    }

    if ((arg2->type == CAP_OBJ_DOUBLE) &&
	(arg1->type == CAP_OBJ_DCX)) {
      dcx_t dcx;

      dcx.re = arg2->obj.d;
      dcx.im = 0.0;
      arg2->obj.dcx = dcx;
      arg2->type = CAP_OBJ_DCX;
      continue;
    }

    if (arg1->type == CAP_OBJ_DOUBLEV) {
      dcxv_t dcxv;
      int i;

      dcxv = Mlib_CNewVec(arg1->size);
      for (i = 0; i < arg1->size; i++) {
	dcxv[i].re = arg1->obj.dv[i];
	dcxv[i].im = 0.0;
      }

      arg1->obj.dcxv = dcxv;
      arg1->type = CAP_OBJ_DCXV;
      continue;
    }

    if (arg2->type == CAP_OBJ_DOUBLEV) {
      dcxv_t dcxv;
      int i;

      dcxv = Mlib_CNewVec(arg2->size);
      for (i = 0; i < arg2->size; i++) {
	dcxv[i].re = arg2->obj.dv[i];
	dcxv[i].im = 0.0;
      }

      arg2->obj.dcxv = dcxv;
      arg2->type = CAP_OBJ_DCXV;
      continue;
    }

    An_CalcError("Unable to promote argument(s) to similar types");

  }

}

