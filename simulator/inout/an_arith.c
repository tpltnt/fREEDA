/*
 * an_arith_op.c
 */

#ifndef ultrix
#include <stdlib.h>
#endif
#include <math.h>
#include <stdio.h>

#include "../compat/stuff.h"
#include "oper.h"
#include "../compat/dsp.h"
#include "../compat/dcx.h"
#include "../compat/ml.h"

/* hypot doesn't seems to work well */
#define hypot(x,y) sqrt((x)*(x) + (y)*(y))

/*
 * An_Op_REAL--
 *	return the real part of a complex scalar or number
 */
int An_Op_REAL(capOutReq_t *arg, capOutReq_t *result)
{
  switch (arg->type) {
  case CAP_OBJ_DCX:
    result->obj.d = arg->obj.dcx.re;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DCXV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DcxvReal(arg->size, result->obj.dv, arg->obj.dcxv);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;

  return 0;
}

/*
 * An_Op_DB20--
 *	return the real part of a complex scalar or number
 */
int An_Op_DB20(capOutReq_t *arg, capOutReq_t *result)
{
  int i;
  switch (arg->type) {
  case CAP_OBJ_DOUBLE:
    result->obj.d = log10(arg->obj.d)*20.;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DoublevDB10(arg->size, result->obj.dv, arg->obj.dv);
    for(i=0; i< arg->size; i++)
      result->obj.dv[i] *= 2.;
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;
  
  return 0;
}

/*
 * An_Op_DB10--
 *	return the real part of a complex scalar or number
 */
int An_Op_DB10(capOutReq_t *arg, capOutReq_t *result)
{
  switch (arg->type) {
  case CAP_OBJ_DOUBLE:
    result->obj.d = log10(arg->obj.d)*10.;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DoublevDB10(arg->size, result->obj.dv, arg->obj.dv);	
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;
  
  return 0;
}



/*
 * An_Op_rad2deg--
 *	convert from radians to degree
 */
int An_Op_rad2deg(capOutReq_t *arg, capOutReq_t *result)
{
  const double r2d = 180. / PI;
  int i;

  switch (arg->type) {
  case CAP_OBJ_DOUBLE:
    result->obj.d = r2d * arg->obj.d;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    for(i=0; i< arg->size; i++)
      result->obj.dv[i] = r2d * arg->obj.dv[i];
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;
  
  return 0;
}





/*
 * An_Op_IMAG--
 *	return the imaginary part of a complex vector.
 */
int An_Op_IMAG(capOutReq_t *arg, capOutReq_t *result)
{
  switch (arg->type) {
  case CAP_OBJ_DCX:
    result->obj.d = arg->obj.dcx.im;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DCXV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DcxvImag(arg->size, result->obj.dv, arg->obj.dcxv);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;

  return 0;
}


/*
 * An_Op_CONJ--
 *	return the conjugate of a complex scalar or vector.
 */
int An_Op_CONJ(capOutReq_t *arg, capOutReq_t *result)
{
  switch (arg->type) {
  case CAP_OBJ_DCX:
    result->obj.dcx.re = arg->obj.dcx.re;
    result->obj.dcx.im = -arg->obj.dcx.im;
    result->type = CAP_OBJ_DCX;
    break;
  case CAP_OBJ_DCXV:
    result->obj.dcxv = Mlib_CNewVec(arg->size);
    Dsp_DcxvConj(arg->size, result->obj.dcxv, arg->obj.dcxv);
    result->type = CAP_OBJ_DCXV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;

  return 0;
}


/*
 * An_Op_ABS--
 *	return the absolute value of a real scalar or vector; returns the
 * magnitude of a complex scalar or vector
 */
int An_Op_ABS(capOutReq_t *arg, capOutReq_t *result)
{
  switch (arg->type) {
  case CAP_OBJ_INT:
    result->obj.i = (arg->obj.i < 0) ? -arg->obj.i : arg->obj.i;
    result->type = CAP_OBJ_INT;
    break;
  case CAP_OBJ_DOUBLE:
    result->obj.d = (arg->obj.d < 0) ? -arg->obj.d : arg->obj.d;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DCX:
    result->obj.d = hypot(arg->obj.dcx.re, arg->obj.dcx.im);
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DoublevAbs(arg->size, result->obj.dv, arg->obj.dv);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  case CAP_OBJ_DCXV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DcxvMag(arg->size, result->obj.dv, arg->obj.dcxv);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;

  return 0;
}



/*
 * An_Op_MAG--
 *	return the magnitude of a complex scalar or vector.
 */
int An_Op_MAG(capOutReq_t *arg, capOutReq_t *result)
{
  switch (arg->type) {
  case CAP_OBJ_DCX:
    result->obj.d = hypot(arg->obj.dcx.re, arg->obj.dcx.im);
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DCXV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DcxvMag(arg->size, result->obj.dv, arg->obj.dcxv);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;

  return 0;
}


/*
 * An_Op_PRINPHASE--
 *	return the principal phase of a complex vector.
 */
int An_Op_PRINPHASE(capOutReq_t *arg, capOutReq_t *result)
{
  switch (arg->type) {
  case CAP_OBJ_DCX:
    if (arg->obj.dcx.re == 0.0)
      if (arg->obj.dcx.im >= 0.0)
	result->obj.d = PI / 2.0;
      else
	result->obj.d = -PI / 2.0;
    else
      result->obj.d = atan2(arg->obj.dcx.im, arg->obj.dcx.re);
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DCXV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DcxvPrinPhase(arg->size, result->obj.dv, arg->obj.dcxv);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;

  return 0;
}


/*
 * An_Op_CONTPHASE--
 *	return the continuous phase of a complex vector.
 */
int An_Op_CONTPHASE(capOutReq_t *arg, capOutReq_t *result)
{
  switch (arg->type) {
  case CAP_OBJ_DCX:
    if (arg->obj.dcx.re == 0.0)
      if (arg->obj.dcx.im >= 0.0)
	result->obj.d = PI / 2.0;
      else
	result->obj.d = -PI / 2.0;
    else
      result->obj.d = atan2(arg->obj.dcx.im, arg->obj.dcx.re);
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DCXV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DcxvContPhase(arg->size, result->obj.dv, arg->obj.dcxv);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;

  return 0;
}


/*
 * An_Op_ADD--
 *	return the sum of two vectors.
 */
int An_Op_ADD(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  if (arg1->type != arg2->type)
    An_Promote(arg1, arg2);

  if (arg1->type != arg2->type)
    An_CalcError("Arguments of ADD are incompatible or can't be \
promoted");

  An_CompareTimebase(arg1, arg2);

  switch (arg1->type) {
  case CAP_OBJ_INT:
    result->obj.i = arg1->obj.i + arg2->obj.i;
    result->type = CAP_OBJ_INT;
    break;
  case CAP_OBJ_DOUBLE:
    result->obj.d = arg1->obj.d + arg2->obj.d;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DCX:
    result->obj.dcx.re = arg1->obj.dcx.re + arg2->obj.dcx.re;
    result->obj.dcx.im = arg1->obj.dcx.im + arg2->obj.dcx.im;
    result->type = CAP_OBJ_DCX;
    break;
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg1->size);
    Dsp_DoublevAdd(arg1->size, result->obj.dv, 
		   arg1->obj.dv, arg2->obj.dv);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg1);
    break;
  case CAP_OBJ_DCXV:
    result->obj.dcxv = Mlib_CNewVec(arg1->size);
    Dsp_DcxvAdd(arg1->size, result->obj.dcxv, 
		arg1->obj.dcxv, arg2->obj.dcxv);
    result->type = CAP_OBJ_DCXV;
    An_CopyTimebase(result, arg1);
    break;
  default:
    An_CalcError("Argument type not allowed for ADD\n");
  }

  result->size = arg1->size;

  return 0;
}


/*
 * An_Op_SUB--
 *	return the difference of two vectors.
 */
int An_Op_SUB(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  if (arg1->type != arg2->type)
    An_Promote(arg1, arg2);

  if (arg1->type != arg2->type)
    An_CalcError("Arguments of SUB are incompatible or can't be \
promoted");

  An_CompareTimebase(arg1, arg2);

  switch (arg1->type) {
  case CAP_OBJ_INT:
    result->obj.i = arg1->obj.i - arg2->obj.i;
    result->type = CAP_OBJ_INT;
    break;
  case CAP_OBJ_DOUBLE:
    result->obj.d = arg1->obj.d - arg2->obj.d;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DCX:
    result->obj.dcx.re = arg1->obj.dcx.re - arg2->obj.dcx.re;
    result->obj.dcx.im = arg1->obj.dcx.im - arg2->obj.dcx.im;
    result->type = CAP_OBJ_DCX;
    break;
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg1->size);
    Dsp_DoublevSub(arg1->size, result->obj.dv, 
		   arg1->obj.dv, arg2->obj.dv);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg1);
    break;
  case CAP_OBJ_DCXV:
    result->obj.dcxv = Mlib_CNewVec(arg1->size);
    Dsp_DcxvSub(arg1->size, result->obj.dcxv, 
		arg1->obj.dcxv, arg2->obj.dcxv);
    result->type = CAP_OBJ_DCXV;
    An_CopyTimebase(result, arg1);
    break;
  default:
    An_CalcError("Argument type not allowed for SUB\n");
  }

  result->size = arg1->size;

  return 0;
}


/*
 * An_Op_MULT--
 *	return the product of two vectors.
 */
int An_Op_MULT(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  if (arg1->type != arg2->type)
    An_Promote(arg1, arg2);

  if (arg1->type != arg2->type)
    An_CalcError("Arguments of MULT are incompatible or can't be \
promoted");

  An_CompareTimebase(arg1, arg2);

  switch (arg1->type) {
  case CAP_OBJ_INT:
    result->obj.i = arg1->obj.i * arg2->obj.i;
    result->type = CAP_OBJ_INT;
    break;
  case CAP_OBJ_DOUBLE:
    result->obj.d = arg1->obj.d * arg2->obj.d;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DCX: 
    result->obj.dcx = Cx_DMult(arg1->obj.dcx, arg2->obj.dcx);
    result->type = CAP_OBJ_DCX;
    break;
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg1->size);
    Dsp_DoublevMult(arg1->size, result->obj.dv, 
		    arg1->obj.dv, arg2->obj.dv);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg1);
    break;
  case CAP_OBJ_DCXV:
    result->obj.dcxv = Mlib_CNewVec(arg1->size);
    Dsp_DcxvMult(arg1->size, result->obj.dcxv, 
		 arg1->obj.dcxv, arg2->obj.dcxv);
    result->type = CAP_OBJ_DCXV;
    An_CopyTimebase(result, arg1);
    break;
  default:
    An_CalcError("Argument type not allowed for MULT\n");
  }

  result->size = arg1->size;

  return 0;
}


/*
 * An_Op_DIV--
 *	return the division of two vectors.
 */
int An_Op_DIV(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  if (arg1->type != arg2->type)
    An_Promote(arg1, arg2);

  if (arg1->type != arg2->type)
    An_CalcError("Arguments of DIV are incompatible or can't be \
promoted");

  An_CompareTimebase(arg1, arg2);

  switch (arg1->type) {
  case CAP_OBJ_INT:
    result->obj.i = arg1->obj.i / arg2->obj.i;
    result->type = CAP_OBJ_INT;
    break;
  case CAP_OBJ_DOUBLE:
    result->obj.d = arg1->obj.d / arg2->obj.d;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DCX:
    result->obj.dcx = Cx_DDiv(arg1->obj.dcx, arg2->obj.dcx);
    result->type = CAP_OBJ_DCX;
    break;
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg1->size);
    Dsp_DoublevDiv(arg1->size, result->obj.dv, 
		   arg1->obj.dv, arg2->obj.dv);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg1);
    break;
  case CAP_OBJ_DCXV:
    result->obj.dcxv = Mlib_CNewVec(arg1->size);
    Dsp_DcxvDiv(arg1->size, result->obj.dcxv, 
		arg1->obj.dcxv, arg2->obj.dcxv);
    result->type = CAP_OBJ_DCXV;
    An_CopyTimebase(result, arg1);
    break;
  default:
    An_CalcError("Argument type not allowed for DIV\n");
  }

  result->size = arg1->size;

  return 0;
}


/*
 * An_Op_RECIP--
 */
int An_Op_RECIP(capOutReq_t *arg, capOutReq_t *result)
{
  double q;

  switch (arg->type) {
  case CAP_OBJ_INT:
    result->obj.d = 1.0 / (double) arg->obj.i;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DOUBLE:
    result->obj.d = 1.0 / arg->obj.d;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DCX:
    q = arg->obj.dcx.re * arg->obj.dcx.re + 
      arg->obj.dcx.im * arg->obj.dcx.im;
    result->obj.dcx.re = arg->obj.dcx.re / q;
    result->obj.dcx.im = -arg->obj.dcx.im / q;
    result->type = CAP_OBJ_DCX;
    break;
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DoublevRecip(arg->size, result->obj.dv, arg->obj.dv);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  case CAP_OBJ_DCXV:
    result->obj.dcxv = Mlib_CNewVec(arg->size);
    Dsp_DcxvRecip(arg->size, result->obj.dcxv, arg->obj.dcxv);
    result->type = CAP_OBJ_DCXV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;

  return 0;
}



/*
 * An_Op_NOOP--
 *	Non funtional operator (No Op)
 */
int An_Op_NOOP(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  printf("Error: Operator under development\n");
  result->type = CAP_OBJ_NONEXISTENT;
  return 0;
}
