#ifndef ultrix
#include <stdlib.h>
#endif
#include <math.h>
#include <stdio.h>

/*
 *  an_math_op.c--
 *
 */

#include "../compat/stuff.h"
#include "../compat/dsp.h"
#include "../compat/dcx.h"
#include "../compat/ml.h"
#include "oper.h"

/*
 * An_Op_MINLMT--
 *	return a real scalar or vector whose minimum value is limited
 */
int An_Op_MINLMT(capOutReq_t *arg, capOutReq_t *lmtarg, capOutReq_t *result)
{
  switch (arg->type) {
  case CAP_OBJ_DOUBLE:
    if (arg->obj.d < lmtarg->obj.d) 
      result->obj.d = lmtarg->obj.d;
    else
      result->obj.d = arg->obj.d;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DoublevMinLmt(arg->size, result->obj.dv,
		      arg->obj.dv, lmtarg->obj.d);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;

  return 0;
}


/*
 * An_Op_MAXLMT--
 *	return a real scalar or vector whose maximum value is limited
 */
int An_Op_MAXLMT(capOutReq_t *arg, capOutReq_t *lmtarg, capOutReq_t *result)
{
  switch (arg->type) {
  case CAP_OBJ_DOUBLE:
    if (arg->obj.d > lmtarg->obj.d) 
      result->obj.d = lmtarg->obj.d;
    else
      result->obj.d = arg->obj.d;
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DoublevMaxLmt(arg->size, result->obj.dv,
		      arg->obj.dv, lmtarg->obj.d);
    result->type = CAP_OBJ_DOUBLEV;
    An_CopyTimebase(result, arg);
    break;
  }
  result->size = arg->size;

  return 0;
}



/*
 * An_Op_SUM--
 *	return the sums (un-normalized integral, sort of) 
 * of elements of a vector
 */
int An_Op_SUM(capOutReq_t *arg, capOutReq_t *result)
{

  result->obj.dv = Mlib_DNewVec(arg->size);
  Dsp_DoublevSum(arg->size, result->obj.dv, arg->obj.dv);

  result->size = arg->size;
  result->type = CAP_OBJ_DOUBLEV;
  An_CopyTimebase(result, arg);

  return 0;
}


