/*
 * an_optbl.c
 *
 * the operator table and lookup routines
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "parser.h"
#include "oper_net.h"

an_OpTbl_t an_OpTbl[] = {
  /*
   * List all operators here.
   * key:
   * name		identifier	argCt		retVal?	func
   * the types of the arguments
   */
  { "system",	CAP_OP_SYSTEM,	ONE_ARG,	FALSE,	(ifunc) An_Op_SYSTEM,
    {STRING_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "crt_out",	CAP_OP_CRT_OUT,	NO_ARG,		FALSE,	NULL,
    {NO_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "write",	CAP_OP_WRITE,	ONE_ARG,	FALSE,	(ifunc) An_Op_WRITE,
    {STRING_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "plot",	CAP_OP_PLOT,	VAR_ARG,	FALSE,	(ifunc) An_Op_PLOT,
    {NO_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "watch",	CAP_OP_PLOT,	VAR_ARG,	FALSE,	(ifunc) An_Op_WATCH,
    {NO_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "read",	CAP_OP_READ,	ONE_ARG,	TRUE,	(ifunc) An_Op_READ,
    {FILE_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  /* Regular Operators */
  { "dup",	CAP_OP_DUP,	ONE_ARG,	TRUE,	NULL,
    {ANY_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "get",	CAP_OP_GET,	TWO_ARG,	TRUE,	NULL,
    {NUM_SV_ARG_TYPE, INT_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "put",	CAP_OP_PUT,	THREE_ARG,	TRUE,	NULL,
    {NUM_SV_ARG_TYPE, INT_ARG_TYPE,	NUM_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "stripx",	CAP_OP_STRIPX, 	ONE_ARG,	TRUE,	(ifunc) An_Op_STRIPX,
    {VECTOR_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "pack",	CAP_OP_PACK, 	VAR_ARG,	TRUE,	(ifunc) An_Op_PACK,
    {NO_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "append",CAP_OP_APPEND, TWO_ARG_RESULT_POINTER,TRUE,	(ifunc) An_Op_APPEND,
    {DCX_SV_ARG_TYPE,STRING_ARG_TYPE,NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "push",	CAP_OP_PUSH, ONE_ARG_RESULT_POINTER,	TRUE,(ifunc) An_Op_PUSH,
    {STRING_ARG_TYPE,NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "cat",	CAP_OP_CAT, 		ONE_ARG,	TRUE,	(ifunc) An_Op_CAT,
    {INT_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},

  /*
   * Include network specific operators unless this is to be a standalone
   * equation version.
   */
#include "an_optbl_net.cx"

  { "real",	CAP_OP_REAL,	ONE_ARG,	TRUE,	(ifunc) An_Op_REAL,
    {DCX_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "imag",	CAP_OP_IMAG,	ONE_ARG,	TRUE,	(ifunc) An_Op_IMAG,
    {DCX_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "mag",	CAP_OP_MAG,	ONE_ARG,	TRUE,	(ifunc) An_Op_MAG,
    {DCX_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "contphase",	CAP_OP_CONTPHASE, ONE_ARG,	TRUE,	(ifunc) An_Op_CONTPHASE,
    {DCX_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "prinphase",	CAP_OP_PRINPHASE, ONE_ARG,	TRUE,	(ifunc) An_Op_PRINPHASE,
    {DCX_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "conj",	CAP_OP_CONJ,	ONE_ARG,	TRUE,	(ifunc) An_Op_CONJ,
    {DCX_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "xy2cx",	CAP_OP_XY2CX,	ONE_ARG,	TRUE,	(ifunc) An_Op_XY2CX,
    {DOUBLEV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "cx2xy",	CAP_OP_CX2XY,	ONE_ARG,	TRUE,	(ifunc) An_Op_CX2XY,
    {DCXV_ARG_TYPE, 	NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE}
    ,1}, 
  { "recip",	CAP_OP_RECIP,	ONE_ARG,	TRUE,	(ifunc) An_Op_RECIP,
    {NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "neg",	CAP_OP_NEG,	ONE_ARG,	TRUE,	NULL,
    {NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "db",		CAP_OP_DB,	ONE_ARG,	TRUE,	(ifunc)An_Op_DB20,
    {NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "db10",	CAP_OP_DB10,	ONE_ARG,	TRUE,  (ifunc)An_Op_DB10,
    {NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "rad2deg",	CAP_OP_RAD2DEG, ONE_ARG,	TRUE,	(ifunc)An_Op_rad2deg,
    {NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "deg2rad",	CAP_OP_DEG2RAD, ONE_ARG,	TRUE,	NULL,
    {NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "abs",	CAP_OP_ABS,	ONE_ARG,	TRUE,	(ifunc) An_Op_ABS,
    {NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "minlmt",	CAP_OP_MINLMT, 	TWO_ARG, 	TRUE,	(ifunc) An_Op_MINLMT,
    {REAL_SV_ARG_TYPE, DOUBLE_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1}, 
  { "maxlmt",	CAP_OP_MAXLMT, 	TWO_ARG, 	TRUE,	(ifunc) An_Op_MAXLMT,
    {REAL_SV_ARG_TYPE, DOUBLE_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1}, 
  { "add",	CAP_OP_ADD,	TWO_ARG,	TRUE,	(ifunc) An_Op_ADD,
    {NUM_SV_ARG_TYPE, NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "sub",	CAP_OP_SUB,	TWO_ARG,	TRUE,	(ifunc) An_Op_SUB,
    {NUM_SV_ARG_TYPE, NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "mult",	CAP_OP_MULT,	TWO_ARG,	TRUE,	(ifunc) An_Op_MULT,
    {NUM_SV_ARG_TYPE, NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "div",	CAP_OP_DIV,	TWO_ARG,	TRUE,	(ifunc) An_Op_DIV,
    {NUM_SV_ARG_TYPE, NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "diff",	CAP_OP_DIFF,	ONE_ARG,	TRUE,	NULL,
    {VECTOR_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "deriv",	CAP_OP_DERIV, 	ONE_ARG,	TRUE,	NULL,
    {VECTOR_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "sum",	CAP_OP_SUM,	ONE_ARG,	TRUE,	(ifunc) An_Op_SUM,
    {VECTOR_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "integ",	CAP_OP_INTEG,	ONE_ARG,	TRUE,	NULL,
    {VECTOR_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "smpltime",	CAP_OP_SMPLTIME, NO_ARG,	TRUE,	(ifunc) An_Op_SMPLTIME,
    {NO_ARG_TYPE, 	NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "sweepfreq",	CAP_OP_SWEEPFREQ, NO_ARG,	TRUE,	(ifunc) An_Op_SWEEPFREQ,
    {NO_ARG_TYPE, 	NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "smplcvt",	CAP_OP_SMPLCVT, TWO_ARG,	TRUE,	(ifunc) An_Op_SMPLCVT,
    {VECTOR_ARG_TYPE, VECTOR_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "sweepcvt",	CAP_OP_SWEEPCVT, TWO_ARG,	TRUE,	(ifunc) An_Op_SWEEPCVT,
    {VECTOR_ARG_TYPE, VECTOR_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "maketime",	CAP_OP_MAKETIME, TWO_ARG,	TRUE,	(ifunc) An_Op_MAKETIME,
    {REAL_ARG_TYPE,	INT_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "scalex", CAP_OP_SCALEX, TWO_ARG, TRUE, (ifunc) An_Op_SCALEX,
    {VECTOR_ARG_TYPE, REAL_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1},
  { "makesweep",	CAP_OP_MAKESWEEP, TWO_ARG,	TRUE,	(ifunc) An_Op_MAKESWEEP,
    {REAL_ARG_TYPE,	INT_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "repeat",	CAP_OP_REPEAT, 	TWO_ARG, 	TRUE,	(ifunc) An_Op_REPEAT,
    {VECTOR_ARG_TYPE, INT_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1},	
  { "getbin",	CAP_OP_GETBIN, 	TWO_ARG, 	TRUE,	(ifunc) An_Op_GETBIN,
    {VECTOR_ARG_TYPE, INT_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1},	
  { "zeropad", CAP_OP_ZEROPAD, ONE_ARG, TRUE, (ifunc) An_Op_ZEROPAD,
    {DOUBLEV_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1},
  { "getlast2n",CAP_OP_GETLAST2N,ONE_ARG, TRUE, (ifunc) An_Op_GETLAST2N,
    {DOUBLEV_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1},
  { "fft",	CAP_OP_FFT,	ONE_ARG,	TRUE,	(ifunc) An_Op_FFT,
    {DOUBLEV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "invfft",	CAP_OP_INVFFT, 	ONE_ARG,	TRUE,	(ifunc) An_Op_INVFFT,
    {DCXV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE,	NO_ARG_TYPE} ,1},
  { "sconv",	CAP_OP_SCONV, 	TWO_ARG,	TRUE,	NULL,
    {DOUBLEV_ARG_TYPE, DOUBLEV_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1},
  { "cconv",	CAP_OP_CCONV, 	TWO_ARG,	TRUE,	(ifunc) An_Op_CCONV,
    {DOUBLEV_ARG_TYPE, DOUBLEV_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1},
  { "upcconv",	CAP_OP_UPCCONV, TWO_ARG,	TRUE,	NULL,
    {DOUBLEV_ARG_TYPE, DOUBLEV_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1},
  { "fconv",	CAP_OP_FCONV, 	TWO_ARG, 	TRUE,	NULL,
    {DOUBLEV_ARG_TYPE, DOUBLEV_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1},
  { "last",	CAP_OP_LAST, 	TWO_ARG, 	TRUE,	(ifunc) An_Op_LAST,
    {INT_ARG_TYPE, VECTOR_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1},
  { "rect2polar",	CAP_OP_RECT2POLAR, 	ONE_ARG, TRUE,(ifunc) An_Op_RECT2POLAR,
    {DCX_SV_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1},
  { "polar2polar",CAP_OP_POLAR2RECT, 	ONE_ARG,TRUE,(ifunc) An_Op_POLAR2RECT,
    {DCX_SV_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE, NO_ARG_TYPE} ,1},
  { "lpbwfrq",	CAP_OP_LPBWFRQ,	THREE_ARG,	TRUE,	An_Op_LPBWFRQ,
    {VECTOR_ARG_TYPE, DOUBLE_ARG_TYPE, INT_ARG_TYPE, NO_ARG_TYPE} ,1},



  /*
   * Place synonym operator names here, after a full unduplicated list of
   * all operators.  Make sure the rest of the data is the same!
   */
  { "+",		CAP_OP_ADD,	TWO_ARG,	TRUE,	(ifunc) An_Op_ADD,
    {NUM_SV_ARG_TYPE, NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "-",		CAP_OP_SUB,	TWO_ARG,	TRUE,	(ifunc) An_Op_SUB,
    {NUM_SV_ARG_TYPE, NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "*",		CAP_OP_MULT,	TWO_ARG,	TRUE,	(ifunc) An_Op_MULT,
    {NUM_SV_ARG_TYPE, NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "/",		CAP_OP_DIV,	TWO_ARG,	TRUE,	(ifunc) An_Op_DIV,
    {NUM_SV_ARG_TYPE, NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE} ,1}, 
  { "^",		CAP_OP_NONE,	TWO_ARG,	TRUE,	(ifunc) An_Op_NOOP,
    {NUM_SV_ARG_TYPE, NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE} ,7}, 
  { "**",		CAP_OP_NONE,	TWO_ARG,	TRUE,	(ifunc) An_Op_NOOP,
    {NUM_SV_ARG_TYPE, NUM_SV_ARG_TYPE, NO_ARG_TYPE,	NO_ARG_TYPE} ,7}, 

  /*
   * Leave a null entry here.
   */
  { "",		0,		0,	FALSE,	NULL,
    {0L,		0L,		0L,		0L} ,0}
};


/*
 * An_LookupOpNum--
 * 	Look up an operator in the opTbl (by its constant).  Returns -1 if
 * not found, else the index in the table.
 */
int An_LookupOpNum(int opConst)
{
  int i;

  for (i = 0; i < sizeof(an_OpTbl)/sizeof(an_OpTbl_t); i++) {
    if (opConst == an_OpTbl[i].opConst)
      return i;
  }
  return -1;
}

/*
 * An_LookupOpName--
 * 	Look up an operator in the opTbl (by its name).  Returns -1 if
 * not found, else the corresponding opConst.
 */
int An_LookupOpName(char *name)
{
  int i;

  for (i = 0; i < sizeof(an_OpTbl)/sizeof(an_OpTbl_t); i++) {
    if (!strcmp(name, an_OpTbl[i].name))
      return an_OpTbl[i].opConst;
  }

  return -1;
}

