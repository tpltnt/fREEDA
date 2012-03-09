#ifndef ultrix
#include <stdlib.h>
#endif
#include <math.h>
#include <stdio.h>


/*
 *  an_misc_op.c--
 *
 */

#include "parser.h"
#include "../compat/ml.h"
#include "../compat/dsp.h"
#include "../compat/dcx.h"

void FreeOutReq(capOutReq_t *outReq);


/*
 * An_Op_STRIPX--
 *	strip timebase or other x data from object
 */
int An_Op_STRIPX(capOutReq_t *arg, capOutReq_t *result)
{
  switch (arg->type) {
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg->size);
    Dsp_DoublevCopy(arg->size, result->obj.dv, arg->obj.dv);
    break;
  case CAP_OBJ_DCXV:
    result->obj.dv = Mlib_CNewVec(arg->size);
    Dsp_DcxvCopy(arg->size, result->obj.dcxv, arg->obj.dcxv);
    break;
  }
  result->size = arg->size;
  result->x = NULL;
  result->type = arg->type;

  return 0;
}


/*
 * An_Op_SCALEX--
 * multiply timebase or other x data
 * arg1 is a vector
 * arg2 is a scalor
 */
int An_Op_SCALEX(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  int count = arg1->size;
  doublev_t x=arg1->x;
  double scale;

  result->x = Mlib_DNewVec(arg1->size);
  doublev_t r=result->x;
  switch (arg2->type) {
    case CAP_OBJ_DOUBLE:
    scale = arg2->obj.d;
    break;
    case CAP_OBJ_INT:
    scale = arg2->obj.i;
    break;
    default:
    An_CalcError("First argument of scalex must be a scalar\n");
  }

  switch (arg1->type) {
    case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(arg1->size);
    Dsp_DoublevCopy(arg1->size, result->obj.dv, arg1->obj.dv);
    break;
    case CAP_OBJ_DCXV:
    result->obj.dcxv = Mlib_CNewVec(arg1->size);
    Dsp_DcxvCopy(arg1->size, result->obj.dcxv, arg1->obj.dcxv);
    break;
    default:
    An_CalcError("Second argument of last must be a vector\n");
  }

  while (count--) *r++ = *x++ * scale;
  result->size = arg1->size;
  result->type = arg1->type;

  return 0;
}



/*
 * An_Op_LAST --
 *       extract the last 2^k points from a given vector
 *       time data is extracted as well
 *
 * Parameters:
 * arg1 is an integer indicating the last number of points needed
 * arg2 contains the vector and the time vector associated with it
 * the new vector is returned with result
 *
 */

int An_Op_LAST(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  int last, noEntries;
  if (arg1->type == CAP_OBJ_INT)
    last = arg1->obj.i;
  else
    An_CalcError("First argument of last must be an integer\n");
  noEntries = arg2->size;
  if(last > noEntries)
    An_CalcError("Number of entries too small\n");

  switch (arg2->type) {
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(last);
    Dsp_DoublevLastCopy(noEntries, last, result->obj.dv, arg2->obj.dv);
    break;
  case CAP_OBJ_DCXV:
    result->obj.dcxv = Mlib_CNewVec(last);
    Dsp_DcxvLastCopy(noEntries, last, result->obj.dcxv, arg2->obj.dcxv);
    break;
  default:
    An_CalcError("Second argument of last must be a vector\n");
  }

  /* Allocate and fill the time vector */
  result->x = Mlib_DNewVec(last);
  Dsp_DoublevLastCopy(noEntries, last, result->x, arg2->x);

  result->size = last;
  result->type = arg2->type;

  return 0;
}

/*
 * An_Op_REPEAT --
 *       Copies data points from a given vector in a file and then
 *       repeats this n amount of times but concatenates the data.
 *
 * Parameters:
 * arg1 contains the vector to be copied over n times in the same file.
 * arg2 is an integer indicating the numer of times to be repeated
 * The new vector is returned with result. 
 * 
 */

int An_Op_REPEAT(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  int repeat, noEntries;
  if (arg2->type == CAP_OBJ_INT)
    repeat = arg2->obj.i;
  else
    An_CalcError("First argument of last must be an integer\n");
  noEntries = arg1->size;
  
  switch (arg1->type) {
  case CAP_OBJ_DOUBLEV:
    result->obj.dv = Mlib_DNewVec(noEntries*repeat);
    Dsp_DoublevRepeatCopy(noEntries, repeat, result->obj.dv, arg1->obj.dv);
    break;
  case CAP_OBJ_DCXV:
    result->obj.dcxv = Mlib_CNewVec(noEntries*repeat);
    Dsp_DcxvRepeatCopy(noEntries, repeat, result->obj.dcxv, arg1->obj.dcxv);
    break;
  default:
    An_CalcError("Second argument of last must be a vector\n");
  }
  /* Allocate and fill the time vector*/ 
  result->x = Mlib_DNewVec(noEntries*repeat);   
  Dsp_TimeRepeatCopy(noEntries, repeat, result->x, arg1->x); 
  result->size = noEntries*repeat;
  result->type = arg1->type;
  return(0);
}


/*
 * An_Op_GETBIN --
 *    Gets the nth bin of any vector
 * 
 *    arg1 is the vector itself 
 *    arg2 is the specified bin number that is needed
 *    result returns the nth bin
 *
 */
int An_Op_GETBIN(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  if (arg2->type != CAP_OBJ_INT)
    An_CalcError("Second argument of getbin must be an integer\n");

  switch (arg1->type) {
  case CAP_OBJ_DOUBLEV:
    result->obj.d = arg1->obj.dv[arg2->obj.i];
    result->type = CAP_OBJ_DOUBLE;
    break;
  case CAP_OBJ_DCXV:  
    result->obj.dcx = arg1->obj.dcxv[arg2->obj.i];
    result->type = CAP_OBJ_DCX;
    break;
  default:
    An_CalcError("First argument of getbin must be a vector\n");
  }
  return(0);
}


/*
 * An_Op_IMPULSE
 *   Takes two integers which correspond to a row and column in the
 * time domain impedance matrix calculated in state variable transient
 * analysis and returns the corresponding impulse response
 *
 *
 */
/* int An_Op_IMPULSE(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result) */
/* { */
/* int t; */
/* doublev_t xv, dv; */
/* if ((arg1->type != CAP_OBJ_INT) | (arg2->type != CAP_OBJ_INT)) */
/*    An_CalcError("Both arguments of impulse must be an integer\n"); */

/* dv = Mlib_DNewVec(NoTimePoints); */
/* xv = Mlib_DNewVec(NoTimePoints); */

/* for (t = 0; t < NoTimePoints; t++) { */
/*   dv[t] = svtrdata->sv_M[arg1->obj.i][arg2->obj.i].Val_v[t]; */
/*   xv[t] = svtrdata->time_P[t]; */
/* } */

/* result->obj.dv = dv;  */
/* result->size = NoTimePoints; */
/* result->x = xv; */
/* result->type = CAP_OBJ_DOUBLEV; */
/* return 0; */

/* } */

/*
 * An_Op_FIMPULSE
 *   Takes two integers which correspond to a row and column in the
 * frequency  domain impedance matrix calculated in state variable transient
 * analysis and returns the corresponding impulse response
 *
 *
 */
/* int An_Op_FIMPULSE(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result) */
/* { */
/* int f; */
/* doublev_t xv; */
/* dcxv_t  dcxv; */
/* if ((arg1->type != CAP_OBJ_INT) | (arg2->type != CAP_OBJ_INT)) */
/*    An_CalcError("Both arguments of impulse must be an integer\n"); */

/* dcxv = Mlib_CNewVec(NoFreqPoints); */
/* xv = Mlib_DNewVec(NoFreqPoints); */

/* for (f = 0; f < NoFreqPoints; f++) { */
/*   dcxv[f].re = svtrdata->svMatrix[f][arg1->obj.i][arg2->obj.i].re; */
/*   dcxv[f].im = svtrdata->svMatrix[f][arg1->obj.i][arg2->obj.i].im; */
/*   xv[f] = FreqV_P[f]; */
/* } */

/* result->obj.dcxv = dcxv; */
/* result->size = NoFreqPoints; */
/* result->x = xv; */
/* result->type = CAP_OBJ_DCXV; */
/* return 0; */

/* } */



/*
 * An_Op_PACK
 * Takes one or more vector quantities and packs them into a single
 * matrix where each column is one of the original vectors.
 * The stack is processed until there are no vectors left or none that matches
 * the original x vector (if any).
 * The x vector of copied across.
 */

int An_Op_PACK(capOutReq_t *result)
{
  int argCount=0,
    i, j,
    xsize,
    rank;
  doublev_t vector;
  doublem_t matrix;
  capOutReq_t *arg;

  /* Find out how many vectors there are */
  for (i = 0; i < stackDepth; i++) {
    arg = Peek(i);
    if(arg->type == CAP_OBJ_DOUBLEV) //was DCXV
      argCount++;
    else
      break;
  }

  if(argCount==0) {
    An_CalcError("No vectors to pack\n");
    return(1);
  }

  /* Get matrix size and copy into matrix */
  arg = Peek(0);
  result->size = xsize = arg->size;
  result->obj.s_dm.rank = rank = argCount;
  result->x = arg->x; /* Fast copy */
  /* result contains the time pts / dc voltages vector */
  /*arg->x = NULL;*/
  result->type = CAP_OBJ_DOUBLEM;
  matrix = result->obj.s_dm.dm = Mlib_DNewMat(xsize, rank);
  for (j = 0 ; j < rank; j++)
  {
    arg = Peek(rank-1-j);  /* Use forward order of arguments */
    if(arg->size != xsize)
    {
      An_CalcError("Vectors are not of equal length\n");
      return(1);
    }
    vector = arg->obj.dv;
    /* Fill up the matrix first */
    for (i = 0; i < xsize; i++)
    {
      matrix[i][j] = vector[i];
    }
  }

  /* Handle the printing onto the terminal */
  if (anType == 0) /* Transient analysis */
    result->xName = "Time (sec)";
  else if (anType == 1) /* DC analysis */
    result->xName = "DC (volts)";
  else
    ParseError("analysis type can't be determined by parser.");

  printf("\n");
  sepLine();
  printf("%s \t", result->xName);
  
  for(i = 0; i < PRINTMAX; i++)
  {
    if(printHeader[i] != NULL)
    {
      if (isVoltage == 1)
        printf("v %s\t\t", printHeader[i]);
      else if (isCurrent == 1)
        printf("i %s\t\t", printHeader[i]);
    }
  }
  isVoltage = 0; isCurrent = 0;
  printf("\n");
  sepLine();

  /* Clear the print labels, in case there is another line */
  countV = 0;
  int n;
  for(n = 0; n < PRINTMAX; n++)
    printHeader[n] = NULL;

  /* Print the vectors row wise (like SPICE) */
  for (i = 0; i < xsize; i++)
  {
    printf("%e\t", result->x[i]);
    for (j = 0; j < rank; j++)
    {
      printf("%f", matrix[i][j]);
      printf("\t");
    }
    printf("\n");
  }
  sepLine();
  printf("\n");

  /* Get rid of the vector arguments */
  for (j = 0; j < rank; j++) FreeOutReq(Pop());

  return(0);
}






/*
 * An_Op_CAT
 * Takes one or more vector quantities and packs them into a single
 * matrix where each column is one of the original vectors.
 * The stack is processed until there are no vectors left or none that matches
 * the original x vector (if any).
 * The x vector of the vectors is copied across.
 */

int An_Op_CAT(capOutReq_t *arg1,  capOutReq_t *result)
{
  int argCount=0,
    i, j,
    xref,
    rank;

  dcxm_t matrix;
  capOutReq_t *arg;
  /* Find out how many vectors there are */
  for (i = 1; i < stackDepth; i++) {
    arg = Peek(i);
    if(arg->type == CAP_OBJ_DCXV)
      argCount++;
    else
      break;
  }

  if(argCount==0) {
    An_CalcError("No vectors to cat\n");
    return(1);
  }

  /* Get matrix size and copy into matrix */ 
  xref = arg1->obj.i;
  arg = Peek(1);
  result->size =1;
  result->obj.s_dm.rank = rank = argCount;
  result->x = Mlib_DNewVec(1);
  (result->x)[0] = (arg->x)[xref]; /* Fast copy  */
  arg->x = NULL;
  result->type = CAP_OBJ_DCXM;
  matrix = result->obj.s_dcxm.dcxm = Mlib_CNewMat(1,rank);
  for (j = 0 ; j < rank; j++)
    {
      arg = Peek(rank-j);  /* Use forward order of arguements */
      matrix[0][j].re =  (arg->obj.dcxv)[xref].re;
      matrix[0][j].im =  (arg->obj.dcxv)[xref].im;
    }

  /* Get rid if the vector arguements */
  for (j = 0; j < rank; j++) FreeOutReq(Pop());

  return(0);
}


/*
 * An_Op_Resize_and_Copy
 * 
 * Increase the size of a vector by one and append a new value
 *
 * x	 - original vector
 * xsize - size of original x
 * y	 - new value to append
 */
int An_Op_Resize_and_Copy(double **x, int oldSize, double y)
{
  double	*temp_d_P;
  int i, newSize;
  newSize = oldSize + 1;
  if( (temp_d_P = *x) && (oldSize) ) {
    *x = Mlib_DNewVec(newSize);
    for(i=0; i < oldSize; i++) (*x)[i]=temp_d_P[i];
    Mlib_DFreeVec(temp_d_P);
    (*x)[oldSize] = y;
  }
  else {
    *x = Mlib_DNewVec(1);
    (*x)[0] = y;
  }
  return 0;
}



/*
 * An_Op_APPEND --
 *    arg1 is a single complex number or a double number 
 *    arg2 is the name of the variable which holds the vector 
 *    result returns the appended vector 
 * 
 *
 */

int An_Op_APPEND(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t **result)
{

  if (arg2->type != CAP_OBJ_STRING) {
    An_CalcError("The second argument of append should be a string\n");
    return(1);
  }

  if (GetOutReqVariable(arg2->obj.s, result) == FALSE) {
    An_CalcError("Error in append: Request not found in the table\n");
    return(1);
  }

  switch(arg1->type){


  case CAP_OBJ_DOUBLE:
    {
      double	*temp_d_P;
      int i, oldSize, newSize;

      oldSize = ((*result)->size)++;
      newSize = oldSize + 1;
      if( (temp_d_P = (*result)->obj.dv) ) {
	(*result)->obj.dv = Mlib_DNewVec((*result)->size);
	for(i=0; i < oldSize; i++)
	  (*result)->obj.dv[i]=temp_d_P[i];
	Mlib_DFreeVec(temp_d_P);
	(*result)->obj.dv[oldSize] = arg1->obj.d;
      }
      else {
	(*result)->obj.dv = Mlib_DNewVec(1);
	(*result)->obj.dv[0] = arg1->obj.d;
      }


      {
	double newX;
	if(sweepList_head)
	  newX=sweepList_head->current;
	else
	  newX=0.;
        An_Op_Resize_and_Copy(&((*result)->x), oldSize, newX);
      }

      (*result)->type = CAP_OBJ_DOUBLEV; 
    }

    break;

  case CAP_OBJ_DCX:
    {
      dcx_t	*temp_dcx_P;
      /*  	double	*temp_d_P; */
      int i, oldSize, newSize;

      oldSize = ((*result)->size)++;
      newSize = oldSize+1;
      if( (temp_dcx_P = (*result)->obj.dcxv) ) {
	(*result)->obj.dcxv = Mlib_CNewVec((*result)->size);
	for(i=0; i < oldSize; i++)
	  (*result)->obj.dcxv[i]=temp_dcx_P[i];
	Mlib_CFreeVec(temp_dcx_P);
      }
      else {
	(*result)->obj.dcxv = Mlib_CNewVec(1);
      }
      (*result)->obj.dcxv[oldSize].re = arg1->obj.dcx.re;
      (*result)->obj.dcxv[oldSize].im = arg1->obj.dcx.im;

      {
	double newX;
	if(sweepList_head)
	  newX=sweepList_head->current;
	else
	  newX=0.;
        An_Op_Resize_and_Copy(&((*result)->x), oldSize, newX);
      }



      (*result)->type = CAP_OBJ_DCXV; 
    }

    break;

  case CAP_OBJ_DCXM:
    {
      dcxm_t temp_dcxm_P;
      int i, j;

      (*result)->obj.s_dcxm.rank = arg1->obj.s_dcxm.rank;

      if(arg1->size !=1)
	{
	  printf("Matrix X size must be 1 for Append\n");
	}

      ((*result)->size)++;
      if( (temp_dcxm_P = (*result)->obj.s_dcxm.dcxm) ) {
	(*result)->obj.s_dcxm.dcxm = Mlib_CNewMat((*result)->size, (*result)->obj.s_dcxm.rank);
	for(i=0; i < (((*result)->size)-1); i++)
	  {
	    for(j=0; j < (((*result)->obj.s_dcxm.rank)); j++)
	      {
		((*result)->obj.s_dcxm.dcxm)[i][j].re =temp_dcxm_P[i][j].re;
		((*result)->obj.s_dcxm.dcxm)[i][j].im =temp_dcxm_P[i][j].im;
	      }
	  }
	Mlib_CFreeMat(temp_dcxm_P);
	 
	for(j=0; j < (((*result)->obj.s_dcxm.rank)); j++)
	  {
	    ((*result)->obj.s_dcxm.dcxm)[((*result)->size)-1][j].re = 
	      (arg1->obj.s_dcxm.dcxm)[0][j].re;
	    ((*result)->obj.s_dcxm.dcxm)[((*result)->size)-1][j].im = 
	      (arg1->obj.s_dcxm.dcxm)[0][j].im;
	  }
	 
      }
      else {
	(*result)->obj.s_dcxm.dcxm = Mlib_CNewMat(1, (*result)->obj.s_dcxm.rank);
	 
	for(j=0; j < (((*result)->obj.s_dcxm.rank)); j++)
	  {
	    (*result)->obj.s_dcxm.dcxm[0][j].re = arg1->obj.s_dcxm.dcxm[0][j].re;
	    (*result)->obj.s_dcxm.dcxm[0][j].im = arg1->obj.s_dcxm.dcxm[0][j].im;
	  }
      }

	
      {
	double newX;
	if(sweepList_head)
	  newX=sweepList_head->current;
	else
	  newX=0.;
        An_Op_Resize_and_Copy(&((*result)->x), (*result)->size - 1, newX);
      }

      (*result)->type = CAP_OBJ_DCXM; 
    }

    break;

  default:
    printf("Can not append type\n");

  }
	
  return(0);
}

/*
 * An_Op_PUSH
 *
 * 
 * 
 * 
 */
int An_Op_PUSH(capOutReq_t *arg1, capOutReq_t **result) 
{
  capOutReq_Pt storedOutReq;

  if (arg1->type != CAP_OBJ_STRING) {
    An_CalcError("The argument of append should be a string\n");
    return(1);
  }

  if (GetOutReqVariable(arg1->obj.s, &storedOutReq) == FALSE) {
    An_CalcError("Error in push: Request not found in the table\n");
    return(1);
  }

  An_COPY(storedOutReq, result);
  return(0);
}






/*
 * An_Op_XY2CX--
 *	put timebase of real vector into real part, and real (Y) into 
 * imaginary part of a complex vector, without timebase.
 */
int An_Op_XY2CX(capOutReq_t *arg, capOutReq_t *result)
{
  int i;

  if (arg->type != CAP_OBJ_DOUBLEV)
    An_CalcError("Operator XY2CX takes a real vector as an argument");
  if (!arg->x)
    An_CalcError("Operator XY2CX requires argument w/ x data");

  result->obj.dcxv = Mlib_CNewVec(arg->size);
  for (i = 0; i < arg->size; i++) {
    result->obj.dcxv[i].re = arg->x[i];
    result->obj.dcxv[i].im = arg->obj.dv[i];
  }

  result->size = arg->size;
  result->type = CAP_OBJ_DCXV;
  result->x = NULL;

  return 0;
}


/*
 * An_Op_CX2XY--
 *	put real part of complex vector into timebase of real vector, and
 * put imaginary part of complex vector into Y (of real vector)
 */
int An_Op_CX2XY(capOutReq_t *arg, capOutReq_t *result)
{
  int i;

  if (arg->type != CAP_OBJ_DCXV)
    An_CalcError("Operator CX2XY takes a complex vector as an argument");

  result->obj.dv = Mlib_DNewVec(arg->size);
  result->x = Mlib_DNewVec(arg->size);
  for (i = 0; i < arg->size; i++) {
    result->x[i] = arg->obj.dcxv[i].re;
    result->obj.dv[i] = arg->obj.dcxv[i].im;
  }

  result->size = arg->size;
  result->type = CAP_OBJ_DOUBLEV;

  return 0;
}



/*
 * An_Op_SYSTEM--
 *
 * Do a System Call
 */
int An_Op_SYSTEM(capOutReq_t *arg, capOutReq_t *result)
{
  if (arg->type != CAP_OBJ_STRING)
    An_CalcError("Operator 'SYSTEM' requires a string argument");

  system(arg->obj.s);

  return 0;
}



