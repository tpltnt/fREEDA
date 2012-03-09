/*
 * ml.c
 *
 * A library of functions for real and complex matrix math.
 *
 * Authors:
 * Mark S. Basel
 * Michael B. Steer
 * 
 * Catagories:  (search for these strings to locate the beginning)
 *				Integer matrix and vector functions.
 *				Double matrix, vector and scalar functions.
 *				Complex matrix and vector functions.
 *				Float matrix and vector functions.
 * 				Index of matrix and vector functions.
 * 				Char matrix and vector functions.
 *                             
 *
 *
 * Index of matrix and vector functions:
 *   Double	 	Complex	  	   Float		 Integer
 *   ------		-------	 	   -----		 -------
 * Mlib_DNewVec		Mlib_CNewVec	   Mlib_FNewVec	        Mlib_INewVec
 * Mlib_DFreeVe		Mlib_CFreeVec	   Mlib_FFreeVec	Mlib_ITypeVec
 * Mlib_DNewMat		Mlib_CNewMat	   Mlib_FNewMat	        Mlib_INewMat
 * Mlib_DFreeMat	Mlib_CFreeMat	   Mlib_FFreeMat	Mlib_IFreeMat
 * Mlib_DInv		Mlib_CInv                               Mlib_INewArray
 * Mlib_DNew3DMat	Mlib_CNew3DMat	   Mlib_FNew3DMat	Mlib_INew3DMat
 * Mlib_DFree3DMat	Mlib_CFree3DMat	   Mlib_FFree3DMat      Mlib_IFree3DMat
 * Mlib_DNew4DMat		           Mlib_FNew4DMat	Mlib_INew4DMat
 * Mlib_DFree4DMat		           Mlib_FFree4DMat      Mlib_IFree4DMat
 *                      Mlib_CNewSparse3DMat 
 *                      Mlib_CFreeSparse3DMat
 * Mlib_DTypeMat	Mlib_CTypeMat
 * Mlib_DTypeVec	Mlib_CTypeVec
 * Mlib_DFprintMat	Mlib_CFprintMat
 * Mlib_DFprintVec	Mlib_CFprintVec
 * Mlib_DSubVec		Mlib_CSubVec
 * Mlib_DMagVec		Mlib_CMagVec
 * Mlib_DCopyVec	Mlib_CCopyVec
 * Mlib_DCopyPartVe
 * 					    Mlib_CCopyMat
 *		 			    Mlib_CAddMat
 * 					    Mlib_CSubMat
 * 					    Mlib_CMaxVec
 *		 			    Mlib_CInfNormVec
 * 					    Mlib_ConjTransp
 * 					    Mlib_CScalMat
 *		 			    Mlib_CMultMat
 * 					    Mlib_CMultMatV
 * 					    Mlib_CMultVec
 *                                          Mlib_CMultVecRec 
 *		  		  	    Mlib_CScalVec
 * 					    Mlib_CInvV
 * 					    Mlib_CExpGL
 * 					    Mlib_CExp
 * 					    Mlib_CElementFlip
 * 					    Mlib_CReduce
 * 					    Mlib_CLUDcomp
 * 					    Mlib_CLUSolve
 *					    Mlib_CType3DMat
 *					    Mlib_CFprint3DMat
 *
 *
 *      Char
 *     -------
 * Mlib_CharNewMat
 * Mlib_CharFreeMat
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "standard.h"
#include "mltypes.h"
#include "dcx.h"
#include "ml.h"

#define MAXIT 10000
#define EULER 0.57721566490153286061
#define FPMIN 1.0e-30
#define MACHEPS 1.0e-15

/*
 *	Declaration of the HEAD and END of the double type, temporary
 *	matrix list.
 */
extern delement 	*DTempHead,
  *DTempEnd;

extern celement 	*CTempHead,
  *CTempEnd;

extern ielement 	*ITempHead,
  *ITempEnd;

/*  static dcx_t	ONE = {1.0,0.0}; */



#define CHECK_TEMP_DOUBLE_P1(n) \
         if((n) > mlTemp_Double_size1) ml_ResizeTempDouble_P1(n)
#define CHECK_TEMP_DOUBLE_P(n) \
         if((n) > mlTemp_Double_size) ml_ResizeTempDouble_P(n)
#define CHECK_TEMP_INT_P(n) \
         if(n > mlTemp_Int_size) ml_ResizeTempInt_P(n)

static  double *mlTemp_Double_P = NULL;
static  int mlTemp_Double_size = 0;
static  double *mlTemp_Double_P1 = NULL;
static  int mlTemp_Double_size1 = 0;
static  int *mlTemp_Int_P = NULL;
static  int mlTemp_Int_size = 0;

void ml_ResizeTempDouble_P(int size);
void ml_ResizeTempDouble_P1(int size);
void ml_ResizeTempInt_P(int size);



/*
 *				Integer matrix and vector functions.
 */


/*
 *  Mlib_INewVec--
 *  allocates an integer vector
 *
 */
int *Mlib_INewVec(int iSize)
{
  return (int *) malloc(sizeof(int) * iSize);
}
/*
 *  Mlib_IFreeVec--
 *  frees an integer vector with range [nl..nh]
 *
 */
void Mlib_IFreeVec(integerv_t integerv)
{
  free((VOIDPTR) integerv);
  return;
}

/*
 *  Mlib_INewMat--
 *  allocates a integer matrix.
 *
 */
integerm_t Mlib_INewMat(int row, int col)
{
  int i;
  integerv_t integerv;
  integerm_t integerm;

  integerm = (integerm_t) malloc(sizeof(integerv_t) * row);
  integerv = (integerv_t) malloc(sizeof(int) * row * col);

  for (i = 0; i < row; i++) {
    integerm[i] = integerv + i * col;
  }

  return integerm;
}

/* added by ece603 aom team 10/10/97
*  Mlib_INewArray
*  allocates an int aray
*/
integerm_t Mlib_INewArray(int row, int *entries_per_row)
{
  int i;
  int no_entries=0;
  int offset=0;
  integerv_t integerv;
  integerm_t integerm;
   
  for(i=0;i<row;i++)
    no_entries+= entries_per_row[i];  /* calculates total # of entries */
   
  integerm = (integerm_t) malloc(sizeof(integerv_t) * row);
  integerv = (integerv_t) malloc(sizeof(int) * no_entries);
   
  for(i=0;i<row;i++) {
    integerm[i]=integerv + offset;
    offset += entries_per_row[i];
  }
  return integerm;
}

/* added by ece603 aom team 10/10/97
 *  Mlib_IFreeArray--
 *  frees an int array m 
 *
 */
void Mlib_IFreeArray(integerm_t integerm)
{
  free((VOIDPTR) *integerm);
  free((VOIDPTR) integerm);
  return;
}

/*
 *  Mlib_IFreeMat--
 *  frees an int matrix m 
 *
 */
void Mlib_IFreeMat(integerm_t integerm)
{
  free((VOIDPTR) *integerm);
  free((VOIDPTR) integerm);
  return;
}


/*
 *  Mlib_ITypeVec--
 *	Type a vector (or portion of a vector) to the screen
 */
void Mlib_ITypeVec(a_P,nrl,nrh,exp_float)                             
                                                                           
     int	           *a_P;           /* the real vector to print out         */ 
     int            nrl,nrh;        /* the lower and upper row dimensions   */ 
     int            exp_float; /* 1: print it out in exponential else float */ 
                                                                           
{
  int         i,endofline;
                                               
  printf ("\n");                             
  if (exp_float){
    endofline = 0;
    for (i=nrl; i<=nrh; i++){
      printf ("(% d) ",a_P[i]);       
      endofline++;
      if(endofline == 6){
	printf("\n");
	endofline = 0;
      }
    }
    printf ("\n");                             
  }
  else{                                       
    endofline = 0;
    for (i=nrl; i<=nrh; i++){
      printf ("(% d) ",a_P[i]);       
      endofline++;
      if(endofline == 6){
	printf("\n");
	endofline = 0;
      }
    }
    printf ("\n");                             
  }
  printf ("\n");                             
                                               
  return;                                       
}                                              


/* ############################################################## */
/*
 *				Complex matrix and vector functions.
 *
 *  Purpose:  These routines perform the following functions:
 *		-allocates memory for dcx_t matrices and vectors
 *              -inverts a complex nxn matrix
 *  		-adds, subtracts, multiplies, raises complex matrices to
 *		  the e power, multiplies a matrix by a vector
 */


/*
 *  Mlib_CNewMat--
 *
 *  allocates a dcx_t matrix 
 */

dcxm_t Mlib_CNewMat(int rows, int cols)
{
  int i, size;
  dcxv_t dcxv;
  dcxm_t dcxm;

  dcxm = (dcxm_t) malloc(sizeof(dcxv_t) * rows);
  size = sizeof(dcx_t) * rows * cols;
  dcxv = (dcxv_t) malloc(size);

  for (i = 0; i < rows; i++) {
    dcxm[i] = dcxv + i * cols;
  }

  return dcxm;
}

/*
 *  Mlib_CNew3DMat--
 *  allocates a dcx_t 3D matrix 
 *
 */

dcxmv_t Mlib_CNew3DMat(int iSize, int jSize, int kSize)
{
  int i;
  dcxv_t dcxv;
  dcxm_t dcxm;
  dcxmv_t dcxmv;

  dcxmv = (dcxmv_t) malloc(sizeof(dcxm_t) * iSize);
  dcxm = (dcxm_t) malloc(sizeof(dcxv_t) * iSize * jSize);
  dcxv = (dcxv_t) malloc(sizeof(dcx_t) * iSize * jSize * kSize);

  for (i = 0; i < iSize; i++) {
    dcxmv[i] = dcxm + i * jSize;
  }

  for (i = 0; i < iSize * jSize; i++) {
    dcxm[i] = dcxv + i * kSize;
  }

  return dcxmv;
}


/*
 *  Mlib_CFree3DMat--
 *  deallocates a dcx_t 3D matrix 
 *
 */

void Mlib_CFree3DMat(dcxmv_t dcxmv)
{
  free((VOIDPTR) **dcxmv);
  free((VOIDPTR) *dcxmv);
  free((VOIDPTR) dcxmv);
  return;
}

/*
 * Mlib_CNewSparse3DMat
 *
 * Allocates a matrix of pointers to dcx_t with range [0..row-1][0..col-1]
 * The matrix is intiallized with NULL pointers.
 *
 */
dcxsmv_t Mlib_CNewSparse3DMat(int row, int col)
{
  int i;
  dcxm_t dcxm;
  dcxsmv_t dcxsmv;

  dcxsmv = (dcxsmv_t) malloc(sizeof(dcxm_t) * row);
  dcxm = (dcxm_t) malloc(sizeof(dcxv_t) * row * col);

  for (i = 0; i < row; i++)
    {
      dcxsmv[i] = dcxm + i * col;
      /* for(j=0; j < col; j++) dcxsmv[j] = NULL; */
    }
  return dcxsmv;
}

/*
 *  Mlib_CFreeSparse3DMat--
 *
 *  frees a dcx_t matrix of pointers
 *
 */
void Mlib_CFreeSparse3DMat(dcxsmv_t dcxsmv, int row, int col)
{
  int i, j;
  for (i = 0; i < row; i++)
    for (j = 0; j < col; j++)
      if((dcxsmv)[i][j] != NULL) free((VOIDPTR) dcxsmv[i][j]);

  free((VOIDPTR) *dcxsmv);
  
  free((VOIDPTR) dcxsmv);
  return;
}


/*
 *  Mlib_CFreeMat--
 *  frees a dx_t matrix 
 *
 */
void Mlib_CFreeMat(dcxm_t dcxm) 
{ 
  free((VOIDPTR) *dcxm); 
  free((VOIDPTR) dcxm); 
  return; 
}

/*
 *  Mlib_CNewVec--
 *  allocates memory for a dcx_t vector 
 *
 */
dcxv_t Mlib_CNewVec(int iSize)
{
  return (dcxv_t) malloc(sizeof(dcx_t) * iSize);
}

/*
 *  Mlib_CFreeVec--
 *  frees the dcx_t vector v_P with range [nl..nh]
 *
 */
void Mlib_CFreeVec(dcxv_t dcxv)
{
  free((VOIDPTR) dcxv);
  return;
}

/*
 *  Mlib_CCopyMat--
 *  copies one complex matrix into another
 *
 */

void Mlib_CCopyMat(c_PP, a_PP, rows, cols)
     dcx_t	**c_PP,		/* result matrix  */
  **a_PP;		/* input matrix	*/
     int	rows,cols;
{
  int 	i,j;

  for (i=0; i<rows; i++) 
    for (j=0; j<cols; j++) 
      c_PP[i][j] = a_PP[i][j];
}

/*
 *  Mlib_CCopyVec--
 *  copies one vector into another
 *
 */

void Mlib_CCopyVec(c_P, a_P, n)
     dcx_t	*c_P,		/* result vector*/
  *a_P;		/* input vector	*/
     int	n;
{
  int 	i;

  for (i=0; i<n; i++) 
    c_P[i] = a_P[i];
}


/*
 *  Mlib_CMult
 *  multiplies two complex numbers
 *
 */


dcx_t Mlib_CMult(dcx_t a, dcx_t b)
{
  dcx_t result;

  result.re = a.re*b.re - a.im*b.im;
  result.im = a.re*b.im + a.im*b.re;
  return(result);

}




/*
 *  Mlib_CAddMat--
 *  adds 2 nxn complex matrices
 *
 */

void Mlib_CAddMat(c_PP, a_PP, b_PP, n)
     dcx_t		**c_PP,		/* result matrix*/
  **a_PP,**b_PP;	/* input complex matrices		*/
     int		n;		
{
  int		i,j;

  for (i=0; i<n; i++) 
    for (j=0; j<n; j++) {
      c_PP[i][j].re = a_PP[i][j].re + b_PP[i][j].re;
      c_PP[i][j].im = a_PP[i][j].im + b_PP[i][j].im;
    }

  return;
}


/*
 *  Mlib_CSubMat--
 *  subtracts 2 nxn complex matrices
 *
 */

void Mlib_CSubMat(c_PP, a_PP, b_PP, n)

     dcx_t		**c_PP,		/* result matrix	*/
  **a_PP,**b_PP;	/* input complex matrices		*/
     int		n;		
{
  int		i,j;

  for (i=0; i<n; i++) 
    for (j=0; j<n; j++) {
      c_PP[i][j].re = a_PP[i][j].re - b_PP[i][j].re;
      c_PP[i][j].im = a_PP[i][j].im - b_PP[i][j].im;
    }

  return;
}

/*
 *  Mlib_CSubVec--
 *  subtracts 2 n complex vectors
 *
 */

void Mlib_CSubVec(c_P, a_P, b_P, n)

     dcx_t		*c_P,		/* result vector */
  *a_P,*b_P;	/* input complex vectors */
     int		n;	
{

  int		i;

  for (i=0; i<n; i++){
    c_P[i].re = a_P[i].re - b_P[i].re;
    c_P[i].im = a_P[i].im - b_P[i].im;
  }

  return;
}

/*
 * Mlib_CScalVec--
 *
 * Multiply a complex vector by a complex number.
 */
void Mlib_CScalVec(dcx_t *a, dcx_t scale, int n)
{
  int i;
  for(i=0; i<n; i++)
    a[i] = Cx_DMult(a[i],scale);

  return;
}

/*
 *  Mlib_CMagVec--
 *  finds the 2 norm of a complex vector
 *
 */
double Mlib_CMagVec(a_P, n)

     dcx_t		*a_P;	/* input complex vectors */
     int		n;	
{
  int		i;
  double	magsqrd,mag,temp;

  magsqrd = 0.0;
  mag     = 0.0;

  for (i=0; i<n; i++){
    temp = Cx_DAbs(a_P[i]);
    magsqrd = magsqrd + temp*temp;
  }
  mag = sqrt(magsqrd);
  return(mag);
}

/*
 *  Mlib_CMaxVec--
 *  returns the absolute value of the largest element in the vector.
 *
 */
double Mlib_CMaxVec(a_P, n)
     dcx_t		*a_P;	
     int		n;	
{
  int		i;
  double	mmax,temp;
  mmax = 0.0;
  for (i=0; i<n; i++){
    temp = Cx_DAbs(a_P[i]);
    if(temp > mmax)
      mmax = temp;
  }
  return(mmax);
}

/*
 * Mlib_CInfNormVec--
 *
 *	Finds the infinity norm of a complex vector.  Returns a complex
 *	number (the member of the vector which has the largest
 *	magnitude).
 */
dcx_t	Mlib_CInfNormVec(x_P,n)
     dcx_t	*x_P;
     int	n;
{
  int	i,
    place;
  double	temp;
  temp = 0.0;
  for(i=0; i<n; i++){
    if(DABS(x_P[i]) > temp){
      temp = DABS(x_P[i]);
      place = i;
    }
  }
  return(x_P[place]);
}

/*
 *  Mlib_ConjTransp--
 *  finds the transpose of a nl x nh complex matrix
 *
 */
void Mlib_ConjTransp(c_PP, a_PP, n)

     dcx_t	**c_PP,			/* result matrix	*/
  **a_PP; 		/* input complex matrices		*/
     int	n;
{

  int		i,j;

  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      c_PP[i][j].re = a_PP[j][i].re;
      c_PP[i][j].im = -a_PP[j][i].im;
    }

  return;
}
/*
 *  Mlib_CScalMat--
 *  Multiplies a matrix by a scalar constant
 *
 *  May overwrite input with output.
 */
void Mlib_CScalMat(c_PP, a_PP, c, n)

     dcx_t	**c_PP,			/* result matrix	*/
  **a_PP, 		/* input complex matrices		*/
  c;			/* scalar constant */
     int	n;
{

  int		i,j;

  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      c_PP[i][j] = Cx_DMult(a_PP[i][j],c);
    }

  return;
}


/*
 *  Mlib_CMultMat--
 *  multiplies 2 nxn complex matrices
 *
 *  Modified: 1/15/91
 *  The first matrix, a_PP can be overwritten but not b_PP.
 *  (i.e. c_PP and a_PP can point to the same matrix)
 */
void Mlib_CMultMat(c_PP, a_PP, b_PP, n)

     dcx_t	**c_PP,			/* result matrix	*/
  **a_PP,**b_PP;		/* input complex matrices		*/
     int	n;
{
  dcx_t	temp;
  int		i,j,k;

  if( a_PP == c_PP ||  b_PP == c_PP){
    printf("ERROR in Mlib_CMultMat.... matrix overwrite problem\n");
    exit(0);
  }
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) 
      {
	c_PP[i][j].re = 0.0;
	c_PP[i][j].im = 0.0;
      }

  for (k=0; k<n; k++) 
    for (i=0; i<n; i++)
      for (j=0; j<n; j++) {

	/*
	 * temp = MULT(a_PP[k][j],b_PP[j][i]);
	 */
	temp.re = a_PP[k][j].re * b_PP[j][i].re - 
	  a_PP[k][j].im * b_PP[j][i].im;

	temp.im = a_PP[k][j].re * b_PP[j][i].im +
	  a_PP[k][j].im * b_PP[j][i].re;

	/*
	 * c_PP[k][i] = ADD(c_PP[k][i],temp);
	 */
	c_PP[k][i].re += temp.re;
	c_PP[k][i].im += temp.im;
      }
  return;
}


/*
 *  Mlib_CMultMatV--
 *  multiplies 2 nxn vectorized complex matrices
 *
 */
void Mlib_CMultMatV(c_PPP, a_PPP, b_PPP, n, lengthOfVector)

     dcx_t	***c_PPP,			/* result matrix	*/
  ***a_PPP,***b_PPP;		/* input complex matrices		*/
     int	n, lengthOfVector;
{
  dcx_t	temp;
  int		h,i,j,k;

  for (h=0; h< lengthOfVector; h++)
    {

      for (i=0; i<n; i++)
	for (j=0; j<n; j++) 
	  {
	    c_PPP[i][j][h].re = 0.0;
            c_PPP[i][j][h].im = 0.0;
	  }

      for (k=0; k<n; k++) 
	for (i=0; i<n; i++)
	  for (j=0; j<n; j++) {

	    /*
	     * temp = MULT(a_PPP[k][j][h],b_PPP[j][i][h]);
	     */
	    temp.re = a_PPP[k][j][h].re * b_PPP[j][i][h].re - 
	      a_PPP[k][j][h].im * b_PPP[j][i][h].im;

	    temp.im = a_PPP[k][j][h].re * b_PPP[j][i][h].im +
	      a_PPP[k][j][h].im * b_PPP[j][i][h].re;

	    /*
	     * c_PPP[k][i][h] = ADD(c_PPP[k][i][h],temp);
	     */
	    c_PPP[k][i][h].re += temp.re;
	    c_PPP[k][i][h].im += temp.im;
	  }
    }
  return;
}


/*
 *  Mlib_CMultVec--
 *  multiplies a nxn complex matrix a_PP with a nx1 complex vector b_P
 *
 */

void Mlib_CMultVec(c_P, a_PP, b_P, n)

     dcx_t		*c_P,		/* result matrix (nx1)	*/
  **a_PP,		/* input complex matrix			*/
  *b_P;		/* input complex vector			*/
     int		n;
{

  dcx_t	temp;
  int 	i,j;

  for (i=0; i<n; i++) { 
    c_P[i].re = 0.0;
    c_P[i].im = 0.0;
  }

  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {    	/* I'm doing it my wayyy!!!! */

      /*
       * temp = MULT(a_PP[i][j],b_P[j]);
       */

      temp.re = a_PP[i][j].re * b_P[j].re - a_PP[i][j].im * b_P[j].im;
      temp.im = a_PP[i][j].re * b_P[j].im + a_PP[i][j].im * b_P[j].re;

      /*
       *  c_P[i] = ADD(c_P[i],temp);
       */

      c_P[i].re += temp.re;
      c_P[i].im += temp.im;

    }

  return;
}


/*
 *  Mlib_CMultVecRec--
 *  multiplies a nxm complex matrix a_PP with a mx1 complex vector b_P
 *
 */

void Mlib_CMultVecRec(c_P, a_PP, b_P, n, m)

     dcx_t		*c_P,		/* result matrix (nx1)	*/
  **a_PP,		/* input complex matrix			*/
  *b_P;		/* input complex vector			*/
     int		n, m;
{

  dcx_t	temp;
  int 	i,j;

  for (i=0; i<n; i++) { 
    c_P[i].re = 0.0;
    c_P[i].im = 0.0;
  }

  for (i=0; i<n; i++)
    for (j=0; j<m; j++) {    	/* I'm doing it my wayyy!!!! */

      /*
       * temp = MULT(a_PP[i][j],b_P[j]);
       */

      temp.re = a_PP[i][j].re * b_P[j].re - a_PP[i][j].im * b_P[j].im;
      temp.im = a_PP[i][j].re * b_P[j].im + a_PP[i][j].im * b_P[j].re;

      /*
       *  c_P[i] = ADD(c_P[i],temp);
       */

      c_P[i].re += temp.re;
      c_P[i].im += temp.im;

    }

  return;
}


/*
 *  Mlib_CMultVecRecTrans--
 *  multiplies the transpose of a nxm complex matrix a_PP 
 *  with a nx1 complex vector b_P
 *
 */

void Mlib_CMultVecRecTrans(c_P, a_PP, b_P, n, m)

     dcx_t		*c_P,		/* result matrix (nx1)	*/
  **a_PP,		/* input complex matrix			*/
  *b_P;		/* input complex vector			*/
     int		n, m;
{

  dcx_t	temp;
  int 	i,j;

  for (i=0; i<m; i++) { 
    c_P[i].re = 0.0;
    c_P[i].im = 0.0;
  }

  for (i=0; i<m; i++)
    for (j=0; j<n; j++) {    	/* I'm doing it my wayyy!!!! */

      /*
       * temp = MULT(a_PP[i][j],b_P[j]);
       */

      temp.re = a_PP[j][i].re * b_P[j].re - a_PP[j][i].im * b_P[j].im;
      temp.im = a_PP[j][i].re * b_P[j].im + a_PP[j][i].im * b_P[j].re;

      /*
       *  c_P[i] = ADD(c_P[i],temp);
       */

      c_P[i].re += temp.re;
      c_P[i].im += temp.im;

    }

  return;
}


/*
 *  Mlib_CInv--
 *  inverts a complex matrix
 *  See "HP-15C Owner's Handbook", pg. 160 for details
 *
 *  Input = Output OK   (output overwrites input)
 */

void Mlib_CInv(c_PP ,a_PP, n)

     dcx_t	**c_PP,		/* result matrix	*/
  **a_PP;		/* the ptr. to the ptr. of the input complex    */
     /* structure/array to invert			*/
     int	n;
{

  int		i,j;
  double	**local1_PP,
    **local2_PP;

  local1_PP = Mlib_DNewMat(2*n,2*n);
  local2_PP = Mlib_DNewMat(2*n,2*n);
  /*
   * form the local1_PP matrix from the a_PP matrix in the following way:
   *
   *    local1_PP     = | X  -Y |
   *       	          | Y   X |
   *
   * where X = real submatrix of a_PP and
   *       Y = imaginary submatrix of a_PP
   *
   */

  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      local1_PP[i][j] = a_PP[i][j].re;

  for (i=n; i<2*n; i++)
    for (j=n; j<2*n; j++)
      local1_PP[i][j] = a_PP[i-n][j-n].re;

  for (i=0; i<n; i++)
    for (j=n; j<2*n; j++)
      local1_PP[i][j] = -a_PP[i][j-n].im;

  for (i=n; i<2*n; i++)
    for (j=0; j<n; j++)
      local1_PP[i][j] = a_PP[i-n][j].im;

  for (i=0; i<2*n; i++){
    if (local1_PP[i][i] == 0.0)
      local1_PP[i][i] = TINY;
  }

  Mlib_DInv(local2_PP,local1_PP,2*n);

  /*
   * put the inverted elements into local2_PP
   */
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      c_PP[i][j].re = local2_PP[i][j];

  for (i=n; i<2*n; i++)
    for (j=0; j<n; j++)
      c_PP[i-n][j].im = local2_PP[i][j];

  Mlib_DFreeMat(local1_PP);
  Mlib_DFreeMat(local2_PP);
  return;
}


/*
 *  Mlib_CInvV--
 *  inverts a vectorized complex matrix
 *  (a 2D matrix of vectors)
 *  See "HP-15C Owner's Handbook", pg. 160 for details
 *
 *  Input = Output OK
 */

void Mlib_CInvV(c_PPP ,a_PPP, n, lengthOfVector)

     dcx_t	***c_PPP,	/* result matrix	*/
  ***a_PPP;	/* the ptr. to the ptr. of the input complex    */
     /* structure/array to invert			*/
     int	n, lengthOfVector;
{

  int		i,j,k;
  double	**local1_PP,
    **local2_PP;

  /* local1_PP = DTempAssignMat(2*n);
   * local2_PP = DTempAssignMat(2*n);
   */

  local1_PP = Mlib_DNewMat(2*n,2*n);
  local2_PP = Mlib_DNewMat(2*n,2*n);
  /*
   * form the local1_PP matrix from the a_PP matrix in the following way:
   *
   *    local1_PP     = | X  -Y |
   *       	          | Y   X |
   *
   * where X = real submatrix of a_PP and
   *       Y = imaginary submatrix of a_PP
   *
   */

  for (k=0; k<lengthOfVector; k++)
    {
      for (i=0; i<n; i++)
        for (j=0; j<n; j++)
	  local1_PP[i][j] = a_PPP[i][j][k].re;

      for (i=n; i<2*n; i++)
        for (j=n; j<2*n; j++)
	  local1_PP[i][j] = a_PPP[i-n][j-n][k].re;

      for (i=0; i<n; i++)
        for (j=n; j<2*n; j++)
	  local1_PP[i][j] = -a_PPP[i][j-n][k].im;

      for (i=n; i<2*n; i++)
	for (j=0; j<n; j++)
	  local1_PP[i][j] = a_PPP[i-n][j][k].im;

      for (i=0; i<2*n; i++)
        {
	  if (local1_PP[i][i] == 0.0)
	    local1_PP[i][i] = TINY;
	}

      Mlib_DInv(local2_PP,local1_PP,2*n);

      /*
       * put the inverted elements into local2_PP
       */
      for (i=0; i<n; i++)
        for (j=0; j<n; j++)
	  c_PPP[i][j][k].re = local2_PP[i][j];

      for (i=n; i<2*n; i++)
	for (j=0; j<n; j++)
	  c_PPP[i-n][j][k].im = local2_PP[i][j];
    }
  /* DTempReleaseMat0(DTempHead,local1_PP);
   * DTempReleaseMat0(DTempHead,local2_PP);
   */
  Mlib_DFreeMat(local1_PP);
  Mlib_DFreeMat(local2_PP);
  return;
}


/*
 *  Mlib_CExpGL--
 *  calculate e raised to a nxn complex matrix
 *
 *  Parameters:
 *	length		the coeff. to mult. by the matrix
 *	eval_P		the eigenvalues of the matrix
 *	evec_PP		the eigenvectors of the matrix
 *	evecinv_PP		the eigenvector inverse matrix
 *	nl,nh			low and high dim. of the matrices
 *
 *  Results:
 *	c_PP			matrix containing e raised to the matrix
 *				 described by the eigensystem
 *
 */

void Mlib_CExpGL(c_PP,L,eval_P,evec_PP,evecinv_PP,n)

     dcx_t	**c_PP,		/* result matrix	*/
  *eval_P,
  **evec_PP,
  **evecinv_PP;
     double	L;
     int	n;
{

  dcx_t	temp, temp2,
    **local1_PP;
  int		i,j;
  local1_PP = Mlib_CNewMat(n,n);

  for (i=0; i<n; i++)
    for (j=0; j<n; j++) 
      c_PP[i][j].re = c_PP[i][j].im = 0.0;

  for (i=0; i<n; i++) {
    temp.re = (L);
    temp.im = 0.0;
    /*
     * (L)*eval_P[i]
     */
    temp2.re = temp.re * eval_P[i].re - temp.im * eval_P[i].im;
    temp2.im = temp.re * eval_P[i].im + temp.im * eval_P[i].re;
    c_PP[i][i] = Cx_DExp(temp2);
  }

  Mlib_CMultMat(local1_PP,evec_PP,c_PP,n);
  Mlib_CMultMat(c_PP,local1_PP,evecinv_PP,n);

  Mlib_CFreeMat(local1_PP);
  return;
}


/*
 *  Cmat_ExpNP--
 *  calculate e raised to a nxn complex matrix
 *
 *  Parameters:
 *	n_P->length		the coeff. to mult. by the matrix
 *	n_P->eval_P		the eigenvalues of the matrix
 *	n_P->evec_PP		the eigenvectors of the matrix
 *	n_P->evecinv_PP		the eigenvector inverse matrix
 *	nl,nh			low and high dim. of the matrices
 *
 *  Results:
 *	c_PP			matrix containing e raised to the matrix
 *				 described by the eigensystem
 *
 */

/*
 * The function Cmat_ExpNP() has been deleted here but may be found
 * in the file CMAT_MATRIX.C which contains the original routines
 */


void Mlib_CExp(c_PP,coeff,eval_P,evec_PP,evecinv_PP,n)
                                                                             
     dcx_t	**c_PP;     	/* result matrix   */   
     float	*coeff;     	/* the coeff. to mult. by the matrix    */   
     dcx_t	*eval_P,		/* the eigenvalues of the matrix  */   
       **evec_PP,      /* the eigenvectors of the matrix       */   
       **evecinv_PP;   /* the eigenvector inverse matrix       */   
     int             n;		/* low and high dim. of matrices        */   
{
  dcx_t	temp, temp2,                                         
    **local1_PP;
  int		i,j;                                                           
  /*
   *	Assign one of the temporary matrices from the linked list.
   */
  local1_PP = Mlib_CNewMat(n,n);
  /*
   *	Zero out the result matrix.
   */
  for (i=0; i<n; i++) 
    for (j=0; j<n; j++) 
      c_PP[i][j].re = c_PP[i][j].im = 0.0;

  for (i=0; i<n; i++) {                                                 
    temp.re = *coeff;                                                    
    temp.im = 0.0;                                                       
    /*                                                                   
     * MULT(temp,eval_P[i]);                                          
     */                                                               
    temp2.re = temp.re * eval_P[i].re - temp.im * eval_P[i].im;       
    temp2.im = temp.re * eval_P[i].im + temp.im * eval_P[i].re;       
    c_PP[i][i] = Cx_DExp(temp2);                                          
  }                                                                     
  Mlib_CMultMat(local1_PP,evec_PP,c_PP,n);                               
  Mlib_CMultMat(c_PP,local1_PP,evecinv_PP,n);                            
  Mlib_CFreeMat(local1_PP);
  return;                                                                  
}


/*
 *	Mlib_CTypeMat--
 *
 *	Types a portion of the matrix to the screen.  This routine
 *	does not assume that the matrix starts at [0][0] so that a
 *	portion of the matrix instead of the whole thing may be 
 *	typed out.
 */
void Mlib_CTypeMat(a_PP,nrl,nrh,ncl,nch,exp_float)

     dcx_t    **a_PP;  /* the complex matrix to print out      */
     int             nrl,nrh, /* the lower and upper row dimensions   */
  ncl,nch;	/* the lower and upper column dimensions*/
     int             exp_float; /* 1: print it out in exponential else float */
{
  int         i,j;
  printf ("\n");
  if (exp_float)
    for (i=nrl; i<=nrh; i++)
      for (j=ncl; j<=nch; j++) {
	printf ("(% 5.2e,% 5.2e) ",a_PP[i][j].re,a_PP[i][j].im);
	if (j==nch)
	  printf ("\n");
      }
  else
    for (i=nrl; i<=nrh; i++)
      for (j=ncl; j<=nch; j++) {
	printf ("(% e, % e) ",a_PP[i][j].re,a_PP[i][j].im);
	if (j==nch)
	  printf ("\n");
      }
  printf ("\n");

  return;
}


/*
 *	Mlib_CType3DMat--
 *
 *  Types a slice of a 3D matrix to the screen.
 */
void Mlib_CType3DMat(a_PPP,nrl,nrh,ncl,nch,k,exp_float)

     dcx_t    ***a_PPP;  /* the complex matrix to print out      */
     int             nrl,nrh, /* the lower and upper row dimensions   */
  ncl,nch,	/* the lower and upper column dimensions*/
  k,
  exp_float; /* 1: print it out in exponential else float */
{
  int         i,j;
  printf ("\n");
  if (exp_float)
    for (i=nrl; i<=nrh; i++)
      for (j=ncl; j<=nch; j++) {
	printf ("(% 5.2e,% 5.2e) ",a_PPP[i][j][k].re,a_PPP[i][j][k].im);
	if (j==nch)
	  printf ("\n");
      }
  else
    for (i=nrl; i<=nrh; i++)
      for (j=ncl; j<=nch; j++) {
	printf ("(% e, % e) ",a_PPP[i][j][k].re,a_PPP[i][j][k].im);
	if (j==nch)
	  printf ("\n");
      }
  printf ("\n");

  return;
}

/*
 *	Mlib_CFprint3DMat--
 *
 *  Prints a slice of a 3D matrix to a file.
 */
void Mlib_CFprint3DMat(out_F,a_PPP,nrl,nrh,ncl,nch,k,exp_float)

     FILE	*out_F;
     dcx_t    ***a_PPP;  /* the complex matrix to print out      */
     int             nrl,nrh, /* the lower and upper row dimensions   */
       ncl,nch,	/* the lower and upper column dimensions*/
       k,
       exp_float; /* 1: print it out in exponential else float */
{
  int         i,j,m,ii;
  fprintf (out_F,"\n");
  switch(exp_float){
  case 0:
    for (i=nrl; i<=nrh; i++){
      for (j=ncl; j<=nch; j++) {
	fprintf (out_F,"(% e, % e) ",
		 a_PPP[i][j][k].re,a_PPP[i][j][k].im);
	if (j==nch)
	  fprintf (out_F,"\n");
      }
      fprintf (out_F,"\n");
    }
    break;
  case 1:
    for (i=nrl; i<=nrh; i++){
      for (j=ncl; j<=nch; j++){
	fprintf (out_F,"(% .2e % .2e)",
		 a_PPP[i][j][k].re,a_PPP[i][j][k].im);
	if (j==nch)
	  fprintf (out_F,"\n");
      }
    }
    break;
  case 2:
    for (i=nrl; i<=nrh; i=i+3){
      fprintf(out_F,"\tCOLUMN %d \t\t\t\tCOLUMN %d \t\t\t\t COLUMN %d\n",
	      i,i+1,i+2);

      m = i+2;
      if(m >= nrh)
	m = nrh;
      for (j=ncl; j<=nch; j++) {
	for (ii=i; ii<=m; ii++){
	  fprintf (out_F,"(% .3e,% .3e) ",
		   a_PPP[j][ii][k].re,a_PPP[j][ii][k].im);
	}
	fprintf(out_F,"\n");
      }
    }
    break;
  }

  return;
}



void Mlib_CFprintMat(out_fp,a_PP,nrl,nrh,ncl,nch,exp_float)
     dcx_t    **a_PP;  /* the complex matrix to print out      */
     int		nrl,nrh, /* the lower and upper row dimensions   */
  ncl,nch;	/* the lower and upper column dimensions*/
     int		exp_float; 
     /*		 0: print it out in exponential format, matrix form
 *		 1: print out in exponential format, pretty matrix form
 *		 2: print out in long exponential format, column form
 */

     FILE	*out_fp;

{
  int         i,j,k,m;

  fprintf (out_fp,"\n");
  switch ( exp_float){
  case 0:
    for (i=nrl; i<=nrh; i++)
      for (j=ncl; j<=nch; j++) {
	fprintf (out_fp,"% .3e % .3e  ",a_PP[i][j].re,a_PP[i][j].im);
	if (j==nch)
	  fprintf (out_fp,"\n");
      }
    break;
  case 1:
    for (i=nrl; i<=nrh; i=i+3){
      fprintf(out_fp,"\tCOLUMN %d \t\tCOLUMN %d \t\t COLUMN %d\n",
	      i,i+1,i+2);

      m = i+2;
      if(m >= nrh)
	m = nrh;
      for (j=ncl; j<=nch; j++) {
	for (k=i; k<=m; k++){
	  fprintf (out_fp,"(% .3e,% .3e) ",
		   a_PP[j][k].re,a_PP[j][k].im);
	}
	fprintf(out_fp,"\n");
      }
    }
    break;

  case 2:
    for (i=nrl; i<=nrh; i++){
      for (j=ncl; j<=nch; j++) {
	fprintf (out_fp,"% .20e  % .20e \n",a_PP[i][j].re,a_PP[i][j].im);
      }
    }
    break;
  }

  fprintf (out_fp,"\n");

  return;
}


void Mlib_CTypeVec(a_P,nrl,nrh,exp_float)

     dcx_t    *a_P;     /* the complex vector to print out      */
     int             nrl,nrh;  /* the lower and upper row dimensions   */
     int             exp_float; /* 1: print it out in exponential else float */
{

  int         i;
  printf ("\n");
  if (exp_float){
    for (i=nrl; i<=nrh; i++){
      printf ("(% 5.2e,% 5.2e)\n",a_P[i].re,a_P[i].im);
    }
  }
  else{
    for (i=nrl; i<=nrh; i++){
      printf ("(% e, % e)\n",a_P[i].re,a_P[i].im);
    }
  }
  printf ("\n");

  return;
}


void Mlib_CFprintVec(out_fp,a_P,nrl,nrh,exp_float)

     dcx_t    *a_P;     /* the complex vector to print out      */
     int             nrl,nrh;  /* the lower and upper row dimensions   */
     int             exp_float; /* 1: print it out in exponential else float */
     FILE	*out_fp;

{
  int         i;

  fprintf (out_fp,"\n");
  switch (exp_float){
  case 0:
    for (i=nrl; i<=nrh; i++)
      fprintf (out_fp,"     % 8.5e  % 8.5e \n",
	       a_P[i].re,a_P[i].im);
    break;

  case 1:
    for (i=nrl; i<=nrh; i++)
      fprintf (out_fp,"     (% .5e,% .5e)\n",
	       a_P[i].re,a_P[i].im);
    break;
  case 2:
    for (i=nrl; i<=nrh; i++)
      fprintf (out_fp,"     % .20e  % .20e \n",
	       a_P[i].re,a_P[i].im);
    break;

  }

  fprintf (out_fp,"\n");

  return;

}



/* Mlib_CElementFlip--
 * flips the individual elements of a matrix.  This is NOT an inversion
 * uses: convert a impedance load [ZL] to an admittance load [YL]
 *
 */

void Mlib_CElementFlip(in_PP,out_PP,n)

     int	n;
     dcx_t	**in_PP,
  **out_PP;
{
  int i,j;
  static dcx_t     one = {1.0,0.0};

  for (i=0; i<n; i++)
    {
      for (j=0; j<n; j++)
	{
	  out_PP[i][j] = Cx_DDiv(one,in_PP[i][j]);  
	}
    }
}

/**********************************************************************
 * Mlib_CReduce--                                                     *
 *        Matrix reduction routine, 2D matrices.                      *
 *                                                                    *
 * Authors:  Mark S. Basel                                            *
 *           C. Chang                                                 *
 *                                                                    *
 *      Inputs:   array   Complex Matrix to be reduced                *
 *                start   Integer indicating the size of the starting *
 *                        matrix. (eg, for a 3x3, start = 3)          *
 *                reduce  Integer indicating the size of the reduced  *
 *                        matrix. (eg, a 3x3 reduced to a 1x1,        *
 *                                 reduce = 1)                        *
 *                                                                    *
 *      return status: 0      ok                                      *
 *                     1      pivot < 1.0e-25,  reduction failure     *
 **********************************************************************/
/* FCABS2 is a float macro for taking the absolute squared of a complex *
 * number                                                               */
#define FCABS2(a) ( (a).re * (a).re + (a).im * (a).im )
#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}

void Mlib_CReduce(start,array, reduce) 
     int start, reduce;
     dcx_t **array;
{
  int i, j,
    row[256], col[256];
  dcx_t pivot, temp,
    *tmp;
  pivot.re = 0.0;
  pivot.im = 0.0;
  /*
   *  Reduce "start" by 1 to match largest extent of [array], i.e. match
   *  matrix description to "C" description.
   */
  start--;
  for(i = reduce; i <=start; ++i)
    if(FCABS2(array[i][i]) > FCABS2(pivot)) {
      pivot = array[i][i];
      row[start] = i;
      col[start] = i;
    }
  if(FCABS2(pivot) < 1.0e-16)
    for(i = reduce; i <=start; ++i) for(j = reduce; j < i; ++j) 
      if(FCABS2(array[i][j]) > FCABS2(pivot)) {
	pivot = array[i][j];
	row[start] = i;
	col[start] = j;
      }
  if(FCABS2(pivot) < 1.0e-16)
    for(i = reduce; i <=start; ++i) for(j = reduce; j <=start; ++j)
      if(FCABS2(array[i][j]) > FCABS2(pivot)) {
	pivot = array[i][j];
	row[start] = i;
	col[start] = j;
      }
  if(FCABS2(pivot) < 1.0e-30) {
    return;
  }
  
  for(i=0; i<=start; ++i) {
    temp = array[i][start];
    array[i][start] = array[i][col[start]];
    array[i][col[start]] = temp;
  }
 
  tmp = array[start];   
  array[start] = array[row[start]];
  array[row[start]] = tmp;

  for(j=0; j<=start; ++j)
    array[start][j] = DDIV(array[start][j], pivot);   
  for(i=0; i<=start; ++i)
    if(FCABS2(array[i][start]) > 1.0e-30)   
      for(j=0; j<=start; ++j)
	array[i][j] = DSUB(array[i][j],
			   Cx_DMult(array[i][start], array[start][j]));    
  return;
}
#undef SWAP
#undef FCABS2

/*
 * Mlib_CLUdcomp--
 *	routine which performs an LU decomposition on the input matrix.
 *
 *	input:	a -- input matrix.
 *		n -- size of matrix.
 *	output: a -- output matrix, contains L and U.
 *		ipivot -- vector containing the pivot information.
 */
void Mlib_CLUdcomp(a,n,ipivot)
     dcx_t	**a;
     int	n,
  *ipivot;
{
  int	i,j,k,m,imax,jj;
  double	mmax;
  dcx_t	temp;

  /*
   *	Setup the initial pivot vector such that ipivot[0] = 0,
   * 	ipivot[1] = 1, etc.
   */
  for(jj=0; jj<n; jj++)
    ipivot[jj] = jj;

  for(k=0; k<n; k++){  /* Beginning of main loop*/
    /*
     *      Debug+++
     */
    int mmm;
    for(mmm=0; mmm<n; mmm++){
      if(a[mmm][mmm].re == 0 && a[mmm][mmm].im == 0){
	printf("\t\t[ml.c] Mlib_CLUdcomp: zero diagonal element\n");
	printf("\t\t\t main loop (k) = %d\n",k);
      }
    }
    /*
     *      Debug---
     */
    for(i=k; i<n; i++){
      temp.re = 0.0;
      temp.im = 0.0;
      for(m=0; m<k; m++){
	temp = DADD(temp,DMUL(a[i][m],a[m][k]));
      }
      a[i][k] = DSUB(a[i][k],temp);
    }
    /*
     *		Find pivot.
     */
    mmax = 0.0;
    imax = -1;
    for(i=k; i<n; i++){
      if(DABS(a[i][k]) >= mmax){
	mmax = DABS(a[i][k]);
	imax = i;
      }
    }
    /*
     *		Rearrange the ipivot elements to reflect the
     *		current positioning of the rows.  This method
     *		allows the values of ipivot[i] to be used in place
     *		of "i" when it comes to performing the forward
     *		substitution.
     */
    jj = ipivot[k];
    ipivot[k] = ipivot[imax];
    ipivot[imax] = jj;
    /*
     *		Perform the pivot if a pivot was found.
     */
    if(imax > 0){
      for(j=0; j<n; j++){
	temp = a[k][j];
	a[k][j] = a[imax][j];
	a[imax][j] = temp;
      }
    }
    /*
     *		Finish the decomposition.
     */
    for(j=k+1; j<n; j++){
      temp.re = 0.0;
      temp.im = 0.0;
      for(m=0; m<k; m++){
	temp = DADD(temp,DMUL(a[k][m],a[m][j]));
      }
      /*
       *          Debug+++
       */
      if(a[k][k].re == 0 && a[k][k].im == 0){
	printf("Going to divide by zero !\n");
      }
      /*
       *          Debug---
       */
      a[k][j] = DDIV(DSUB(a[k][j],temp),a[k][k]);
    }

  }      /* End of Main Loop */
}

/*
 *  Mlib_CLUsolve--
 *  performs forward and back substitution
 *
 *  input:	a -- input matrix containing LU decomposed matrix.
 *		n -- size of matrix.
 *		ipivot -- vector containing the pivot information.
 *		b -- vector containing the R.H.S. information, i.e.
 *			Ax = b.
 *		x -- outpout vector which contains the "answer".
 */
void Mlib_CLUsolve(a,n,ipivot,b,x)
     dcx_t 	**a,
  *b,
  *x;
     int 	n,*ipivot;
{
  int	i,j;
  dcx_t	temp;
  /*
   *	Setup first element in "y" vector, preparing for forward
   *	substitution.  Note that "x" is actually used since a separate
   *	vector "y" is unnecessary,  hence "x" does double duty.
   */
  x[0] = DDIV(b[ipivot[0]],a[0][0]);
  for(i=1; i<n; i++){
    temp.re = 0.0;
    temp.im = 0.0;
    for(j=0; j<i; j++){
      temp = DADD(temp,DMUL(a[i][j],x[j]));
    }
    x[i] = DDIV(DSUB(b[ipivot[i]],temp),a[i][i]);
  }
  /*
   *	Perform Backsubstitution.
   */
  for(i=n-2; i>=0; i--){
    temp.re = 0.0;
    temp.im = 0.0;
    for(j=i+1; j<n; j++){
      temp = DADD(temp,DMUL(a[i][j],x[j]));
    }
    x[i] = DSUB(x[i],temp);
  }
}
/* ############################################################## */
/*
 *				Double matrix, vector and functions.
 */
/*
 *  Mlib_DNewMat--
 *  allocates a double matrix with range [0..(row-1)][0..(col-1)]
 *
 */

doublem_t Mlib_DNewMat(int row, int col)
{
  int i;
  doublev_t doublev;
  doublem_t doublem;

  doublem = (doublem_t) malloc(sizeof(doublev_t) * row);
  doublev = (doublev_t) malloc(sizeof(double) * row * col);

  for (i = 0; i < row; i++) {
    doublem[i] = doublev + i * col;
  }

  return doublem;
}



/*
 *  Mlib_DNew3DMat--
 *  allocates a double 3D matrix 
 *
 */
doublemv_t Mlib_DNew3DMat(int iSize, int jSize, int kSize)
{
  int i;
  doublev_t doublev;
  doublem_t doublem;
  doublemv_t doublemv;

  doublemv = (doublemv_t) malloc(sizeof(doublem_t) * iSize);
  doublem = (doublem_t) malloc(sizeof(doublev_t) * iSize * jSize);
  doublev = (doublev_t) malloc(sizeof(double) * iSize * jSize * kSize);

  for (i = 0; i < iSize; i++) {
    doublemv[i] = doublem + i * jSize;
  }

  for (i = 0; i < iSize * jSize; i++) {
    doublem[i] = doublev + i * kSize;
  }

  return doublemv;
}


/*
 *  Mlib_DNew4DMat--
 *  allocates a double 4D matrix
 *
 */
doublemm_t Mlib_DNew4DMat(int iSize, int jSize, int kSize, int lSize)
{
  int i;
  doublev_t doublev;
  doublem_t doublem;
  doublemv_t doublemv;
  doublemm_t doublemm;

  doublemm = (doublemm_t) malloc(sizeof(doublemv_t) * iSize);
  doublemv = (doublemv_t) malloc(sizeof(doublem_t) * iSize * jSize);
  doublem = (doublem_t) malloc(sizeof(doublev_t) * iSize * jSize * kSize);
  doublev = (doublev_t) malloc(sizeof(double) * iSize * jSize * kSize * lSize);

  for (i = 0; i < iSize; i++){
    doublemm[i] = doublemv + i * jSize ;
  }

  for (i = 0; i < iSize * jSize; i++) {
    doublemv[i] = doublem + i * kSize;

  }

  for (i = 0; i < iSize * jSize * kSize; i++) {
    doublem[i] = doublev + i * lSize;
  }
	
  return doublemm ;

}



/*
 *  Mlib_DNew4DMat--
 *  allocates a double 4D matrix
 *
 */
dcxmm_t Mlib_CNew4DMat(int iSize, int jSize, int kSize, int lSize)
{
  int i;
  dcxv_t cx_doublev;
  dcxm_t cx_doublem;
  dcxmv_t cx_doublemv;
  dcxmm_t cx_doublemm;

  cx_doublemm = (dcxmm_t) malloc(sizeof(dcxmv_t) * iSize);
  cx_doublemv = (dcxmv_t) malloc(sizeof(dcxm_t) * iSize * jSize);
  cx_doublem = (dcxm_t) malloc(sizeof(dcxv_t) * iSize * jSize * kSize);
  cx_doublev = (dcxv_t) malloc(sizeof(dcx_t) * iSize * jSize * kSize * lSize);

  for (i = 0; i < iSize; i++){
    cx_doublemm[i] = cx_doublemv + i * jSize ;
  }

  for (i = 0; i < iSize * jSize; i++) {
    cx_doublemv[i] = cx_doublem + i * kSize;

  }

  for (i = 0; i < iSize * jSize * kSize; i++) {
    cx_doublem[i] = cx_doublev + i * lSize;
  }
	
  return cx_doublemm ;

}


/*
 *  Mlib_DNewVec--
 *  allocates memory for a double vector 
 *
 */
doublev_t Mlib_DNewVec(int iSize)
{
  return (doublev_t) malloc(sizeof(double) * iSize);
}

/*
 *  Mlib_DFreeMat--
 *  frees a double matrix m 
 *
 */
void Mlib_DFreeMat(doublem_t doublem)
{
  free((VOIDPTR) *doublem);
  free((VOIDPTR) doublem);
  return;
}

/*
 *  Mlib_DFree3DMatrix--
 *  deallocates a double 3D matrix 
 *
 */
void Mlib_DFree3DMat(doublemv_t doublemv)
{
  free((VOIDPTR) **doublemv);
  free((VOIDPTR) *doublemv);
  free((VOIDPTR) doublemv);
  return;
}

/*
 *  Mlib_DFree4DMatrix--
 *  deallocates a double 4D matrix 
 *
 */
void Mlib_DFree4DMat(doublemm_t doublemm)
{
  free((VOIDPTR) ***doublemm);
  free((VOIDPTR) **doublemm);
  free((VOIDPTR) *doublemm);
  free((VOIDPTR) doublemm);
  return;
}


/*
 *  Mlib_DFree4DMatrix--
 *  deallocates a float 4D matrix 
 *
 */
void Mlib_FFree4DMat(floatmm_t floatmm)
{
  free((VOIDPTR) ***floatmm);
  free((VOIDPTR) **floatmm);
  free((VOIDPTR) *floatmm);
  free((VOIDPTR) floatmm);
  return;
}

/*
 *  Mlib_DFreeVec--
 *  frees the double vector v
 *
 */
void Mlib_DFreeVec(doublev_t doublev)
{
  free((VOIDPTR) doublev);
  return;
}



/*
 *  Mlib_DTypeMat--
 *	Types a matrix (or portion of a matrix) to the screen
 */
void Mlib_DTypeMat(a_PP,nrl,nrh,ncl,nch,exp_float)                     
                                                                            
     double           **a_PP;         /* the real matrix to print out         */  
     int             nrl,nrh,        /* the lower and upper row dimensions   */  
  ncl,nch;        /* the lower and upper column dimensions*/  
     int             exp_float; /* 1: print it out in exponential else float */  
{
  int         i,j;                                                        
                                                                            
  printf ("\n");                                                          
  if (exp_float){
    for (i=nrl; i<=nrh; i++){
      for (j=ncl; j<=nch; j++) {                                      
	printf ("(% .2e) ",a_PP[i][j]);                             
	if (j==nch)                                                 
	  printf ("\n");                                          
      }
    }
  }
  else{
    for (i=nrl; i<=nrh; i++){
      for (j=ncl; j<=nch; j++) {                                     
	printf ("(% .2f) ",a_PP[i][j]);                            
	if (j==nch)                                                
	  printf ("\n");                                         
      }                                                              
    }
  }
  printf ("\n");                                                         
  return;                                                                   
}
                                                                           
/*
 *  Mlib_DTypeVec--
 *	Type a vector (or portion of a vector) to the screen
 */
void Mlib_DTypeVec(a_P,nrl,nrh,exp_float)                             
                                                                           
     double           *a_P;           /* the real vector to print out         */ 
     int             nrl,nrh;        /* the lower and upper row dimensions   */ 
     int             exp_float; /* 1: print it out in exponential else float */ 
                                                                           
{
  int         i,endofline;
                                               
  printf ("\n");                             
  if (exp_float){
    endofline = 0;
    for (i=nrl; i<=nrh; i++){
      printf ("(% .2e) ",a_P[i]);       
      endofline++;
      if(endofline == 6){
	printf("\n");
	endofline = 0;
      }
    }
    printf ("\n");                             
  }
  else{                                       
    endofline = 0;
    for (i=nrl; i<=nrh; i++){
      printf ("(% .2f) ",a_P[i]);       
      endofline++;
      if(endofline == 6){
	printf("\n");
	endofline = 0;
      }
    }
    printf ("\n");                             
  }
  printf ("\n");                             
                                               
  return;                                       
}                                              


/*
 *  Mlib_DFprintMat--
 *	Print a matrix (or portion of one) to a file.
 */
void Mlib_DFprintMat(out_fp,a_PP,nrl,nrh,ncl,nch,exp_float)

     double    **a_PP;  /* the matrix to print out      */
     int             nrl,nrh, /* the lower and upper row dimensions   */
  ncl,nch;	/* the lower and upper column dimensions*/
     int             exp_float; /* 1: print it out in exponential else float */
     FILE	*out_fp;

{
  int         i,j;

  fprintf (out_fp,"\n");
  switch (exp_float){
  case 0:
    for (i=nrl; i<=nrh; i++){
      for (j=ncl; j<=nch; j++){
	fprintf (out_fp,"% 8.5e ",a_PP[i][j]);
      }
      fprintf (out_fp,"\n");
    }
    break;

  case 1:
    for (i=nrl; i<=nrh; i++){
      for (j=ncl; j<=nch; j++) {
	fprintf (out_fp,"(% 8.5e) ",a_PP[i][j]);
      }
      fprintf (out_fp,"\n");
    }
    break;

  case 2:
    for (i=nrl; i<=nrh; i++){
      for (j=ncl; j<=nch; j++) {
	fprintf (out_fp,"% .15e \n",a_PP[i][j]);
      }
    }
    break;
  }
  fprintf (out_fp,"\n");

  return;
}


/*
 *  Mlib_DFprintVec--
 *	Print a double vector (or part of one) to a file
 */
void Mlib_DFprintVec(out_fp,a_P,nrl,nrh,exp_float)

     double *a_P;     /* the vector to print out      */
     int nrl,nrh;  /* the lower and upper row dimensions   */
     int exp_float; /* 1: print it out in exponential else float */
     FILE *out_fp;

{
  int         i;

  fprintf (out_fp,"\n");
  switch (exp_float){
  case 0:
    for (i=nrl; i<=nrh; i++)
      fprintf(out_fp,"% 5.2e \n",a_P[i]);
    break;

  case 1:
    for (i=nrl; i<=nrh; i++)
      fprintf (out_fp,"(% 5.2e) \n",a_P[i]);
    break;

  case 2:
    for (i=nrl; i<=nrh; i++)
      fprintf (out_fp,"% .15e \n",a_P[i]);
    break;


  }
  fprintf (out_fp,"\n");

  return;
}

/*
 *  Mlib_DVecxVec--
 *  Multiple two double vectors (dot product)
 */
double Mlib_DVecxVec(a_P,b_P,n)
     double	*a_P,
  *b_P;
     int	n;
{
  int i;
  double	c;
  c = 0.0;
  for(i=0; i<n; i++)
    c = c + a_P[i]*b_P[i];

  return(c);
}


/*
 *  Mlib_DMagVec--
 *  finds the 2 norm of a double vector
 *
 */
double Mlib_DMagVec(double *a_P, int n)
{
  int		i;
  double	mag;

  mag     = 0.0;
  for (i=0; i<n; i++){
    mag = mag + a_P[i]*a_P[i];
  }
  mag = sqrt(mag);
  return(mag);
}

/*
 *  Mlib_DSubVec--
 *  Subtract one double vector from another:  c_P = a_P - b_P
 */
void Mlib_DSubVec(c_P,a_P,b_P,n)
     double	*a_P,
  *b_P,
  *c_P;
     int	n;
{
  int	i;
  for(i=0; i<n; i++)
    c_P[i] = a_P[i] - b_P[i];
}

/*
 *  Mlib_DCopyPartVec--
 *  Copy part of a double vector to another double vector.
 *
 *  Copies elements "start" to "end" INCLUSIVE.
 */
void Mlib_DCopyPartVec(double *b_P,double *a_P,int start, int end, int bn)
{
  int	i;

  if((end - start) >= bn){
    printf("ERROR COPYING PARTIAL VECTOR..... DESTINATION TOO SMALL\n");
    return;
  }
	
  for(i=start; i<=end; i++)
    b_P[i-start] = a_P[i];
}

/*
 *  Mlib_DCopyVec--
 *  Copy a double vector to another double vector.
 */
void Mlib_DCopyVec(b_P,a_P,n)
     double	*a_P,
  *b_P;
     int	n;
{
  int	i;
  for(i=0; i<n; i++)
    b_P[i] = a_P[i];
}
/*
 *  Mlib_DScalxVec--
 *  Multiple a double vector by a scalar constant.
 */
void Mlib_DScalxVec(x_P,scal,n)
     double	*x_P,
  scal;
     int	n;
{
  int i;
  for(i=0; i<n; i++)
    x_P[i] = scal*x_P[i];

}


/*
 *  Mlib_DMatxVec--
 *  Multiplies a double matrix by a double vector:  A*x
 */
void Mlib_DMatxVec(y_P,a_PP,x_P,n)
     double	*y_P,
  **a_PP,
  *x_P;
     int	n;
{
  int	i,j;

  for(i=0; i<n; i++){
    y_P[i] = 0.0;
    for(j=0; j<n; j++){
      y_P[i] = y_P[i] + x_P[j]*a_PP[i][j];
    }
  }

}


/*
 *  Mlib_DInv--
 *  inverts a nxn real matrix
 *
 */
void Mlib_DInv(y_PP,a_PP,n)

     double 	**y_PP,		/* a pointer to a pointer of the output matrix	*/
  **a_PP;		/*      "    "  "    "    "   "  input  "       */
     int 	n;

{
  double 	d;
  double	*col_P;
  int 	i,j,
    *indx_P;
  /*indx_P[40];*/
  void	LUDcmp(),LUBackSubst();

  CHECK_TEMP_DOUBLE_P(n);
  CHECK_TEMP_INT_P(n);
  col_P = mlTemp_Double_P;
  indx_P = mlTemp_Int_P;
  LUDcmp(a_PP,0,n-1,indx_P,&d);

  for(j=0; j<n; j++) {
    for(i=0;i<n;i++) 
      col_P[i]=0.0;
    col_P[j]=1.0;
    LUBackSubst(a_PP,0,n-1,indx_P,col_P);
    for (i=0; i<n; i++) 
      y_PP[i][j]=col_P[i];
  }
  return;
}
 
 
 

/*
 *  LUDcmp--
 *  performs a LU decomposition on the a matrix
 *  See "Numerical Recipes in C", pg. 43 for details
 *
 */
void LUDcmp(a_PP,nl,nh,indx_P,d_P)

     double 	**a_PP;
     int 	nl,nh,
  *indx_P;
     double	*d_P;

{
  int 	i,imax,j,k;
  double 	big,dum,sum,temp;
  double 	*vv_P;

  CHECK_TEMP_DOUBLE_P1(nh+1);
  vv_P =  mlTemp_Double_P1;

  *d_P=1.0;
  for (i=nl; i<=nh; i++) {
    big=0.0;
    for (j=nl; j<=nh; j++)
      if((temp=fabs(a_PP[i][j])) > big) 
	big =temp;
    if (big == 0.0) 
      printf("Singular matrix in LUDcmp\n");
    vv_P[i]=1.0/big;
  }
  for (j=nl; j <=nh; j++) {
    for (i=nl; i<j; i++) {
      sum=a_PP[i][j];
      for (k=nl; k<i; k++) 
	sum -= a_PP[i][k]*a_PP[k][j];
      a_PP[i][j]=sum;
    }
    big=0.0;
    for (i=j; i<=nh; i++) {
      sum=a_PP[i][j];
      for (k=nl; k<j; k++)
	sum -= a_PP[i][k]*a_PP[k][j];
      a_PP[i][j]=sum;
      if ((dum=vv_P[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=nl; k<=nh; k++) {
	dum=a_PP[imax][k];
	a_PP[imax][k]=a_PP[j][k];
	a_PP[j][k]=dum;
      }
      *d_P= -(*d_P);
      vv_P[imax]=vv_P[j];
    }
    indx_P[j]=imax;
    if (a_PP[j][j] == 0.0) 
      a_PP[j][j] = TINY;
    if (j != nh) {
      dum=1.0/(a_PP[j][j]);
      for (i=j+1; i<=nh; i++) 
	a_PP[i][j] *= dum;
    }
  }
  return;
}




/*
 *  LUBackSubst--
 *  performs forward and back substitution
 *
 */
 
void LUBackSubst(a_PP,nl,nh,indx_P,b_P)

     double 	**a_PP,		/* pointer to the pointer of the LU decomposed
			   matrix a					*/
  *b_P;
     int 	nl,nh,*indx_P;
{

  int 	i,ii,ip,j;
  double 	sum;
  ii=nl-1;
 
  for (i=nl; i<=nh; i++) {
    ip=indx_P[i];
    sum=b_P[ip];
    b_P[ip]=b_P[i];
    if (ii != nl- 1)
      for (j=ii; j<=i-1; j++) 
	sum -= a_PP[i][j]*b_P[j];
    else if (sum) 
      ii=i;
    b_P[i]=sum;
  }
  for (i=nh; i>=nl; i--) {
    sum=b_P[i];
    for (j=i+1; j<=nh; j++) 
      sum -= a_PP[i][j]*b_P[j];
    b_P[i]=sum/a_PP[i][i];
  }
}

/*
  * Mlib_DLog10--
  */
double Mlib_DLog10(double x)
{
  return log(x)/2.302585093;
}

/* ############################################################## */
/*
 *				Float matrix and vector functions.
 *
 */
 
/*
 *  Mlib_FNewVec--
 *  allocates memory for a float vector
 *
 */
floatv_t Mlib_FNewVec(int iSize)
{
  return (floatv_t) malloc(sizeof(float) * iSize);
}
 

/*
 *  Mlib_FFreeVec--
 *  frees the float vector v_P 
 *
 */
void Mlib_FFreeVec(floatv_t floatv)
{
  free((VOIDPTR) floatv);
  return;
}


/*
 *  Mlib_FNewMat--
 *  allocates a float matrix with range [0..(row-1)][0..(col-1)]
 *
 */

floatm_t Mlib_FNewMat(int row, int col)
{
  int i;
  floatv_t floatv;
  floatm_t floatm;

  floatm = (floatm_t) malloc(sizeof(floatv_t) * row);
  floatv = (floatv_t) malloc(sizeof(float) * row * col);

  for (i = 0; i < row; i++) {
    floatm[i] = floatv + i * col;
  }

  return floatm;
}

 
/*
 *  Mlib_FFreeMat--
 *  frees a float matrix m 
 *
 */
void Mlib_FFreeMat(floatm_t floatm)
{
  free((VOIDPTR) *floatm);
  free((VOIDPTR) floatm);
  return;
}



/*
 *  Mlib_FNew3DMat--
 *  allocates a float 3D matrix 
 *
 */
floatmv_t Mlib_FNew3DMat(int iSize, int jSize, int kSize)
{
  int i;
  floatv_t floatv;
  floatm_t floatm;
  floatmv_t floatmv;

  floatmv = (floatmv_t) malloc(sizeof(floatm_t) * iSize);
  floatm = (floatm_t) malloc(sizeof(floatv_t) * iSize * jSize);
  floatv = (floatv_t) malloc(sizeof(float) * iSize * jSize * kSize);

  for (i = 0; i < iSize; i++) {
    floatmv[i] = floatm + i * jSize;
  }

  for (i = 0; i < iSize * jSize; i++) {
    floatm[i] = floatv + i * kSize;
  }

  return floatmv;
}

/*
 *  Mlib_FFree3DMatrix--
 *  deallocates a float 3D matrix 
 *
 */
void Mlib_FFree3DMat(floatmv_t floatmv)
{
  free((VOIDPTR) **floatmv);
  free((VOIDPTR) *floatmv);
  free((VOIDPTR) floatmv);
  return;
}


/*
 *  Mlib_FNew4DMat--
 *  allocates a float 4D matrix
 *
 */
floatmm_t Mlib_FNew4DMat(int iSize, int jSize, int kSize, int lSize)
{
  int i;
  floatv_t floatv;
  floatm_t floatm;
  floatmv_t floatmv;
  floatmm_t floatmm;

  floatmm = (floatmm_t) malloc(sizeof(floatmv_t) * iSize);
  floatmv = (floatmv_t) malloc(sizeof(floatm_t) * iSize * jSize);
  floatm = (floatm_t) malloc(sizeof(floatv_t) * iSize * jSize * kSize);
  floatv = (floatv_t) malloc(sizeof(float) * iSize * jSize * kSize * lSize);

  for (i = 0; i < iSize; i++){
    floatmm[i] = floatmv + i * jSize ;
  }

  for (i = 0; i < iSize * jSize; i++) {
    floatmv[i] = floatm + i * kSize;

  }

  for (i = 0; i < iSize * jSize * kSize; i++) {
    floatm[i] = floatv + i * lSize;
  }
	
  return floatmm ;

}


/* ##############################################################
 *
 *				Char matrix and vector functions.
 */

/*
 *  Mlib_CharNewMat--
 *  allocates a char matrix with range [0..(row-1)][0..(col-1)]
 *
 */

charm_t Mlib_CharNewMat(int row, int col)
{
  int i;
  charv_t charv;
  charm_t charm;

  charm = (charm_t) malloc(sizeof(charv_t) * row);
  charv = (charv_t) malloc(sizeof(char) * row * col);

  for (i = 0; i < row; i++) {
    charm[i] = charv + i * col;
  }

  return charm;
}



/*
 *  Mlib_CharFreeMat--
 *  frees a char matrix m 
 *
 */
void Mlib_CharFreeMat(charm_t charm)
{
  free((VOIDPTR) *charm);
  free((VOIDPTR) charm);
  return;
}


/*
 * Mlib_INew3DMat
 *
 * Allocate an int 3D matrix
 */ 
integermv_t Mlib_INew3DMat(int iSize, int jSize, int kSize )
{
  int i;
  int *intv;
  int **intm;
  int ***intmv;

  intmv = (int ***) malloc(sizeof(int **) * iSize );
  intm  = (int **) malloc(sizeof(int *) * iSize * jSize);
  intv  = (int *) malloc(sizeof(int) * iSize * jSize * kSize );

  for (i = 0; i < iSize; i++ ) {
    intmv[i] = intm + i * jSize ;
  }

  for (i = 0; i < iSize * jSize ; i++) {
    intm[i] = intv + i * kSize ;
  }

  return intmv;
}



/*
 * Mlib_INew4DMat
 *
 * Allocate an int 4D matrix
 */ 
integermm_t Mlib_INew4DMat(int iSize, int jSize, int kSize, int lSize )
{
  int i;
  int *intv;
  int **intm;
  int ***intmv;
  integermm_t intmm;

  intmm = (integermm_t) malloc(sizeof(int ***) * iSize );
  intmv  = (int ***) malloc(sizeof(int **) * iSize * jSize);
  intm  = (int **) malloc(sizeof(int *) * iSize * jSize * kSize );
  intv = (int *) malloc(sizeof(int) * iSize * jSize * kSize *
			lSize );

  for (i = 0; i < iSize; i++ ) {
    intmm[i] = intmv + i * jSize ;
  }

  for (i = 0; i < iSize * jSize; i++ ) {
    intmv[i] = intm + i * kSize ;
  }

  for (i = 0; i < iSize * jSize * kSize ; i++) {
    intm[i] = intv + i * lSize ;
  }

  return(intmm);
}

/*
 *  Mlib_IFree3DMatrix--
 *  deallocates a double 3D matrix 
 *
 */
void Mlib_IFree3DMat(integermv_t integermv)
{
  free((VOIDPTR) **integermv);
  free((VOIDPTR) *integermv);
  free((VOIDPTR) integermv);
  return;
}

/*
 *  Mlib_IFree4DMatrix--
 *  deallocates a double 4D matrix 
 *
 */
void Mlib_IFree4DMat(integermm_t integermm)
{
  free((VOIDPTR) ***integermm);
  free((VOIDPTR) **integermm);
  free((VOIDPTR) *integermm);
  free((VOIDPTR) integermm);
  return;
}


double expint(int n, double x)
{
  int i,ii,nm1;
  double a,b,c,d,del,fact,h,psi,ans;

  nm1=n-1;
  if (n < 0 || x < 0.0 || (x==0.0 && (n==0 || n==1)))
    printf("bad arguments in expint");
  else 
    {
      if (n == 0)
	ans=exp(-x)/x;
      else
	{
	  if (x == 0.0)
	    ans=1.0/nm1;
	  else
	    {
	      if (x > 1.0)
		{
		  b=x+n;
		  c=1.0/FPMIN;
		  d=1.0/b;
		  h=d;
		  for (i=1;i<=MAXIT;i++)
		    {
		      a = -i*(nm1+i);
		      b += 2.0;
		      d=1.0/(a*d+b);
		      c=b+a/c;
		      del=c*d;
		      h *= del;
		      if (fabs(del-1.0) < MACHEPS)
			{
			  ans=h*exp(-x);
			  return ans;
			}
		    }
		  printf("continued fraction failed in expint");
		}
	      else 
		{
		  ans = (nm1!=0 ? 1.0/nm1 : -log(x)-EULER);
		  fact=1.0;
		  for (i=1;i<=MAXIT;i++)
		    {
		      fact *= -x/i;
		      if (i != nm1)
			del = -fact/(i-nm1);
		      else
			{
			  psi = -EULER;
			  for (ii=1;ii<=nm1;ii++)
			    psi += 1.0/ii;
			  del=fact*(-log(x)+psi);
			}
		      ans += del;
		      if (fabs(del) < fabs(ans)*MACHEPS)
			return ans;
		    }
		  printf("series failed in expint");
		}
	    }
	}
    }
  return ans;
}
#undef MAXIT
#undef MACHEPS
#undef FPMIN
#undef EULER
/* (C) Copr. 1986-92 Numerical Recipes Software -)#21|{Y". */

/*
 * Sinple integer linked liust handling routines.
 *
 * IAddToList	Add an integer to then end of an integer linked list and
 *		create the linked list of required.
 * IFreeList	Free the integer linked list
 *
 */
void IAddToList(int x, int_link_list_Pt *list)
{
  int_link_list_Pt entry_P, newEntry_P;
  newEntry_P = (int_link_list_Pt) malloc(sizeof(int_link_list_t));
  /* If pointer is null then it will need to be created. */
  if( (entry_P=*list) )
    {
      for(; entry_P->next; entry_P=entry_P->next);
      entry_P->next = newEntry_P;
    }
  else
    {
      *list = newEntry_P;
    }
  newEntry_P->x = x; 
  newEntry_P->next = NULL; 
  newEntry_P->prior = entry_P;
  /* Find end of list */
  return;
}

void IFreeList(int_link_list_Pt *list)
{
  int_link_list_Pt entry_P, nextEntry_P;
  for(entry_P = *list; entry_P; entry_P=nextEntry_P)
    {
      nextEntry_P = entry_P->next;
      free(entry_P);
    }
  *list = NULL;
  return;
}


/*
 * ml_Init
 *
 * Create space for vectors and other small quantities which ndo not need
 * to be deleted all of the time.
 */
void ml_Init()
{
}
 
/*
 * ml_CleanUp()
 *
 * Free space for vectors and other small quantities which ndo not need
 * to be deleted all of the time.
 */
void ml_CleanUp()
{
  if(mlTemp_Double_P) free(mlTemp_Double_P);
  if(mlTemp_Double_P1) free(mlTemp_Double_P1);
  if(mlTemp_Int_P) free(mlTemp_Int_P);
}

void ml_ResizeTempDouble_P(int size)
{
  if(mlTemp_Double_P) free(mlTemp_Double_P);
  mlTemp_Double_P = (double *) malloc(size*sizeof(double));
  mlTemp_Double_size = size;
}

void ml_ResizeTempDouble_P1(int size)
{
  if(mlTemp_Double_P1) free(mlTemp_Double_P1);
  mlTemp_Double_P1 = (double *) malloc(size*sizeof(double));
  mlTemp_Double_size1 = size;
}
void ml_ResizeTempInt_P(int size)
{
  if(mlTemp_Int_P) free(mlTemp_Int_P);
  mlTemp_Int_P = (int *) malloc(size*sizeof(int));
  mlTemp_Int_size = size;
}
