/*
 * ml.h
 *
 * header file for math and matrix library
 *
 * Authors:
 * Joseph Nathan Hall
 * Mark S. Basel
 */

#ifndef Mlib_H
#define Mlib_H
/*
 * Deleted from version 1.3.  typedef is made in cap.h.  The following
 * definition of VOIDPTR does not always work.
 *
 *  #ifndef VOIDPTR
 *  #ifdef mips
 *  #define VOIDPTR char *
 *  #else
 *  #define VOIDPTR void *
 *  #endif
 *  #endif
 */

#define	TINY	1.0e-20
/* extern dcx_t	ONE; */
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

extern	double **DTempAssignMat(int size);
extern	int DTempReleaseMat(delement *start,double **mat_ptr);
extern	int DTempReleaseMat0(delement *start,double **mat_ptr);
extern	double *DTempAssignVec(int size);
extern	int DTempReleaseVec(delement *start,double *vec_ptr);
extern	int DTempReleaseVec0(delement *start,double *vec_ptr);
extern	void DTempCreateMat(int size);
extern	void DTempDisplayBySize(delement *start, int size);
extern	void DTempDisplayMat(delement *start);
extern	void DTempDisplayVec(delement *start);
extern	void DTempAddToList(delement *el_new,
			    delement **first,
			    delement **last);

extern	dcx_t **CTempAssignMat(int size);
extern	int CTempReleaseMat(celement *start,dcx_t **mat_ptr);
extern	int CTempReleaseMat0(celement *start, dcx_t **mat_ptr);
extern	dcx_t *CTempAssignVec(int size);
extern	int CTempReleaseVec(celement *start, dcx_t *vec_ptr);
extern	int CTempReleaseVec0(celement *start,  dcx_t *vec_ptr);
extern	void CTempCreateMat(int size);
extern	void CTempDisplayBySize(celement *start, int size);
extern	void CTempDisplayMat(celement *start);
extern	void CTempDisplayVec(celement *start);
extern	void CTempAddToList(celement *el_new,
			    celement **first,
			    celement **last);
extern	void CTempStatus(celement *start);

/*
 * Declarations for integer matrice's
 */

extern integerv_t Mlib_INewVec(int iSize);
extern void Mlib_IFreeVec(integerv_t intv);
extern integermv_t Mlib_INew3DMat(int iSize, int jSize, int kSize);
extern integermm_t Mlib_INew4DMat(int iSize, int jSize, int kSize, int lSize);
extern void Mlib_IFree3DMat(integermv_t integermv);
extern void Mlib_IFree4DMat(integermm_t integermm);
extern integerm_t Mlib_INewMat(int isize, int jsize);
extern void Mlib_IFreeMat(integerm_t integerm);

extern	int **ITempAssignMat(int size);
extern	int ITempReleaseMat(ielement *start,int **mat_ptr);
extern	int ITempReleaseMat0(ielement *start, int **mat_ptr);
extern	int *ITempAssignVec(int size);
extern	int ITempReleaseVec(ielement *start, int *vec_ptr);
extern	int ITempReleaseVec0(ielement *start,  int *vec_ptr);
extern	void ITempCreateMat(int size);
extern	void ITempDisplayBySize(ielement *start, int size);
extern	void ITempDisplayMat(ielement *start);
extern	void ITempDisplayVec(ielement *start);
extern	void ITempAddToList(ielement *el_new,
			    ielement **first,
			    ielement **last);

/*
 * Declarations for double matrice's
 */
extern doublev_t Mlib_DNewVec(int size);
extern void Mlib_DFreeVec(doublev_t doublev);
extern doublem_t Mlib_DNewMat(int isize, int jsize);
extern void Mlib_DFreeMat(doublem_t doublem);
extern void Mlib_DInv(double **y_PP,double **a_PP,int n);
extern doublemv_t Mlib_DNew3DMat(int iSize, int jSize, int kSize);
extern void Mlib_DFree3DMat(doublemv_t doublemv);
extern doublemm_t Mlib_DNew4DMat(int iSize, int jSize, int kSize, int lSize);
extern void Mlib_DFree4DMat(doublemm_t doublemm);
extern floatmm_t Mlib_FNew4DMat(int iSize, int jSize, int kSize, int lSize);
extern void Mlib_FFree4DMat(floatmm_t floatmm);
extern void Mlib_DTypeMat(double **a_PP, int nrl, int nrh, int ncl,
			  int nch, int exp_float);
extern void Mlib_DTypeVec(double *a_P,int nrl, int nrh, int exp_float);
extern void Mlib_DFprintMat(FILE *out_fp, double **a_PP, int nrl, int nrh,
			    int ncl, int nch, int exp_float);
extern void Mlib_DFprintVec(FILE *out_fp, double *a_P, int nrl, 
			    int nrh, int exp_float);
extern void Mlib_DSubVec(double *c_P,double *a_P,double *b_P,int n);
extern double Mlib_DMagVec(double *a_P, int n);
extern void Mlib_DCopyPartVec(double *b_P,double *a_P,int start,int end,int
			      bn);
extern void Mlib_DCopyVec(double *b_P,double *a_P,int n);


extern dcxsmv_t Mlib_CNewSparse3DMat(int row, int col);
extern void Mlib_CFreeSparse3DMat(dcxsmv_t dcxsmv, int row, int col);
extern dcxm_t Mlib_CNewMat(int rows, int cols);
extern dcxmv_t Mlib_CNew3DMat(int iSize, int jSize, int kSize);
extern dcxmm_t Mlib_CNew4DMat(int iSize, int jSize, int kSize, int lSize);
extern void Mlib_CFree3DMat(dcxmv_t dcxmv);
extern void Mlib_CFreeMat(dcxm_t dcxm); 
extern dcxv_t Mlib_CNewVec(int iSize);
extern void Mlib_CFreeVec(dcxv_t dcxv);
extern void Mlib_CCopyMat(dcxm_t c_PP,dcxm_t  a_PP,int rows,int cols);
extern void Mlib_CCopyVec(dcxv_t c_P,dcxv_t  a_P,int n);
extern dcx_t Mlib_CMult(dcx_t,dcx_t);
extern void Mlib_CAddMat(dcxm_t c_PP,dcxm_t  a_PP,dcxm_t  b_PP,int n);
extern void Mlib_CSubMat(dcxm_t c_PP,dcxm_t a_PP,dcxm_t b_PP,int n);
extern void Mlib_CSubVec(dcxv_t c_P,dcxv_t a_P,dcxv_t b_P,int n);
extern double Mlib_CMagVec(dcxv_t a_P,int n);
extern double Mlib_CMaxVec(dcxv_t a_P,int n);
extern dcx_t Mlib_CInfNormVec(dcxv_t a_P,int n);
extern void Mlib_ConjTransp(dcxm_t c_PP,dcxm_t a_PP,int n);
extern void Mlib_CScalMat(dcxm_t c_PP,dcxm_t a_PP,dcx_t c,int n);
extern void Mlib_CMultMat(dcxm_t c_PP,dcxm_t a_PP,dcxm_t b_PP,int n);
extern void Mlib_CMultMatV(dcxmv_t c_PPP,dcxmv_t a_PPP,dcxmv_t b_PPP,int n,
			   int lengthOfVector);
extern void Mlib_CMultVec(dcxv_t c_P,dcxm_t a_PP,dcxv_t b_P,int n);
extern void Mlib_CMultVecRec(dcxv_t c_P, dcxm_t a_PP, dcxv_t b_P,int n,int m);
extern void Mlib_CMultVecRecTrans(dcxv_t c_P, dcxm_t a_PP, dcxv_t b_P,int n,int m);
extern void Mlib_CScalVec(dcxv_t a_P, dcx_t scale, int n);
extern void Mlib_CInv(dcxm_t c_PP, dcxm_t a_PP, int n);
extern void Mlib_CInvV(dcxmv_t c_PPP, dcxmv_t a_PPP, int n, int lengthOfVector);
extern void Mlib_CExpGL(dcxm_t c_PP, double L, dcxv_t eval_P,
			dcxm_t evec_PP, dcxm_t evecinv_PP, int n);
extern void Mlib_CExp(dcxm_t c_PP, float *coeff, dcxv_t eval_P, 
		      dcxm_t evec_PP, dcxm_t evecinv_PP, int n);
extern void Mlib_CTypeMat(dcxm_t a_PP, int nrl, int nrh, int ncl,
			  int nch, int exp_float);
extern void Mlib_CFprintMat(FILE *out_fp,dcxm_t a_PP,
			    int nrl,int nrh,int ncl,int nch,int exp_float);
extern void Mlib_CTypeVec(dcxv_t a_P,int nrl,int nrh,int exp_float);
extern void Mlib_CFprintVec(FILE *out_fp,dcxv_t a_P,
			    int nrl,int nrh,int exp_float);
extern void Mlib_CElementFlip(dcxm_t in_PP,dcxm_t out_PP,int n);
extern void Mlib_CReduce(int rank,dcxm_t array,int reduce); 
extern void Mlib_CLUDcomp(dcx_t **a, int n, int *ipivot);
extern void Mlib_CLUSolve(dcx_t **a, int n, int *ipivot, dcx_t *b, 
			  dcx_t *x);
extern void Mlib_CFprint3DMat(FILE *out_F,dcxmv_t a_PPP,int nrl,int nrh,
			      int ncl,int nch,int k,int exp_float);
extern void Mlib_CType3DMat(dcxmv_t a_PPP,int nrl,int nrh,
			    int ncl,int nch,int k,int exp_float);

/*
 * Declarations for float vectors
 */
extern floatv_t Mlib_FNewVec(int iSize);
extern void Mlib_FFreeVec(floatv_t floatv);
extern floatmv_t Mlib_FNew3DMat(int iSize, int jSize, int kSize);
extern void Mlib_FFree3DMat(floatmv_t floatmv);
extern floatm_t Mlib_FNewMat(int isize, int jsize);
extern void Mlib_FFreeMat(floatm_t floatm);

/*
 * Declarations for Char matrices
 */
extern char **Mlib_CharNewMat(int row, int col);
extern void Mlib_CharFreeMat(char **charm);

/* Integer Linked List Handling */
 
extern void IAddToList(int x, int_link_list_Pt *list);
extern void IFreeList(int_link_list_Pt *list);

#endif
