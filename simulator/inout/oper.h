/*
 * Operator analysis routines
 */

typedef struct{
   int rank;
   dcxmv_t   dcxmv;
   } s_dcxmv_t, *s_dcxmv_t_Pt;

 typedef struct{
   int rank;
   dcxm_t   dcxm;
   } s_dcxm_t, *s_dcxm_t_Pt;

typedef struct{
   int rank;
   doublem_t   dm;
   } s_dm_t, *s_dm_t_Pt;

typedef union {
	gr_Id_t		tid;
	gr_Id_t		nid;
	gr_Id_t		eid;
	short		opType;
	char		*s;
	int		i;
	double		d;
	dcx_t		dcx;
	doublev_t	dv;
	s_dm_t		s_dm;
	dcxv_t		dcxv;
	s_dcxm_t        s_dcxm;
	s_dcxmv_t       s_dcxmv;
	int		(*func)();
} capOutObj_t, *capOutObj_Pt;
	
/*
 * type for representing output requests
 */
typedef struct capOutReqStruct {
	int		type;		/* type of object		*/
	int		protect;	/* if TRUE, do not delete during stack
					   operations */
	capOutObj_t	obj;		/* object itself		*/
	int		size;		/* size parameter if appropriate */
	doublev_t	x;		/* timebase or other x data	*/
					/* 	if appropriate		*/
	char 		*xName,		/* x name if available		*/
			*objName;	/* Object name if available	*/

	struct capOutReqStruct *next;	/* next in chain		*/
} capOutReq_t, *capOutReq_Pt;

extern void DumpTable(FILE *output_F, st_Table_Pt st_P);

extern int stackDepth;
extern capOutReq_t *Peek(int depth);
extern char *PeekName(int depth);
extern capOutReq_t *Pop(void);
extern int GetOutReqVariable(char *variableName, capOutReq_Pt *outReq_P);
extern int UpdateOutReqVariable(char *variableName, capOutReq_Pt *outReq_P);

/*
 * Utilities
 */
double An_LinInterpPoint(int n, doublev_t yv, doublev_t xv, double x);
void An_LinInterpVec(int n, doublev_t yv, doublev_t xv,
				int newn, doublev_t newyv, doublev_t newxv);
dcx_t An_LinInterpCPoint(int n, dcxv_t yv, doublev_t xv, double x);
void An_LinInterpCVec(int n, dcxv_t yv, doublev_t xv,
				int newn, dcxv_t newyv, doublev_t newxv);
void An_CopyTimebase(capOutReq_t *dest, capOutReq_t *src);
int An_CompareTimebase(capOutReq_t *arg1, capOutReq_t *arg2);
void An_Promote(capOutReq_t *arg1, capOutReq_t *arg2);


/*
 * Operators and arguments from an_oper.c
 */
int doc_tr(FILE *f_P);	/* don't include this with DOS version */

/*
 * File I/O
 */
int An_Op_WRITE(capOutReq_t *arg, capOutReq_t *result);
int An_Op_PLOT(capOutReq_t *result);
int An_Op_WATCH(capOutReq_t *result);
int An_Op_READ(capOutReq_t *arg, capOutReq_t *result);

/*
 * Manipulating complex vectors
 */
int An_Op_REAL(capOutReq_t *arg, capOutReq_t *result);
int An_Op_DB20(capOutReq_t *arg, capOutReq_t *result);
int An_Op_DB10(capOutReq_t *arg, capOutReq_t *result);
int An_Op_rad2deg(capOutReq_t *arg, capOutReq_t *result);
int An_Op_IMAG(capOutReq_t *arg, capOutReq_t *result);
int An_Op_MAG(capOutReq_t *arg, capOutReq_t *result);
int An_Op_PRINPHASE(capOutReq_t *arg, capOutReq_t *result);
int An_Op_CONTPHASE(capOutReq_t *arg, capOutReq_t *result);
int An_Op_CONJ(capOutReq_t *arg, capOutReq_t *result);
int An_Op_RECIP(capOutReq_t *arg, capOutReq_t *result);
int An_Op_ABS(capOutReq_t *arg, capOutReq_t *result);
int An_Op_XY2CX(capOutReq_t *arg, capOutReq_t *result);
int An_Op_CX2XY(capOutReq_t *arg, capOutReq_t *result);

int An_Op_DUP(capOutReq_t *arg, capOutReq_t *result);
int An_Op_GET(capOutReq_t *arg, capOutReq_t *index, 
						capOutReq_t *result);
int An_Op_PUT(capOutReq_t *arg, capOutReq_t *index,
				capOutReq_t *val, capOutReq_t *result);
int An_Op_STRIPX(capOutReq_t *arg, capOutReq_t *result);
int An_Op_PUSH(capOutReq_t *arg, capOutReq_t **result);
int An_Op_PACK(capOutReq_t *result);
int An_Op_CAT(capOutReq_t *arg1, capOutReq_t *result);
int An_Op_APPEND(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t **result);
int An_Op_SYSTEM(capOutReq_t *arg, capOutReq_t *result);

/*
 * Arithmetic vector ops -- unary
 */
int An_Op_NEG(capOutReq_t *arg, capOutReq_t *result);
int An_Op_RECIP(capOutReq_t *arg, capOutReq_t *result);

/*
 * Arithmetic vector ops -- binary
 */
int An_Op_ADD(capOutReq_t *arg1, capOutReq_t *arg2,
						capOutReq_t *result);
int An_Op_SUB(capOutReq_t *arg1, capOutReq_t *arg2,
						capOutReq_t *result);
int An_Op_MULT(capOutReq_t *arg1, capOutReq_t *arg2,
						capOutReq_t *result);
int An_Op_DIV(capOutReq_t *arg1, capOutReq_t *arg2,
						capOutReq_t *result);
int An_Op_MINLMT(capOutReq_t *arg, capOutReq_t *lmtarg,
						capOutReq_t *result);
int An_Op_MAXLMT(capOutReq_t *arg, capOutReq_t *lmtarg,
						capOutReq_t *result);
int An_Op_NOOP(capOutReq_t *arg1, capOutReq_t *arg2,
						capOutReq_t *result);
/*
 * Other vector ops
 */
int An_Op_DIFF(capOutReq_t *arg, capOutReq_t *result);
int An_Op_DERIV(capOutReq_t *arg, capOutReq_t *result);
int An_Op_SUM(capOutReq_t *arg, capOutReq_t *result);
int An_Op_INTEG(capOutReq_t *arg, capOutReq_t *result);
int An_Op_REPEAT(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result);
int An_Op_GETBIN(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result);
int An_Op_IMPULSE(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result);
int An_Op_SCALEX(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result);
int An_Op_FIMPULSE(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result);
/*
 * Misc. signal processing
 */
int An_Op_SMPLTIME(capOutReq_t *result);
int An_Op_SWEEPFREQ(capOutReq_t *result);
int An_Op_SMPLCVT(capOutReq_t *srcarg, capOutReq_t *tbarg,
						capOutReq_t *result);
int An_Op_SWEEPCVT(capOutReq_t *srcarg, capOutReq_t *sfarg,
						capOutReq_t *result);
int An_Op_MAKETIME(capOutReq_t *lenarg, capOutReq_t *ptsarg,
						capOutReq_t *result);
int An_Op_MAKESWEEP(capOutReq_t *frqarg, capOutReq_t *ptsarg,
						capOutReq_t *result);
/*
 * FFT and convolution
 */
int An_Op_FFT(capOutReq_t *arg, capOutReq_t *result);
int An_Op_INVFFT(capOutReq_t *arg, capOutReq_t *result);
int An_Op_SCONV(capOutReq_t *arg1, capOutReq_t *arg2,
						capOutReq_t *result);
int An_Op_CCONV(capOutReq_t *arg1, capOutReq_t *arg2,
						capOutReq_t *result);
int An_Op_UPCCONV(capOutReq_t *arg1, capOutReq_t *arg2,
						capOutReq_t *result);
int An_Op_FCONV(capOutReq_t *arg1, capOutReq_t *arg2,
						capOutReq_t *result);
int An_Op_LAST(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result);
int An_Op_RECT2POLAR(capOutReq_t *arg1, capOutReq_t *result);
int An_Op_POLAR2RECT(capOutReq_t *arg1, capOutReq_t *result);
int An_Op_LPBWFRQ(capOutReq_t *tbarg, capOutReq_t *corner,
				capOutReq_t *order, capOutReq_t *result);
int An_Op_ZEROPAD(capOutReq_t *arg, capOutReq_t *result);
int An_Op_GETLAST2N(capOutReq_t *arg, capOutReq_t *result);


/* From an_output.c */
void An_CalcError(const char *s);
int An_COPY(capOutReq_t *outReq, capOutReq_t **newReq);

/* From an_optbl.c */
int An_LookupOpNum(int opConst);
int An_LookupOpName(char *opName);

/*
 * Type for optbl
 */
typedef int (*ifunc)();
#define MAX_NO_OF_ARGS	4
typedef struct {
	char	*name;
	int	opConst,
		argCt,
		retVal;
	ifunc	func;
	int	argMask[MAX_NO_OF_ARGS];
        int     opPrecedence;  /* Operator evaluation precedence */

} an_OpTbl_t;

extern an_OpTbl_t an_OpTbl[];

#define GET_OPERATOR_NAME(outReq)      an_OpTbl[outReq->obj.opType].name
#define OPERATOR_PRECEDENCE(outReq)    an_OpTbl[outReq->obj.opType].opPrecedence


#define NO_ARG		0
#define ONE_ARG		1
#define ONE_ARG_RESULT_POINTER		-10000
#define ONE_OR_TWO_ARG 	-12
#define TWO_ARG		2
#define TWO_ARG_RESULT_POINTER		-20000
#define THREE_ARG	3
#define FOUR_ARG	4
#define	VAR_ARG		-1000

/*
 * Argument type masks
 */
#define INT_ARG_TYPE		0x00000001L
#define DOUBLE_ARG_TYPE		0x00000002L
#define	DCX_ARG_TYPE		0x00000004L
#define REAL_ARG_TYPE		0x00000003L
#define NUM_ARG_TYPE		0x00000007L

#define STRING_ARG_TYPE		0x00000008L

#define DOUBLEV_ARG_TYPE	0x00000010L
#define DCXV_ARG_TYPE		0x00000020L
#define VECTOR_ARG_TYPE		0x00000030L
#define	REAL_SV_ARG_TYPE	0x00000013L /* real "scalar or vector" type */
#define NUM_SV_ARG_TYPE		0x00000037L /* numeric "scalar or vector" */
#define	DCX_SV_ARG_TYPE		0x00000024L
#define	DOUBLEM_ARG_TYPE	0x00000040L
#define DOUBLEMV_ARG_TYPE	0x00000080L
#define DCXM_ARG_TYPE		0x00000100L
#define DCXMV_ARG_TYPE		0x00000200L

#define TERM_ARG_TYPE		0x00001000L
#define NODE_ARG_TYPE		0x00002000L
#define EDGE_ARG_TYPE		0x00004000L

#define FILE_ARG_TYPE		0x00010000L
#define	VAR_ARG_TYPE		0x00020000L
#define FILENAME_ARG_TYPE	0x00030000L

#define ANY_ARG_TYPE		0xffffffffL
#define	NO_ARG_TYPE		0x00000000L


#define GEN_TYPE_OUTREQ_LIST 1001	/* type of output request list	*/
					/*   in symbol table		*/

#define CAP_OBJ_NONEXISTENT 0
#define CAP_OBJ_TERM 1
#define CAP_OBJ_NODE 2
#define CAP_OBJ_EDGE 3

#define CAP_OBJ_DATAFILE 11
#define CAP_OBJ_FILENAME 12

#define CAP_OBJ_VAR 21

#define CAP_OBJ_OPER 101

#define CAP_OBJ_INT 201
#define CAP_OBJ_DOUBLE 202
#define CAP_OBJ_DCX 203
#define CAP_OBJ_DOUBLEV 204
#define CAP_OBJ_DCXV 205
#define CAP_OBJ_STRING 206

#define CAP_OBJ_DOUBLEM 211
#define CAP_OBJ_DOUBLEMV 212
#define CAP_OBJ_DCXM 213
#define CAP_OBJ_DCXMV 214

#define CAP_OBJ_FUNC 301

/*
 * Constant numbers for operators in the output section
 *
 * Try to keep these in the same order as in an_optbl.c, just for the
 * sake of neatness.
 */
#define CAP_OP_NONE	0
#define CAP_OP_DUP	11
#define CAP_OP_GET	12
#define CAP_OP_PUT	13
#define CAP_OP_PULL     14
#define CAP_OP_PUSH     15
#define CAP_OP_CAT     16

#define CAP_OP_STRIPX	21
#define CAP_OP_PACK	22
#define CAP_OP_APPEND	23


#define CAP_OP_SYSTEM	51
#define CAP_OP_CRT_OUT	52
#define CAP_OP_WRITE	53
#define CAP_OP_PLOT	54
#define CAP_OP_READ	55

#define CAP_OP_VT 101
#define CAP_OP_IT 102
#define CAP_OP_UT 103
#define CAP_OP_XT 104

#define CAP_OP_VF 105
#define CAP_OP_IF 106
#define CAP_OP_XF 107

#define CAP_OP_REAL 201
#define CAP_OP_IMAG 202
#define CAP_OP_MAG 203
#define CAP_OP_CONTPHASE 204
#define CAP_OP_PRINPHASE 205
#define CAP_OP_CONJ 206
#define CAP_OP_XY2CX 211
#define CAP_OP_CX2XY 212
#define CAP_OP_SCALEX 213

#define CAP_OP_RECIP 301
#define CAP_OP_NEG 302
#define CAP_OP_DB 303
#define CAP_OP_DB10 304
#define CAP_OP_RAD2DEG 305
#define CAP_OP_DEG2RAD 306
#define CAP_OP_ABS 307
#define CAP_OP_MINLMT 308
#define CAP_OP_MAXLMT 309

#define CAP_OP_ADD 401
#define CAP_OP_SUB 402
#define CAP_OP_MULT 403
#define CAP_OP_DIV 404

#define CAP_OP_DIFF 501
#define CAP_OP_DERIV 502
#define CAP_OP_SUM 503
#define CAP_OP_INTEG 504
#define CAP_OP_SMPLTIME 601
#define CAP_OP_SWEEPFREQ 602
#define CAP_OP_SMPLCVT 603
#define CAP_OP_SWEEPCVT 604
#define CAP_OP_MAKETIME 611
#define CAP_OP_MAKESWEEP 612
#define CAP_OP_REPEAT 613
#define CAP_OP_GETBIN 614
#define CAP_OP_IMPULSE 615
#define CAP_OP_FIMPULSE 616

#define CAP_OP_FFT 701
#define CAP_OP_INVFFT 702
#define CAP_OP_SCONV 703
#define CAP_OP_CCONV 704
#define CAP_OP_UPCCONV 705
#define CAP_OP_FCONV 706
#define CAP_OP_LAST 707
#define CAP_OP_RECT2POLAR 708
#define CAP_OP_POLAR2RECT 709
#define CAP_OP_ZEROPAD 710
#define CAP_OP_GETLAST2N 711

#define CAP_OP_LPBWFRQ 801
#define CAP_OP_IPD 802
#define CAP_OP_BIPD 803
#define CAP_OP_IPDINFO 804
#define CAP_OP_IPD_TM 805

