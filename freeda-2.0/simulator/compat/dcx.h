/*
 * dcx.h
 *
 * declarations for double-precision complex math routines
 *
 * Authors:
 * Joseph Nathan Hall
 * Mark S. Basel
 *
 */

extern	dcx_t	Cx_DComplex(double a, double b);
extern  dcx_t   Cx_DConjugate(dcx_t a);
extern	dcx_t	Cx_DNegative(dcx_t x);
extern	dcx_t	Cx_DAdd(dcx_t x, dcx_t y);
extern	dcx_t	Cx_DMult(dcx_t x, dcx_t y);
extern	dcx_t	Cx_DSub(dcx_t x, dcx_t y);
extern	dcx_t	Cx_DDiv(dcx_t x, dcx_t y);
extern	dcx_t	Cx_DExp(dcx_t x);
extern	dcx_t	Cx_DTanh(dcx_t x);
extern	dcx_t	Cx_DSinh(dcx_t x);
extern	dcx_t	Cx_DCosh(dcx_t x);
extern	dcx_t	Cx_DSqrt(dcx_t x);
extern	double	Cx_DAbs(dcx_t x);
extern	dcx_t	Cx_DInv(dcx_t x);
extern	double	Cx_DAng(dcx_t x);

/*
 * more convenient declarations
 */

#define	DCMPLX(x, y) Cx_DComplex(x, y)
#define	DCON(x)      Cx_DConjugate(x)
#define	DNEG(x)      Cx_DNegative(x)
#define	DADD(x, y)   Cx_DAdd(x, y)
#define	DMULT(x, y)  Cx_DMult(x, y)
#define	DMUL(x, y)   Cx_DMult(x, y)
#define	DSUB(x, y)   Cx_DSub(x, y)
#define DDIV(x, y)   Cx_DDiv(x, y)
#define DEXP(x)      Cx_DExp(x)
#define DTANH(x)     Cx_DTanh(x)
#define DSINH(x)     Cx_DSinh(x)
#define DCOSH(x)     Cx_DCosh(x)
#define DSQRT(x)     Cx_DSqrt(x)
#define DABS(x)      Cx_DAbs(x)
#define DINV(x)      Cx_DInv(x)
#define DANG(x)      Cx_DAng(x)

#define DLOG10(x)    Mlib_DLog10(x)


/*
 * for debugging
 */

#define DCXPRINT(f, x, str) {\
	fprintf(f, "%s ", str); \
	fprintf(f, "%.4e %c %.4ej\n", x.re, (x.im < 0 ? '-' : '+'), \
	fabs(x.im)); \
	}
