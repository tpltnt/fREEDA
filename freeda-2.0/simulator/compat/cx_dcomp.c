/*
 * cx_dcomplex.c
 *
 * A library of functions for complex math.
 *
 * Author:
 * Mark S. Basel
 *
 */


#include <stdio.h>
#include <math.h>
#include "stuff.h"
#include "dcx.h"

/*
 *  Create a complex number from two doubles.
 */
dcx_t Cx_DComplex(a,b)
     double	a,b;
{
  dcx_t	x;

  x.re = a;
  x.im = b;

  return(x);
}
/*
 *  Find the complex conjugate.
 */
dcx_t Cx_DConjugate(dcx_t a)
{
  dcx_t x;
  x.re = a.re;
  x.im = -a.im;

  return(x);
}
/*
 * Find the negative of a complex number.
 */
dcx_t Cx_DNegative(dcx_t x)
{
  dcx_t y;
  y.re = -x.re;
  y.im = -x.im;
  return(y);
}
/*
 * Find the inverse of a complex number.
 *
 *     x.re - jx.im
 * y = ------------
 *         |x|
 */
dcx_t Cx_DInv(dcx_t x)
{
  return(DDIV(DCMPLX(x.re,-x.im),DCMPLX(x.re*x.re+x.im*x.im,0.0)));
}

/*
 * Calculate the absolute value of x, x complex
 *
 * Method:
 *		If x.re and x.im is small, use a taylor series
 *		expansion to find the absolute value of x.  Otherwise
 *		use sqrt(x.re^2 + x.im^2).
 */
double Cx_DAbs(x)
     dcx_t     x;
{
  double	a,b,u,r,z,p;
  static	double	minimum = 1.0e-35;

  a = fabs(x.re);
  b = fabs(x.im);
  if(a>b){
    u = a;
    z = b;
  }
  else{
    u = b;
    z = a;
  }
  if(u < minimum && z < minimum) {
    if( u != 0.0 ){
      r = z/u;
      r = r*r;
      p = u*(1 + r/2 - r*r/8 + r*r*r/16 - 5*r*r*r*r/128
	     + 7*r*r*r*r*r/256 - 21*r*r*r*r*r*r/1024);
      return(p);
    }
    else{
      p = 0.0;
      return(p);
    }
  }
  else {
    p = sqrt(a*a + b*b);
    return(p);
  }
}

/*
 * Calculate x + y, x and y complex
 *
 * x = x1 + jx2
 * y = y1 + jy2
 * z = (x1 + y1) + j(x2 + y2)
 */ 

dcx_t Cx_DAdd(x,y)
     dcx_t x,y;
{

  dcx_t	z;

  z.re = x.re + y.re;
  z.im = x.im + y.im;
  return(z);
}

/*
 * Subtract: x - y, x and y complex
 *
 * x = x1 + jx2
 * y = y1 + jy2
 * z = (x1 - y1) + j(x2 - y2)
 */

dcx_t Cx_DSub(x,y)
     dcx_t x,y;
{

  dcx_t	z;

  z.re = x.re - y.re;
  z.im = x.im - y.im;
  return(z);
}


/*
 * Calculate x * y, x and y complex
 *
 * x = x1 + jx2
 * y = y1 + jy2
 * z = (x1 * y1 - x2 * y2) + j(x1 * y2 + x2 * y1)
 */

dcx_t Cx_DMult(x,y)
     dcx_t x,y;
{

  dcx_t	z;

  z.re = x.re * y.re - x.im * y.im;
  z.im = x.re * y.im + x.im * y.re;
  return(z);
}

/*       divide two complex numbers:   z = x/y   
 *       From MIX math library.
 */
/*dcx_t Cx_DDiv(dcx_t x, dcx_t y)
 *{
 *dcx_t z;
 *double temp,denom;
 * 
 *       if (fabs(y.im) > fabs(y.re)) {
 *               temp=y.re/y.im;
 *               denom=1.0/(temp*y.re + y.im);
 *               z.im= (double) denom*(x.im*temp-x.re);
 *               z.re= (double) denom*(x.re*temp+x.im);
 *       }
 *      else {
 *                temp=y.im/y.re;
 *               denom=1.0/(temp*y.im + y.re);
 *               z.im= (double) denom*(x.im-temp*x.re);
 *               z.re= (double) denom*(x.re+temp*x.im);
 *       }
 *       return z;
 *}
 */

/*
 * Calculate x / y, x and y complex
 *
 * x = x1 + jx2
 * y = y1 + jy2
 * z = [(x1 * y1 + x2 * y2)/(y1 * y1 + y2 * y2)] +
 *    j[(x2 * y1 - x1 * y2)/(y1 * y1 + y2 * y2)]     
 */
dcx_t Cx_DDiv(x, y)
     dcx_t x,y;
{
 
  dcx_t	z;
  register double	tmp;
 
  tmp = y.re * y.re + y.im * y.im;
  if (tmp == 0.0) {
    fprintf(output_F, "ERROR[cx_dcomplex.c] divide by 0.0 in Cx_DDiv\n");
    z.re = 1.0e10;
    z.im = 0.0;
  }
  else {
    z.re = (x.re * y.re + x.im * y.im) / tmp;
    z.im = (x.im * y.re - x.re * y.im) / tmp;
  }
  return(z);
}
 

/*
 * Calculate e raised to the complex power x
 *
 * x = x1 + jx2
 * z = [exp(x1) * cos(x2)] + j[exp(x1) * sin(x2)] 
 */

dcx_t Cx_DExp(x)
     dcx_t x;
{

  dcx_t	z;

  z.re = exp(x.re) * cos(x.im);
  z.im = exp(x.re) * sin(x.im);
  return(z);
}


/*
 *  cosh(x)
 *     exp(x) + exp(-x)
 *   = ----------------
 *            2
 */
dcx_t Cx_DCosh(dcx_t x)
{
  dcx_t  nx;

  nx.re = -x.re;
  nx.im = -x.im;

  return(DDIV(DADD(DEXP(x),DEXP(nx)),DCMPLX(2,0)));
}

/*
 *  sinh(x)
 *     exp(x) - exp(-x)
 *   = ----------------
 *            2
 */
dcx_t Cx_DSinh(dcx_t x)
{
  dcx_t  nx;

  nx.re = -x.re;
  nx.im = -x.im;

  return(DDIV(DSUB(DEXP(x),DEXP(nx)),DCMPLX(2,0)));
}


/*
 * tanh(x)
 *	  exp(x) - exp(-x)
 *	= ----------------
 *	  exp(x) + exp(-x)
 */
dcx_t Cx_DTanh(x)
     dcx_t x;
{
  dcx_t	nx,num,den;

  nx = x;
  nx.re = -nx.re;
  nx.im = -nx.im;
  num = DSUB(DEXP(x),DEXP(nx));
  den = DADD(DEXP(x),DEXP(nx));
  nx = DDIV(num,den);

  return(nx);

}

/*
 * Calculate a square root
 *
 * convert to polar, take the square root of the magnitude and halve theta,
 * then convert back to rectangular
 *
 * Now also check for the case where the point lies on one of the quad.
 * Where the real or the imaginary part is zero.
 *
 *              Im
 *               | 
 *           3   2   1
 *               |
 *        ---4-------8---Re
 *               |
 *           5   6   7
 *               |
 */
dcx_t Cx_DSqrt(x)
     dcx_t x;
{

  double	m, theta;
  int quad;
  static double pi = 3.1415926535897932384626433;
  quad = 0;
  if(x.re > 0){
    if(x.im > 0.0)
      quad = 1;
    if(x.im == 0.0)
      quad = 8;
    if(x.im < 0.0)
      quad = 7;
  }
  if(x.re == 0.0){
    if(x.im > 0.0)
      quad = 2;
    if(x.im == 0.0)
      quad = 9;
    if(x.im < 0.0)
      quad = 6;
  }
  if(x.re < 0.0){
    if(x.im > 0.0)
      quad = 3;
    if(x.im == 0.0)
      quad = 4;
    if(x.im < 0.0)
      quad = 5;
  }
  /*	m = hypot(x.re, x.im);  This doesn't work with a 387sx and POWERC */
  m = sqrt(x.re*x.re + x.im*x.im);
  m = sqrt(m);
  switch(quad){
  case 1:
    theta = atan(x.im/x.re)/2.0;
    x.re = m * cos(theta);
    x.im = m * sin(theta);
    break;
  case 2:
    x.re = m/sqrt(2.0);
    x.im = m/sqrt(2.0);
    break;
  case 3:
    theta = (atan(x.im/x.re) + pi)/2.0;
    x.re = m * cos(theta);
    x.im = m * sin(theta);
    break;
  case 4:
    x.re = 0.0;
    x.im = m;
    break;
  case 5:
    theta = (atan(x.im/x.re) - pi)/2.0;
    x.re = m * cos(theta);
    x.im = m * sin(theta);
    break;
  case 6:
    x.re = m/sqrt(2.0);
    x.im = -m/sqrt(2.0);
    break;
  case 7:
    theta = atan(x.im/x.re)/2.0;
    x.re = m * cos(theta);
    x.im = m * sin(theta);
    break;
  case 8:
    x.re = m;
    x.im = 0.0;
    break;
  case 9:
    x.re = 0.0;
    x.im = 0.0;
    break;
  default:
    theta = atan(x.im/x.re)/2.0;
    x.re = m * cos(theta);
    x.im = m * sin(theta);
    fprintf(output_F,"ERROR IN Cx_DSqrt: Can't find quadrant.  ");
    fprintf(output_F,"Using default value.\n");
    break;
  }
  return(x);
}

/*
 * Calculate Phase Angle
 *
 *              Im
 *               | 
 *           3   2   1
 *               |
 *        ---4-------8---Re
 *               |
 *           5   6   7
 *               |
 */
double Cx_DAng(dcx_t x)
{

  double	theta;
  int quad;
  static double pi = 3.1415926535897932384626433;
  quad = 0;
  if(x.re > 0){
    if(x.im > 0.0)
      quad = 1;
    if(x.im == 0.0)
      quad = 8;
    if(x.im < 0.0)
      quad = 7;
  }
  if(x.re == 0.0){
    if(x.im > 0.0)
      quad = 2;
    if(x.im == 0.0)
      quad = 9;
    if(x.im < 0.0)
      quad = 6;
  }
  if(x.re < 0.0){
    if(x.im > 0.0)
      quad = 3;
    if(x.im == 0.0)
      quad = 4;
    if(x.im < 0.0)
      quad = 5;
  }
  switch(quad){
  case 1:
    theta = atan(x.im/x.re);
    break;
  case 2:
    theta = pi/2.0;
    break;
  case 3:
    theta = (atan(x.im/x.re) + pi);
    break;
  case 4:
    theta = pi;
    break;
  case 5:
    theta = (atan(x.im/x.re) + pi);
    break;
  case 6:
    theta = 3*pi/2.0;
    break;
  case 7:
    theta = 2*pi + atan(x.im/x.re);
    break;
  case 8:
    theta = 0.0;
    break;
  case 9:
    theta = 0.0;
    break;
  default:
    theta = atan(x.im/x.re);
    fprintf(output_F,"ERROR IN Cx_DSqrt: Can't find quadrant. ");
    fprintf(output_F,"Using default value.\n");
    break;
  }
  theta = 180.0/pi * theta;
  return(theta);
}
