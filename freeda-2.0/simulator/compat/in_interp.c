#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "standard.h"
#include "in_interp.h"

/*
 *
 * In_Interpolation.c
 *                                                                        *
 * Collection of Interpolation Routines
 *
 *
 * Contents:
 *
 * Routine     Description
 *
 * In_Ipoltn
 * In_Cspline
 * In_Csplint
 *
 *
 */

int In_Ipoltn(float f, float h[], float z[], float g[],
	      float *gre, float *gim, float *zre, float *zim)  
{
  int  status;
  static  float zret[74], gret[74]; 
  static int first_time = TRUE;
  /* Until the cubic spline coefficients are stored we need to do a
   * first_time trip always. */
  first_time = TRUE;
  if (f<1.0 || f>5000000.0)
    {
      status = -9;
      printf("error in input no interpolation.\n");
    }
  else
    { 
      if (first_time)
	{
	  In_Cspline (h,z, zret);
	  In_Cspline (h,g, gret);
	}
      else
	{
	}
      In_Csplint (h,z, zret,f, zre,zim); 
      In_Csplint (h,g, gret,f, gre,gim);
      status = 0;
      first_time = FALSE;
      return(status);
    }
  return 0;
}

void In_Cspline(x, y,y2)
     float x[], y[], y2[];
{
  int k,i,dum;
  float u[37], p, qn, un, sig,first,last;
  for(dum = 0;dum<=1;dum++)
    {
      first = (y[dum+2]-y[dum])/(x[1]-x[0]);
      last = (y[72 + dum]-y[70 + dum])/(x[36]-x[35]);
      if (first>0.99e30) 
	{
	  y2[dum] = 0.0;
	  u[0] = 0.0;
	}
      else
	{
	  y2[dum] = -0.5;
	  u[0] = (3.0/(x[1]-x[0]))*((y[dum+2]-y[dum])/(x[1]-x[0])-first);
	}

      for  (i = 1; i<=35; i++ )
	{
	  sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
	  p = sig*y2[2*(i-1) + dum]+2.0;
	  y2[2*i + dum ] = (sig-1.0)/p;
	  u[i] = (6.0*((y[2*(i+1)+ dum]-y[2*i+ dum])/(x[i+1]-x[i])-
		       (y[2*i + dum]-y[2*(i-1) + dum])/(x[i]-x[i-1]))
		  /(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}

      if (last>0.99e30) 
	{
	  qn = 0.0;
	  un = 0.0;
	}
      else
	{
	  qn = 0.5;
	  un = (3.0/(x[36]-x[35]))*(last-(y[72 + dum]-y[70 + dum])/(x[36]-x[35]));
	}

      y2[72+dum] = (un-qn*u[35])/(qn*y2[70+dum]+1.0);

      for ( k=35; k>=0; k--)
	{
	  y2[2*k + dum] = y2[2*k +dum]*y2[2*(k+1) + dum]+u[k];
	}
    }
  return;
}


/*  Subroutine for searching an ordered table*/
void In_Csplint(float x[37], float y[74], float y2[74], float xa,
		float *yap1, float *yap2)       
{
  float a, b, h1,ya;
  int k, klo, khi;

  klo=0;
  khi=36;

  while(khi-klo > 1) 
    {
      k=(khi+klo) >> 1;
      if (x[k] > xa) 
	khi=k;
      else  
	klo=k;
    }                      

  h1=x[khi]-x[klo];
  if (fabs(h1) <= 1.0e-15) 
    {
      fprintf(stderr, "bad x value\n");
      exit(0);
    }

  a = (x[khi]-xa)/h1;
  b = (xa-x[klo])/h1;
  ya = a*y[2*klo ]+b*y[2*khi]+((a*a*a-a)*y2[2*klo]
			       +(b*b*b-b)*y2[2*khi])*(h1*h1)/6.0;
  *yap1 = ya;
  ya = a*y[2*klo + 1]+b*y[2*khi + 1]+((a*a*a-a)*y2[2*klo
						  +1]+(b*b*b-b)*y2[2*khi + 1])*(h1*h1)/6.0;
  *yap2 = ya;
  return;

}


/*
 * floatHunt
 *
 * Smart table look-up routine.
 *
 * From numerical recipes, modified to work from zero indexed xx 
 */
void floatHunt(float *xx, int n, float x, int *jlo)
{
  int jm,jhi,inc,ascnd;

  ascnd=(xx[n-1] > xx[0]);
  if (*jlo <= 0 || *jlo > n) {
    *jlo=0;
    jhi=n+1;
  } else {
    inc=1;
    if (x >= xx[*jlo-1] == ascnd) {
      if (*jlo == n) return;
      jhi=(*jlo)+1;
      while (x >= xx[jhi-1] == ascnd) {
	*jlo=jhi;
	inc += inc;
	jhi=(*jlo)+inc;
	if (jhi > n) {
	  jhi=n+1;
	  break;
	}
      }
    } else {
      if (*jlo == 1) {
	*jlo=0;
	return;
      }
      jhi=(*jlo);
      *jlo -= 1;
      while (x < xx[*jlo-1] == ascnd) {
	jhi=(*jlo);
	inc += inc;
	*jlo=jhi-inc;
	if (*jlo < 1) {
	  *jlo=0;
	  break;
	}
      }
    }
  }
  while (jhi-(*jlo) != 1) {
    jm=(jhi+(*jlo)) >> 1;
    if (x > xx[jm-1] == ascnd)
      *jlo=jm;
    else
      jhi=jm;
  }
}




/*
 * doubleHunt
 *
 * Smart table look-up routine.
 *
 * From numerical recipes, modified to work from zero indexed xx 
 */
void doubleHunt(double *xx, int n, double x, int *jlo)
{
  int jm,jhi,inc,ascnd;

  ascnd=(xx[n-1] > xx[0]);
  if (*jlo <= 0 || *jlo > n) {
    *jlo=0;
    jhi=n+1;
  } else {
    inc=1;
    if (x >= xx[*jlo-1] == ascnd) {
      if (*jlo == n) return;
      jhi=(*jlo)+1;
      while (x >= xx[jhi-1] == ascnd) {
	*jlo=jhi;
	inc += inc;
	jhi=(*jlo)+inc;
	if (jhi > n) {
	  jhi=n+1;
	  break;
	}
      }
    } else {
      if (*jlo == 1) {
	*jlo=0;
	return;
      }
      jhi=(*jlo);
      *jlo -= 1;
      while (x < xx[*jlo-1] == ascnd) {
	jhi=(*jlo);
	inc += inc;
	*jlo=jhi-inc;
	if (*jlo < 1) {
	  *jlo=0;
	  break;
	}
      }
    }
  }
  while (jhi-(*jlo) != 1) {
    jm=(jhi+(*jlo)) >> 1;
    if (x > xx[jm-1] == ascnd)
      *jlo=jm;
    else
      jhi=jm;
  }
}

