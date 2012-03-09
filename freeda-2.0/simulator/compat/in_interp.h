/*
 * Header file for the interpolation routines in in_interp.c
 */
int In_Ipoltn(float f, float h[], float z[], float g[],
	      float *gre, float *gim, float *zre, float *zim);
void In_Cspline(float x[], float y[], float y2[]);
void In_Csplint(float x[], float y[], float y2[], float xa,
		float *yap1, float *yap2);
void floatHunt(float *xx, int n, float x, int *jlo);

void doubleHunt(double *xx, int n, double x, int *jlo);


