/*************************************************************************

  These Global variables are needed by the output routines. They
  must be set by the analysis code.

  USE THE PROVIDED FUNCTIONS TO ALLOCATE OR FREE THE VECTORS!

  Functions and variables defined in init_clean.cc

**************************************************************************/

extern int NoFreqPoints;
extern int NoTimePoints;

extern double* FreqV_P;
extern double* TimeV_P;


/* Allocate memory for FreqV_P and set NoFreqPoints */
double* allocFreqV_P(int nfp);

/* Free memory */
void freeFreqV_P();

/* Allocate memory for TimeV_P and set NoTimePoints */
double* allocTimeV_P(int ntp);

/* Free memory */
void freeTimeV_P();

