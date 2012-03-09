/*
 * This file is the interface with gnuplot.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>


/* System-dependent command line */
#ifdef __unix__
#define GNUPLOT_COMMAND "xterm -e gnuplot %s &"
#else
/* This is known to work on windows NT using egcs */
#define GNUPLOT_COMMAND "gnuplot %s &"
#endif


void gnuplot(char *plotfile)
{ 
  char buff[1024];

  /* returns the  termination  status  of  the  command */
  sprintf(buff, GNUPLOT_COMMAND, plotfile);
  system(buff);
}




