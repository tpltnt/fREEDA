/**************************************************************************
 *
 * This function is used to send output to both stdout and the output file.
 *
 **************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "../compat/stuff.h"
#include "report.h"

void report(enum ReportType type, const char *message)
{
  switch (type) {

  case LOGFILE:
    fprintf(output_F," %s\n", message);
    break;

  case MESSAGE:
    fprintf(stdout,"   %s\n",  message);
    fprintf(output_F,"   %s\n", message);
    break;

  case WARNING:
    fprintf(stderr,"\n   Warning: %s\n",  message);
    fprintf(output_F,"\n   Warning: %s\n", message);
    break;

  case FATAL:
  default:
    fprintf(stderr,"\n\n   *** Fatal: %s\n",  message);
    fprintf(output_F,"\n\n   *** Fatal: %s\n", message);
    exit(1);
  }
}

void sepLine()
{
    fprintf(output_F, 
	    "-----------------------------------------------------------------------------\n");
    fprintf(stdout, 
	    "-----------------------------------------------------------------------------\n");
}

void newLine()
{
  fprintf(output_F, "\n");
  fprintf(stdout, "\n");
}
