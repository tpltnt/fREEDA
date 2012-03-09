// This module contains the top level dispatcher for network analysis.
// This routines will be replaced later by the Simulator class object
// with enhaced functionality.

#include <fstream>
#include <cassert>

extern "C"
{
#include <stdio.h>
#include <string.h>
#include "../compat/stuff.h"
#include "../inout/report.h"
#include "../inout/variables.h"
#include "../inout/environment.h"

// Defined in an_output.c
int An_DoOutput(void);
int An_DoEndOutput(void);
}

#include "Analysis.h"

// Change this definition if a new analysis type is added
#define ANALYSIS_TYPES 6
#include "SVTr.h"
#include "AC.h"
#include "DC.h"
#include "SVHB.h"
#include "SVTran.h"
#include "SVTran2.h"

// By now keep a global pointer with the analysis to perform.
Analysis *analysis = NULL;

void anMain(Circuit* cir)
{
  report(MESSAGE, "*** Starting analysis ...\n");
  // Now do analysis using old-transim style by now.
  do
  {
    sepLine();
    fflush(output_F);
    v_evaluateExpressions();

    if (analysis)
    {
      analysis->run(cir);
    }
    else
    {
      report(WARNING, "No recognized analysis specified.\n");
      break;
    }

    An_DoOutput();
    minorCleanUp();
  }
  while (v_setupNextSweep());
  An_DoEndOutput();

  // Delete analysis which is no longer used
  if (analysis)
    delete analysis;
}

// This routine will be replaced later by the Simulator class object.
Analysis* createAnalysis(const string& name)
{
  // Only one analysis at a time is supported by now.
  assert(!analysis);

  // Add cases for each analysis type.
  if (name == "svtr")
    analysis = new SVTr();
  else if (name == "ac")
    analysis = new AC();
  else if (name == "dc")
    analysis = new DC();
  else if (name == "svhb")
    analysis = new SVHB();
  else if (name == "svtran")
    analysis = new SVTran();
  else if (name == "svtran2")
    analysis = new SVTran2();
  else 
  {
    cerr << "Fatal: " << name << ": analysis type unknown." << endl;
    exit(1);
  }
  return analysis;
}

// The following routine generates the analysis catalog file. By now, everything
// in this routine has to be modified by hand.
void printACatalog()
{
  Analysis* a_vector[ANALYSIS_TYPES];
  char htmlFilename[FILENAME_MAX];
  char analysisName[FILENAME_MAX];

  a_vector[0] = new AC();
  a_vector[1] = new DC();
  a_vector[2] = new SVTr();
  a_vector[3] = new SVHB();
  a_vector[4] = new SVTran();
  a_vector[5] = new SVTran2();

  // open output file
  // Use documentation environment variable
  int length = strlen(env_freeda_documentation) + 18;
  if(length >= FILENAME_MAX)
  {
    fprintf(stderr," ** Filename too long.\n");
    return;
  }
  htmlFilename[0] = 0;
  strcpy(htmlFilename,env_freeda_documentation);
  strcat(htmlFilename,"/");
  strcat(htmlFilename,"fr_analysis.html");
  FILE *out_f = fopen(htmlFilename, "w");
  if (out_f == NULL)
  {
    fprintf(stderr, "Could not open %s for writing.\n", htmlFilename);
    return;
  }
  else
  {
    printf("** Creating analysis catalog in %s\n",htmlFilename);
  }

  fprintf(out_f,
    "<HTML> \n <HEAD><TITLE>fREEDA(TM) Analysis Catalog</TITLE></HEAD> \n");
  fprintf(out_f, "<HR SIZE=4>\n<h3 ALIGN=CENTER>");
  fprintf(out_f, "fREEDA(TM) Analysis Catalog</h3>\n<HR SIZE=4>\n");

  for (int i = 0; i < ANALYSIS_TYPES; i++)
  {
    Analysis* an = a_vector[i];
    strcpy(analysisName,an->getName().c_str());

    fprintf(out_f, "<H3>%s</H3>\n", an->getDescription().c_str());
    fprintf(out_f, "Author(s): %s<br>\n", an->getAuthor().c_str());
    fprintf(out_f, "Analysis engine: %s</B><br>", an->getName().c_str());
    fprintf(out_f, "<B>Usage:</B> <br>\n");


    //If a new analysis type is added the following must be addressed
    // Changing the following requires that tr2_parse.y be studied.
    if(!strcmp("ac",analysisName))
    {
      fprintf(out_f, "<B>.ac</B>");
      fprintf(out_f, " &lt;parameter list&gt; </ul><br>\n");
    }
    else if(!strcmp("dc",analysisName))
    {
      fprintf(out_f, "<B>.dc</B>");
      fprintf(out_f, " &lt;parameter list&gt; </ul><br>\n");
    }
    else if(!strcmp("svtr",analysisName))
    {
      fprintf(out_f, "<B>.svtr</B>");
      fprintf(out_f, " &lt;parameter list&gt; </ul><br>\n");
    }
    else if(!strcmp("svhb",analysisName))
    {
      fprintf(out_f, "<B>.svhb</B>");
      fprintf(out_f, " &lt;parameter list&gt; </ul><br>\n");
    }
    else if(!strcmp("SVTran",analysisName))
    {
      fprintf(out_f, "<B>.tran</B>");
      fprintf(out_f, " &lt;parameter list&gt; </ul><br>\n");
    }
    else if(!strcmp("SVTran2",analysisName))
    {
      fprintf(out_f, "<B>.tran2</B>");
      fprintf(out_f, " &lt;parameter list&gt; </ul><br>\n");
    }
    else // this is the best we can do, may not be right
    {
      fprintf(out_f, "<B>.*</B>",analysisName);
      fprintf(out_f, " &lt;parameter list&gt; </ul><br>\n");
    }

    if (an->getNumberOfParams())
    {
      fprintf(out_f, "<br><table border cellpadding=5 cellspacing=0> \n\n");
      fprintf(out_f, "<tr> <th>Parameter</th> <th>Type</th>");
      fprintf(out_f, "<th>Default value</th> <th>Required?</th> </tr> \n");

      for (unsigned j = 0; j < an->getNumberOfParams(); j++)
      {
	       string parname, comment, type, dflt_val, required;
	       an->getParamDesc(j, parname, comment, type, dflt_val, required);
	       fprintf(out_f, "<tr><td><B>%s:</B> %s</td>\n",
		     parname.c_str(), comment.c_str());
	       fprintf(out_f, "<td> %s</td>\n", type.c_str());
	       fprintf(out_f, "<td> %s</td>\n", dflt_val.c_str());
	       fprintf(out_f, "<td> %s</td></tr>\n\n", required.c_str());
      }
      fprintf(out_f, "</table>\n");
    }
    else
      fprintf(out_f, "</ul><br><br>\n");

    fprintf(out_f, "<br><HR size=2>\n");
  }
  fprintf(out_f, "</BODY>\n</HTML>\n");
  fclose(out_f);

  {
    char s[FILENAME_MAX];
    printf("** htmlFilename = %s\n",htmlFilename);
    length = strlen(env_freeda_browser) +1 + strlen(htmlFilename) +1;
    if(length >= FILENAME_MAX)
    {
      printf(" ** Filename too long.\n");
      return;
    }
    s[0] = 0;
    strcpy(s,env_freeda_browser);
    strcat(s," ");
    strcat(s,htmlFilename);
    system(s);
  }

  // Free the dummy analysis
  for (int i=0; i<ANALYSIS_TYPES; i++)
    delete (a_vector[i]);
}

