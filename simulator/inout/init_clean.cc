/*
 *  This file contains initialization code used directly by freeda.cc
 */

extern "C"
{
	#include <stdio.h>
	#include <string.h>
	#include "parser.h"
	#include "ftvec.h"
	#include "../compat/ml.h"
}

#include <cassert>
#include <cstdlib>
#include <iostream>
#include "../help.h"
#include "../network/Circuit.h"
#include "../network/ADInterface.h"
#include "../elements/element_headers.h"
#include "../network/ElementManager.h"
#include "environment.h"

/* Symbol table pointers */
st_Table_Pt capOptionsT_P = NULL,
	    capModelT_P = NULL,
	    capOutputT_P = NULL,
	    capEndOutputT_P = NULL,
	    capOutputVariablesT_P = NULL;

extern char *freeda_version;
extern void str2lower(char *);
void printACatalog(); // Defined in anMain.cc

// Global pointers to environment variables (in environment.h)
char *env_freeda_home;
char *env_freeda_library;
char *env_freeda_projects;
char *env_freeda_path;
char *env_freeda_bin;
char *env_freeda_simulator;
char *env_freeda_elements;
char *env_freeda_documentation;
char *env_freeda_web_documentation;
char *env_freeda_browser;

char spiceTitle[SPICE_TITLE_LEN]; // Title, taken as first line of input
char inputFilename[MAX_TOK_LEN];
char outputFilename[MAX_TOK_LEN];
char plotFilename[MAX_TOK_LEN];
/* Input/output file pointers */
FILE* input_F;
FILE* output_F;
FILE* plotlist_F;
/* Other global variables */
spice_t spiceTable;


// Frequency and time stuff
int NoFreqPoints = 0;
int NoTimePoints = 0;
doublev_t FreqV_P = NULL;
doublev_t TimeV_P = NULL;

void Init()
{

  /* Initialize Tables */
  capOptionsT_P = St_NewTable("OPTIONS", HASH_SIZE);
  capModelT_P   = St_NewTable("MODEL", HASH_SIZE);
  capOutputT_P  = St_NewTable("OUTPUT", HASH_SIZE);
  capOutputVariablesT_P  = St_NewTable("OUTPUT_VARIABLES", HASH_SIZE);
  capEndOutputT_P  = St_NewTable("END_OUTPUT", HASH_SIZE);

  setUpEnvironmentVariables();

  /* Initialize Globals */
  input_F = fopen(inputFilename, "r");
  if(!input_F) {
    fprintf(stderr,"Error: couldn't open inputfile.\n");
    exit(1);
  }

  /* Set up output file */
  output_F = fopen(outputFilename, "w");
  if(!output_F) {
    fprintf(stderr, "Error: couldn't open output file.\n");
    exit(1);
  }

  /* Open plot list output file */
  strcpy(plotFilename, outputFilename);
  strcat(plotFilename,".plot");
  plotlist_F = fopen(plotFilename, "w");
  if(!plotlist_F) {
    fprintf(stderr, "Error: couldn't open plot list output file.\n");
    exit(1);
  }

  SpiceDefault();      /* Now set spice defaults */
  InitScan(input_F);   /* Init Scan */


  // Print starting time
  time_t curtime;
  time(&curtime);
  char str[90];
  strftime(str, 80, "%c", localtime(&curtime));
  printf ("\n**********  fREEDA %s running on %s  **********\n",
		freeda_version, str);
  fprintf(output_F,
		"\n**********  fREEDA %s running on %s  **********\n\n",
		freeda_version, str);

  /* Output environment variables */
  fprintf(output_F, "**  Environment variables: **\n");
  fprintf(output_F,"FREEDA_HOME = %s\n",env_freeda_home);
  fprintf(output_F,"FREEDA_LIBRARY = %s\n",env_freeda_library);
  fprintf(output_F,"FREEDA_PROJECTS = %s\n",env_freeda_projects);
  fprintf(output_F,"FREEDA_PATH = %s\n",env_freeda_path);
  fprintf(output_F,"FREEDA_BIN = %s\n",env_freeda_bin);
  fprintf(output_F,"FREEDA_SIMULATOR = %s\n",env_freeda_simulator);
  fprintf(output_F,"FREEDA_ELEMENTS = %s\n",env_freeda_elements);
  fprintf(output_F,"FREEDA_DOCUMENTATION = %s\n",env_freeda_documentation);
  fprintf(output_F,"FREEDA_WEB_DOCUMENTATION = %s\n",
    env_freeda_web_documentation);
  fprintf(output_F,"FREEDA_BROWSER = %s\n",env_freeda_browser);
  fprintf(output_F,"\n");
}



// Allocate memory for FreqV_P and set NoFreqPoints
double* allocFreqV_P(int nfp)
{
  assert(!NoFreqPoints);
  NoFreqPoints = nfp;
  return FreqV_P = Mlib_DNewVec(NoFreqPoints);
}

// Free memory
void freeFreqV_P()
{
  if (NoFreqPoints) {
    NoFreqPoints = 0;
    Mlib_DFreeVec(FreqV_P);
  }
}

// Allocate memory for TimeV_P and set NoTimePoints
double* allocTimeV_P(int ntp)
{
  assert(!NoTimePoints);
  NoTimePoints = ntp;
  return TimeV_P = Mlib_DNewVec(NoTimePoints);
}

// Free memory
void freeTimeV_P()
{
  if (NoTimePoints) {
    NoTimePoints = 0;
    Mlib_DFreeVec(TimeV_P);
  }
}


void minorCleanUp()
{
  /* Free frequency and time vectors if they where allocated */
  freeFreqV_P();
  freeTimeV_P();
}


void majorCleanUp()
{
  /* Deletes and frees up tables which are initialized by Init()    */
  St_DelTable(capOptionsT_P); capOptionsT_P = NULL;
  St_DelTable(capModelT_P); capModelT_P = NULL;
  St_DelTable(capOutputT_P); capOutputT_P = NULL;
  St_DelTable(capOutputVariablesT_P); capOutputVariablesT_P = NULL;
  St_DelTable(capEndOutputT_P); capEndOutputT_P = NULL;

  /* Free equations in these tables: SUBPARAMSymT_P,PARAMSymT_P     */
  /* MEMTODO : v_freeExpressions() and v_freeSweeps() do not appear */
  /* MEMTODO : to free.. needs further checking                     */
  v_freeExpressions();
  v_freeSweeps();

  /* Free equations in these tables: SUBPARAMSymT_P,PARAMSymT_P     */
  /* MEMTODO : free_equation() does not appear to free.. needs      */
  /* MEMTODO : further checking 					  */
  /* free_equation(); */

  /* Cleanup ml.c */
  ml_CleanUp();

  /* Clean up buffer in parser which is allocated in lex.yy.c */
  PCleanUpFLEX();

  // Close files
  fclose(input_F);
  fclose(output_F);
  fclose(plotlist_F);

}


/*
 * Spice Default
 *
 * Set up spice defaults
 */
void SpiceDefault()
{
  generic_t gval;
  gval.d = OPTIONS_DEFAULT_DEFL ;
  St_DefSym(capOptionsT_P, "defl"	, GEN_TYPE_DOUBLE	, gval);
  gval.d = OPTIONS_DEFAULT_DEFW ;
  St_DefSym(capOptionsT_P, "defw"	, GEN_TYPE_DOUBLE	, gval);
  gval.d = OPTIONS_DEFAULT_DEFAD ;
  St_DefSym(capOptionsT_P, "defad"	, GEN_TYPE_DOUBLE	, gval);
  gval.d = OPTIONS_DEFAULT_DEFAS ;
  St_DefSym(capOptionsT_P, "defas"	, GEN_TYPE_DOUBLE	, gval);
  gval.d = OPTIONS_DEFAULT_TNOM ;
  St_DefSym(capOptionsT_P, "tnom"	, GEN_TYPE_DOUBLE	, gval);
  gval.i = OPTIONS_DEFAULT_NUMDGT ;
  St_DefSym(capOptionsT_P, "numdgt"	, GEN_TYPE_INT		, gval);
  gval.d = OPTIONS_DEFAULT_CPTIME ;
  St_DefSym(capOptionsT_P, "cptime"	, GEN_TYPE_DOUBLE	, gval);
  gval.d = OPTIONS_DEFAULT_LIMPTS ;
  St_DefSym(capOptionsT_P, "limpts"	, GEN_TYPE_DOUBLE	, gval);
  gval.i = OPTIONS_DEFAULT_ITL1 ;
  St_DefSym(capOptionsT_P, "itl1"	, GEN_TYPE_INT		, gval);
  gval.i = OPTIONS_DEFAULT_ITL2 ;
  St_DefSym(capOptionsT_P, "itl2"	, GEN_TYPE_INT		, gval);
  gval.i = OPTIONS_DEFAULT_ITL4 ;
  St_DefSym(capOptionsT_P, "itl4"	, GEN_TYPE_INT		, gval);
  gval.i = OPTIONS_DEFAULT_ITL5 ;
  St_DefSym(capOptionsT_P, "itl5"	, GEN_TYPE_INT		, gval);
  gval.d = OPTIONS_DEFAULT_RELTOL ;
  St_DefSym(capOptionsT_P, "reltol"	, GEN_TYPE_DOUBLE	, gval);
  gval.d = OPTIONS_DEFAULT_TRTOL ;
  St_DefSym(capOptionsT_P, "trtol"	, GEN_TYPE_DOUBLE	, gval);
  gval.d = OPTIONS_DEFAULT_ABSTOL ;
  St_DefSym(capOptionsT_P, "abstol"	, GEN_TYPE_DOUBLE	, gval);
  gval.d = OPTIONS_DEFAULT_CHGTOL ;
  St_DefSym(capOptionsT_P, "chgtol"	, GEN_TYPE_DOUBLE	, gval);
  gval.d = OPTIONS_DEFAULT_VNTOL ;
  St_DefSym(capOptionsT_P, "vntol"	, GEN_TYPE_DOUBLE	, gval);
  gval.d = OPTIONS_DEFAULT_PIVREL ;
  St_DefSym(capOptionsT_P, "pivrel"	, GEN_TYPE_DOUBLE	, gval);
  gval.d = spiceTable.gmin = OPTIONS_DEFAULT_GMIN ;
  St_DefSym(capOptionsT_P, "gmin"	, GEN_TYPE_DOUBLE	, gval);

  spiceTable.echo = 1;

}


/*
 * Common code for setting environmental variables
 */
void setUpEnvironmentVariablesCommon(const char *environmentName, char **envVariable,
  const char *firstPart, const char *secondPart, const char *thirdPart)
  {
  char *s;
  if((s = getenv(environmentName)) == NULL) {
    *envVariable = (char *)
      malloc(strlen(firstPart) + strlen(secondPart) + strlen(thirdPart) + 1);
    strcpy(*envVariable,firstPart);
    strcat(*envVariable,secondPart);
    strcat(*envVariable,thirdPart);
    // printf("%s\n",*envVariable);
    }
  else {
    *envVariable = (char *) malloc(strlen(s) + 1);
    strcpy(*envVariable,s);
    }
  }


/*
 * Get environmental variables
 */
int setUpEnvironmentVariables()
  {
  char *s;

  s = getenv("USER");
  setUpEnvironmentVariablesCommon("FREEDA_HOME",
    &env_freeda_home, "/", s, "/freeda");
  setUpEnvironmentVariablesCommon("FREEDA_LIBRARY",
    &env_freeda_library,env_freeda_home,"/library","");
  setUpEnvironmentVariablesCommon("FREEDA_PROJECTS",
    &env_freeda_projects,env_freeda_home,"/projects","");
  setUpEnvironmentVariablesCommon("FREEDA_PATH",
    &env_freeda_path,env_freeda_home,"/freeda","");
  setUpEnvironmentVariablesCommon("FREEDA_BIN",
    &env_freeda_bin,env_freeda_path,"/bin","");
  setUpEnvironmentVariablesCommon("FREEDA_SIMULATOR",
    &env_freeda_simulator,env_freeda_path,"/simulator","");
  setUpEnvironmentVariablesCommon("FREEDA_ELEMENTS",
    &env_freeda_elements,env_freeda_simulator,"/elements","");
  setUpEnvironmentVariablesCommon("FREEDA_DOCUMENTATION",
    &env_freeda_documentation,"/tmp","","");
  setUpEnvironmentVariablesCommon("FREEDA_WEB_DOCUMENTATION",
    &env_freeda_web_documentation,"http://www.freeda.org/doc","","");


  // Set up environment to launch browser
#ifdef __CYGWIN__
  setUpEnvironmentVariablesCommon("FREEDA_BROWSER",
      &env_freeda_browser,"cygstart","","");
#else
  // We are dealing with unix, assume firefox as browser
  setUpEnvironmentVariablesCommon("FREEDA_BROWSER",
      &env_freeda_browser,"firefox","","");
#endif

  }

/*
 * Do a simple sanity check.
 * Of most use for developers.
 * Not much here for now.
 */
void sanityCheck()
{
  int sane = 1;

  printf("Conducting sanity check\n");

  // Go through all of the elements.
  {
    char c;
    char elementName[MAX_STRING_LEN];
    printf("Checking elements ");
    // Vector to hold the original object for each element type.
    Element** elem_vector = NULL;
    // Allocate memory for the vector
    elem_vector = new Element*[ELEM_TYPES];
    // Allocate one object each element type
    // The code for this is automatically generated
    #include "../elements/create_dummy_elem.cc"

    for (int i = 0; i < ELEM_TYPES; i++)
    {
      printf(".");
      Element* elem = elem_vector[i];
      strcpy(elementName,elem->getName().c_str());

      // Check that the element have lower case names.
      {
        int j;
        for (j=0; j<strlen(elementName); j++)
        {
          c = elementName[j];
          tolower(c);
          if(c != elementName[j])
          {
            printf("\nelement %s must be lower case.\n",elementName);
            sane = 0;
            break;
          }
        }
      }
    }
    delete [] elem_vector;
    printf("\n");
  }

  // Final report/
  if(sane)
    printf("\nSimple sanity found no development errors.\n");
  else
    printf("\nSanity error found, generally this is fatal.\n");

}


/*
 * Process command line arguments
 */
void processCommandLine(int argc, char **argv)
  {
  char s[MAX_TOK_LEN];
  char* dot;

// Element Catalog
// -c		produces a catalog listing of elements
// -c element	produces the documentation for the element "element"
// -a		produces analysis documentation
// -s		do a simple sanity check
  if(argc >= 2 && (!strcmp(argv[1], "-c") || !strcmp(argv[1], "--catalog")))
    {
    setUpEnvironmentVariables();
    if(argc == 3) // So the name of an element is specified.
      {
      pCatalogElement(argv[2]);
      }
    else // Produce a listing of the elements.
      {
      pCatalog();
      }
    exit(0);
    }

// Analysis Catalog
// -a		produces analysis documentation
  if(argc >= 2 && (!strcmp(argv[1], "-a") || !strcmp(argv[1], "--analysis")))
    {
    setUpEnvironmentVariables();
    printACatalog();
    exit(0);
    }

  else if (argc == 2 &&	(!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")))
   {
   pHelpMessage();
   exit(0);
   }
  else if (argc == 2 &&	(!strcmp(argv[1], "-v") || !strcmp(argv[1],
		"--version")))
   {
   pVersion();
   exit(0);
   }
  else if (argc == 2 &&	(!strcmp(argv[1], "-s") || !strcmp(argv[1],
		"--sanity")))
   {
   sanityCheck();
   exit(0);
   }
  else if (argc == 2 &&	(!strcmp(argv[1], "-l") || !strcmp(argv[1],
		"--licence")))
   {
   pLicence();
   exit(0);
   }

  // If none of the previous match, get input/output file names and run
  // the netlist.
  // get input file name
  if(argc > 1)
    strcpy(inputFilename, argv[1]);
  else
    {
    printf("Input file: ");
    scanf("%s", inputFilename);
    }
  // get output file name
  if(argc > 2)
    {
    strcpy(outputFilename, argv[2]); // Use second argument as name
    }
  else
    {
    strcpy(s,inputFilename);
    dot = strstr(s, ".net");
    if(dot)
      dot[0] = 0;
    strcpy(outputFilename, s); // Append .out to input file
    strcat(outputFilename,".out");
    }
  }
