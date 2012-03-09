/* 
 * Declarations from init_clean.cc
 */

extern st_Table_Pt
  capOptionsT_P,	/* symbol table for options */
  capModelT_P,		/* symbol table for models */
  capOutputT_P,		/* symbol table for output */
  capOutputVariablesT_P,/* symbol table for output variables*/
  capEndOutputT_P;	/* symbol table for final output */
 
extern char spiceTitle[];
extern char inputFilename[];
extern char outputFilename[];
extern char plotFilename[];

extern FILE *input_F;			/* file pointer for input */
extern FILE *output_F;			/* file pointer for output */
extern FILE *plotlist_F;                /* plot list file */

int setUpEnvironmentVariables();
void processCommandLine(int argc, char **argv);
void sanityCheck();

// Process command line and init input and output files.
void Init();

// minorCleanUp() -- Used to free memory between frequency sweeps inorder to
// 		     optimize the memory space
void minorCleanUp();

// majorCleanUp() -- Used to free memory after fREEDA has completed its
// 		     simulation.
void majorCleanUp();
