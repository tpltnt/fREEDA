/****************************************************************************
 * 
 * This header provides all the definitions required by the parser
 * and the parser functions.
 *
 ****************************************************************************/


#include "../compat/stuff.h"
#include "spice.h"
#include "oper.h"
#include "variables.h"
#include "parserGlobals.h"

#define SPICE_TITLE_LEN 80        	/* max length of title          */
#define MAX_TOK_LEN MAX_STRING_LEN	/* maximum acceptable token or	*/
					/* string const. length		*/
#define TERM_UNDEFINED -1
#define GR_MAX_CONN	2000		/* Maximum connectivity.        */

#define YYMAXDEPTH 40
#define MAX_SPICE_COMMAND_TOKENS 100
#define INPUT_BUF_SIZE 2048
#define ST_COULDNT_COERCE -3
#define ST_SYM_READ_FMT_ERR -2
#define BUFSIZE1 1024
/* #define TOKSIZE 256 */
#define EOFTOKEN -1

# define PRINTMAX 32			/* Max number of arguments 	*/
					/* in a .print line		*/

/********************
 * Global variables *
 ********************/

extern FILE *input_F;			/* file pointer for input */
extern FILE *output_F;			/* file pointer for output */

/*
 * Globals used by the tokenizer (lexical analyzer)
 */
extern char *lastNumberAsText;
extern int tokType;		/* type of current token */
extern char buf[];
extern FILE *yyin;
extern char   spiceTitle[];	/* Spice title */


/*
 * Globals used by the parser
 */

extern char thisId[],		/* current identifier name */
	thisElementType[],	/* current element type	   */
	thisElementName[],	/* current element name    */
	thisModelName[],		/* current model name	   */
	thisVariableName[],	/* current variable name   */
        thisModelType[],

	*thisOutputName,
	*sweepId;				/* pointer to string */

extern sweep_Pt sweepList_head; 

extern int	usageType,thisElementLevel;			/* Identifies type of line for
						   expressions and sweeps */
extern int	usageVal;			/* Identifier value */

extern char	*usageTableName;		/* Identifies table used */
extern char    thisParent[];   /* current parent */
extern char    thisSubCkt[];   /* current subcircuit name */
extern char    *lastText;
extern char    *printHeader[32];
extern int     countV;
extern int   anType;
extern int   spSyn;
extern int isCurrent, isVoltage;

extern generic_t globalGval;
extern sweep_Pt  sweepList_head;
extern int	  globalType;
extern double	  sweepInitial_temp, sweepFinal_temp, sweepStep_temp;

extern gr_Id_t thisTid;		/* current terminal id		    */
extern gr_Id_t thisNid;		/* current node id		    */

extern gr_Id_t thisTermList[];  /* current list of terminals	    */
extern capOutReq_t *thisOutReqList;	    /* current output request list	    */

extern int thisListCount,		    /* current # of terminals	    */
  nameIndex,			    /* current # of names in C group    */
  thisLevel;			    /* layer number of current element  */

extern double doubleVector[]; /* Used as temporary storage
						   * for vectors in input */
extern char *stringVector_P[]; 
                                            /* Used as temporary storage for
					     * pointers to strings in an input
					     * string vector*/
extern int vectorType;


typedef struct {
  char s[MAX_NAME_LEN +1];
} parameter_string_t;

extern parameter_string_t   parameterList[];
extern int thisParameter;

extern int noVectorElements;    	/* number of vector elements entered
					 * so far in doubleVector */

extern st_Table_t *thisSymT_P;		/* temporary symbol table pointer
					 * changes value from time to time */

extern char inputBuf[];
					/* input line buffer		    */
extern char *inputBufp;         	/* position in input line buffer    */
extern int lineNo;			/* number of current line	    */
extern char inputFilename[];	        /* name of input file		    */

extern int maimed;			/* TRUE if fatal errors in netlist  */

extern char tc;		        	/* tokenizer character		    */
extern int sequence;		        /* sequence index used to construct
					 * unique names by concatenating with
					 * nonunique string.		    */

/*
 * Variables and arrays to store a single POLY specification before it is
 * digested by an element.
 */
extern gr_Id_t polyList[];	        /* polynomial list of terminals or
					 * edges */
extern int 	listType;		/* type of element in list
					 *  =  GR_ID_TERM_TYPE if terminal id's
					 *  =  GR_ID_Edge_TYPE if edge id's
					 */
extern double poly_coeff[];
extern int number_poly_coeffs,		/* max. is ST_MAX_VECTOR_LENGTH      */
	polydimension;			/* The maximum dimension is
					 * GR_MAX_CONN / 2          	     */
extern int noElementsInList;		/* When elements are listed as in
					 * a .COUPLE statement or in an
					 * F element card this is the number of
					 * elements in the list		     */
					/* When polynomial information is
					 * stored this is the number of edges
					 * or terminals in list.	     */     


/***************************
 * Function prototypes 
 ***************************/

/* defined in misc.c */
void SpiceDefault();
void majorCleanUp();
void minorCleanUp();

/* defined in parser_functions.cc */
int PCleanUpFLEX();
int yyparse();
int yylex();
int yyerror(char* s);
void DumpTable(FILE *output_F, st_Table_Pt st_P);
void InitScan(FILE *f);

void InitOutReqList();
void ClearOutReqList();
void AppendOutReqList(capOutReq_t newOutReq);
void AppendOutReqString(char *inputString, int outReqType);
void AppendOutReqOperator(char *operatorName);

void AddParameterUsage(const char *variableName, gr_Id_t id,
		       const char *parameterName);
void AddSymbolTableUsage(const char *variableName, char *tableName,
			 const char *parameterName);
void AddAnalysisUsage(const char *variableName, char *tableName,
			 const char *parameterName);
void AddVariableUsage(char *variableName, char *prevVariableName);

void initElement(char *elementName, char *elementType);
void checkLevel();
void checkModel();
void closeElement();
int AddSubCkt();
void SaveOutReqList(st_Table_Pt table_P, char *name);
char *MakeSequence(char *name);
void SymAssgnInt(st_Table_t *t_P, char *name, int i);
void SymAssgnDouble(st_Table_t *t_P, char *name, double d);
void SymAssgnString(st_Table_t *t_P, char *name, char *s);

void ParseWarning(const char *message);
int ParseError(const char *message);

gr_Id_t getTermNamed(char *name);
gr_Id_t getElemNamed(char *name);
gr_Id_t newTerm(char* name);
gr_Id_t newElem(char* elemtype, char* name);
int popCircuit();
int pushCircuit(char* name);
void setupCircuits();
void setReference(gr_Id_t tid);

void setTermX(gr_Id_t id, double x);
void setTermY(gr_Id_t id, double y);

void replaceElemParam(gr_Id_t id, int type, generic_t gval, char* pname);
void initAnalysis(char* analysisName);
void closeAnalysis();
void replaceAnalysisParam(char* atype, generic_t gval, char* pname);
void insertPolyParam(char* elementName) ;
