// The functions in this file are used by the parser to have access
// to the C++ network model from C.
// The code is mainly a hack of the old Transim code. It is a mix of C++
// and C, so do not expect elegance here.

#include <cstring>
#include <cctype>
#include <cstdlib>
#include <cstdio>
#include <cmath>
using namespace std;

extern "C"
{
	#include "parser.h"
	#include "report.h"
}

#include "../network/CircuitManager.h"
#include "../elements/x/Xsubckt/src/Xsubckt.h"
#include "../analysis/Analysis.h"

extern Analysis* analysis;

// Defined in anMain.cc
Analysis* createAnalysis(const string& name);

// Globals used by the tokenizer (lexical analyzer)
char yylvalStr[YYMAXDEPTH][MAX_TOK_LEN + 2],	/* text of yylval.id */
*lastNumberAsText;
char buf[BUFSIZE1];

// Globals used by the parser

int tokType;	                         // type of current token
char thisId[MAX_NAME_LEN + 2],		 // current identifier name
thisElementType[MAX_NAME_LEN + 2],	 // current element type
thisElementName[MAX_NAME_LEN + 2],	 // current element name
thisModelName[MAX_NAME_LEN + 2],	 // current model name
thisVariableName[MAX_NAME_LEN + 2], // current variable name
*thisOutputName,
thisModelType[MAX_NAME_LEN+2],
*sweepId;				 // pointer to string
int usageType,thisElementLevel;	         // Identifies type of line for expressions and sweeps
int usageVal;                            // Identifier value
char *usageTableName;		         // Identifies table used
char thisParent[MAX_NAME_LEN + 2];       // current parent
char thisSubCkt[MAX_NAME_LEN + 2];       // current subcircuit name
char *lastText;

char *printHeader[32];                   // labels for .print
int countV;
int anType, spSyn, isCurrent, isVoltage;

generic_t globalGval;

// First entry in sweep linked list.
sweep_Pt sweepList_head = 0;
int	  globalType;
double	  sweepInitial_temp, sweepFinal_temp, sweepStep_temp;

gr_Id_t thisTid;	// current terminal id
gr_Id_t thisNid;	// current node id

gr_Id_t thisTermList[GR_MAX_CONN];  // current list of terminals
capOutReq_t *thisOutReqList=NULL;  // current output request list

int thisTermCount, // current # of terminals
thisListCount, // current # of items in OPTIONS list
nameIndex, // current # of names in C group
thisLevel; // layer number of current element

double doubleVector[ST_MAX_VECTOR_LENGTH]; // Used as temporary storage for vectors in input
char *stringVector_P[ST_MAX_VECTOR_LENGTH]; // Used as temporary storage for
// pointers to strings in an input
// string vector
int vectorType = GEN_TYPE_NONE;

parameter_string_t   parameterList[MAX_SPICE_COMMAND_TOKENS];
int thisParameter;

int noVectorElements = 0; // number of vector elements entered so far in doubleVector

st_Table_t *thisSymT_P; // temporary symbol table pointer changes value from time to time

char inputBuf[INPUT_BUF_SIZE] = { 0 };
// input line buffer
char *inputBufp = inputBuf;	// position in input line buffer
int lineNo = 0;			// number of current line

char tc = ' ';	// tokenizer character
int sequence = 0; // sequence index used to construct unique names by concatenating with
// nonunique string

// Variables and arrays to store a single POLY specification before it is
// digested by an element.
gr_Id_t polyList[GR_MAX_CONN]; 		// nomial list of terminals or
// edges
int 	listType;			// type of element in list
//  =  GR_ID_TERM_TYPE if terminal id's
//  =  GR_ID_Edge_TYPE if edge id's
//
double poly_coeff[ST_MAX_VECTOR_LENGTH];
int	number_poly_coeffs,		// max. is ST_MAX_VECTOR_LENGTH
polydimension;			// The maximum dimension is
// GR_MAX_CONN / 2
int	noElementsInList;		// When elements are listed as in
// a .COUPLE statement or in an
// F element card this is the number of
// elements in the list
// When polynomial information is
// stored this is the number of edges
// or terminals in a list
// Some function prototypes
int Coerce(int sType, ParamType pType, generic_t *gval_P);
void EditParams(Element* elem, bool overwrite = true);
void EditItemParams(NetListItem* item, bool overwrite = true);
Element* ConstrPrim();

//-------------------------------------------------------------------------
// Function definitions begin
//-------------------------------------------------------------------------
void InitScan(FILE *f)
{
  yyin = f;
  thisSymT_P = NULL;
  maimed = FALSE;
  polydimension = 0;  /* zero initial polynomial dimension */
}

int yyerror(char* s)
{
  ParseError(s);
  return 0;
}

void SymAssgnInt(st_Table_t *t_P, char *name, int i)
{
  generic_t gval;

  gval.i = i;
  St_ReplSym(t_P, name, GEN_TYPE_INT, gval);
}

void SymAssgnDouble(st_Table_t *t_P, char *name, double d)
{
  generic_t gval;

  gval.d = d;
  St_ReplSym(t_P, name, GEN_TYPE_DOUBLE, gval);
}

void SymAssgnString(st_Table_t *t_P, char *name, char *s)
{
  generic_t gval;

  gval.s = (char *) malloc(strlen(s) + 1);
  strcpy(gval.s, s);
  St_ReplSym(t_P, name, GEN_TYPE_STRING, gval);
}

void ClearOutReqList()
{
  capOutReq_t *outReq = thisOutReqList,
	*tmp;

  while (outReq)
  {
    tmp = outReq;
    if (outReq->type == CAP_OBJ_DOUBLEV)
      free(outReq->obj.dv);
    else if (outReq->type == CAP_OBJ_DCXV)
      free(outReq->obj.dcxv);
    if (outReq->x)
      free(outReq->x);
    if (outReq->xName)
      free(outReq->xName);
    free(outReq);
    outReq = tmp->next;
  }

  thisOutReqList = NULL;
}

void InitOutReqList()
{
  if (thisOutReqList)
    ParseError("Internal error in InitOutReqList");
  thisOutReqList = (capOutReq_t *) calloc(1, sizeof(capOutReq_t));
  thisOutReqList->next = NULL;
}

void AppendOutReqList(capOutReq_t newOutReq)
{
  capOutReq_t *outReq;
  for (outReq = thisOutReqList; outReq->next; outReq = outReq->next);
  outReq->next = (capOutReq_t *) calloc(1, sizeof(capOutReq_t));
  *outReq->next = newOutReq;
  outReq->next->protect = TRUE;
  outReq->next->next = NULL;
}

// AppendOutReqOperator
// Similar to AppendOutReqList but specialized for operators
void AppendOutReqOperator(char *operatorName)
{
  int n;
  capOutReq_t outReq;
  outReq.type = CAP_OBJ_OPER;
  n = An_LookupOpName(operatorName);
  if(n < 0) ParseError("Unknown Operator");
  outReq.obj.opType = n;
  outReq.size = 0;
  outReq.x = NULL;
  AppendOutReqList(outReq);
}

// AppendOutReqString
// Append string to OutReqList
// Similar to AppendOutReqList but special treatment of strings.
//  1.  A copy of the string is made.
//  2.  Length of the string is determined
// Parameters
// outReqType	The type of the output request must be supplied.
// inputString	String to be placed in output request list (A permanent copy
//   		is made first.
void AppendOutReqString(char *inputString, int outReqType)
{
  capOutReq_t outReq;
  int length;
  char *s;
  length  = strlen(inputString);
  s = (char *) malloc(length+1);
  strcpy(s, inputString);
  outReq.type = outReqType;
  outReq.obj.s = s;
  outReq.size = length;
  outReq.x = NULL;
  AppendOutReqList(outReq);
}

void DumpOutReqList(capOutReq_Pt outReq)
{
  for (; outReq->next; outReq = outReq->next)
  {
    fprintf(output_F,"\t\ttype = %d\n",outReq->type);
    switch(outReq->type)
    {
      case CAP_OBJ_NONEXISTENT:
      fprintf(output_F,"\t termination\n");
      break;
      case CAP_OBJ_TERM:
      fprintf(output_F,"\t terminal id(val,type) = (%d,%d)\n",
        outReq->obj.tid.val,outReq->obj.tid.type);
      break;
      case CAP_OBJ_NODE:
      fprintf(output_F,"\t node id(val,type) = (%d,%d)\n",
        outReq->obj.tid.val,outReq->obj.nid.type);
      break;
      case CAP_OBJ_EDGE:
      fprintf(output_F,"\t edge id(val,type) = (%d,%d)\n",
        outReq->obj.tid.val,outReq->obj.nid.type);
      break;
      case CAP_OBJ_DATAFILE:
      fprintf(output_F, "\t datafile = %s\n",outReq->obj.s);
      break;
      case CAP_OBJ_FILENAME:
      fprintf(output_F, "\t filename = %s\n",outReq->obj.s);
      break;
      case CAP_OBJ_INT:
      fprintf(output_F, "\t int = %d\n",outReq->obj.i);
      break;
      case CAP_OBJ_DOUBLE:
      fprintf(output_F, "\t double = %g\n",outReq->obj.d);
      break;
      case CAP_OBJ_DCX:
      fprintf(output_F, "\t double complex\n");
      break;
      case CAP_OBJ_DOUBLEV:
      fprintf(output_F, "\t double vector\n");
      break;
      case CAP_OBJ_DCXV:
      fprintf(output_F, "\t double complex vector\n");
      break;
      case CAP_OBJ_STRING:
      fprintf(output_F, "\t string = %s\n",outReq->obj.s);
      break;
      case CAP_OBJ_DOUBLEM:
      fprintf(output_F, "\t double matrix\n");
      break;
      case CAP_OBJ_DOUBLEMV:
      fprintf(output_F, "\t double matrix of vectors\n");
      break;
      case CAP_OBJ_DCXM:
      fprintf(output_F, "\t double complex matrix\n");
      break;
      case CAP_OBJ_DCXMV:
      fprintf(output_F, "\t double complex matrix of vectors\n");
      break;
      case CAP_OBJ_OPER:
      fprintf(output_F, "\t operator\n");
      break;
      default:
      fprintf(output_F,"\n");
    }
  }
}

void SaveOutReqList(st_Table_Pt table_P, char *name)
{
  generic_t gval;
  gval.v = (VOIDPTR) thisOutReqList;
  if (St_DefSym(table_P, name, GEN_TYPE_OUTREQ_LIST, gval) == ST_SYM_FOUND)
    ParseError("Duplicate output request");
  thisOutReqList = NULL;
}

char *MakeSequence(char *name)
{
  // Allocate memory here. Who knows where will it be deallocated.
  char* ss = (char*)(malloc(sizeof(char) * (MAX_NAME_LEN+10)));
  if(strlen(name) > MAX_NAME_LEN)
    name = 0; /* Prevent Overflows */
  sprintf(ss,"%s%d",name,sequence++);
  return(ss);
}

void ErrMsg(const char *s)
{
  int i;

  fprintf(stderr, "Error: %s at line %d\n", s, lineNo);
  fprintf(stderr, "%s", inputBuf);
  for (i = 0; i < inputBufp - inputBuf; i++) {
    fputc((inputBuf[i] == '\t' ? '\t' : ' '), stderr);
  }
  fprintf(stderr, "^ before here\n\n");
}


// /*
//  * CleanUpTable--
//  *	Clean up temporary table thisSymT_P
//  * NOT USED ANY MORE
//  */
// void CleanUpTable()
// {
	//   char **sym_A = NULL;
	//   char **sym_B = NULL;
	//   char *chars_P = NULL;

	//   /*
	//    * Get rid of the temporary symbol table.  We have to de-allocate
	//    * storage for any strings in the table.
	//    */
	//   assert(thisSymT_P);

	//   generic_t gval;
	//   int type;

	//   St_ListTable(thisSymT_P, &sym_A, &chars_P);
	//   if (sym_A) {
		//     sym_B = sym_A;
		//     for (; *sym_A; sym_A++) {

			//       St_GetSym(thisSymT_P, *sym_A, &type, &gval);
			//       // Strings are freed by symbol table function
			// //       if (type == GEN_TYPE_STRING)
			// // 	free(gval.s);
			//       if (type == GEN_TYPE_DOUBLE_VECTOR)
			// 	free(gval.dv.v);
			//       else if (type == GEN_TYPE_DOUBLE_VECTOR_OF_VECTORS) {
				// 	for(int i = gval.dvv.size; i; i--) free(gval.dvv.dv_P[i].v);
				// 	free(gval.dvv.dv_P);
			//       }
		//     }
		//     St_ListTableFree(&sym_B, &chars_P);
	//   }
	//   St_DelTable(thisSymT_P);
	//   thisSymT_P = NULL;
// }


// Remains old CleanUpParse()
//   /* Clear or reset the parser globals.  */

//   if(vectorType == GEN_TYPE_STRING) {
	//     int i;
	//     for(i=0 ; i < noVectorElements ; i++) {
		//       free(stringVector_P[i]);
		//       stringVector_P[i] = NULL;
	//     }
//   }

//   thisId[0] = thisElementType[0] = thisModelName[0] = inputBuf[0] =
//     inputFilename[0] = noVectorElements = 0;
//   vectorType = GEN_TYPE_NONE;
//   thisNid.val = 0;
//   thisTermCount = thisListCount = 0;
//   thisLevel = -1;
//   inputBufp = inputBuf;
//   lineNo = 0;
//   tc = ' ';

// AddParameterUsage
// Add a parameter name to the usage table for a variable.
// This tables indicates where the variable is used.
void AddParameterUsage(const char *variableName, gr_Id_t id, 
  const char *parameterName)
{
  usageId_Pt	prev_P, new_P, next_P;
  st_Entry_Pt	st_P;
  /* Get pointer to variable info */
  st_P = LookupSym(capOptionsT_P, variableName);
  new_P = (usageId_Pt) malloc(sizeof(usageId_t));
  /* Get the end of the variable usage chain */
  prev_P = st_P->usage_P;
  if(prev_P == NULL)
	{
		st_P->usage_P = new_P;
	}
  else
	{
		next_P = prev_P->next_P;
		while(next_P != NULL)
		{
			prev_P = next_P;
			next_P = next_P->next_P;
		}
		prev_P->next_P = new_P;
	}

  strcpyCareFul(&(new_P->name), parameterName);
  new_P->val = id.val;
  new_P->type = id.type;
  new_P->next_P = NULL;

  return;
}

// AddSymbolTableUsage
// Add a parameter name to the usage table for a variable.
// This is used for symbol tables.
// This tables indicates where the variable is used.
void AddSymbolTableUsage(const char *variableName, char *tableName,
const char *parameterName)
{
  usageId_Pt	prev_P, new_P, next_P;
  st_Entry_Pt	st_P;
  /* Get pointer to variable info */
  st_P = LookupSym(capOptionsT_P, variableName);
  new_P = (usageId_Pt) malloc(sizeof(usageId_t));
  /* Get the end of the variable usage chain */
  prev_P = st_P->usage_P;
  if(prev_P == NULL)
	{
		st_P->usage_P = new_P;
	}
  else
	{
		next_P = prev_P->next_P;
		while(next_P != NULL)
		{
			prev_P = next_P;
			next_P = next_P->next_P;
		}
		prev_P->next_P = new_P;
	}

  strcpyCareFul(&(new_P->name), parameterName);
  strcpyCareFul(&(new_P->tableName), tableName);
  new_P->val = GR_ID_NO_VAL;
  new_P->type = USAGE_ID_SYMBOL_TABLE;
  new_P->next_P = NULL;

  return;
}

void AddAnalysisUsage(const char *variableName, char *tableName,
const char *parameterName)
{
  usageId_Pt	prev_P, new_P, next_P;
  st_Entry_Pt	st_P;
  /* Get pointer to variable info */
  st_P = LookupSym(capOptionsT_P, variableName);
  new_P = (usageId_Pt) malloc(sizeof(usageId_t));
  /* Get the end of the variable usage chain */
  prev_P = st_P->usage_P;
  if(prev_P == NULL)
	{
		st_P->usage_P = new_P;
	}
  else
	{
		next_P = prev_P->next_P;
		while(next_P != NULL)
		{
			prev_P = next_P;
			next_P = next_P->next_P;
		}
		prev_P->next_P = new_P;
	}

  strcpyCareFul(&(new_P->name), parameterName);
  strcpyCareFul(&(new_P->tableName), tableName);
  new_P->val = GR_ID_NO_VAL;
  new_P->type = USAGE_ID_ANALYSIS_TYPE;
  new_P->next_P = NULL;

  return;

}

// AddVariableUsage
// Add a variable name to the usage table for a variable.
// This tables indicates where the variable is used.
void AddVariableUsage(char *variableName, char *prevVariableName)
{
  usageId_Pt	prev_P, new_P, next_P;
  st_Entry_Pt	st_P;
  /* Get pointer to variable info */
  st_P = LookupSym(capOptionsT_P, prevVariableName);
  new_P = (usageId_Pt) malloc(sizeof(usageId_t));
  /* Get the end of the variable usage chain */
  prev_P = st_P->usage_P;
  if(prev_P == NULL) {
    st_P->usage_P = new_P;
  }
  else {
    next_P = prev_P->next_P;
    while(next_P != NULL) {
      prev_P = next_P;
      next_P = next_P->next_P;
    }
    prev_P->next_P = new_P;
  }

  strcpyCareFul(&(new_P->name), variableName);
  new_P->val = 0;
  new_P->type = USAGE_ID_OPTIONS_TYPE;
  new_P->next_P = NULL;

  return;
}

// DumpTable--
// Dump a table (for debugging or utility purposes)
void DumpTable(FILE *output_F, st_Table_Pt st_P)
{
  st_Entry_Pt	entry_P;
  int bucket, i, j, size, bigSize;

  if (!st_P)
  {
    fprintf(output_F, "Table ptr is NULL.\n");
    return;
  }

  fprintf(output_F, "\n'%s' table, %d entries\n", st_P->name, st_P->entryCt);

  for (bucket = 0; bucket < st_P->size; bucket++)
  {
    entry_P = st_P->tab_A[bucket];
    while (entry_P) {
      fprintf(output_F, "%d: \t", bucket);

      switch(entry_P->type)
			{

				case GEN_TYPE_INT:
				fprintf(output_F, " '%s' \t= %d \t(int)\n", entry_P->name,
				entry_P->gval.i);
				break;

				case GEN_TYPE_LONG:
				fprintf(output_F, " '%s' \t= %ld \t(long)\n", entry_P->name,
				entry_P->gval.l);
				break;

				case GEN_TYPE_FLOAT:
				fprintf(output_F, " '%s' \t= %g \t(float)\n", entry_P->name,
				entry_P->gval.f);
				break;

				case GEN_TYPE_DOUBLE:
				fprintf(output_F, " '%s' \t= %g \t(double)\n", entry_P->name,
				entry_P->gval.d);
				break;

				case GEN_TYPE_DOUBLE_VECTOR:
				size = entry_P->gval.dv.size;
				fprintf(output_F, " '%s' \t \t(double vector) size = %d \n",
				entry_P->name, size);
				for(i = 0; i < size; i++)
					fprintf(output_F,"\t\t %g\n", entry_P->gval.dv.v[i]);
				break;

				case GEN_TYPE_DOUBLE_VECTOR_OF_VECTORS:
				bigSize = entry_P->gval.dvv.size;
				fprintf(output_F,
				" '%s' \t \t(vector of double vectors) size = %d \n",
				entry_P->name, bigSize);
				for(j = 0; j < bigSize; j++) {
					size = entry_P->gval.dvv.dv_P[j].size;
					for(j = 0; i < size; i++)
						fprintf(output_F,"\t\t %g\n",
		      entry_P->gval.dvv.dv_P[j].v[i]);
				}
				break;

				case GEN_TYPE_CHAR:
				fprintf(output_F, " '%s' \t= '%c'/%d \t(char)\n", entry_P->name,
				entry_P->gval.c, entry_P->gval.c);
				break;

				case GEN_TYPE_STRING:
				if (entry_P->gval.s)
					fprintf(output_F, " '%s' \t= '%s' \t(string)\n",
		    entry_P->name, entry_P->gval.s);
				else
					fprintf(output_F, " '%s' \t= NULL \t(string)\n",
		    entry_P->name);
				break;

				case GEN_TYPE_STRING_A:

				if (!entry_P->gval.s_A) {
					fprintf(output_F, " '%s' \t= NULL \t(string array)\n",
					entry_P->name);
				} else {
					char	**s_A;
					fprintf(output_F, " '%s' is string array:\n",
					entry_P->name);
					for (s_A = entry_P->gval.s_A; *s_A; s_A++)
						fprintf(output_F, "\t\t '%s'\n", *s_A);
				}
				break;

				case GEN_TYPE_VOID_PTR:
				if (entry_P->gval.v)
					fprintf(output_F, " '%s' \tis user type \t(void ptr)\n",
		    entry_P->name);
				else
					fprintf(output_F, " '%s' \tis NULL user type \t(void ptr)\n",
		    entry_P->name);
				break;

				case GEN_TYPE_VARIABLE_USAGE:
				if (entry_P->gval.v)
					fprintf(output_F, " '%s' \t= '%s' \t(variable usage)\n",
		    entry_P->name, ((st_Entry_Pt) entry_P->gval.v)->name);
				else
					fprintf(output_F, " '%s' \t= NULL \t(variable usage)\n",
		    entry_P->name);
				break;

				case GEN_TYPE_OUTREQ_LIST:
				if (entry_P->gval.v)
				{
					fprintf(output_F, " '%s' \t \t(output request) =\n",
		      entry_P->name);
					DumpOutReqList( (capOutReq_Pt) (entry_P->gval.v) );
				}
				else
					fprintf(output_F, " '%s' \t \t(output request) = \t\t NULL\n",
		    entry_P->name);
				break;

				default:
				fprintf(output_F, " '%s' \tis unknown type\n", entry_P->name);
				break;
			}
      dumpUsage(output_F,entry_P->usage_P);
      entry_P = entry_P->next_P;
    }
  }
  return;
}

// VariableIncorrect
void VariableIncorrect(char *variable)
{
  char ss[MAX_STRING_LEN];
  sprintf(ss,
	"Variable `%s' incorrectly specified or not specified\n",variable);
  ParseError(ss);
  maimed = TRUE;
  return;
}

void initElement( char *elementName, char *elementType)
{
  strncpy(thisElementName, elementName, MAX_NAME_LEN);
  strncpy(thisElementType, elementType, MAX_NAME_LEN);
  thisSymT_P = St_NewTable("temp",HASH_SIZE);
  thisTermCount = 0;
  thisModelName [0] = 0;
  usageType = USAGE_ID_ELEMENT_TYPE;
  usageVal  = GR_ID_NO_VAL;
}

void checkLevel()
{
  char s[MAX_NAME_LEN +2], modelName[MAX_NAME_LEN+2];
  st_Table_Pt temp;
  temp = thisSymT_P;
	if (St_GetSymAsInt(thisSymT_P, "level", &thisElementLevel) == ST_SYM_FOUND)
	{
		//  strcpy(thisElementLevel,s);
		St_MarkSym(thisSymT_P, "level");
	}
  else if(St_GetSymAsString(thisSymT_P, "model", modelName) == ST_SYM_FOUND)
	{
		thisSymT_P = St_GetTable(modelName);
		if(thisSymT_P)
	  {
	    if(St_GetSymAsInt(thisSymT_P, "level",&thisElementLevel)==ST_SYM_FOUND)
	    {
	      //strcpy(thisElementLevel,s);
	      St_MarkSym(thisSymT_P,"level");
	    }
	    else
	    {
	      fprintf(stderr, "Model level not found in model %s\n",modelName);
	      fprintf(stderr, "Continuing\n");
	    }
	  }
	}
	else
	{
		fprintf(stderr, "Model type not found in model %s\n",s);
		fprintf(stderr, "Continuing\n");
	}
	thisSymT_P = temp;
}

void checkModel()
{
  char s[MAX_NAME_LEN +2], modelName[MAX_NAME_LEN+2];
  st_Table_Pt temp;
  temp = thisSymT_P;
  //  if (St_GetSymAsString(thisSymT_P, "modeltype", s) == ST_SYM_FOUND)
  //{
  // strcpy(thisModelType,s);
  // St_MarkSym(thisSymT_P,"modeltype");
  //	 return;
  // }
  if(St_GetSymAsString(thisSymT_P, "model", modelName) == ST_SYM_FOUND)
    {
    thisSymT_P = St_GetTable(modelName);
    if(thisSymT_P)
      {
      if (St_GetSymAsString(thisSymT_P, "TYPE", s) ==  ST_SYM_FOUND)
        {
        if (!strcmp(s, thisElementType))
          fprintf(stderr,"freeda syntax of element");
          //fprintf(stderr,"Attempt to use model %s of type %s in element %s\n",
          //  modelName,s, thisElementType);
        else
          //SPICE syntax
          //if(St_GetSymAsString(thisSymT_P, "modeltype",s)==ST_SYM_FOUND)
          {
          strcpy(thisModelType,s);
          St_MarkSym(thisSymT_P,"TYPE");
          }
        }
      }
    else {
      fprintf(stderr, "Model %s not found for element %s\n",s,thisElementType);
      fprintf(stderr, "Continuing\n");
      }
    }
  thisSymT_P= temp;
}

void closeElement()
{
  // Create element and connect terminals.
  char modelName[MAX_NAME_LEN+2], s[MAX_NAME_LEN+2];
  st_Table_Pt temp;
  temp = thisSymT_P;
  generic_t gval;
	if(St_GetSymAsString(thisSymT_P, "model", modelName) == ST_SYM_FOUND)
	{
		thisSymT_P = St_GetTable(modelName);
		if(thisSymT_P)
	  {
	    if (St_GetSymAsString(thisSymT_P, "TYPE", s) ==  ST_SYM_FOUND)
			{
				gval.s = (char *) malloc(strlen(thisElementType)+1);
				strcpy(gval.s,thisElementType);
				St_ReplSym(thisSymT_P,"TYPE",GEN_TYPE_STRING,gval);
			}
	  }
	}
	thisSymT_P = temp;
  Element* elem = ConstrPrim();

  if (elem)
  {
    /* We need to update the expression and sweep usage information as
		* we did not have the element id when we did the parsing */
    gr_Id_t elementId;
    elementId.val = elem->getId();
    elementId.type = GR_ID_NODE_TYPE;
    v_updateExpressionIds(elementId.val, elementId.type);
    v_updateSweepIds(elementId.val, elementId.type);
  }
  // free memory used by thisSymT_P
  St_DelTable(thisSymT_P);
  thisSymT_P = NULL;
}

void insertPolyParam(char* elementName)
{
  //locate the element in the current ckt.
  Circuit* cir = the_CM->getCurrent() ;
  Element* elem = cir->getElement(elementName) ;

  //next 3 lines insert polynomial related parameters.
  DenseDoubleVector coeff_temp(Teuchos::View, &poly_coeff[0], number_poly_coeffs);
  elem->setParam("poly_coeff",&coeff_temp,TR_DOUBLE_VECTOR) ;
  elem->setParam("polydimension",&polydimension,TR_INT) ;

	//connects the controlling terminals.
  for(int i=0;i<noElementsInList;i++)
    cir->connect(elem->getId(),polyList[i].val) ;
}

void initAnalysis(char* analysisName)
{
  thisSymT_P = St_NewTable("temp", HASH_SIZE);
  usageType = USAGE_ID_ANALYSIS_TYPE;
  usageTableName = analysisName;
};

void closeAnalysis()
{
  // Create the analysis object.
  // This routine will assign the global variable "analysis".
  createAnalysis(usageTableName);

  // Edit analysis parameters
  EditItemParams(analysis);

  // free memory used by thisSymT_P
  St_DelTable(thisSymT_P);
  thisSymT_P = NULL;
}

// SetUpModels--
// Goes through list of elements.  If a model is specified, the
// parameter values are set if not already specified.
int setupModels()
{
	// Add in symbols from model, if any, that don't
	// conflict with symbols defined in the element def.  St_DefSym()
	// doesn't re-define a symbol if it already exists.
  // Go through all elements in all circuits
  unsigned nc = the_CM->getNumberOfCircuits();
  Circuit* cir;
  for (unsigned j=0; j < nc; j++)
  {
    cir = the_CM->getCircuit(j);
    cir->setFirstElement(0);
    Element* elem = cir->nextElement();
    while (elem)
    {
      // This is to avoid the contravariance violation
      char modelName[MAX_NAME_LEN];
      strcpy(modelName, elem->getModelName());
      if(modelName[0])
      {
				thisSymT_P = St_GetTable(modelName);
				if(thisSymT_P)
				{
					char s[MAX_NAME_LEN + 2];
					if (St_GetSymAsString(thisSymT_P, "TYPE", s) == ST_SYM_FOUND)
					{
						//  strcpy(s ,elem->getName().c_str());
						if (!strcmp(s, elem->getName().c_str()))
							// Add parameters but do not overwrite
						  EditParams(elem, false);
						else
						{
							fprintf(stderr,
							"Attempt to use model %s of type %s in element %s\n",
							modelName, s,
							(elem->getInstanceName()).c_str());
							fprintf(stderr, "Continuing\n");
						}
					}
					else
					{
						fprintf(stderr,
						"Mode type not found in model %s\n",
						modelName);
						fprintf(stderr, "Continuing\n");
					}
				}
				else
				{
					fprintf(stderr, "Model %s not available for element %s\n",
					modelName, (elem->getInstanceName()).c_str());
					fprintf(stderr, "Continuing\n");
				}
      }
      elem = cir->nextElement();
    }
  }

  return(0);
}

// Set parameter value and convert between gval and new formats.
void setItemParam(NetListItem* item, const int& pIndex,
const ParamType& p_type, generic_t gvalStore)
{
  switch (p_type)
  {
		case TR_INT:
		case TR_LONG:
		case TR_FLOAT:
		case TR_DOUBLE:
		case TR_CHAR:
    {
      item->setParam(pIndex, &gvalStore, p_type);
    }
    break;

		case TR_BOOLEAN:
    {
      // Convert here
      bool var;
      if (gvalStore.i)
				var = true;
      else
				var = false;
      item->setParam(pIndex, &var, p_type);
    }
    break;

		case TR_COMPLEX:
    {
      double_complex x1(gvalStore.x.re, gvalStore.x.im);
      item->setParam(pIndex, &x1, p_type);
    }
    break;

		case TR_STRING:
    {
      // Build a string
      string tmpstr(gvalStore.s);
      item->setParam(pIndex, &tmpstr, p_type);
    }
    break;

		case TR_INT_VECTOR:
    {
      // Build an integer vector
      // The vector is parsed as double, convert here
      DenseIntVector tmpiv(gvalStore.dv.size);
      for (int i=0; i<gvalStore.dv.size; i++)
				tmpiv[i] = int(gvalStore.dv.v[i]);
      item->setParam(pIndex, &tmpiv, p_type);
    }
    break;

		case TR_DOUBLE_VECTOR:
    {
      // Build a double vector
      DenseDoubleVector tmpdv(gvalStore.dv.size);
			for (int i = 0; i < gvalStore.dv.size; i++)
				tmpdv[i] = gvalStore.dv.v[i];
      item->setParam(pIndex, &tmpdv, p_type);
    }
    break;

		case TR_DOUBLE_MATRIX:
    {
      // Build a double matrix
      // First find matrix dimensions (number of columns)
      int cols = 0;
      for (int i=0; i<gvalStore.dvv.size; i++)
				if (gvalStore.dvv.dv_P[i].size > cols)
					cols = gvalStore.dvv.dv_P[i].size;
      // Create the matrix initialized to zero
      DoubleDenseMatrix tmpdm(gvalStore.dvv.size, cols);
      tmpdm.putScalar(zero);
      // Assign nonzero elements
      for (int i=0; i<gvalStore.dvv.size; i++)
				for (int j=0; j<gvalStore.dvv.dv_P[i].size; j++)
					tmpdm(i,j) = gvalStore.dvv.dv_P[i].v[j];

      item->setParam(pIndex, &tmpdm, p_type);
    }
    break;

		case TR_INT_MATRIX:
    {
      // Build an integer matrix
      // First find matrix dimensions (number of columns)
      // The matrix comes as a double_vector_of_doubles, convert here.
      int cols = 0;
      for (int i=0; i<gvalStore.dvv.size; i++)
				if (gvalStore.dvv.dv_P[i].size > cols)
					cols = gvalStore.dvv.dv_P[i].size;
      // Create the matrix initialized to zero
      IntDenseMatrix tmpim(gvalStore.dvv.size, cols);
      // Assign nonzero elements
      for (int i=0; i<gvalStore.dvv.size; i++)
				for (int j=0; j<gvalStore.dvv.dv_P[i].size; j++)
					tmpim(i,j) = int(gvalStore.dvv.dv_P[i].v[j]);

      item->setParam(pIndex, &tmpim, p_type);
    }
  }
}

// EditParams --
// Edit element parameters.
// If overwrite is false, we assume that thisSymT_P is a model symbol table
// and therefore can not contain a symbol called "model".
void EditParams(Element* elem, bool overwrite)
{
  generic_t gval, gvalStore;
  int pIndex, sType, sTypeStore, pCount;
  ParamType p_type;
  const char *variableName;
 	const char *parameterName;
  gr_Id_t id;

  // Set a valid value for id
  id.type = GR_ID_NODE_TYPE;
  id.val = elem->getId();

  // Set parameters using values in temporary symbol table.
  pCount = elem->getNumberOfParams();
  for (pIndex = 0; pIndex < pCount; pIndex++)
  {
    elem->getParamSpec(pIndex, p_type, parameterName);
    if(parameterName == NULL)
    {
      ParseError("Internal Error in specification of parameters.");
      maimed = TRUE;
      return;
    }

    // Attempt to get values from temporary symbol table. (The element line
    // parameters).
    if (St_GetSym(thisSymT_P, parameterName, &sType, &gval) ==
			ST_SYM_FOUND)
    {
      // Mark symbol as used
      St_MarkSym(thisSymT_P, parameterName);
      // Set the parameter only if needed
      if (!(elem->isSet(unsigned(pIndex))) || overwrite)
      {
				switch(sType)
				{
					case GEN_TYPE_VARIABLE_USAGE:
					/* Record variable usage */
					variableName = (char *) gval.v;
					St_GetSym(capModelT_P, variableName, &sType, &gval);
					AddParameterUsage(variableName, id, parameterName);
					default:
					;
				}

				gvalStore = gval;
				sTypeStore = sType;

				if (Coerce(sTypeStore, p_type, &gvalStore))
					fprintf(stderr, "^error defining '%s'\n\n", parameterName);
				else
				{
					// Set element parameter
					setItemParam(elem, pIndex, p_type, gvalStore);
				}
      }
    }
  }

  // There should not be any unresolved parameters left in the
  // temporary symbol table. If there is, this indicates an
  // error. The unresolved symbols will be listed. Otherwise this
  // condition will be ignored.
  if(!St_CheckAllMarked(stderr, thisSymT_P))
  {
		char s[MAX_NAME_LEN];
		sprintf(s,"    Warning: unused parameters for element %s\n",
		(elem->getInstanceName()).c_str());
		ParseWarning(s);
  }

  return;
}

// Replace the parameter value of an element
void replaceElemParam(gr_Id_t id, int type, generic_t gval, char* pname)
{
  // The current variable code seems to work only for type double
  assert(type == GEN_TYPE_DOUBLE);
  Element* elem = the_CM->getCurrent()->getElement(id.val);
  // I guess this code will only work with variables outside
  // subcircuits
  assert(elem);
  try
  {
    ParamType p_type = elem->askParamType(pname);
    assert(p_type == TR_DOUBLE);
    elem->setParam(pname, &gval, p_type);
  }
  catch (string&)
  {
    ParseError("Attempt to use a variable inside a subcircuit");
  }
}

// Edit Item parameters from symbol table
void EditItemParams(NetListItem* item, bool overwrite)
{
  generic_t gval, gvalStore;
  int pIndex, sType, sTypeStore, pCount;
  ParamType p_type;
  const char *parameterName;

  // Set parameters using values in temporary symbol table.
  pCount = item->getNumberOfParams();
  for (pIndex = 0; pIndex < pCount; pIndex++)
  {
    item->getParamSpec(pIndex, p_type, parameterName);

    if(parameterName == NULL)
    {
      ParseError("Internal Error in specification of parameters.");
      maimed = TRUE;
      return;
    }

    // Attempt to get values from temporary symbol table. (The item line
    // parameters).
    if (St_GetSym(thisSymT_P, parameterName, &sType, &gval) == ST_SYM_FOUND)
    {
      // Mark symbol as used
      St_MarkSym(thisSymT_P, parameterName);

      // Set the parameter only if needed
      if (!(item->isSet(unsigned(pIndex))) || overwrite)
      {
				switch(sType)
				{
					case GEN_TYPE_VARIABLE_USAGE:
					{
						/* Record variable usage */
						char* variableName = (char *) gval.v;
						St_GetSym(capModelT_P, variableName, &sType, &gval);
						AddAnalysisUsage(variableName, usageTableName, parameterName);
					}
					default:
					;
				}

				gvalStore = gval;
				sTypeStore = sType;

				if (Coerce(sTypeStore, p_type, &gvalStore))
					fprintf(stderr, "^error defining '%s'\n\n", parameterName);
				else
				{
					// Set item parameter
					setItemParam(item, pIndex, p_type, gvalStore);
				}
      }
    }
  }

  // There should not be any unresolved parameters left in the
  // temporary symbol table. If there is, this indicates an
  // error. The unresolved symbols will be listed. Otherwise this
  // condition will be ignored.
  if(!St_CheckAllMarked(stderr, thisSymT_P))
  {
		char s[MAX_NAME_LEN];
		sprintf(s,"    Warning: unused parameters for netlist item %s\n",
		(item->getName()).c_str());
		ParseWarning(s);
  }
  return;
}

// This code will have to be replaced if we support more than
// one analysis in the same run.
void replaceAnalysisParam(char* atype, generic_t gval, char* pname)
{
  assert(analysis->getName() == string(atype));
  // I guess this code will only work with variables outside
  // subcircuits
  try
  {
    ParamType p_type = analysis->askParamType(pname);
    assert(p_type == TR_DOUBLE || p_type == TR_INT);
    analysis->setParam(pname, &gval, p_type);
  }
  catch (string&)
  {
    ParseError("Something is wrong with the variable in the analysis");
  }
}

// ConstrPrim--
//	Construct an element and add it to the current network.
Element* ConstrPrim()
{
  int i;
  unsigned elem_id;
  Element* elem;

  // The element type is stored as thisElementType.
  //
  // e.g. if name = rload      thisElementType is r
  //      if name = mstrip:120 thisElementType is mstrip
  //
  // The actual name is stored in thisElementName.
  try
  {
    // Create element (and add it to circuit).
    elem_id =
		the_CM->getCurrent()->addElement(thisElementType, thisElementName);

    // Connect terminals
    // number of terminals in netlist: thisTermCount
    for (i = 0; i < thisTermCount; i++)
    {
      // thisTermList contains the ids of the terminals
      the_CM->getCurrent()->connect(elem_id, thisTermList[i].val);
    }

    elem = the_CM->getCurrent()->getElement(elem_id);
    assert(elem);

    // Special treatment for subcircuit instances.
    if(!strcmp(elem->getName().c_str(),"x"))
    {
      // We know that we have a pointer to a Xsubckt element, so
      // try to convert the element pointer (hope it works).
      Xsubckt* subckt = (Xsubckt*)(elem);
      // Set additional information required by the instance.
      subckt->attachDefinition(the_CM->getSubCircuit(thisModelName));
      // Now add this instance to the main circuit's list
      the_CM->getCurrent()->addXinstance(subckt);
    }
    else
    {
      // Set model name if specified for other elements
      if (thisModelName[0])
				elem->setModelName(thisModelName);
      else {
				char s[MAX_NAME_LEN + 2];
				// Get modelname from symbol table
				if (St_GetSymAsString(thisSymT_P, "model", s) ==
				ST_SYM_FOUND) {
					elem->setModelName(s);
					// Erase modelname from symbol table
					St_DelSym(thisSymT_P, "model");
				}
      }
    }


    // set up parameters
    EditParams(elem);

  }
  catch (string& error)
  {
    ParseError(error.c_str());
    elem = NULL;
    maimed = TRUE;
  }

  return(elem);
}

// AddSubCkt--
//     This will change after adding the new subcircuit handler.
int AddSubCkt()
{
  // thisId: "subckt".
  // subcircuit name: char* thisSubCkt (already given in push()).
  // parameters in: symbol table "thisSubCkt".
  // This function must return 0 if everything Ok.
  // Parameters and terminals are stored in the way of any other element.

  // Set external subcircuit terminals
  // Terminal ids are in thisTermList[]
  for (int i = 0; i < thisTermCount; i++)
    // thisTermList contains the ids of the terminals
	the_CM->getCurrent()->setConnection(thisTermList[i].val);

  // Later add parameter processing, if any.
  return 0;
}

// The following are Transim2 brand new functions ...
//
// Note: id structs are no longer used in the network code.
// Nevertheless, each element and terminal still have a unique
// id number.

// Returns a ficticious id struct from a terminal
gr_Id_t getTermNamed(char *name)
{
  gr_Id_t tid;
  tid.type = GR_ID_TERM_TYPE;
  Terminal* term = the_CM->getCurrent()->getTerminal(name);
  if (term)
    tid.val = term->getId();
  else
    tid.val = 0;

  return tid;
}


// Returns a ficticious id struct from an element
gr_Id_t getElemNamed(char *name)
{
  gr_Id_t nid;
  nid.type = GR_ID_NODE_TYPE;
  Element* elem = the_CM->getCurrent()->getElement(name);
  if (elem)
    nid.val = elem->getId();
  else
    nid.val = 0;
  return nid;
}

// Returns a ficticious id struct from a new terminal
gr_Id_t newTerm(char* name)
{
  gr_Id_t tid;
  tid.type = GR_ID_TERM_TYPE;
  tid.val = the_CM->getCurrent()->addTerminal(name);
  // The terminals named "0" or "gnd" are considered to be references
  if (!strcmp(name, "0") || !strcmp(name, "gnd"))
    the_CM->getCurrent()->setRefTerm(tid.val);
  return tid;
}

// Returns a ficticious id struct from a new element
gr_Id_t newElem(char* elemtype, char* name)
{
  gr_Id_t nid;
  nid.type = GR_ID_NODE_TYPE;
  // Add the element to current circuit.
  nid.val = the_CM->getCurrent()->addElement(elemtype,name);
  return nid;
}

// This is just a front-end function
int popCircuit()
{
  return int(the_CM->pop());
}

// This is just a front-end function
int pushCircuit(char* name)
{
  the_CM->push(name);
  return 0;
}

// This function will probably not be needed
//
// This is just a front-end function
// void getParentCircuitName(char* name)
// {
	//   name = (const char*)(the_CM->getParent()->getInstanceName());
// }

// This function is called at the end of the parsing
void setupCircuits()
{
  setupModels();
  //if (!maimed) {
		// try {
			//  the_CM->check();
		// }
		//catch (string& error) {
			//  report(FATAL, error.c_str());
		// }
  // }
}

// Set the x coodinate of a terminal. 
void setTermX(gr_Id_t id, double x)
{
  //ParseWarning("Coordinates implemented");
  the_CM->getCurrent()->getTerminal(id.val)->setX(x);
}

// Same thing for y.
void setTermY(gr_Id_t id, double y)
{
  //ParseWarning("Coordinates implemented");
  the_CM->getCurrent()->getTerminal(id.val)->setY(y);
}

// Same thing for z. OR calculate Z from layer.
void setTermZ(gr_Id_t id, double z)
{
  the_CM->getCurrent()->getTerminal(id.val)->setZ(z);
}

// Set one terminal as (local or global) reference.
void setReference(gr_Id_t tid)
{
  the_CM->getCurrent()->setRefTerm(tid.val);
}


// Coerce--
// Coerce a symbol to the appropriate parameter type.
//
// sType -- original type
// pType -- final type
int Coerce(int sType, ParamType pType, generic_t *gval_P)
{
  generic_t gtmp;
  int tmpType;
  int status = 0;

  switch (sType)
  {

		case GEN_TYPE_NONE:
    ParseError("Couldn't coerce TYPE_NONE\n");
    status = 1;
    break;

		case GEN_TYPE_INT:
    switch (pType)
    {
			case TR_FLOAT:
      gtmp.f = gval_P->i;
      *gval_P = gtmp;
      break;
			case TR_DOUBLE:
      gtmp.d = gval_P->i;
      *gval_P = gtmp;
      break;
			case TR_BOOLEAN:
			case TR_INT:
      break;
			case TR_CHAR:
      gtmp.c = gval_P->i;
      *gval_P = gtmp;
      break;
			case TR_STRING:
      ParseError("Couldn't coerce INT to STRING\n");
      status = 1;
      break;
			default:
      ParseError("Couldn't coerce INT to unknown\n");
      status = 1;
      break;
    }
    break;

		case GEN_TYPE_LONG:
    switch (pType)
    {
			case TR_FLOAT:
      gtmp.f = gval_P->l;
      *gval_P = gtmp;
      break;
			case TR_DOUBLE:
      gtmp.d = gval_P->l;
      *gval_P = gtmp;
      break;
			case TR_INT:
			case TR_BOOLEAN:
      gtmp.i = gval_P->l;
      *gval_P = gtmp;
      break;
			case TR_CHAR:
      gtmp.c = gval_P->l;
      *gval_P = gtmp;
      break;
			case TR_STRING:
      ParseError("Couldn't coerce LONG to STRING\n");
      status = 1;
      break;
			default:
      ParseError("Couldn't coerce LONG to unknown\n");
      status = 1;
      break;
    }
    break;

		case GEN_TYPE_FLOAT:
    switch (pType)
    {
			case TR_FLOAT:
      break;
			case TR_DOUBLE:
      gtmp.d = gval_P->f;
      *gval_P = gtmp;
      break;
			case TR_INT:
      gtmp.i = int(gval_P->f);
      *gval_P = gtmp;
      break;
			case TR_CHAR:
			case TR_BOOLEAN:
      ParseError("Couldn't coerce FLOAT to CHAR, or BOOLEAN\n");
      status = 1;
      break;
			case TR_STRING:
      ParseError("Couldn't coerce FLOAT to STRING\n");
      status = 1;
      break;
			default:
      ParseError("Couldn't coerce FLOAT to unknown\n");
      status = 1;
      break;
    }
    break;

		case GEN_TYPE_DOUBLE:
    switch (pType)
    {
			case TR_FLOAT:
      gtmp.f = gval_P->d;
      *gval_P = gtmp;
      break;
			case TR_DOUBLE:
      break;
			case TR_INT:
      gtmp.i = int(gval_P->d);
      *gval_P = gtmp;
      break;
			case TR_CHAR:
			case TR_BOOLEAN:
      ParseError("Couldn't coerce DOUBLE to CHAR, or BOOLEAN\n");
      status = 1;
      break;
			case TR_STRING:
      ParseError("Couldn't coerce DOUBLE to STRING\n");
      status = 1;
      break;
			default:
      ParseError("Couldn't coerce DOUBLE to unknown\n");
      status = 1;
      break;
    }
    break;

		case GEN_TYPE_CHAR:
    switch (pType)
    {
			case TR_FLOAT:
      gtmp.f = gval_P->c;
      *gval_P = gtmp;
      break;
			case TR_DOUBLE:
      gtmp.d = gval_P->c;
      *gval_P = gtmp;
      break;
			case TR_BOOLEAN:
			case TR_INT:
      gtmp.i = gval_P->c;
      *gval_P = gtmp;
      break;
			case TR_CHAR:
      break;
			case TR_STRING:
      ParseError("Couldn't coerce CHAR to STRING\n");
      status = 1;
      break;
			default:
      ParseError("Couldn't coerce CHAR to unknown\n");
      status = 1;
      break;
    }
    break;

		case GEN_TYPE_STRING:
    switch (pType)
    {
			case TR_STRING:
      break;
			default:
      ParseError("Couldn't coerce STRING to numeric quantity\n");
      status = 1;
      break;
    }
    break;

		case GEN_TYPE_INT_VECTOR:
    // This case will never happen because the parser parse all
    // vectors as double vectors
    ParseError("Couldn't coerce vector\n");
    status = 1;
    break;

		case GEN_TYPE_DOUBLE_VECTOR:
    switch (pType)
    {
			case TR_INT_VECTOR:
			case TR_DOUBLE_VECTOR:
      // Leave the vectors as they are and convert later.
      gtmp.dv = gval_P->dv;
      *gval_P = gtmp;
      break;
			default:
      ParseError("Couldn't coerce vector\n");
      status = 1;
      break;
    }
    break;

		case GEN_TYPE_DOUBLE_VECTOR_OF_VECTORS:
    switch (pType)
    {
			case TR_DOUBLE_MATRIX:
			case TR_INT_MATRIX:
      // Leave the vectors of vectors as they are and convert later.
      gtmp.dvv = gval_P->dvv;
      *gval_P = gtmp;
      break;
			default:
      ParseError("Couldn't coerce vector of vectors\n");
      status = 1;
      break;
    }
    break;

		case GEN_TYPE_STRING_VECTOR:
    ParseError("Couldn't coerce vector\n");
    status = 1;
    break;

		case GEN_TYPE_BOOLEAN:
    switch (pType)
    {
			case TR_INT:
			case TR_BOOLEAN:
      break;
			default:
      ParseError("Couldn't coerce boolean\n");
      status = 1;
      break;
    }
    break;

		case GEN_TYPE_ENUM:
    switch (pType)
    {
			case TR_INT:
      break;
			default:
      ParseError("Couldn't coerce enum\n");
      status = 1;
      break;
    }
    break;

		case GEN_TYPE_VARIABLE_USAGE:
    St_GetSym(capOptionsT_P, gval_P->s, &tmpType, &gtmp);
    switch(pType)
		{
      case TR_DOUBLE:
			switch(tmpType)
			{
				case GEN_TYPE_FLOAT:
				gtmp.d = gtmp.f;
				*gval_P = gtmp;
				break;
				case GEN_TYPE_DOUBLE:
				*gval_P = gtmp;
				break;
				case GEN_TYPE_INT:
				gtmp.d = gtmp.i;
				*gval_P = gtmp;
				break;
				default:
				ParseError("Couldn't coerce variable.\n");
				status = 1;
				break;
			}
			break;

      case TR_INT:
			switch(tmpType)
			{
				case GEN_TYPE_FLOAT:
				gtmp.i = int(gtmp.f);
				*gval_P = gtmp;
				break;
				case GEN_TYPE_DOUBLE:
				gtmp.i = int(gtmp.d);
				*gval_P = gtmp;
				break;
				case GEN_TYPE_INT:
				*gval_P = gtmp;
				break;
				default:
				ParseError("Couldn't coerce variable.\n");
				status = 1;
				break;
			}
			break;

      default:
			ParseError("Couldn't coerce variable.\n");
		}
    break;

		default:
    ParseError("Couldn't coerce quantity\n");
    status = 1;
    break;

  }

  return(status);
}
