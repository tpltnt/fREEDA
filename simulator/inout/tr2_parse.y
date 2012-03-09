/***********************************************************************
* tr2_parse.y--
*	YACC grammar and actions for parsing netlists.
*
***********************************************************************/

%{
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parser.h"

int noVectors = 0;
int modelflag=0,pCount=0;
/* number of vectors entered
* so far in doubleVector_P */

genericDoubleV_t doubleVector_P[ST_MAX_VECTOR_LENGTH];
/* Used as temporary storage for
* pointers to vectors */

%}

/*
* TOKEN DEFINITIONS
*/

/*
* The parser's stack can hold ints, doubles, and strings.
*/
%union {
int ival;
double dval;
char *id;
}

/*
* Reserved words in input netlist
*/
%token NOECHO WRITE PLOT WATCH IN_ ELEMENT TERM ELEM FILE_ SYSTEM POLY_
%token VARIABLE_ END HB_ TRAN_ PARAMS_ NOISE_

/*
* SPICE DOTTED COMMANDS
*/
%token DOT_AC	DOT_DC
%token DOT_END DOT_ENDS DOT_FOUR DOT_IC
%token DOT_INC DOT_KEEP
%token DOT_LIB DOT_MODEL DOT_NOISE
%token DOT_NODESET DOT_OP DOT_OPTIONS	DOT_PLOT
%token DOT_PRINT DOT_PROBE	DOT_SENS DOT_SUBCKT
%token DOT_TF DOT_TRAN DOT_WIDTH DOT_TEMP
%token DOT_OUT DOT_WATCH

/*
* SPECIAL KEY WORDS
*/
%token EOL	WHITESPACE	VARIABLE	PARAMETER     EXPRESSION

/*
* SPECIAL KEY WORDS FOR SPICE COMMAND ARGUEMENTS
*/
%token DC	TRAN	AC	PULSE	SIN	EXP	PWL	SFFM	NOISE

/*
* ADDITIONAL SPICE-LIKE DOTTED COMMANDS
*/
%token DOT_LOCATE DOT_REF DOT_SVHB
%token DOT_TRAN2 DOT_TRAN3 DOT_SVTR DOT_WAVTRAN DOT_TRAN4 DOT_TRANT

/*
* TOKENS RETURNED FOR .PRINT FOR BOTH V( AND I(
*/
%token PRINT_V PRINT_I


/*
* SPICE element tokens
*/
%token A_ELEMENT	B_ELEMENT	C_ELEMENT	D_ELEMENT
%token E_ELEMENT	F_ELEMENT	G_ELEMENT	H_ELEMENT
%token I_ELEMENT	J_ELEMENT	K_ELEMENT	L_ELEMENT
%token M_ELEMENT	N_ELEMENT	O_ELEMENT	P_ELEMENT
%token Q_ELEMENT	R_ELEMENT	S_ELEMENT	T_ELEMENT
%token U_ELEMENT	V_ELEMENT	W_ELEMENT	X_ELEMENT
%token Y_ELEMENT	Z_ELEMENT

/*
* Entities
*/

%token <id>   IDENTIFIER	/* identifier, not a keyword		    */
%token <dval> FLOAT_LIT		/* floating-pt literal			    */
%token <ival> INT_LIT		/* integer literal			    */
%token <id>   STRING_LIT	/* string literal			    */
%token <ival> BOOL_LIT		/* boolean literal			    */

/*
* ASCII chars are their own tokens
*/



%start	netlist
/*
* THE PARSER
*/

%%

netlist		: /* empty */
| netlist dotted_command
| netlist additional_element
| netlist spice_element
| netlist EOL /* Ignore leading blank lines */
| error
{ ErrMsg("Error in netlist");
maimed = TRUE;
yyerrok;
}
;

maybe_comma     :
| ','
;
dotted_command	: DOT_AC
{
initAnalysis("ac");
}
sym_assgn_list_1 EOL
{
closeAnalysis();
}
| DOT_DC
{
initAnalysis("dc");
}
sym_assgn_list_1 EOL
{
closeAnalysis();
}
| DOT_END
{ if (popCircuit())
ParseError(".ENDS card missing.");
setupCircuits();
}
EOL
| DOT_ENDS maybe_subname
{ popCircuit(); }
EOL
| DOT_FOUR
{ ErrMsg("Command Not Implemented");
maimed = TRUE;
yyerrok;
}
EOL
| DOT_IC
{ ErrMsg("Command Not Implemented");
maimed = TRUE;
yyerrok;
}
EOL
| DOT_INC
{ 
printf("XXXX Just parsed dot_inc \n");
}
| DOT_KEEP
{ ErrMsg("Command Not Implemented");
maimed = TRUE;
yyerrok;
}
EOL
| DOT_LIB
{ ErrMsg("Command Not Implemented");
maimed = TRUE;
yyerrok;
}
EOL
| DOT_LOCATE any_literal
{
thisTid = getTermNamed(lastText);
if (!thisTid.val) {
thisTid = newTerm(lastText);
}
}
number_as_double
{ setTermX(thisTid, globalGval.d); }
number_as_double
{ setTermY(thisTid, globalGval.d); }
EOL
| DOT_MODEL
any_identifier
{
strncpy(&thisId[0], lastText, MAX_NAME_LEN);
}
any_identifier
{ generic_t gval;
int type;

gval.v = NULL;
if (St_DefSym(capModelT_P, thisId,
GEN_TYPE_VOID_PTR, gval) ==
ST_SYM_FOUND) {
St_GetSym(capModelT_P,thisId,&type,&gval);
thisSymT_P = (st_Table_Pt) gval.v;
} else {
gval.v = (VOIDPTR) (thisSymT_P =
St_NewTable(thisId, HASH_SIZE));
St_ReplSym(capModelT_P, thisId,
GEN_TYPE_VOID_PTR, gval);
}
strncpy(thisElementType,lastText,
MAX_NAME_LEN);
gval.s = (char *)
malloc(strlen(thisElementType)+1);
strcpy(gval.s, thisElementType);
/* Store the element type of the model */
St_ReplSym(thisSymT_P, "TYPE",
GEN_TYPE_STRING, gval);
/* Mark symbol as used */
St_MarkSym(thisSymT_P, "TYPE");

}

maybe_bracketed_sym_assgn_list
EOL


| DOT_NODESET
{ ErrMsg("Command Not Implemented");
maimed = TRUE;
yyerrok;
}
EOL
| DOT_OP
{ ErrMsg("Command Not Implemented");
maimed = TRUE;
yyerrok;
}
EOL
| DOT_OPTIONS
{ thisSymT_P = capOptionsT_P;
usageTableName = "OPTIONS";
usageType =  USAGE_ID_OPTIONS_TYPE; }
sym_assgn_list_1
EOL
| DOT_PLOT
{ InitOutReqList(); }
plotPrintLine EOL
{ AppendOutReqOperator("plot");
SaveOutReqList(capEndOutputT_P, MakeSequence("out"));
}
| DOT_PRINT
{ thisSymT_P = capOutputT_P; }
{ InitOutReqList(); }
plotPrintLine EOL
{ AppendOutReqOperator("pack");
SaveOutReqList(thisSymT_P, MakeSequence("out"));
}
| DOT_PROBE
{ ErrMsg("Command Not Implemented");
maimed = TRUE;
yyerrok;
}
EOL
| DOT_REF any_literal
{
thisTid = getTermNamed(lastText);
if (!thisTid.val) {
thisTid = newTerm(lastText);
}
setReference(thisTid);
}
EOL
| DOT_SENS
{ ErrMsg("Command Not Implemented");
maimed = TRUE;
yyerrok;
}
EOL

| DOT_SUBCKT any_identifier
{
thisTermCount=0;
thisModelName[0]=0;
strncpy(thisId,"subckt",MAX_NAME_LEN);
strncpy(thisSubCkt,lastText,MAX_NAME_LEN);
pushCircuit(thisSubCkt);
}
connections
{
gr_Id_t tid;
char parent_plus_tok[MAX_NAME_LEN+2];
strcpy(parent_plus_tok,"0");
tid=getTermNamed(parent_plus_tok);
if(!tid.val)
{
tid= newTerm(parent_plus_tok);
}
thisTermList[thisTermCount++]=tid;
}
maybe_subckt_params EOL
{
if(AddSubCkt()) {
maimed = TRUE;
yyerrok;
}
}

| DOT_TF
{ ErrMsg("Command Not Implemented");
maimed = TRUE;
yyerrok;
}
EOL

| DOT_TRAN number_as_double
{
initAnalysis("svtran");
St_ReplSym(capOptionsT_P, "tstep", GEN_TYPE_DOUBLE, globalGval);
}
number_as_double
{ St_ReplSym(capOptionsT_P, "tstop", GEN_TYPE_DOUBLE, globalGval); }
DOT_TRAN_REST
 EOL

| DOT_WATCH
{ InitOutReqList(); }
plotPrintLine EOL
{ thisOutputName = MakeSequence("watch");
AppendOutReqString(thisOutputName, CAP_OBJ_FILENAME);
AppendOutReqOperator("watch");
SaveOutReqList(capOutputT_P, thisOutputName);
}
| DOT_WIDTH
{ ErrMsg("Command Not Implemented");
maimed = TRUE;
yyerrok;
}
EOL
| DOT_TEMP
{ ErrMsg("Command Not Implemented");
maimed = TRUE;
yyerrok;
}
EOL

| DOT_SVHB
{ initAnalysis("svhb"); }
sym_assgn_list_1 EOL
{ closeAnalysis(); }

| DOT_TRAN
{ initAnalysis("svtran"); }
sym_assgn_list_1 EOL
{ closeAnalysis(); }

| DOT_NOISE
{ initAnalysis("noise"); }
sym_assgn_list_1 EOL
{ closeAnalysis(); }

| DOT_TRAN2
{ initAnalysis("svtran2"); }
sym_assgn_list_1 EOL
{ closeAnalysis(); }

| DOT_TRAN3
{ initAnalysis("svtran3"); }
sym_assgn_list_1 EOL
{ closeAnalysis(); }

| DOT_TRAN4
{ initAnalysis("svtran4"); }
sym_assgn_list_1 EOL
{ closeAnalysis(); }

| DOT_TRANT
{ initAnalysis("svtrant"); }
sym_assgn_list_1 EOL
{ closeAnalysis(); }

| DOT_SVTR
{ initAnalysis("svtr"); }
sym_assgn_list_1 EOL
{ closeAnalysis(); }

| DOT_WAVTRAN
{ initAnalysis("wavtran"); }
sym_assgn_list_1 EOL
{ closeAnalysis(); }

| DOT_OUT
{ thisSymT_P = capOutputT_P; }
DOT_OUT_ANY

| DOT_OUT END
{ thisSymT_P = capEndOutputT_P; }

DOT_OUT_ANY
| DOT_OUT HB_
{ thisSymT_P = capOutputT_P; }

DOT_OUT_ANY
| DOT_OUT TRAN_
{ thisSymT_P = capOutputT_P; }

DOT_OUT_ANY
| DOT_OUT VARIABLE_ any_out_identifier
{
generic_t gval;
capOutReq_Pt newReq_P;
newReq_P = (capOutReq_t *)
calloc(1, sizeof(capOutReq_t));
newReq_P->protect = TRUE;
gval.v = (VOIDPTR) newReq_P;
St_DefSym(capOutputVariablesT_P, lastText,
GEN_TYPE_VOID_PTR, gval);
/* ???
St_ReplSym(capOutputVariablesT_P, lastText,
GEN_TYPE_VOID_PTR, gval);
*/
/* 			printf("dot_out variable lastText = %s\n",lastText); */
}
;

/*
* Subsequences for .OUT
*/
DOT_OUT_ANY	:	{ InitOutReqList(); }
outreq_list EOL
{ SaveOutReqList(thisSymT_P, MakeSequence("crt_out"));}
|	DOT_OUT_WRITE
|	DOT_OUT_PLOT
|	DOT_OUT_SYSTEM
|	DOT_OUT_WATCH
;

any_out_identifier :    IDENTIFIER
{ tokType=IDENTIFIER;
globalGval.s = lastText = $<id>1; }
|       VARIABLE
{ tokType=IDENTIFIER;
globalGval.s = lastText = $<id>1; }
|	STRING_LIT
{ tokType=IDENTIFIER;
globalGval.s = lastText = $<id>1; }
;



DOT_OUT_WRITE	:	WRITE
{ InitOutReqList();
thisOutputName = MakeSequence("write");}
outreq_list file_dest EOL
{ AppendOutReqOperator("write");
SaveOutReqList(thisSymT_P, MakeSequence("out"));}
;
DOT_OUT_PLOT 	:	PLOT
{ InitOutReqList();
thisOutputName = MakeSequence("plot");}
outreq_list file_dest EOL
{ AppendOutReqOperator("plot");
SaveOutReqList(thisSymT_P, MakeSequence("out"));}
;
DOT_OUT_SYSTEM	:	SYSTEM
{ InitOutReqList(); }
outreq_list EOL
{ AppendOutReqOperator("system");
SaveOutReqList(thisSymT_P, MakeSequence("out"));}
;
DOT_OUT_WATCH	:	WATCH
{ InitOutReqList();
thisOutputName = "watch"; }
outreq_list file_dest EOL
{ AppendOutReqOperator(thisOutputName);
SaveOutReqList(thisSymT_P, MakeSequence("watch"));}
;

/*
 * Subsequences for .TRAN
 */
DOT_TRAN_REST	:
| number_as_double
{ St_ReplSym(capOptionsT_P, "tstart", GEN_TYPE_DOUBLE, globalGval); }
number_as_double EOL
{ St_ReplSym(capOptionsT_P, "tmax", GEN_TYPE_DOUBLE, globalGval); }
;


/*
 * Subsequences for .OUT
 */

outreq_list	:	outreq_list outreq
|	outreq
;

outreq		:	TERM INT_LIT
{ capOutReq_t outReq;
thisTid = getTermNamed(lastNumberAsText);
if (!thisTid.val)
ParseError("No such terminal.\n");
outReq.type = CAP_OBJ_TERM;
outReq.size = 0;
outReq.obj.tid = thisTid;
outReq.x = NULL;
AppendOutReqList(outReq);
}
|	TERM STRING_LIT
{ capOutReq_t outReq;
thisTid = getTermNamed($<id>2);
if (!thisTid.val)
ParseError("No such terminal.\n");
outReq.type = CAP_OBJ_TERM;
outReq.obj.tid = thisTid;
outReq.size = 0;
outReq.x = NULL;
AppendOutReqList(outReq);
}
|	ELEM STRING_LIT
{ capOutReq_t outReq;
thisNid = getElemNamed($<id>2);
if (!thisNid.val) {
fprintf(stderr,
"No such element. %s\n", $<id>2);
exit(0);
}
outReq.type = CAP_OBJ_NODE;
outReq.obj.nid = thisNid;
outReq.size = 0;
outReq.x = NULL;
AppendOutReqList(outReq);
}
|	outreq_int
|	FLOAT_LIT
{ capOutReq_t outReq;
outReq.type = CAP_OBJ_DOUBLE;
outReq.obj.d = $<dval>1;
outReq.size = 0;
outReq.x = NULL;
AppendOutReqList(outReq);
}
|	FILE_ STRING_LIT
{ AppendOutReqString($<id>2, CAP_OBJ_DATAFILE);
AppendOutReqOperator("read"); }
|	outreq_string
|	ELEMENT STRING_LIT
{ capOutReq_t outReq;
thisNid = getElemNamed($<id>2);
if (!thisNid.val)
ParseError("No such element.\n");
outReq.type = CAP_OBJ_NODE;
outReq.size = 0;
outReq.obj.nid = thisNid;
outReq.x = NULL;
AppendOutReqList(outReq);
}
|	IDENTIFIER
{ AppendOutReqOperator($<id>1); }
|	VARIABLE
{ 
capOutReq_t outReq;
double value;
int type;
generic_t gval;
char string[MAX_STRING_LEN];
/* Check type of variable */
if(St_GetSym(capOptionsT_P,$<id>1,&type, &gval) == ST_SYM_FOUND) {
  if(type == GEN_TYPE_STRING) {
    // If variable is a string, keep it that way
    AppendOutReqString(gval.s, CAP_OBJ_STRING);
  }
  else if(type == GEN_TYPE_DOUBLE) {
    outReq.type = CAP_OBJ_DOUBLE;
    outReq.size = 0;
    outReq.obj.d = value;
    outReq.x = NULL;
    AppendOutReqList(outReq);
  }
  else if(St_GetSymAsDouble(capOptionsT_P,$<id>1,&value) == ST_SYM_FOUND) {
    outReq.type = CAP_OBJ_DOUBLE;
    outReq.size = 0;
    outReq.obj.d = value;
    outReq.x = NULL;
    AppendOutReqList(outReq);
  }
  else
    ParseError("Variable Undefined.");
}
else
ParseError("Variable Undefined.");
}
|	arith_op
{ AppendOutReqOperator(lastText); }
;

arith_op	:	'+'	{ lastText = "+"; }
|	'-'	{ lastText = "-"; }
|	'*'	{ lastText = "*"; }
|	'/'	{ lastText = "/"; }
;

/*
outreq_literal	:	outreq_string
|	outreq_int
;
*/

outreq_string	:	STRING_LIT
{ AppendOutReqString($<id>1, CAP_OBJ_STRING); }
;

outreq_int	:	INT_LIT
{ capOutReq_t outReq;
outReq.type = CAP_OBJ_INT;
outReq.obj.i = $<ival>1;
outReq.size = 0;
outReq.x = NULL;
AppendOutReqList(outReq);
}
;

file_dest	: 	IN_ STRING_LIT
{ AppendOutReqString($<id>2, CAP_OBJ_FILENAME);}
|
{ AppendOutReqString(thisOutputName, CAP_OBJ_FILENAME);}
;



plotPrintLine	:	TRAN_
{
anType = 0;
countV = 0;
spSyn = 1;
int n;
for(n = 0; n < PRINTMAX; n++)
printHeader[n] = NULL;
}
plotPrintSpecsV
|	TRAN_
{
anType = 0;
countV = 0;
spSyn = 1;
int n;
for(n = 0; n < PRINTMAX; n++)
printHeader[n] = NULL;
}
plotPrintSpecsI
|	DC
{
anType = 1;
countV = 0;
spSyn = 1;
int n;
for(n = 0; n < PRINTMAX; n++)
printHeader[n] = NULL;
}
plotPrintSpecsV
|	DC
{
anType = 1;
countV = 0;
spSyn = 1;
int n;
for(n = 0; n < PRINTMAX; n++)
printHeader[n] = NULL;
}
plotPrintSpecsI
;


plotPrintSpecsV	:	plotPrintSpecV
|	plotPrintSpecsV plotPrintSpecV
;


plotPrintSpecV	:	PRINT_V
print_outreq_literal_v ')'
{ lastText = "vt";
AppendOutReqOperator(lastText);
}
;


print_outreq_literal_v  : print_outreq_literal_v print_outreq_v
|  print_outreq_v
;


print_outreq_v	:   INT_LIT
{
capOutReq_t outReq;
thisTid = getTermNamed(lastNumberAsText);
if (!thisTid.val)
ParseError("No such terminal.\n");
outReq.type = CAP_OBJ_TERM;
outReq.size = 0;
outReq.obj.tid = thisTid;
outReq.x = NULL;
outReq.xName = lastNumberAsText;
AppendOutReqList(outReq);
}
|  IDENTIFIER
{ capOutReq_t outReq;
thisTid = getTermNamed($<id>1);
if (!thisTid.val)
ParseError("No such terminal.\n");
outReq.type = CAP_OBJ_TERM;
outReq.obj.tid = thisTid;
outReq.size = 0;
outReq.x = NULL;
outReq.xName = $<id>1;
AppendOutReqList(outReq);
}
;


plotPrintSpecsI	:	plotPrintSpecI
|	plotPrintSpecsI plotPrintSpecI
;


plotPrintSpecI	:	PRINT_I
print_outreq_literal_i ')'
{ lastText = "it";
AppendOutReqOperator(lastText);
}
;


print_outreq_literal_i  : print_outreq_literal_i print_outreq_i
|  print_outreq_i
;


print_outreq_i	:    IDENTIFIER
{ capOutReq_t outReq, terminal_number;
thisNid = getElemNamed($<id>1);
if (!thisNid.val)
ParseError("No such element.\n");
outReq.type = CAP_OBJ_NODE;
outReq.obj.nid = thisNid;
outReq.size = 0;
outReq.x = NULL;
outReq.xName = $<id>1;

terminal_number.type = CAP_OBJ_INT;
terminal_number.size = 0;
terminal_number.obj.tid = thisTid;
terminal_number.x = NULL;

AppendOutReqList(outReq);
AppendOutReqList(terminal_number);
}
;

any_literal	:	FLOAT_LIT
{ tokType=FLOAT_LIT;
globalGval.d = $<dval>1;
lastText = lastNumberAsText; }
|	INT_LIT
{ tokType=INT_LIT;
globalGval.i = $<ival>1;
lastText = lastNumberAsText; }
|	BOOL_LIT
{ tokType=BOOL_LIT;
globalGval.i = $<ival>1;
lastText = lastNumberAsText; }
|	STRING_LIT
{ tokType=STRING_LIT;
globalGval.s = lastText = $<id>1; }
|	IDENTIFIER
{ tokType=IDENTIFIER;
globalGval.s = lastText = $<id>1; }
|	VARIABLE
{ tokType=VARIABLE;
globalGval.s = lastText = $<id>1; }
;

any_literal_but_identifier :	FLOAT_LIT
{ tokType=FLOAT_LIT;
globalGval.d = $<dval>1;
lastText = lastNumberAsText; }
|	INT_LIT
{ tokType=INT_LIT;
globalGval.i = $<ival>1;
lastText = lastNumberAsText;}
|	BOOL_LIT
{ tokType=BOOL_LIT;
globalGval.i = $<ival>1;
lastText = lastNumberAsText;}
|	STRING_LIT
{ tokType=STRING_LIT;
globalGval.s = lastText = $<id>1; }
|	VARIABLE
{ tokType=VARIABLE;
globalGval.s = lastText = $<id>1; }
;

/*
any_literal_but_string :	FLOAT_LIT
{ tokType=FLOAT_LIT;
globalGval.d = $<dval>1;
lastText = lastNumberAsText; }
|	INT_LIT
{ tokType=INT_LIT;
globalGval.i = $<ival>1;
lastText = lastNumberAsText;}
|	BOOL_LIT
{ tokType=BOOL_LIT;
globalGval.i = $<ival>1;
lastText = lastNumberAsText;}
|	IDENTIFIER
{ tokType=IDENTIFIER;
globalGval.s = lastText = $<id>1; }
|	VARIABLE
{ tokType=VARIABLE;
globalGval.s = lastText = $<id>1; }
;
*/


any_number	:	FLOAT_LIT
{ tokType=FLOAT_LIT;
globalGval.d = $<dval>1;
lastText = lastNumberAsText; }
|	INT_LIT
{ tokType=INT_LIT;
globalGval.i = $<ival>1;
lastText = lastNumberAsText; }
|	BOOL_LIT
{ tokType=BOOL_LIT;
globalGval.i = $<ival>1;
lastText = lastNumberAsText; }
|	VARIABLE
{ tokType=VARIABLE;
globalGval.s = lastText = $<id>1; }
;


literal		:	FLOAT_LIT
{ tokType=FLOAT_LIT;
globalGval.d = $<dval>1;
lastText = lastNumberAsText; }
|	INT_LIT
{ tokType=INT_LIT;
globalGval.i = $<ival>1;
lastText = lastNumberAsText; }
|	BOOL_LIT
{ tokType=BOOL_LIT;
globalGval.i = $<ival>1;
lastText = lastNumberAsText; }
|	STRING_LIT
{ tokType=STRING_LIT;
globalGval.s = lastText = $<id>1; }
;

any_identifier	:	IDENTIFIER
{ tokType=IDENTIFIER;
globalGval.s = lastText = $<id>1; }
;


/*
* spice element parsing
*/
spice_element	:	A_ELEMENT
{ ParseError("Unknown element"); }
|	B_ELEMENT
{ ParseError("Unsupported element"); }
|	C_ELEMENT
{ initElement($<id>1,"capacitor");}
connection connection maybe_model
{ strncpy(thisId,"c", MAX_NAME_LEN); }
sym_assgn_value sym_assgn_list EOL
{ closeElement(); }
|	D_ELEMENT
{
initElement($<id>1,"d");}
connection connection maybe_model
maybe_area_sym_assgn_list EOL
{ closeElement(); }
|	E_ELEMENT
{ initElement($<id>1,"vcvs");}
connection connection
e_or_g_arguments EOL
{ closeElement();
insertPolyParam($<id>1) ;
}
|	F_ELEMENT
{ initElement($<id>1,"f");}
connection connection f_or_h_arguments EOL
{ closeElement(); }
|	G_ELEMENT
{ initElement($<id>1,"g");}
connection connection e_or_g_arguments EOL
{ closeElement(); }
|	H_ELEMENT
{ initElement($<id>1,"h");}
connection connection f_or_h_arguments EOL
{ closeElement(); }
|	I_ELEMENT
{ initElement($<id>1,"isource");}
connection connection
maybe_idc maybe_iac maybe_itransient EOL
{ closeElement(); }
| 	J_ELEMENT
{initElement($<id>1,"j");}
connection connection connection need_model maybe_area EOL
{
checkModel();
if(!strcmp(thisModelType,"njf"))
  strcpy(thisElementType,"jfetn");
else if(!strcmp(thisModelType,"jfetn"))
  strcpy(thisElementType,"jfetn");
else if(!strcmp(thisModelType,"pjf"))
  strcpy(thisElementType,"jfetp");
else if(!strcmp(thisModelType,"jfetp"))
  strcpy(thisElementType,"jfetp");
else
  ParseError("Model and Level not supported.");

closeElement();
}
|	K_ELEMENT
{ ParseError("Not supported in Spice format yet"); }
|	L_ELEMENT
{ initElement($<id>1,"inductor");}
connection connection maybe_model
{ strncpy(thisId,"l", MAX_NAME_LEN); }
sym_assgn_value sym_assgn_list EOL
{ closeElement(); }
|	M_ELEMENT
{ initElement($<id>1,"mos");}
connection connection connection connection
need_model
maybe_area_sym_assgn_list EOL
{

checkModel();
checkLevel();
if(!strcmp(thisModelType,"nmos") && (thisElementLevel==1))
  strcpy(thisElementType,"mosn1");
else if(!strcmp(thisModelType,"pmos") && (thisElementLevel==1))
  strcpy(thisElementType,"mosp1");
else if(!strcmp(thisModelType,"nmos") && (thisElementLevel==2))
  strcpy(thisElementType,"mosn2");
else if(!strcmp(thisModelType,"pmos") && (thisElementLevel==2))
  strcpy(thisElementType,"mosp2");
else if(!strcmp(thisModelType,"nmos") && (thisElementLevel==3))
  strcpy(thisElementType,"mosn3");
else if(!strcmp(thisModelType,"pmos") && (thisElementLevel==3))
  strcpy(thisElementType,"mosp3");
else if(!strcmp(thisModelType,"nmos") &&(thisElementLevel==8))
  strcpy(thisElementType,"mosnbsim3");
else if(!strcmp(thisModelType,"pmos") &&(thisElementLevel==8))
  strcpy(thisElementType,"mospbsim3");
else if(!strcmp(thisModelType,"nmos") &&(thisElementLevel==9))
  strcpy(thisElementType,"mosn9");
else if(!strcmp(thisModelType,"pmos") &&(thisElementLevel==9))
  strcpy(thisElementType,"mosp9");
else if(!strcmp(thisModelType,"nmos") &&(thisElementLevel==14))
  strcpy(thisElementType,"mosnbsim4");
/* There is no PMOS for BSIM4
else if(!strcmp(thisModelType,"pmos") &&(thisElementLevel==14))
  strcpy(thisElementType,"pmos14"); */
else
  ParseError("Model and Level not supported.");

closeElement();

}

|	N_ELEMENT
{ ParseError("Unknown element"); }
|	O_ELEMENT
{ ParseError("Unknown element"); }
|	P_ELEMENT
{ ParseError("Unknown element"); }
|	Q_ELEMENT
{initElement($<id>1,"q");}
connection connection connection connection need_model
need_area EOL
{
checkModel();
if((!strcmp(thisModelType,"npn"))||(!strcmp(thisModelType,"bjtnpn")))
  strcpy(thisElementType,"bjtnpn");
else if((!strcmp(thisModelType,"pnp"))||(!strcmp(thisModelType,"bjtpnp")))
  strcpy(thisElementType,"bjtpnp");
else
  ParseError("Model and Level not supported.");
closeElement();
}

|	R_ELEMENT
{ initElement($<id>1,"res");}
connection connection maybe_model
{ strncpy(thisId,"r", MAX_NAME_LEN); }
sym_assgn_value sym_assgn_list EOL
{ closeElement(); }
|	S_ELEMENT
{ ParseError("Unsupported element"); }
|	T_ELEMENT
{ ParseError("Unsupported element"); }
|	U_ELEMENT
{ ParseError("Unknown element"); }
|	V_ELEMENT
{ initElement($<id>1,"v");}
connection connection
maybe_vdc maybe_vac maybe_vtransient EOL
{ closeElement(); }
|	W_ELEMENT
{ ParseError("Unsupported element"); }
|	X_ELEMENT
{ initElement($<id>1,"x"); }
connections subname
{
gr_Id_t tid;
char parent_plus_tok[MAX_NAME_LEN+2];
strcpy(parent_plus_tok,"0");
tid=getTermNamed(parent_plus_tok);
if(!tid.val)
{
tid= newTerm(parent_plus_tok);
}
thisTermList[thisTermCount++]=tid;
}
sym_assgn_list EOL
{ closeElement(); }
|	Y_ELEMENT
{ ParseError("Unknown element"); }
|	Z_ELEMENT
{ initElement($<id>1,"z");}
connection connection connection connection
need_model need_area EOL
{
checkModel();
checkLevel();
if(!strcmp(thisModelType,"nmf") && (thisElementLevel==1))
  strcpy(thisElementType,"mesfetc");
else if(!strcmp(thisModelType,"nmf") && (thisElementLevel==2))
  strcpy(thisElementType,"mesfetm");
else if(!strcmp(thisModelType,"nmf") && (thisElementLevel==3))
  strcpy(thisElementType,"tmesfetc");
else
  ParseError("Model and Level not supported.");

closeElement();
}
;


/*
* element argument subsequences
*/
e_or_g_arguments:	poly_1
| 	connection connection
{ strncpy(thisId, "gain", MAX_NAME_LEN); }
sym_assgn_list
;



f_or_h_arguments:	poly_2
|	 connection connection
{ strncpy(thisId, "gain", MAX_NAME_LEN); }
sym_assgn_list
;

subname         :       IDENTIFIER
{
strncpy(thisModelName,$<id>1,MAX_NAME_LEN);
}
;

/*
* I_ELEMENT and V_ELEMENT subsequences
*/

maybe_vdc	:
|	DC
{
strcpy(thisElementType,"vsource");
strncpy(thisId, "vdc", MAX_NAME_LEN); }
sym_assgn_single_value
|	{strcpy(thisElementType,"vsource");
strncpy(thisId, "vdc", MAX_NAME_LEN);
}
sym_assgn_single_value

;

maybe_vac	:
|	AC
{
strcpy(thisElementType,"vsource");
thisParameter = 0;
strncpy(parameterList[thisParameter++].s, "vac",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "f",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "phase",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "delay",
MAX_NAME_LEN);
pCount=thisParameter;
thisParameter=0;
}
sym_assgn_list_2
;

maybe_vtransient	:
|	PULSE 
{ generic_t gval;
gval.i = SPICE_SOURCE_PULSE;
St_ReplSym(thisSymT_P,"transientType",GEN_TYPE_INT,gval);
strcpy(thisElementType,"vpulse");
thisParameter = 0;
strncpy(parameterList[thisParameter++].s, "v1",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "v2",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "td",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "tr",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "tf",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "pw",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "per",
MAX_NAME_LEN);
pCount=thisParameter;
thisParameter=0;
}
maybe_bracketed_sym_assgn_list_2
|	EXP 
{ generic_t gval;
gval.i = SPICE_SOURCE_EXP;
St_ReplSym(thisSymT_P,"transientType",GEN_TYPE_INT,gval);
strcpy(thisElementType,"vexp");
thisParameter = 0;
strncpy(parameterList[thisParameter++].s, "v1",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "v2",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "trd",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "trc",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "tfd",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "tfc",
MAX_NAME_LEN);
pCount=0;
thisParameter=0;
}
maybe_bracketed_sym_assgn_list_2
|	SIN 
{ generic_t gval;
gval.i = SPICE_SOURCE_SIN;
St_ReplSym(thisSymT_P,"transientType",GEN_TYPE_INT, gval);
strcpy(thisElementType,"vsin");
thisParameter = 0;
strncpy(parameterList[thisParameter++].s, "vo",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "va",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "freq",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "td",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "alpha",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "theta",
MAX_NAME_LEN);
pCount=thisParameter;
thisParameter=0;
}
maybe_bracketed_sym_assgn_list_2
|	PWL
{ generic_t gval;
gval.i = SPICE_SOURCE_PWL;
St_ReplSym(thisSymT_P,"transientType",GEN_TYPE_INT,gval);
strcpy(thisElementType,"vpwl");
thisParameter = 0;
strncpy(parameterList[thisParameter++].s, "waveform", MAX_NAME_LEN);
pCount=thisParameter;
thisParameter=0;
}
maybe_bracketed_sym_assgn_list_2
|	SFFM 
{ generic_t gval;
gval.i = SPICE_SOURCE_SFFM;
St_ReplSym(thisSymT_P,"transientType",GEN_TYPE_INT,gval);
strcpy(thisElementType,"vsfFm");
thisParameter = 0;
strncpy(parameterList[thisParameter++].s, "vo",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "va",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "fcarrier",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "mdi",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "fsignal",
MAX_NAME_LEN);
pCount=thisParameter;
thisParameter=0;
}
maybe_bracketed_sym_assgn_list_2
;



maybe_idc	:
|	DC
{
strcpy(thisElementType,"isource");
strncpy(thisId, "idc", MAX_NAME_LEN); }
sym_assgn_single_value
|	{strcpy(thisElementType,"isource");
strncpy(thisId, "idc", MAX_NAME_LEN);
}
sym_assgn_single_value

;

maybe_iac	:
|	AC
{
strcpy(thisElementType,"isource");
thisParameter = 0;
strncpy(parameterList[thisParameter++].s, "iac",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "f",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "phase",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "delay",
MAX_NAME_LEN);
pCount=thisParameter;
thisParameter=0;
}
sym_assgn_list_2
;

maybe_itransient	:
|	PULSE
{ generic_t gval;
gval.i = SPICE_SOURCE_PULSE;
St_ReplSym(thisSymT_P,"transientType",GEN_TYPE_INT,gval);
strcpy(thisElementType,"ipulse");
thisParameter = 0;
strncpy(parameterList[thisParameter++].s, "v1",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "v2",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "td",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "tr",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "tf",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "pw",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "per",
MAX_NAME_LEN);
pCount=thisParameter;
thisParameter=0;
}
maybe_bracketed_sym_assgn_list_2
|	EXP
{ generic_t gval;
gval.i = SPICE_SOURCE_EXP;
St_ReplSym(thisSymT_P,"transientType",GEN_TYPE_INT,gval);
strcpy(thisElementType,"iexp");
thisParameter = 0;
strncpy(parameterList[thisParameter++].s, "v1",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "v2",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "trd",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "trc",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "tfd",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "tfc",
MAX_NAME_LEN);
pCount=0;
thisParameter=0;
}
maybe_bracketed_sym_assgn_list_2
|	SIN
{ generic_t gval;
gval.i = SPICE_SOURCE_SIN;
St_ReplSym(thisSymT_P,"transientType",GEN_TYPE_INT, gval);
strcpy(thisElementType,"isin");
thisParameter = 0;
strncpy(parameterList[thisParameter++].s, "io",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "ia",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "freq",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "td",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "alpha",
MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "theta",
MAX_NAME_LEN);
pCount=thisParameter;
thisParameter=0;
}
maybe_bracketed_sym_assgn_list_2
|	PWL
{ generic_t gval;
gval.i = SPICE_SOURCE_PWL;
St_ReplSym(thisSymT_P,"transientType",GEN_TYPE_INT,gval);
strcpy(thisElementType,"ipwl");
thisParameter = 0;
strncpy(parameterList[thisParameter++].s, "waveform", MAX_NAME_LEN);
pCount=thisParameter;
thisParameter=0;
}
maybe_bracketed_sym_assgn_list_2
|	SFFM
{ generic_t gval;
gval.i = SPICE_SOURCE_SFFM;
St_ReplSym(thisSymT_P,"transientType",GEN_TYPE_INT,gval);
strcpy(thisElementType,"isfFm");
thisParameter = 0;
strncpy(parameterList[thisParameter++].s, "io", MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "ia", MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "fcarrier", MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "mdi", MAX_NAME_LEN);
strncpy(parameterList[thisParameter++].s, "fsignal", MAX_NAME_LEN);
pCount=thisParameter;
thisParameter=0;
}
maybe_bracketed_sym_assgn_list_2
;


/*
* table name for MOSFET device type

maybe_table	:
file_dest
{ strncpy(thisId, $<id>1, MAX_NAME_LEN); }
sym_assgn_single_value
{ strncpy(thisId, $<id>1, MAX_NAME_LEN); }
sym_assgn_single_value
;
*/



/*
* subcircuit paramenters
*/

maybe_subckt_params	:
|	PARAMS_
sym_assgn_list
;

/*
* non-traditional spice element parsing
*/
additional_element :	special_element
{ initElement( $<id>1, $<id>1);
}
':'
{ int length;
length = strlen(thisElementName);
thisElementName[length] = ':';
thisElementName[length+1] = 0;
}
any_literal
{ int length, i, length1;
length = strlen(thisElementName);
length1 = strlen(lastText);
for(i = 0; i < length1; i++)
{
thisElementName[length+i] = lastText[i];
}
thisElementName[length+length1] = 0;
thisTermCount = 0;
thisModelName[0] = 0;
}
connections
sym_assgn_list maybe_poly_1 EOL
{ closeElement(); }
;

special_element	:	A_ELEMENT
|	B_ELEMENT
|	C_ELEMENT
|	D_ELEMENT
|	E_ELEMENT
|	F_ELEMENT
|	G_ELEMENT
|	H_ELEMENT
|	I_ELEMENT
|	J_ELEMENT
|	K_ELEMENT
|	L_ELEMENT
|	M_ELEMENT
|	N_ELEMENT
|	O_ELEMENT
|	P_ELEMENT
|	Q_ELEMENT
|	R_ELEMENT
|	S_ELEMENT
|	T_ELEMENT
|	U_ELEMENT
|	V_ELEMENT
|	W_ELEMENT
|	X_ELEMENT
|	Y_ELEMENT
|	Z_ELEMENT
;

// Allow an element without a connection
//connections	: 	connection
connections   :
|	connections connection
;

connection	:	any_literal
{
gr_Id_t tid;
char parent_plus_tok[MAX_NAME_LEN+2];
strncpy(parent_plus_tok,lastText, MAX_NAME_LEN);
/* Use the old way to add terminals mainly to avoid changing the parser. */
tid = getTermNamed(parent_plus_tok);
if (!tid.val) {
tid = newTerm(parent_plus_tok);
}
thisTermList[thisTermCount++] = tid;
}
;

maybe_model	:
|       IDENTIFIER
{  generic_t gval ;
globalGval.s = lastText = $<id>1;
strncpy(thisModelName, lastText, MAX_NAME_LEN);
gval.s = (char *) malloc(strlen(thisModelName)+1);
strcpy(gval.s,thisModelName);
St_DefSym(thisSymT_P,"model",GEN_TYPE_STRING,gval);
St_MarkSym(thisSymT_P,"model");
}
;

need_model	:      IDENTIFIER
{ generic_t gval;
globalGval.s = lastText = $<id>1;
strncpy(thisModelName, lastText, MAX_NAME_LEN);
gval.s = (char *) malloc(strlen(thisModelName)+1);
strcpy(gval.s,thisModelName);
St_DefSym(thisSymT_P,"model",GEN_TYPE_STRING,gval);
St_MarkSym(thisSymT_P,"model");
}
;

maybe_subname   :       /* could be empty... */
{
strncpy(thisSubCkt,"\0",MAX_NAME_LEN);
}
|       IDENTIFIER
{
strncpy(thisSubCkt,"\0",MAX_NAME_LEN);
}

/* By now, we'll not check nesting. Later we'll check it in the
* network classes. */
/*                         { */
/*                          strncpy(thisSubCkt,$<id>1,MAX_NAME_LEN); */
/*                          getParentCircuitName(thisParent); */
/*                          if  (strncmp(thisSubCkt,thisParent,MAX_NAME_LEN)) */
/*                                ParseError("SubCkt nesting error"); */
/*                         } */


/*
* As used by .OPTIONS to allow for an identifier such as NOPAGE to appear
St_GetSymAsFloat(capOptionsT_P, thisVariableName, globalGval);
it_GetSymAsFloat(capOptionsT_P, thisVariableName, globalGval);
* on its own
*/
;

sym_assgn_list_1:	/* empty */
|	sym_assgn_list_1 sym_assgn_1
;

sym_assgn_1	:	PARAMETER
{ strncpy(thisId, $<id>1, MAX_NAME_LEN);
thisId[MAX_NAME_LEN] = 0; }
sym_assgn_vector
|	PARAMETER
{ strncpy(thisId, $<id>1, MAX_NAME_LEN);
thisId[MAX_NAME_LEN] = 0;
}
sym_assgn_single_value
|	PARAMETER
{ strncpy(thisId, $<id>1, MAX_NAME_LEN);
thisId[MAX_NAME_LEN] = 0;
}
sym_assgn_expression
|	PARAMETER
{ strncpy(thisId, $<id>1, MAX_NAME_LEN);
thisId[MAX_NAME_LEN] = 0; }
sym_assgn_sweep
|	IDENTIFIER
{ generic_t gval;
gval.i = 1;
strncpy(thisId, $<id>1, MAX_NAME_LEN);
thisId[MAX_NAME_LEN+1] = 0;
St_ReplSym(thisSymT_P, thisId, GEN_TYPE_INT, gval); }
|	NOECHO
{ generic_t gval;
gval.i = 1;
spiceTable.echo = 0;
St_ReplSym(thisSymT_P, "noecho", GEN_TYPE_INT, gval); }
;




/* assign values to symbols in dotted_commands */
maybe_bracketed_sym_assgn_list_2:
    '(' sym_assgn_list_2 ')'
|	sym_assgn_list_2 
;

sym_assgn_list_2:	/* empty */
|	sym_assgn_list_2 sym_assgn_2
;

sym_assgn_2	:	{ if(pCount|| (thisParameter<pCount))
strncpy(thisId, parameterList[thisParameter++].s,
MAX_NAME_LEN);
else
ParseError("Unrecognized entry");
}
sym_assgn_single_value
;




/*
* symbol assignment sequences
*/

maybe_bracketed_sym_assgn_list: '(' sym_assgn_list ')'
|	'(' sym_assgn_list
|	sym_assgn_list ')'
|	sym_assgn_list 
;

maybe_area_sym_assgn_list :
	{strncpy(thisId,"area", MAX_NAME_LEN); }
	sym_assgn_value
|	sym_assgn_list sym_assgn
;

maybe_area :
|	{strncpy(thisId,"area", MAX_NAME_LEN); }
	sym_assgn_single_value
;

need_area :
	{strncpy(thisId,"area", MAX_NAME_LEN); }
	sym_assgn_single_value
;

sym_assgn_list	:
|	sym_assgn_list sym_assgn
;

sym_assgn	:	PARAMETER
{ strncpy(thisId, $<id>1, MAX_NAME_LEN);
thisId[MAX_NAME_LEN+1] = 0; }
	sym_assgn_value
;


sym_assgn_value :
|       sym_assgn_vector
|	sym_assgn_single_value
|	sym_assgn_sweep
|	sym_assgn_expression
|	sym_assgn_vector_of_vectors
;


sym_assgn_expression : 	EXPRESSION
{ generic_t gval;
/* Add expression to list, evaluate it, and set the current of the symbol. */
gval.d = v_addExpression(($<id>1), usageVal, usageType, thisId);
St_ReplSym(thisSymT_P,thisId,GEN_TYPE_DOUBLE,gval);
}
;
sym_assgn_single_value
: 	any_literal_but_identifier
{ generic_t gval; char *s;
switch(tokType)
{
case INT_LIT:
case BOOL_LIT:
gval = globalGval;
St_ReplSym(thisSymT_P,thisId,GEN_TYPE_INT,gval);

break;
case FLOAT_LIT:
gval = globalGval;
St_ReplSym(thisSymT_P,thisId,GEN_TYPE_DOUBLE,gval);
break;
case STRING_LIT:
{
int length;
length = strlen(lastText);
s = (char *) malloc(length + 1);
strcpy(s, lastText);
gval.s = s;
St_ReplSym(thisSymT_P,thisId,GEN_TYPE_STRING,gval);

break;
}
case VARIABLE:
strncpy(thisVariableName,globalGval.s,MAX_NAME_LEN);
if((gval.v = (VOIDPTR) LookupSym(capOptionsT_P,
thisVariableName))== ST_SYM_FOUND)
ParseError("Variable Undefined.");
if(thisSymT_P == capOptionsT_P)
{
AddVariableUsage(thisId, thisVariableName);
St_GetSymAsDouble(capOptionsT_P,
thisVariableName, &(gval.d));
St_ReplSym(thisSymT_P, thisId,
GEN_TYPE_DOUBLE, gval);
}
else if(usageType == USAGE_ID_SYMBOL_TABLE)
{
St_GetSymAsDouble(capOptionsT_P,
thisVariableName, &(gval.d));
St_ReplSym(thisSymT_P, thisId,
GEN_TYPE_DOUBLE, gval);
AddSymbolTableUsage(thisVariableName,
usageTableName, thisId);
}
else if(usageType == USAGE_ID_ANALYSIS_TYPE)
{
St_GetSymAsDouble(capOptionsT_P,
thisVariableName, &(gval.d));
St_ReplSym(thisSymT_P, thisId,
GEN_TYPE_DOUBLE, gval);
AddAnalysisUsage(thisVariableName,
usageTableName, thisId);
}
else
St_ReplSym(thisSymT_P, thisId,
GEN_TYPE_VARIABLE_USAGE, gval);
break;
default:
ParseError("Syntax Error.");
}
}
;

sym_assgn_sweep :	'\\'
sym_assgn_sweep_value
{ sweepInitial_temp = globalGval.d; }
sym_assgn_sweep_value
{ sweepFinal_temp = globalGval.d; }
sym_assgn_sweep_value
{ sweepStep_temp = globalGval.d; }
'\\'
{ v_addSweep(sweepInitial_temp, sweepFinal_temp,
sweepStep_temp, usageVal, usageType, thisId);
globalGval.d = sweepInitial_temp;
St_ReplSym(thisSymT_P, thisId, GEN_TYPE_DOUBLE,
globalGval);
}
;


sym_assgn_sweep_value
: 	any_number
{ generic_t gval;
switch(tokType)
{
case INT_LIT:
case BOOL_LIT:
globalGval.d = (double) $<ival>1;
/* Before was:
Coerce(GEN_TYPE_INT,GEN_TYPE_DOUBLE,&globalGval);
*/
break;
case FLOAT_LIT:
globalGval.d = $<dval>1;
break;
case VARIABLE:
strncpy(thisVariableName,$<id>1,MAX_NAME_LEN);
if((gval.v = (VOIDPTR) LookupSym(capOptionsT_P,
thisVariableName))== ST_SYM_FOUND)
ParseError("Variable Undefined.");
else
St_ReplSym(thisSymT_P, sweepId,
GEN_TYPE_VARIABLE_USAGE, gval);
break;
default:
ParseError("Syntax Error.");
}
}
;
number_as_double	:	any_number
{
double tmp;
int i;
generic_t gval;
switch(tokType)
{
case INT_LIT:
case BOOL_LIT:
i = $<ival>1;
tmp = (double) i;
globalGval.d = tmp;
globalType = GEN_TYPE_DOUBLE;
break;
case FLOAT_LIT:
globalGval.d = $<dval>1;
break;
case VARIABLE: /* Note variable usage not recorded */
strncpy(thisVariableName,$<id>1,MAX_NAME_LEN);
if((gval.v = (VOIDPTR) LookupSym(capOptionsT_P,
thisVariableName))== ST_SYM_FOUND)
ParseError("Variable Undefined.");
St_GetSymAsDouble(capOptionsT_P, thisVariableName,
&(globalGval.d));
break;
}
}
;

sym_assgn_vector_of_vectors:	'<'
{
noVectors = 0;
/* 			  printf("A");  */
}
sym_assgn_vectors '>'
{ int i, type;
generic_t gval;
/* 			  printf("B"); */
gval.dvv.dv_P = (genericDoubleV_Pt) malloc(noVectors
* sizeof(genericDoubleV_t));
gval.dvv.size = noVectors;
for(i = 0; i < noVectors; i++)
gval.dvv.dv_P[i] = doubleVector_P[i];
type = GEN_TYPE_DOUBLE_VECTOR_OF_VECTORS;
St_ReplSym(thisSymT_P, thisId, type, gval);
noVectors = 0;
}
;

sym_assgn_vectors:	sym_assgn_vector_part
|	sym_assgn_vectors maybe_comma
sym_assgn_vector_part
;

sym_assgn_vector:	'<' sym_assgn_vector_body '>'
{ int type = GEN_TYPE_DOUBLE_VECTOR;
St_ReplSym(thisSymT_P, thisId, type, globalGval);
noVectorElements = 0;
}
;

sym_assgn_vector_part:	'<' sym_assgn_vector_body '>'
{
doubleVector_P[noVectors++] = globalGval.dv;
/* 			  printf("D"); */
noVectorElements = 0;
}
;


sym_assgn_vector_body:	{ noVectorElements = 0; }
sym_assgn_vector_elements
{ int i;
globalGval.dv.v = (double *) malloc(noVectorElements
* sizeof(double));
globalGval.dv.size = noVectorElements;
for(i = 0; i < noVectorElements; i++)
globalGval.dv.v[i] = doubleVector[i];
}
;

sym_assgn_vector_elements : sym_assgn_vector_element
|   	sym_assgn_vector_elements
/* maybe_eol ',' maybe_eol */
maybe_comma
sym_assgn_vector_element
;

sym_assgn_vector_element
:	number_as_double
{ doubleVector[noVectorElements++] = globalGval.d; }

|	STRING_LIT
{
if(noVectorElements && vectorType != GEN_TYPE_STRING)
ParseError("Vector of mixed type");
else
{
char *tmps;
int length;
vectorType = GEN_TYPE_STRING;
length = strlen($<id>1);
stringVector_P[noVectorElements++] = tmps =
(char *) malloc((length+1) * sizeof(char));
strcpy(tmps,$<id>1);
tmps[length] = 0;
}
}
;



/*
* POLY_1 sequences
*
* polynomial sequences for voltage controlled elements
*/

maybe_poly_1	:  /* empty */
|  poly_1
;

poly_1		:	POLY_ '(' INT_LIT ')'
{ number_poly_coeffs = 0;
listType = GR_ID_TERM_TYPE;
noElementsInList = 0;
polydimension = $<ival>3;
/* printf("polydimension = %d\n", polydimension); */
if(polydimension > (GR_MAX_CONN /2))
ParseError("Too many Dimensions.\n");
number_poly_coeffs = -2 * polydimension; }
poly_lit poly_lit poly_list

;


poly_list	:	/* empty */
|	poly_list poly_lit
;

poly_lit	:	literal
{
if (number_poly_coeffs >= 0)
{
/* So we must be dealing with polynomial
* coefficients */
if(number_poly_coeffs > ST_MAX_VECTOR_LENGTH)
{
ParseError("Too many polynomial coefficients.\n");
break;
}
switch(tokType) {
case INT_LIT:
case BOOL_LIT:
poly_coeff[number_poly_coeffs++]
= (double) $<ival>1;
break;
case FLOAT_LIT:
poly_coeff[number_poly_coeffs++]
= $<dval>1;
break;
case STRING_LIT:
ParseError("Cannot handle strings as polynomial coefficient \n");
break;
}
}
else
{
/* So we must be dealing with terminals */
switch(tokType)
{
case INT_LIT:
case BOOL_LIT:
case STRING_LIT:
{ gr_Id_t tid;
tid = getTermNamed(lastText);
if (!tid.val) {
tid = newTerm(lastText);
}
polyList[noElementsInList++] = tid;
number_poly_coeffs++;
}
break;
default:
ParseError("Error in reading polynomial");
}
}
}
;


/*
* POLY_2 sequences
*
* polynomial sequences for current controlled elements
*/
poly_2		:	POLY_ INT_LIT
{ number_poly_coeffs = 0;
listType = GR_ID_NODE_TYPE;
noElementsInList = 0;
polydimension = $<ival>2;
/* printf("polydimension = %d\n", polydimension); */
if(polydimension > (GR_MAX_CONN /2))
ParseError("Too many Dimensions.\n");
number_poly_coeffs = polydimension; }
poly_2_lit poly_2_list
;


poly_2_list	:	/* empty */
|	poly_list poly_lit
;

poly_2_lit	:	literal
{
if (number_poly_coeffs >= 0)
{
/* So we must be dealing with polynomial
* coefficients */
if(number_poly_coeffs > ST_MAX_VECTOR_LENGTH)
{
ParseError("Too many polynomial coefficients.\n");
break;
}
switch(tokType) {
case INT_LIT:
case BOOL_LIT:
poly_coeff[number_poly_coeffs++]
= (double) $<ival>1;
break;
case FLOAT_LIT:
poly_coeff[number_poly_coeffs++]
= $<dval>1;
break;
case STRING_LIT:
ParseError( "Cannot handle strings as polynomial coefficient \n");
break;
}
}
else
{
/* So we must be dealing with nodes */
/* It only makes sense if we are talking about
* two terminal nodes.  We cannot check this at
* this stage. Perhaps checking will be done
* someday. */
switch(tokType)
{
case INT_LIT:
case BOOL_LIT:
case STRING_LIT:
{ ErrMsg("List of nodes not supported yet.");
maimed = TRUE;
yyerrok;
}
/* 			        { gr_Id_t eid; */
/* 			          eid = getElemNamed(lastText); */
/* 				  if (!eid.val) { */
/* 				    eid = newElem(lastText); */
/* 				  } */
/* 			          polyList[noElementsInList++] = eid; */
/* 			          number_poly_coeffs++; */
/* 			          break; */
/* 				} */
}
}
}
;

%%


