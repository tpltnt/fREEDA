/*
 * flex definitions for netlist scanner
 *
 * Some Syntax notes:
 *
 * string
 *    A string is defined as any printable character plus newline between
 *    quotes (").  If a quote is to be included it must be escaped eg (\").
 */

/* Some versions of flex do not recognize %option
   If so use the flex command line options -1 -Ce -Cm */
/* %option caseless */
/* %option ecs */
/* %option meta-ecs */

%{

#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parser.h"
#include "tr2_parse.tab.h"

#include "environment.h"

#define	MAX_FILENAME_LEN MAX_NAME_LEN /* maximum file name length	*/
#define IMPL_NAME_LEN 20	/* max length of implementation names	*/
#define MAX_SPICE_STRING_LEN 64	/* max length of a string		*/
#define MAX_SPICE_LINE_LEN 132	/* max length of a line			*/
#define SPICE_TITLE_LEN 80	/* max length of title			*/

#define MAX_INCLUDE_DEPTH 128
YY_BUFFER_STATE include_stack[MAX_INCLUDE_DEPTH];
FILE* include_files[MAX_INCLUDE_DEPTH];
int include_stack_mylineno[MAX_INCLUDE_DEPTH];
int include_stack_ptr = 0;
int no_print_line = 0; // skip printing a line

#define RET(token)	return(token);

#define MIN(x,y) ((x) < (y) ? (x) : (y))

#ifdef DEBUG
#undef RET
#define RET(token)	 if(token!=EOL)\
            fprintf(stderr, "mylineno = %d\t \ttoken = %d\n", \
              mylineno,token);\
		  else\
            fprintf(stderr,"\nmylineno = %d\teol \ttoken = %d\n",mylineno - 1,token);\
		return(token);
#endif

  extern char 	yylvalStr[YYMAXDEPTH][MAX_TOK_LEN + 1], /* text of yylval.id */
    *lastNumberAsText;
  /** lastNumberAsText[MAX_TOK_LEN + 1]; **/
  static int	thisDepth = 0, 	/* token index */
    mylineno = 0;
  int     	idIndex = 0;	/* used to index tokens */
  void copy_yytext(void);
  void copy_yytext_to_white(void);
  double copy_yytext_return_base(void);
  void copy_yytext_string(void);
  void str2lower(char *s);
  int PrintCurrentLine(IN FILE *file_F, OUT char **b);
%}

/************************************************************************
Note: the following lines go between the lines

    YY_FATAL_ERROR( "read() in flex scanner failed" ); \

        and

  for(i = result, b = buf; i ; i--)\

in #define YY_INPUT below.

#ifdef LEX_DEBUG
  printf("result = %d \n",result);\
  printf("buf    =\n");\
  for(i = result, b = buf; i ; i--) putchar(*(b++));;\
  printf("\n");\ 
#endif

************************************************************************/


%{
/*
 * In debug mode redefine YY_INPUT to print out what is read
 */
#ifdef LEX_DEBUG
#undef YY_INPUT
#define YY_INPUT(buf,result,max_size) \
  { \
  int i; char *b;\
  if ( (result = read( fileno(yyin), (char *) buf, max_size )) < 0 ) \
    YY_FATAL_ERROR( "read() in flex scanner failed" ); \
  printf("result = %d \n",result);\
  printf("buf    =\n");\
  for(i = result, b = buf; i ; i--) putchar(*(b++));;\
  printf("\n");\
  }
#endif

#define PRINT_LINE \
  { \
  if (spiceTable.echo) \
    if(no_print_line) \
      {no_print_line = 0;} \
    else \
      { \
      char *b; \
      PrintCurrentLine(output_F, &b); \
      } \
  }


int yywrap(void)
  {
  return(1);
  }
%}

lf		[\x0d]*[\n]
alpha		[A-Za-z_]
alpha_or_dollar_or_	[A-Za-z$_]
alphanumeric		[0-9A-Za-z_.]
white		[\x09\x0B\x0C\x20]
digits		([0-9])+
float_core	(([0-9])*\.([0-9])+)|(([0-9])+\.([0-9])*)
predigit	"-"|"+"
expon		[Ee]({predigit}?){digits}
scale_core	{predigit}?({float_core}|{digits})

%x SPICE CONTINUATION INCLUDE LIBRARY

%%

<INITIAL>([^\n])*	{ /* Save first line */
			strncpy(spiceTitle,yytext,SPICE_TITLE_LEN);
			++mylineno;
			PRINT_LINE BEGIN(SPICE);
			}


<SPICE>"true"		{ RET(TRUE) }
<SPICE>"false"		{ RET(FALSE) }
<SPICE>"noecho"		{ RET(NOECHO) }
<SPICE>"write"		{ RET(WRITE) }
<SPICE>"in"		{ RET(IN_) }
<SPICE>"element"	{ RET(ELEMENT) }
<SPICE>"term"		{ RET(TERM) }
<SPICE>"plot"		{ RET(PLOT) }
<SPICE>"watch"		{ RET(WATCH) }
<SPICE>"file"		{ RET(FILE_) }
<SPICE>"system"		{ RET(SYSTEM) }
<SPICE>"poly"		{ RET(POLY_) }
<SPICE>"elem"		{ RET(ELEM) }
<SPICE>"end"		{ RET(END) }
<SPICE>"variable"	{ RET(VARIABLE_) }
<SPICE>"hb"		{ RET(HB_) }
<SPICE>"tran"		{ RET(TRAN_) }
<SPICE>"noise"		{ RET(NOISE_) }
<SPICE>"dc"		{ RET(DC) }
<SPICE>"ac"		{ RET(AC) }
<SPICE>"exp"		{ RET(EXP) }
<SPICE>"pwl"		{ RET(PWL) }
<SPICE>"sffm"		{ RET(SFFM) }
<SPICE>"sin"		{ RET(SIN) }
<SPICE>"pulse"		{ RET(PULSE) }
<SPICE>"params:"	{ RET(PARAMS_) }

<SPICE>"v("		{ RET(PRINT_V) }
<SPICE>"i("		{ RET(PRINT_I) }

<SPICE>^"a"({alphanumeric})*	{ copy_yytext(); RET(A_ELEMENT) }
<SPICE>^"b"{alphanumeric}*	{ copy_yytext(); RET(B_ELEMENT) }
<SPICE>^"c"{alphanumeric}*	{ copy_yytext(); RET(C_ELEMENT) }
<SPICE>^"d"{alphanumeric}*	{ copy_yytext(); RET(D_ELEMENT) }
<SPICE>^"e"{alphanumeric}*	{ copy_yytext(); RET(E_ELEMENT) }
<SPICE>^"f"{alphanumeric}*	{ copy_yytext(); RET(F_ELEMENT) }
<SPICE>^"g"{alphanumeric}*	{ copy_yytext(); RET(G_ELEMENT) }
<SPICE>^"h"{alphanumeric}*	{ copy_yytext(); RET(H_ELEMENT) }
<SPICE>^"i"{alphanumeric}*	{ copy_yytext(); RET(I_ELEMENT) }
<SPICE>^"j"{alphanumeric}*	{ copy_yytext(); RET(J_ELEMENT) }
<SPICE>^"k"{alphanumeric}*	{ copy_yytext(); RET(K_ELEMENT) }
<SPICE>^"l"{alphanumeric}*	{ copy_yytext(); RET(L_ELEMENT) }
<SPICE>^"m"{alphanumeric}*	{ copy_yytext(); RET(M_ELEMENT) }
<SPICE>^"n"{alphanumeric}*	{ copy_yytext(); RET(N_ELEMENT) }
<SPICE>^"o"{alphanumeric}*	{ copy_yytext(); RET(O_ELEMENT) }
<SPICE>^"p"{alphanumeric}*	{ copy_yytext(); RET(P_ELEMENT) }
<SPICE>^"q"{alphanumeric}*	{ copy_yytext(); RET(Q_ELEMENT) }
<SPICE>^"r"{alphanumeric}*	{ copy_yytext(); RET(R_ELEMENT) }
<SPICE>^"s"{alphanumeric}*	{ copy_yytext(); RET(S_ELEMENT) }
<SPICE>^"t"{alphanumeric}*	{ copy_yytext(); RET(T_ELEMENT) }
<SPICE>^"u"{alphanumeric}*	{ copy_yytext(); RET(U_ELEMENT) }
<SPICE>^"v"{alphanumeric}*	{ copy_yytext(); RET(V_ELEMENT) }
<SPICE>^"w"{alphanumeric}*	{ copy_yytext(); RET(W_ELEMENT) }
<SPICE>^"x"{alphanumeric}*	{ copy_yytext(); RET(X_ELEMENT) }
<SPICE>^"y"{alphanumeric}*	{ copy_yytext(); RET(Y_ELEMENT) }
<SPICE>^"z"{alphanumeric}*	{ copy_yytext(); RET(Z_ELEMENT) }


<*>{scale_core}"f"([a-z]*) {
	        yylval.dval = 1.e-15 * copy_yytext_return_base();
	        lastNumberAsText = yylvalStr[idIndex];
		RET(FLOAT_LIT)
		}

<*>{scale_core}"p"([a-z]*) {
	        yylval.dval = 1.e-12 * copy_yytext_return_base();
	        lastNumberAsText = yylvalStr[idIndex];
		RET(FLOAT_LIT)
		}

<*>{scale_core}"n"([a-z]*) {
	        yylval.dval = 1.e-9 * copy_yytext_return_base();
	        lastNumberAsText = yylvalStr[idIndex];
		RET(FLOAT_LIT)
		}

<*>{scale_core}"u"([a-z]*) {
	        yylval.dval = 1.e-6 * copy_yytext_return_base();
	        lastNumberAsText = yylvalStr[idIndex];
		RET(FLOAT_LIT)
		}

<*>{scale_core}("m"|("m"[a-df-z][a-z]*)) {
	        yylval.dval = 1.e-3 * copy_yytext_return_base();
	        lastNumberAsText = yylvalStr[idIndex];
		RET(FLOAT_LIT)
		}

<*>{scale_core}"k"([a-z]*) {
	        yylval.dval = 1.e3 * copy_yytext_return_base();
	        lastNumberAsText = yylvalStr[idIndex];
		RET(FLOAT_LIT)
		}

<*>{scale_core}"meg"([a-z]*) {
	        yylval.dval = 1.e6 * copy_yytext_return_base();
	        lastNumberAsText = yylvalStr[idIndex];
		RET(FLOAT_LIT)
		}

<*>{scale_core}"g"([a-z]*) {
	        yylval.dval = 1.e9 * copy_yytext_return_base();
	        lastNumberAsText = yylvalStr[idIndex];
		RET(FLOAT_LIT)
		}

<*>{scale_core}"t"([a-z]*) {
	        yylval.dval = 1.e12 * copy_yytext_return_base();
	        lastNumberAsText = yylvalStr[idIndex];
		RET(FLOAT_LIT)
		}

<*>{scale_core}[av] { /* Ignore voltage, current. */
	        yylval.dval = copy_yytext_return_base();
	        lastNumberAsText = yylvalStr[idIndex];
		RET(FLOAT_LIT)
		}

<*>{predigit}?{digits} {
		copy_yytext();
		lastNumberAsText = yylvalStr[idIndex];
		yylval.ival = atoi(yytext);
		RET(INT_LIT)
		}

<*>{predigit}?({float_core}{expon}?)|({digits}{expon})	{
		copy_yytext();
		lastNumberAsText = yylvalStr[idIndex];
		yylval.dval = atof(yytext);
		RET(FLOAT_LIT)
		}

<SPICE>{alpha_or_dollar_or_}{alphanumeric}*{white}*"=" {
	copy_yytext_to_white(); /* leave out {white}*"=" */
	RET(PARAMETER) }

<SPICE>{alpha_or_dollar_or_}{alphanumeric}*            {
	int pType;
	generic_t gval;
	copy_yytext();
	if(St_GetSym(capOptionsT_P, yytext, &pType, &gval) == ST_SYM_FOUND)
           { RET(VARIABLE) }
        else { RET(IDENTIFIER) }
	}

<SPICE>^".ac"		{ RET(DOT_AC) }
<SPICE>^".dc"		{ RET(DOT_DC) }
<SPICE>^".end"		{ RET(DOT_END) }
<SPICE>^".ends"		{ RET(DOT_ENDS) }
<SPICE>^".four"		{ RET(DOT_FOUR) }
<SPICE>^".ic"		{ RET(DOT_IC) }
<SPICE>^".inc"		{ PRINT_LINE BEGIN(INCLUDE); }
<SPICE>^".keep" { RET(DOT_KEEP) }
<SPICE>^".lib"		{ PRINT_LINE BEGIN(LIBRARY); }
<SPICE>^".model"	{ RET(DOT_MODEL) }
<SPICE>^".noise"	{ RET(DOT_NOISE) }
<SPICE>^".nodeset"	{ RET(DOT_NODESET) }
<SPICE>^".op"		{ RET(DOT_OP) }
<SPICE>^".options"	{ RET(DOT_OPTIONS) }
<SPICE>^".plot"		{ RET(DOT_PLOT) }
<SPICE>^".print"	{ RET(DOT_PRINT) }
<SPICE>^".probe"	{ RET(DOT_PROBE) }
<SPICE>^".sens"		{ RET(DOT_SENS) }
<SPICE>^".subckt"	{ RET(DOT_SUBCKT) }
<SPICE>^".tf"		{ RET(DOT_TF) }
<SPICE>^".tran"		{ RET(DOT_TRAN) }
<SPICE>^".watch"	{ RET(DOT_WATCH) }
<SPICE>^".width"	{ RET(DOT_WIDTH) }
<SPICE>^".temp"		{ RET(DOT_TEMP) }
<SPICE>^".locate"	{ RET(DOT_LOCATE) }
<SPICE>^".ref"		{ RET(DOT_REF) }
<SPICE>^".out"		{ RET(DOT_OUT) }
<SPICE>^".svhb"		{ RET(DOT_SVHB) }
<SPICE>^".tran2"	{ RET(DOT_TRAN2) }
<SPICE>^".tran3"	{ RET(DOT_TRAN3) }
<SPICE>^".tran4"	{ RET(DOT_TRAN4) }
<SPICE>^".trant"	{ RET(DOT_TRANT) }
<SPICE>^".svtr"		{ RET(DOT_SVTR) }
<SPICE>^".wavtran"	{ RET(DOT_WAVTRAN) }


<SPICE>"\""[^"\""]*"\"" {
		copy_yytext_string();
 		RET(STRING_LIT)
		}

<SPICE>"{"[^"}"]*"}" {
		copy_yytext_string();
 		RET(EXPRESSION)
		}

<SPICE>"("	{ RET('(') }
<SPICE>")"	{ RET(')') }
<SPICE>","	{ RET(',') }
<SPICE>"="	{ RET('=') }
<SPICE>";"	{ RET(';') }
<SPICE>":"	{ RET(':') }
<SPICE>"*"	{ RET('*') }
<SPICE>"+"	{ RET('+') }
<SPICE>"-"	{ RET('-') }
<SPICE>"/"	{ RET('/') }
<SPICE>"\\"	{ RET('\\') }
<SPICE>"<"	{ RET('<') }
<SPICE>">"	{ RET('>') }

<*>{white}* /* eat whitespace */

<*>{lf}"+"	{ ++mylineno; PRINT_LINE}


<*>{lf}"*"	{ ++mylineno; PRINT_LINE BEGIN(CONTINUATION); }

<CONTINUATION>([^\n]*){lf}"+"  {
		  ++mylineno; PRINT_LINE BEGIN(SPICE); }

<CONTINUATION>([^\n]*){lf}"*"  {
		  ++mylineno; PRINT_LINE }

<CONTINUATION>([^\n]*){lf}  {
		  ++mylineno; PRINT_LINE  BEGIN(SPICE); RET(EOL) }


<*>{lf}		{ ++mylineno; PRINT_LINE RET(EOL) }

<SPICE>^"*"([^\n]*){lf}     { ++mylineno; PRINT_LINE }


<INCLUDE>[^ \t\n]+ { /* Support for .INC, got the include filename */
		   if (include_stack_ptr >= MAX_INCLUDE_DEPTH )
		     ParseError(".INC and/or .LIB nested too deeply.\n");
		   include_stack[include_stack_ptr] = YY_CURRENT_BUFFER;
		   include_stack_mylineno[include_stack_ptr] = mylineno;
		   yyin = fopen(yytext, "r");
		   if(!yyin) 
		     ParseError(".INC file not found.\n");
		   include_files[include_stack_ptr++] = yyin;
		   yy_switch_to_buffer(yy_create_buffer(yyin,YY_BUF_SIZE));
		   mylineno = 1;
		   BEGIN(SPICE);
		   }

<LIBRARY>[^ \t\n]+ { /* Support for .LIB, got the include filename */
		   if (include_stack_ptr >= MAX_INCLUDE_DEPTH )
		     ParseError(".LIB and/or .INC nested too deeply.\n");
		   include_stack[include_stack_ptr] = YY_CURRENT_BUFFER;
		   include_stack_mylineno[include_stack_ptr] = mylineno;
		   /* Add library path prefix */
		   {
		   char* libstring;
		   libstring = (char *) malloc(strlen(env_freeda_library) + 1
		     + strlen(yytext) + 1);
		   strcpy(libstring,env_freeda_library);
		   strcat(libstring,"/");
		   strcat(libstring,yytext);
		   yyin = fopen(libstring, "r");
		   if(!yyin) 
		     ParseError(".LIB file not found.\n");
		   fprintf(output_F,"** Opening library: %s\n",libstring);
		   free(libstring);
                   }
		   include_files[include_stack_ptr++] = yyin;
		   yy_switch_to_buffer(yy_create_buffer(yyin,YY_BUF_SIZE));
		   mylineno = 1;
		   BEGIN(SPICE);
		   }

<<EOF>>		{
		if(--include_stack_ptr < 0)
		  return(EOF);
		else
		  {
		  yy_delete_buffer(YY_CURRENT_BUFFER);
		  mylineno=include_stack_mylineno[include_stack_ptr];
		  yy_switch_to_buffer(include_stack[include_stack_ptr]);
		  fclose(include_files[include_stack_ptr]) ;
		  no_print_line =1; // Avoids printing .INC, .LIB line twice
		  }
		}

%%
/*
 * copy_yytext
 *
 * Copy token in yytext to save string.
 *
 */
void copy_yytext(void)
{
  int maxLength;
  str2lower(yytext);
  maxLength = MIN(MAX_TOK_LEN,yyleng);
  idIndex = (thisDepth++) % YYMAXDEPTH;
  strncpy(yylvalStr[idIndex], yytext, maxLength);
  yylvalStr[idIndex][maxLength] = 0;
  yylval.id=yylvalStr[idIndex];
}

/*
 * copy_yytext_string
 *
 * Copy token in yytext to save string.
 *
 * String delimited by "   "
 *
 * Remove - delimeters
 *        - remove escape sequences from string
 *        - remove "lf"
 *        - remove the sequence "lf""+"
 */
void copy_yytext_string(void)
{
  int i, maxLength;
  char *outStr, *inStr;
  str2lower(yytext);
  inStr = yytext;
  idIndex = (thisDepth++) % YYMAXDEPTH;
  outStr = yylvalStr[idIndex];
  maxLength = MIN(MAX_TOK_LEN,yyleng-2); /* skip first/last charcter, it is " */
  inStr++; /* skip first character, it is a " */
  for(i = 0; i < maxLength; i++)
    {
      if(*inStr == '\n')
	{
	  inStr++;
	  if(*inStr == '+' ) inStr++;
	}
      else if(*inStr == '\\' && *(inStr+1) == 'n' )
	{
	  *outStr++ = '\n';
	  inStr = inStr + 2;
	}
      else if(*inStr == '\\' && *(inStr+1) == '\\' )
	{
	  *outStr++ = *inStr++;
	  inStr++;
	}
      else if(*inStr != '\\')
	{ *outStr++ = *inStr++; }
    }
  *outStr = 0;
  yylval.id=yylvalStr[idIndex];
  return;
}


/*
 * copy_yytext_to_white
 *
 * Copy token in yytext up to but not including first white space or =
 * to save string.
 *
 */
void copy_yytext_to_white(void)
{
  int i,
    maxLength;
  char *outStr, *inStr;
  str2lower(yytext);
  inStr = yytext;
  idIndex = (thisDepth++) % YYMAXDEPTH;
  outStr = yylvalStr[idIndex];
  maxLength = MIN(MAX_TOK_LEN,yyleng);
  for(i = 0; i < maxLength; i++)
    {
      *outStr = *inStr++;
      if(   (*outStr == '\x00') || (*outStr == '\x09') || (*outStr == '\x0B')
	    || (*outStr == '\x0C') || (*outStr == '\x20') || (*outStr == '=')
	    ) break;
      outStr++;
    }
  *outStr = 0;
  yylval.id=yylvalStr[idIndex];
  return;
}


/*
 * copy_yytext_return_base
 *
 * Copy token in yytext and return base number obtained without scaling
 * factor.
 *
 */
double copy_yytext_return_base(void)
{
  int i, maxLength;
  char *outStr, *inStr, *startStr;
  double number;
  str2lower(yytext);
  inStr = yytext;
  idIndex = (thisDepth++) % YYMAXDEPTH;
  startStr = yylvalStr[idIndex];
  maxLength = MIN(MAX_TOK_LEN,yyleng);
  for(i = 0, outStr = startStr; i < maxLength; i++) /* Copy to scale factor */
    {
      *outStr = *inStr++;
      if( ( *outStr >= 'a') && (*outStr <= 'z')  ) break;
      outStr++;
    }
  *outStr = 0;
  number = atof(startStr);
  inStr--;
  while(*outStr++ = *inStr++); /* Copy rest of string */
  yylval.id = startStr;
  return(number);
}


/*
 * PrintCurrentLine
 *
 * Print current line
 *
 * Parameters:
 *	file_F	Pointer to file to write to
 *	b	Ppointer to character pointer to start of line.
 *
 * Use following flex generated scanner
 *
 */
int PrintCurrentLine(IN FILE *file_F, OUT char **b)
{
  char *c;
  int i, holdCharFound=0;
  /* back up to last newline */
  for (*b = yytext - 1;
       (*b >= YY_CURRENT_BUFFER->yy_ch_buf && (**b != '\n')); *b = *b-1);
  /* now print line */
  for (c = *b+1,i=0; i < MAX_SPICE_LINE_LEN; i++, c++)
    {
      if (!(*c))
	{
	  if (holdCharFound) break;
	  putc(yy_hold_char, file_F); /* Nulls are used for holding a position.
				       * Only do this once. */
	  holdCharFound=1;
	}
      else
	if (*c == '\n') break;
      /* else if (*c == universalCommentIndicator)
	 putc(commentChar, file_F); */
	else
	  putc(*c, file_F);
    }
  putc('\n', file_F);
  return 0;
}



/*
 * PrintTextPointer
 *
 * Parameters:
 *	file_F	Pointer to file to write to
 *	b	Character pointer to start of line.
 *
 * Use following flex generated scanner
 *
 */
int PrintTextPointer(IN FILE *file_F, IN char *b)
{
  char *c;
  int i;
  /* print out spaces until we get to offending text */
  /* for ( c = b+1, i=0; strncmp( c,yytext,yyleng ) && *c && (*c != '\n') &&
     i < MAX_SPICE_LINE_LEN; i++, c++ )  */
  for ( c = b+1, i=0; c != yytext && *c && (*c != '\n') &&
	  i < MAX_SPICE_LINE_LEN; i++, c++ )
    {
      if (*c == '\t')
	putc('\t', file_F);
      else
	putc(' ', file_F);
    }
  putc('^', file_F);
  putc('\n', file_F);
  return 0;
}




/*
 * ParseMessage-
 *	Print message message.
 */
int ParseMessage(FILE *output_F, const char *message, char *warningType)
{
  char *b;
/*    extern char *yytext; */
  PrintCurrentLine(output_F, &b);
  PrintTextPointer(output_F, b);
  fprintf(output_F,"  Parse %s near line %d\n", warningType, mylineno);
  fprintf(output_F,"  %s\n", message);
  return 0;
}



/*
 * ParseWarning
 *	Print parse warning.
 */
void ParseWarning(const char *s)
{
/*    extern char *yytext; */
  ParseMessage(stderr, s, "Warning");
  ParseMessage(output_F, s, "Warning");
}

/*
 * ParseError--
 *	Print parse error message.
 */
int ParseError(const char *s)
{
/*    extern char *yytext; */
  ParseMessage(stderr, s, "Error");
  ParseMessage(output_F, s, "Error");
  exit(1);
}


/*
 * CleanUpFLEX
 *	Do something (print it for now) with an error message.
 *      This routine must follow the inclusion of the flex generated scanner.
 */
int PCleanUpFLEX()
{
  yy_delete_buffer( YY_CURRENT_BUFFER );
  return 0;
}

/* Convert string to lowercase */
void str2lower(char *s)
{
  int i;
  for (i=0; i<strlen(s); i++)
    s[i] = tolower(s[i]);
}
