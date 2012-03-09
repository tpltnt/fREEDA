#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
/*
*  an_output.c--  
*   
*/

#include "parser.h"
#include "../compat/dcx.h"
#include "report.h"

/* This is a Cygwin bug workaround. It should not be needed after the
compiler bug is fixed. */
#ifndef __CYGWIN__
#define OUTNUMFORM "%16.9e"
#else
#define OUTNUMFORM "%16.10g"
#endif

static capOutReq_t *stack[50];
char *thisRequestName;
int stackDepth;
long ObjType2ArgType(int objType);

/* From gnuplot.c */
void gnuplot(char *plotfile);

void An_CalcError(const char *s)
{
  int i;
	
  report(WARNING, s);
  report(WARNING, "See log file");
	
  fprintf(output_F, "Stack :\n");
	
  for (i = 0; i < stackDepth; i++) 
	{
    fprintf(output_F, "[%1d]\t", stackDepth - i);
    switch(stack[i]->type) {
			case CAP_OBJ_STRING:
      fputs(stack[i]->obj.s, output_F);
      putc('\n', output_F);
      break;
			case CAP_OBJ_INT:
      fprintf(output_F, "%1d\n", stack[i]->obj.i);
      break;
			case CAP_OBJ_OPER:
      fprintf(output_F,"(operator) %s\n",GET_OPERATOR_NAME(stack[i]));
      break;
			case CAP_OBJ_DOUBLE:
      fprintf(output_F, OUTNUMFORM"\n", stack[i]->obj.d);
      break;
			case CAP_OBJ_DCX:
      fprintf(output_F, OUTNUMFORM" "OUTNUMFORM"\n", 
			stack[i]->obj.dcx.re, 
			stack[i]->obj.dcx.im);
      break;
			case CAP_OBJ_DOUBLEV:
      report(WARNING, "double vector (x,y)\n");
      break;
			case CAP_OBJ_DOUBLEM:
      report(WARNING,"double matrix (x,y1,y2,...)\n");
      break;
			case CAP_OBJ_DCXV:
      report(WARNING,"double vector of complex numbers\n");
      break;
			case CAP_OBJ_TERM:
      report(WARNING,"terminal\n");
      break;
			case CAP_OBJ_NODE:
      report(WARNING,"element\n");
      break;
			case CAP_OBJ_DATAFILE:
      report(WARNING,"data filename\n");
      break;
			case CAP_OBJ_FILENAME:
      report(WARNING,"filename\n");
      break;
			default:
      report(WARNING,"unknown operator type\n");
      break;
    }
  }
  exit(1);
}

/*
* TypecheckOneArg
*
* Return errFlag
*   TRUE if fatal error
*   FALSE if okay
*/
int TypecheckOneArg(int argType, int argMask)
{
  int okay, errFlag = FALSE;
  okay = (ObjType2ArgType(argType) & argMask) != 0L;
  if (!okay) errFlag = TRUE;
  return(errFlag);
}

/* Get the contents of the stack and decrement the stackPointer */
capOutReq_t *Pop(void)
{
  if (stackDepth == 0)
	{
		capOutReq_Pt result;
		result = (capOutReq_t *) malloc(sizeof(capOutReq_t));
		result->protect = FALSE;
		result->type = CAP_OBJ_NONEXISTENT;
		return result;
	}
  stackDepth--;
  return stack[stackDepth];
}

/* The stack operates as a FILO */
/* Get the contents of the stack at a specified position
* but don't decrement the stack position */
capOutReq_t *Peek(int depth)
{
  if (stackDepth - depth == 0)
    An_CalcError("Stack underflow\n");
  return stack[stackDepth - 1 - depth];
}

/*
* GetOutReqVariable
*
* Get Output Request Variable from output variable symbol table
*
* The pointer is returned and so the request remains on the symbol table.
*
* Returns 0 id symbol not found.
*
*/
int GetOutReqVariable(char *variableName, capOutReq_Pt *outReq_P)
{
	int type;
	generic_t gval;
	if(St_GetSym(capOutputVariablesT_P, variableName, &type, &gval) 
		== ST_SYM_FOUND)
	{
		if(type != GEN_TYPE_VOID_PTR) return(FALSE);
		*outReq_P = (capOutReq_Pt) gval.v;
		return(TRUE);
	}
	return(FALSE);
}



int UpdateOutReqVariable(char *variableName, capOutReq_Pt *outReq_P)
{
	int type;
	generic_t gval;
	
	type = GEN_TYPE_VOID_PTR;
	gval.v = outReq_P;
	if (St_ReplSym(capOutputVariablesT_P, variableName, type, gval) != ST_SYM_FOUND) 
		return(1);
	return(0);
}



void FreeOutReq(capOutReq_t *outReq)
{
	if(!outReq) return;
	if(outReq->protect == TRUE) return; /* There is a special reason not to
	delete it */
	switch(outReq->type)
	{
		case CAP_OBJ_DOUBLEV:
		if (outReq->obj.dv) free(outReq->obj.dv);
		if (outReq->x) free(outReq->x);
		break;
		case CAP_OBJ_DOUBLEM:
		if (outReq->obj.s_dm.dm) free(outReq->obj.s_dm.dm);
		if (outReq->x) free(outReq->x);
		break;
		case CAP_OBJ_DCXV:
		if (outReq->obj.dcxv) free(outReq->obj.dcxv);
		if (outReq->x) free(outReq->x);
		break;
		
		case CAP_OBJ_DCXM:
		if (outReq->obj.s_dcxm.dcxm) free(outReq->obj.s_dcxm.dcxm);
		if (outReq->x) free(outReq->x);
		break;
		
		case CAP_OBJ_STRING:
		if (outReq->obj.s) free(outReq->obj.s);
		break;
	}
	free(outReq);
}

/*
* ReadFile--
*	Read an object from a file.
*/
capOutReq_t *ReadFile(char *filename)
{
	FILE *f;
	char s[256], s1[80], s2[80], s3[80];
	double *xv, *yv;
	dcx_t *dcxv;
	static int blockSize = 1024;
	int n, blockOff, bufi;
	capOutReq_t *it;
	
	f = fopen(filename, "r");
	if (!f) An_CalcError("Could not open data file");
	
	/*
	* Examine the first line of the file.  If it contains
	* one word (we'll assume it's a number) then we read y data
	* only.  If it contains two words then we read x,y pairs.
	* If it contains three words then we read x + complex y.
	* No error checking or anything subtle at work here for now.
	* Skip empty lines.
	*/
	fgets(s, 256, f);
	while (!(n = sscanf(s, "%s%s%s", s1, s2, s3)) && !feof(f) && 
		!strcmp(s1, "--"))
	fgets(s, 256, f);
	
	if (feof(f))
		An_CalcError("Unexpected end of data file");
	
	bufi = blockOff = 0;
	
	/*
	* Allocate a return object
	*/
	it = (capOutReq_t *) malloc(sizeof(capOutReq_t));
	
	/*
	* Read data from the file, expanding arrays dynamically as we go,
	* if necessary.  Contract the arrays to the proper size after
	* we're done reading.
	*/
	if (n == 1) 
	{
		yv = (double *) malloc(sizeof(double) * blockSize);
		for (;;) 
		{
			yv[bufi + blockOff] = atof(s1);
			for (;;) 
			{
				if (feof(f))
					break;
				fgets(s, 256, f);
				if (sscanf(s, "%s", s1) && strcmp(s1, "--"))
					break;
			}
			bufi++;
			if (feof(f))
				break;
			if (bufi == blockSize) 
			{
				blockOff += blockSize;
				yv = (double *) realloc((VOIDPTR) yv,
				sizeof(double) * (blockOff * blockSize));
				bufi = 0;
			}
		}
		fclose(f);
		yv = (double *) realloc((VOIDPTR) yv, sizeof(double) *
		(blockOff + bufi));
		it->type = CAP_OBJ_DOUBLEV;
		it->size = blockOff + bufi;
		it->obj.dv = yv;
		it->x = NULL;
	} 
	else if (n == 2) 
	{
		xv = (double *) malloc(sizeof(double) * blockSize);
		yv = (double *) malloc(sizeof(double) * blockSize);
		for (;;) 
		{
			xv[bufi + blockOff] = atof(s1);
			yv[bufi + blockOff] = atof(s2);
			for (;;) 
			{
				if (feof(f))
					break;
				fgets(s, 256, f);
				if (sscanf(s, "%s%s", s1, s2) && 
					strcmp(s1, "--"))
				break;
			}
			bufi++;
			if (feof(f))
				break;
			if (bufi == blockSize) 
			{
				blockOff += blockSize;
				xv = (double *) realloc((VOIDPTR) xv,
				sizeof(double) * (blockOff * blockSize));
				yv = (double *) realloc((VOIDPTR) yv,
				sizeof(double) * (blockOff * blockSize));
				bufi = 0;
			}
		}
		fclose(f);
		xv = (double *) realloc((VOIDPTR) xv, sizeof(double) *
		(blockOff + bufi));
		yv = (double *) realloc((VOIDPTR) yv, sizeof(double) *
		(blockOff + bufi));
		it->type = CAP_OBJ_DOUBLEV;
		it->size = blockOff + bufi;
		it->obj.dv = yv;
		it->x = xv;
	} 
	else 
	{
		xv = (double *) malloc(sizeof(double) * blockSize);
		dcxv = (dcx_t *) malloc(sizeof(dcx_t) * blockSize);
		for (;;) 
		{
			xv[bufi + blockOff] = atof(s1);
			dcxv[bufi + blockOff].re = atof(s2);
			dcxv[bufi + blockOff].im = atof(s3);
			for (;;) 
			{
				if (feof(f))
					break;
				fgets(s, 256, f);
				if (sscanf(s, "%s%s%s", s1, s2, s3) &&
					strcmp(s1, "--"))
				break;
			}
			bufi++;
			if (feof(f))
				break;
			if (bufi == blockSize) 
			{
				blockOff += blockSize;
				xv = (double *) realloc((VOIDPTR) xv,
				sizeof(double) * (blockOff * blockSize));
				dcxv = (dcx_t *) realloc((VOIDPTR) dcxv,
				sizeof(dcx_t) * (blockOff * blockSize));
				bufi = 0;
			}
		}
		fclose(f);
		xv = (double *) realloc((VOIDPTR) xv, sizeof(double) *
		(blockOff + bufi));
		dcxv = (dcx_t *) realloc((VOIDPTR) dcxv, sizeof(dcx_t) *
		(blockOff + bufi));
		it->type = CAP_OBJ_DCXV;
		it->size = blockOff + bufi;
		it->obj.dcxv = dcxv;
		it->x = xv;
	}
	return it;
}

int WriteObj(FILE *f, capOutReq_t *arg)
{
	int i;
	switch(arg->type) 
	{
		case CAP_OBJ_STRING:
		fputs(arg->obj.s, f);
		putc('\n', f);
		break;
		case CAP_OBJ_INT:
		fprintf(f, "%1d\n", arg->obj.i);
		break;
		case CAP_OBJ_DOUBLE:
		fprintf(f, OUTNUMFORM"\n", arg->obj.d);
		break;
		case CAP_OBJ_DCX:
		fprintf(f, "(."OUTNUMFORM", "OUTNUMFORM")\n", arg->obj.dcx.re, arg->obj.dcx.im);
		break;
		case CAP_OBJ_DOUBLEV:
		if (arg->x) 
		{
			for (i = 0; i < arg->size; i++) 
			{
				fprintf(f, OUTNUMFORM"\t", arg->x[i]);
				fprintf(f, OUTNUMFORM"\n", arg->obj.dv[i]);
			}
		}
		else 
		{
			for (i = 0; i < arg->size; i++)
				fprintf(f, OUTNUMFORM"\n", arg->obj.dv[i]);
		}
		break;
		case CAP_OBJ_DOUBLEM:
		{ int isize = arg->size;
			int jsize = arg->obj.s_dm.rank;
			doublev_t  x = arg->x;
			doublem_t  dm = arg->obj.s_dm.dm ;
			int j;
			if (x) 
			{
				for (i = 0; i < isize; i++) 
				{
					fprintf(f,OUTNUMFORM"", x[i]);
					for (j = 0; j < jsize; j++) fprintf(f,"\t"OUTNUMFORM"",dm[i][j]);
					fprintf(f,"\n");
				}
			}
			else 
			{
				for (i = 0; i < isize; i++) 
				{
					if(jsize) fprintf(f, OUTNUMFORM"", dm[i][0]);
					for (j = 1; j < jsize; j++) fprintf(f, "\t"OUTNUMFORM"", dm[i][j]);
					fprintf(f,"\n");
				}
			}
		}
		break;
		case CAP_OBJ_DCXV:
		for (i = 0; i < arg->size; i++) 
		{
			if (arg->x)
				fprintf(f, OUTNUMFORM"\t", arg->x[i]);
			fprintf(f, OUTNUMFORM" "OUTNUMFORM"\n",
			arg->obj.dcxv[i].re, arg->obj.dcxv[i].im);
		}
		break;
		case CAP_OBJ_DCXM:
		
		{ int isize = arg->size;
			int jsize = arg->obj.s_dcxm.rank;
			doublev_t  x = arg->x;
			dcxm_t  dxm = arg->obj.s_dcxm.dcxm;
			int j;
			if (x) 
			{
				for (i = 0; i < isize; i++) 
				{
					fprintf(f,OUTNUMFORM"", x[i]);
					for (j = 0; j < jsize; j++)
					{
						fprintf(f,"\t"OUTNUMFORM" \t"OUTNUMFORM"",dxm[i][j].re,dxm[i][j].im);
					}
					fprintf(f,"\n");
				}
			}
			else 
			{
				for (i = 0; i < isize; i++) 
				{
					if(jsize) fprintf(f, OUTNUMFORM" \t"OUTNUMFORM"", dxm[i][0].re, dxm[i][0].im);
					for (j = 1; j < jsize; j++)
					{
						fprintf(f, "\t"OUTNUMFORM" \t"OUTNUMFORM"", dxm[i][j].re, dxm[i][j].im);
					}
					fprintf(f,"\n");
				}
			}
		}
		break;
		case CAP_OBJ_FUNC:  /* function pointer */
		arg->obj.func(f);
		break;
		default:
		An_CalcError("Can't write object of result type");
	}
	return 0;
}

/* Returns next argument from stack if string
* Used in an_op_plot() */
void checkarg(capOutReq_t **arg)
{
	int errFlag;
	*arg = Peek(0);
	errFlag = TypecheckOneArg((*arg)->type, STRING_ARG_TYPE);
	if (errFlag) 
		*arg = NULL;
	else 
		Pop();
}

/*
* An_Op_WATCH--
*/
int An_Op_WATCH(capOutReq_t *result)
{
	/* Note: this routine is almost the same as An_Op_Plot().
	* The reason for the existence of this is unknown to me. */
	FILE *f; 
	FILE *command;
	
	char plotfile[1024];
	
	capOutReq_t *arg1 = NULL;
	capOutReq_t *arg2 = NULL;
	capOutReq_t *arg3 = NULL;
	capOutReq_t *arg4 = NULL;
	int errFlag;
	
	char *currentOutputFileName;
	char *gnuplotPreambleScript, *gnuplotPostambleScript;
	
	/* Get filename */
	arg1 = Pop();
	errFlag = TypecheckOneArg((arg1)->type, STRING_ARG_TYPE);
	if (errFlag) 
	{
		An_CalcError("Missing output file name.");
		return 1;
	}
	currentOutputFileName = (arg1)->obj.s;
	
	gnuplotPreambleScript = "# Preamble";
	gnuplotPostambleScript = "# Postamble";
	
	checkarg(&arg2);
	if (arg2) 
	{
		gnuplotPreambleScript = arg2->obj.s;
		checkarg(&arg3);
		if (arg3) 
		{
			gnuplotPostambleScript = arg3->obj.s;
		}
	}
	
	/* Write command file */
	sprintf(plotfile, "%s.cmd", currentOutputFileName);
	command = fopen(plotfile,"w");
	if(!command) 
	{
		printf("Error: In gnuplot");
		exit(1);
	}
	/*  Write gnuplot script file. */
	fprintf(command,"set title '%s';\n",currentOutputFileName);  
	fprintf(command,"set zero 1e-15; set grid; set style data lines\n");
	fprintf(command, "%s \n", gnuplotPreambleScript);
  fprintf(command,"plot '%s' \\\n", currentOutputFileName);
	fprintf(command, "%s \n", gnuplotPostambleScript);
	fprintf(command,"pause -1 'Type <Enter> to remove the plot %s'\n"
    ,currentOutputFileName);
	fclose(command);
	
	f = fopen(currentOutputFileName, "w");
	if (!f) {
		An_CalcError("Could not open output file.");
		return 1;
	}
	/* Get output data from stack. */
	while ((arg4=Pop())->type != CAP_OBJ_NONEXISTENT) 
	{
		WriteObj(f, arg4);
		FreeOutReq(arg4);
	}
	fclose(f);
	
	fprintf(output_F," Plotting output file: %s.\n",\
	currentOutputFileName);    
	printf(" Plotting output file: %s.\n",\
	currentOutputFileName);
	/* Call gnuplot program */
	gnuplot(plotfile);
	
	/* Free output requests */
	FreeOutReq(arg1);
	FreeOutReq(arg2); 
	FreeOutReq(arg3);
	FreeOutReq(arg4);
	
	return 0;
}

/*
* An_Op_PLOT-- This routine writes files for gnuplot
*
* Assumed order of arguments in netlist:
* *_vector           [required]
* gnuplotPreambleScript   [optional] 
* gnuplotPostambleScript  [optional] Required if above optional strings used
*/
int An_Op_PLOT(capOutReq_t *result)
{
	FILE *f; 
	FILE *command;
	
	char plotfile[1024];
	
	capOutReq_t *arg1 = NULL;
	capOutReq_t *arg2 = NULL;
	capOutReq_t *arg3 = NULL;
	capOutReq_t *arg4 = NULL;
	int errFlag;
	
	char *currentOutputFileName;
	char *gnuplotPreambleScript, *gnuplotPostambleScript;
	
	/* Get filename */
	arg1 = Pop();
	errFlag = TypecheckOneArg((arg1)->type, STRING_ARG_TYPE);
	if (errFlag) 
	{
		An_CalcError("Missing output file name.");
		return 1;
	}
	currentOutputFileName = (arg1)->obj.s;
	
	gnuplotPreambleScript = "# Preamble";
	gnuplotPostambleScript = "# Postamble";
	
	checkarg(&arg2);
	if (arg2) 
	{
		gnuplotPreambleScript = arg2->obj.s;
		checkarg(&arg3);
		if (arg3) 
		{
			gnuplotPostambleScript = arg3->obj.s;
		}
	}
	
	/* Write command file */
	sprintf(plotfile, "%s.cmd", currentOutputFileName);
	command = fopen(plotfile,"w");
	if(!command) 
	{
		printf("Error: In gnuplot");
		exit(1);
	}
	/*  Write gnuplot script file. */
	fprintf(command,"set title '%s';\n",currentOutputFileName);  
	fprintf(command,"set zero 1e-15; set grid; set style data lines\n");
	fprintf(command, "%s \n", gnuplotPreambleScript);
  fprintf(command,"plot '%s' \\\n", currentOutputFileName);
	fprintf(command, "%s \n", gnuplotPostambleScript);
	fprintf(command,"pause -1 'Type <Enter> to remove the plot %s'\n"
    ,currentOutputFileName);
	fclose(command);
	
	f = fopen(currentOutputFileName, "w");
	if (!f) 
	{
		An_CalcError("Could not open output file.");
		return 1;
	}
	/* Get output data from stack. */
	while ((arg4=Pop())->type != CAP_OBJ_NONEXISTENT) 
	{
		WriteObj(f, arg4);
		FreeOutReq(arg4);
	}
	fclose(f);
	
	fprintf(output_F," Plotting output file: %s.\n",\
	currentOutputFileName);    
	printf(" Plotting output file: %s.\n",\
	currentOutputFileName);
	
	/* Write name of output file for user interface */
	fprintf(plotlist_F, "%s\n", currentOutputFileName);
	
	if( (LookupSym(capOptionsT_P, "gnuplot")) ) 
		/* Call gnuplot program */
	gnuplot(plotfile);
	
	/* Free output requests */
	FreeOutReq(arg1);
	FreeOutReq(arg2); 
	FreeOutReq(arg3);
	FreeOutReq(arg4);
	
	return 0;
}

/*
* An_Op_READ--
*
*/
int An_Op_READ(capOutReq_t *arg, capOutReq_t *result)
{
	report(WARNING,"An_Op_READ: not implemented.");
	return 0;
}


/*
* An_COPY--
* 
*/ 
int An_COPY(capOutReq_t *outReq, capOutReq_t **newReq)
{
	(*newReq) = (capOutReq_t *) calloc(1, sizeof(capOutReq_t));
	switch (outReq->type) 
	{
		case CAP_OBJ_NONEXISTENT:
		break;
		case CAP_OBJ_TERM:
		(*newReq)->obj.tid = outReq->obj.tid; 
		break;
		case CAP_OBJ_EDGE:
		(*newReq)->obj.eid = outReq->obj.eid; 
		break;
		case CAP_OBJ_NODE:
		(*newReq)->obj.nid = outReq->obj.nid; 
		break;
		case CAP_OBJ_DCX:
		(*newReq)->obj.dcx = outReq->obj.dcx; 
		break;
		case CAP_OBJ_DCXV:
		/* ERROR Not really copied */
		(*newReq)->obj.dcxv = outReq->obj.dcxv; 
		break;
		case CAP_OBJ_STRING:
		(*newReq)->obj.s = 
		(char *) malloc(sizeof(char) * (strlen(outReq->obj.s) + 1));
		strcpy((*newReq)->obj.s, outReq->obj.s);
		break;
		case CAP_OBJ_DOUBLEM:
		/* ERROR Not really copied */
		(*newReq)->obj.s_dm = outReq->obj.s_dm; 
		break;
		case CAP_OBJ_DCXMV:
		/* ERROR Not really copied */
		(*newReq)->obj.s_dcxmv = outReq->obj.s_dcxmv; 
		break;
		default:
		(*newReq)->obj = outReq->obj;
	}
	(*newReq)->type = outReq->type;
	(*newReq)->size = outReq->size;
	/* ERROR Not really copied */
	(*newReq)->x = outReq->x;
	/* ERROR Not really copied */
	(*newReq)->xName = outReq->xName;
	/* ERROR Not really copied */
	(*newReq)->objName = outReq->objName;
	(*newReq)->next = outReq->next;
	return(0);
}

static void Push(capOutReq_t *outReq)
{
	/*
	printf("pushing ");
	WriteObj(stdout, outReq);
	*/
	stack[stackDepth] = outReq;
	stackDepth++;
}


long ObjType2ArgType(int objType)
{
	switch (objType) 
	{
		case CAP_OBJ_NONEXISTENT:
		/* Possible stack underflow */
		return NO_ARG_TYPE;
		case CAP_OBJ_TERM:
		return TERM_ARG_TYPE;
		case CAP_OBJ_NODE:
		return NODE_ARG_TYPE;
		case CAP_OBJ_EDGE:
		return EDGE_ARG_TYPE;
		case CAP_OBJ_DATAFILE:
		return FILE_ARG_TYPE;
		case CAP_OBJ_FILENAME:
		return STRING_ARG_TYPE;
		case CAP_OBJ_VAR:
		return VAR_ARG_TYPE;
		case CAP_OBJ_INT:
		return INT_ARG_TYPE;
		case CAP_OBJ_DOUBLE:
		return DOUBLE_ARG_TYPE;
		case CAP_OBJ_DCX:
		return DCX_ARG_TYPE;
		case CAP_OBJ_DOUBLEV:
		return DOUBLEV_ARG_TYPE;
		case CAP_OBJ_DCXV:
		return DCXV_ARG_TYPE;
		case CAP_OBJ_STRING:
		return STRING_ARG_TYPE;
		case CAP_OBJ_DOUBLEM:
		return DOUBLEM_ARG_TYPE;
		case CAP_OBJ_DOUBLEMV:
		return DOUBLEMV_ARG_TYPE;
		case CAP_OBJ_DCXM:
		return DCXM_ARG_TYPE;
		case CAP_OBJ_DCXMV:
		return DCXMV_ARG_TYPE;
		case CAP_OBJ_OPER:
		default:
		An_CalcError("Internal error in ObjType2ArgType");
	}
	return(0);
}


static void PrintType(long argMask)
{
	switch(argMask) 
	{
		case ANY_ARG_TYPE:
		fprintf(output_F, "any type");
		break;
		case NUM_SV_ARG_TYPE:
		fprintf(output_F, "a numeric vector or scalar type");
		break;
		case REAL_SV_ARG_TYPE:
		fprintf(output_F, "a real vector or scalar type");
		break;
		case DCX_SV_ARG_TYPE:
		fprintf(output_F, "a complex vector or scalar type");
		break;
		case VECTOR_ARG_TYPE:
		fprintf(output_F, "a numeric vector type");
		break;
		case NUM_ARG_TYPE:
		fprintf(output_F, "a numeric scalar type");
		break;
		case REAL_ARG_TYPE:
		fprintf(output_F, "a real scalar type");
		break;
		case INT_ARG_TYPE:
		fprintf(output_F, "an integer");
		break;
		case DOUBLE_ARG_TYPE:
		fprintf(output_F, "a floating-point number");
		break;
		case DCX_ARG_TYPE:
		fprintf(output_F, "a complex number");
		break;
		case STRING_ARG_TYPE:
		fprintf(output_F, "a string");
		break;
		case DOUBLEV_ARG_TYPE:
		fprintf(output_F, "a floating-point vector");
		break;
		case DCXV_ARG_TYPE:
		fprintf(output_F, "a complex vector");
		break;
		case DOUBLEM_ARG_TYPE:
		fprintf(output_F, "a floating-point matrix");
		break;
		case DOUBLEMV_ARG_TYPE:
		fprintf(output_F, "a vector of floating-point matrices");
		break;
		case DCXM_ARG_TYPE:
		fprintf(output_F, "a complex matrix");
		break;
		case DCXMV_ARG_TYPE:
		fprintf(output_F, "a vector of complex matrices");
		break;
		case TERM_ARG_TYPE:
		fprintf(output_F, "a terminal reference");
		break;
		case NODE_ARG_TYPE:
		fprintf(output_F, "a junction reference");
		break;
		case EDGE_ARG_TYPE:
		fprintf(output_F, "a line reference");
		break;
		case FILE_ARG_TYPE:
		fprintf(output_F, "a data file");
		break;
		case VAR_ARG_TYPE:
		fprintf(output_F, "a variable name");
		break;
		default:
		fprintf(output_F, "a weird type (internal error in PrintType)");
		break;
	}
}


static void Typecheck(int opIndex)
{
	int i, errFlag, argCt;
	capOutReq_t *arg;
	argCt = an_OpTbl[opIndex].argCt;
	if(argCt == VAR_ARG) return;
	for (i = 0; i < argCt; i++)
	{
		arg = Peek(argCt - i - 1);
		errFlag = TypecheckOneArg(arg->type, an_OpTbl[opIndex].argMask[i]);
		if (errFlag)
		{
			report(WARNING,"Type error(s)");
			
			fprintf(output_F, "Argument #%1d for op '%s' must be ",
			i, an_OpTbl[opIndex].name);
			PrintType(an_OpTbl[opIndex].argMask[i]);
			fprintf(output_F, "\n");
			report(WARNING,"Can't continue");
			An_CalcError("Can't continue");
			return;
		}
	}
}


static int DoOneRequest(capOutReq_t *outReq)
{
	capOutReq_t *tmp, *newReq, *result, *arg1, *arg2, *arg3, *arg4;
	int i, status;
	stack[0] = NULL;
	stackDepth = 0;
	
	for (tmp = outReq->next; tmp; tmp = tmp->next) 
	{
		if (tmp->type != CAP_OBJ_OPER) 
		{
			An_COPY(tmp, &newReq);
			Push(newReq);
			continue;
		}
		result = (capOutReq_t *) calloc(1, sizeof(capOutReq_t));
		i = An_LookupOpNum(tmp->obj.opType);
		if (i < 0) An_CalcError("Unrecognized operator");
		if (an_OpTbl[i].func == NULL) An_CalcError("Operator not yet supported");
		Typecheck(i);
		
		switch(an_OpTbl[i].argCt) 
		{
			case NO_ARG:
			status = (*an_OpTbl[i].func)(result);
			break;
			case ONE_ARG_RESULT_POINTER:
			arg1 = Peek(0);
			FreeOutReq(result);
			status = (*an_OpTbl[i].func)(arg1, &result);
			FreeOutReq(Pop());
			break;
			case ONE_ARG: /* This case is called as many times as necessary */
			arg1 = Peek(0);
			status = (*an_OpTbl[i].func)(arg1, result);
			if (spSyn)
			{
				isVoltage = 1;
				printHeader[countV] = arg1->xName;
				countV++;
			}
			FreeOutReq(Pop());
			break;
			case TWO_ARG:
			arg2 = Peek(0);
			arg1 = Peek(1);
			status = (*an_OpTbl[i].func)(arg1, arg2, result);
			if (spSyn)
			{
				isCurrent = 1;
				printHeader[countV] = arg1->xName;
				countV++;
			}
			FreeOutReq(Pop());
			FreeOutReq(Pop());
			break;
			case TWO_ARG_RESULT_POINTER:
			arg2 = Peek(0);
			arg1 = Peek(1);
			FreeOutReq(result);
			status = (*an_OpTbl[i].func)(arg1, arg2, &result);
			FreeOutReq(Pop());
			FreeOutReq(Pop());
			break;
			case THREE_ARG:
			arg3 = Peek(0);
			arg2 = Peek(1);
			arg1 = Peek(2);
			status = (*an_OpTbl[i].func)(arg1, arg2, arg3, result);
			FreeOutReq(Pop());
			FreeOutReq(Pop());
			FreeOutReq(Pop());
			break;
			case FOUR_ARG:
			arg4 = Peek(0);
			arg3 = Peek(1);
			arg2 = Peek(2);
			arg1 = Peek(3);
			status = (*an_OpTbl[i].func)(arg1, arg2, arg3, arg4, result);
			FreeOutReq(Pop());
			FreeOutReq(Pop());
			FreeOutReq(Pop());
			FreeOutReq(Pop());
			break;
			case VAR_ARG:
			status = (*an_OpTbl[i].func)(result);
			break;
			default:
			An_CalcError("funny internal error in DoOneRequest");
		}
		if (an_OpTbl[i].retVal) Push(result);
	}
	return 0; 
}

int An_initWatch(void)
{
	char **sym_A;
	char *chars_P;
	char **Next_sym_A;
	int type;
	generic_t gval;
	St_SListTable(capOutputT_P, &sym_A, &chars_P);
	if(sym_A) 
	{
		for (Next_sym_A = sym_A; *Next_sym_A; Next_sym_A++) 
		{
			St_GetSym(capOutputT_P, *Next_sym_A, &type, &gval);
			thisRequestName = *Next_sym_A;
			/*       if(strncmp(requestName,"watch",5)
			{
				DoOneInitWatch((capOutReq_t *) gval.v);
				set up up gnuplot pipes
			} */
		}
		St_SListTableFree( &sym_A, &chars_P); /* Free temprary memory usage */
	}
	return 0;
}

int An_updateWatch(void)
{
	/* WriteObj */
	return 0;
}

int An_closeWatch(void)
{
	return 0;
}



int An_DoOutput(void)
{
	char **sym_A;
	char **Next_sym_A;
	char *chars_P;
	int type;
	generic_t gval;
	
	St_SListTable(capOutputT_P, &sym_A, &chars_P);
	if(sym_A) {
		for (Next_sym_A = sym_A; *Next_sym_A; Next_sym_A++) 
		{
			St_GetSym(capOutputT_P, *Next_sym_A, &type, &gval);
			thisRequestName = *Next_sym_A;
			DoOneRequest((capOutReq_t *) gval.v);
		}
		St_SListTableFree(&sym_A, &chars_P); /* Free temporary memory usage */
	}
	return 0;
}




int An_DoEndOutput(void)
{
	char **sym_A;
	char **Next_sym_A;
	char *chars_P;
	int type;
	generic_t gval;
	
	St_SListTable(capEndOutputT_P, &sym_A, &chars_P);
	if(sym_A) 
	{
		for (Next_sym_A = sym_A; *Next_sym_A; Next_sym_A++) 
		{
			St_GetSym(capEndOutputT_P, *Next_sym_A, &type, &gval);
			thisRequestName = *Next_sym_A;
			DoOneRequest((capOutReq_t *) gval.v);
		}
		St_SListTableFree(&sym_A, &chars_P); /* Free temporary memory usage */
	}
	return 0;
	
}


/*
* An_Op_WRITE--
*
* Dump what is on the stack to the output.
*/
int An_Op_WRITE(capOutReq_t *arg, capOutReq_t *result)
{
	char *cwd;
	FILE *f;
	static char *currentOutputFileName;
	capOutReq_t *arg1;
	
	currentOutputFileName = arg->obj.s;
	if ((cwd = getcwd(NULL, 164)) == NULL) 
	{
		perror("pwd");
		exit(1);
	}
	f = fopen(currentOutputFileName, "w");
	if (!f) 
	{
		An_CalcError("Could not open data file");
		return 0;
	}
	Pop(); 
	
	/* Get output data from stack. */
	while ((arg1=Pop())->type != CAP_OBJ_NONEXISTENT) 
	{
		WriteObj(f, arg1);
		FreeOutReq(arg1);
	}
	
	fprintf(output_F,"writing data file %s.\n",currentOutputFileName);
	printf("writing data file %s.\n",currentOutputFileName);
	
	fclose(f);
	
	FreeOutReq(arg1);
	FreeOutReq(arg);
	return(0);
}

