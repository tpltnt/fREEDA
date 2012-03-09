/*
 * variables.c
 *
 * routines for variable quantities
 * sweep   	vaiable usage		expressions
 *
 * Author:
 * Michael B. Steer
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "parser.h"
#include "report.h"
#include "../compat/equation.h"

int	variableSequence;  /* This is the sequence number for sweeps and
			    * expressions */

/*
 * v_addSweep
 *
 * Add n new sweep to the bottom of the linked list and set sequence.
 *
 * Parameters:
 * idVal	see st.h
 * idType	type (where it is used)
 */
void v_addSweep(double initial, double final, double step, int idVal,
		int idType, char *name)
{
  /*    generic_t gval; */
  char *s;
  int length;
  sweep_Pt entry_P, lastEntry_P; /*  , previousEntry_P; */
  
  /* Find end of expression linked list */
  for (lastEntry_P = NULL, entry_P = sweepList_head;
       entry_P; lastEntry_P = entry_P, entry_P = entry_P->next_sweep);
  entry_P = malloc(sizeof(sweep_t));
  if (lastEntry_P == NULL)
    sweepList_head = entry_P;
  else
    lastEntry_P->next_sweep = entry_P;
  entry_P->next_sweep = NULL;
  entry_P->sequence = variableSequence++;
  entry_P->initial = initial;
  entry_P->final = final;
  entry_P->step = step;
  entry_P->current = initial;
  entry_P->id.val = (short) idVal;
  entry_P->id.type = (short) idType;

  /* Store name */
  length = strlen(name);
  s = (char *) malloc(length + 1);
  strcpy(s, name);
  entry_P->id.name = s;
  entry_P->id.next_P = NULL;
  
  /* Now update the sweep variable with its initial value */
  /*gval.d = initial; 
    St_ReplSym(capOptionsT_P, s, entry_P->id.type, gval);*/
  return;
  free(s);
}

/*
 * v_addExpression
 *
 * Add an new expression to the bottom of the linked list and set sequence.
 *
 * If idType is null then the id can be updated subsequently by 
 * v_udateEpressionId
 * This is often the case when setting up an element.
 *
 * Parameters:
 * expression	string containing expression
 * idVal	see st.h
 * idType	type (where it is used)
 * 
 * Return current value of expression
 */
double v_addExpression(char *expression, int idVal, int idType, char *name)
{
  double value;
  char *s1, *s2;
  int length;
  expression_Pt entry_P, lastEntry_P;
  
  /* Find end of expression linked list */
  for (lastEntry_P = NULL, entry_P = expressionList_head;
       entry_P;
       lastEntry_P = entry_P, entry_P = entry_P->next_expression);
  
  entry_P = (expression_Pt) malloc(sizeof(expression_t));
  if (lastEntry_P == NULL)
    expressionList_head = entry_P;
  else
    lastEntry_P->next_expression = entry_P;
  entry_P->next_expression = NULL;
  entry_P->sequence = variableSequence++;
  
  /* Store expression */
  length = strlen(expression);
  s1 = (char *) malloc(length + 1);
  strcpy(s1, expression);
  entry_P->expression = s1;
  entry_P->id.val = (short) idVal;
  entry_P->id.type = (short) idType;
  
  /* Store name */
  length = strlen(name);
  s2 = (char *) malloc(length + 1);
  strcpy(s2, name);
  entry_P->id.name = s2;
  entry_P->id.next_P = NULL;

  /* Evaluate expression */
  v_evaluateOneExpression(entry_P);
  value = entry_P->current;
  return(value);
}

/*
 * v_updateExpressionIds
 *
 * Goes through all expressions and if the usage type is USAGE_ID_NO_TYPE
 * the type av val of usage is updated.
 * This is often required in setting uping up an element.
 *
 * Parameters:
 * idVal	see st.h
 * idType	type (where it is used)
 */
void v_updateExpressionIds(int idVal, int idType)
{
  /*    generic_t gval; */
  expression_Pt entry_P;
  usageId_Pt usage_P;
  
  /* Go through all expressions checking the usage entry */
  for (entry_P=expressionList_head; 
       entry_P; entry_P=entry_P->next_expression) {
    for (usage_P=&(entry_P->id); usage_P; usage_P=usage_P->next_P) {
      if(usage_P->type == USAGE_ID_NO_TYPE)
	usage_P->type = (short) idType;
      usage_P->val = (short) idVal;
    }
  }
}

/*
 * evaluateExpressions
 *  Goes through all of the expressions and evaluates them by calling
 *  expr in equation.c
 */
void v_evaluateExpressions()
{
  // generic_t gval;
  expression_Pt entry_P;
  /*    usageId_Pt usage_P; */
  /*    int type; */
  // int charsRead;
  
  for (entry_P=expressionList_head; 
       entry_P; entry_P=entry_P->next_expression) {
    v_evaluateOneExpression(entry_P);
    // expr(entry_P->expression,&entry_P->current,&charsRead);
    // entry_P->used = 1; /* To indicate that it has been evaluated */
    /* Now update the evaluated expression in the table */
    // gval.d = entry_P->current;
    // St_ReplSym(capOptionsT_P, entry_P->id.name, GEN_TYPE_DOUBLE, gval); 
  }
}


/*
 * evaluateOneExpression
 *  Evaluate a single expression by calling expr in equation.c
 */
int v_evaluateOneExpression(expression_Pt entry_P)
{
  generic_t gval;
  int charsRead;
  //printf("++++evaluateOneExpression enter\n");

  expr(entry_P->expression,&entry_P->current,&charsRead);
  //printf("++++evaluateOneExpression after expr\n");
  entry_P->used = 1; /* To indicate that it has been evaluated */
  /* Now update the evaluated expression in the table */
  gval.d = entry_P->current;
  St_ReplSym(capOptionsT_P, entry_P->id.name, GEN_TYPE_DOUBLE, gval); 
  return(0);
}


/*
 * v_updateSweepIds
 *
 * Goes through all sweeps and if the usage type is USAGE_ID_NO_TYPE
 * the type av val of usage is updated.
 * This is often required in setting uping up an element.
 *
 * Parameters:
 * idVal	see st.h
 * idType	type (where it is used)
 */
void v_updateSweepIds(int idVal, int idType)
{
  /*    generic_t gval; */
  sweep_Pt entry_P;
  usageId_Pt usage_P;
  
  /* Go through all sweeps checking the usage entry */
  for (entry_P=sweepList_head; entry_P; entry_P=entry_P->next_sweep) {
    for (usage_P=&(entry_P->id); usage_P; usage_P=usage_P->next_P) {
      if(usage_P->type == USAGE_ID_NO_TYPE)
	usage_P->type = (short) idType;
      usage_P->val = (short) idVal;
    }
  }
}


void v_dumpSweeps()
{
  sweep_Pt entry_P;
  if(sweepList_head)
    {
      fprintf(output_F, "List of sweeps:\n");
      fprintf(output_F, "sequence current initial final step\n");
      for (entry_P=sweepList_head; entry_P; entry_P=entry_P->next_sweep) {
	fprintf(output_F, 
		"%d\t%g\t%g\t%g\t%g\n",entry_P->sequence,entry_P->current,
		entry_P->initial, entry_P->final, entry_P->step);
	dumpUsage(output_F, &(entry_P->id));
      }
    }
  else
    fprintf(output_F, "No sweeps\n");
  return;
}


void v_dumpExpressions()
{
  expression_Pt entry_P;
  if(expressionList_head)
  {
    fprintf(output_F, "\n***\nLIST OF EXPRESSIONS:\n");
    fprintf(output_F, "sequence current name\t expression\n");
    for (entry_P=expressionList_head; 
        entry_P; entry_P=entry_P->next_expression)
    {
      fprintf(output_F, "%d\t %g\t %s\t%s\n",
          entry_P->sequence, entry_P->current, entry_P->id.name,
          entry_P->expression);
      dumpUsage(output_F, &(entry_P->id));
    }
  }
  else
    fprintf(output_F, "No expressions\n");

  fprintf(output_F, "***\n");
  return;
}


/* Set up sweeps as though they have never been used */
void v_resetSweeps() 
{
  sweep_Pt entry_P;
  for (entry_P=sweepList_head; entry_P; entry_P=entry_P->next_sweep)
    entry_P->used = 0;
}


/*
 * v_freeExpressions
 *
 * Free expressions and any related information
 */
void v_freeExpressions()
{
  expression_Pt nextEntry_P, entry_P;
  for (entry_P=expressionList_head; entry_P; entry_P=nextEntry_P)
    {
      free(entry_P->expression);
      free(entry_P->id.name);
      nextEntry_P=entry_P->next_expression;
      free(entry_P);
    }
  return;
}

/*
 * v_freeSweeps
 *
 * Free sweeps and any related information
 */
void v_freeSweeps()
{
  sweep_Pt nextEntry_P, entry_P;
  for (entry_P=sweepList_head; entry_P; entry_P=nextEntry_P)
    {
      free(entry_P->id.name);
      nextEntry_P=entry_P->next_sweep;
      free(entry_P);
    }
  return;
}



/*
 * dumpUsage
 *	Dump information about usage.
 */
void dumpUsage(FILE *output_F, usageId_Pt usage_P)
{
  /*    gr_Id_t		id; */
  /*
if(usage_P != NULL)
  {
  fprintf(output_F, "\t\tUsed in the following places.\n");
  while(usage_P != NULL)
    {
    switch(usage_P->type)
      {
      case USAGE_ID_OPTIONS_TYPE:
        fprintf(output_F,"\tVARIABLE:\t%s\n",usage_P->name);
      break;
      case USAGE_ID_SYMBOL_TABLE:
        fprintf(output_F,"\tSYMBOL TABLE VARIABLE:\t%s\tTABLE: %s\n",
	  usage_P->name, usage_P->tableName);
      break;
      case USAGE_ID_TERM_TYPE:
      case USAGE_ID_ELEMENT_TYPE:
        {
        char name[MAX_NAME_LEN+1];
        id.val = usage_P->val;
        id.type = usage_P->type;
        Gr_GetName(id, name);
        fprintf(output_F,
          "\tELEMENT:\t%s\tPARAMETER:\t%s\n",name,usage_P->name);
        }
      break;
      case USAGE_ID_NO_TYPE:
      case USAGE_ID_SYMBOL_TABLE:
      default:
        REPORT(WARNING_FATAL_INTERNAL," ");
      }
    usage_P=usage_P->next_P;
    }
  }
  */
  return;
}



/*
 * updateUsage
 *
 * Take the current value for a variable, sweep or expression everywhere it is
 * used.
 */
/*void updateUsage(currentValue, usageId_t usage)
{
usageId_Pt usage_P;
generic_t gval;
for (usage_P = &usage; usage_P; usage_P = usage_P->next_P)
  {
  switch(usage_P->type)
    {
    case USAGE_ID_OPTIONS_TYPE:
      St_ReplSym(capOptionsT_P, usage_P->name, GEN_TYPE_DOUBLE, gval);
    break;

    case USAGE_ID_TERM_TYPE:
    case USAGE_ID_ELEMENT_TYPE:
    default:
        {
        REPORT(ERROR_NETLIST,"Usage of sweep, expression, cahngeable variables not supported\ for Nodes, Terms, Edge, ECG or NCG");
        exit(1);
    }
  }
return;
}*/

/*
 * v_setupNextSweep returns 1 if sweep if sweep is to continue 
 *                  returns 0 if sweep is finished
 *
 */
int v_setupNextSweep() 
{
  double final, step, current_value;
  char msg[80];
  
  if (sweepList_head == NULL) return(0);
  final = sweepList_head->final;
  step = sweepList_head->step;
  current_value = sweepList_head->current;

  sepLine();
  sprintf(msg, "--- Sweep variable %s = %g", 
	  sweepList_head->id.name, current_value + step);
  report(MESSAGE, msg);
  
  if (current_value < final)
    { 
      sweepList_head->current = current_value + step; 
      v_updateVariables(sweepList_head->current, &sweepList_head->id); 
      return(1);
    }
  else return(0);
}

/* 
 * v_updateVariables
 * 
 * Update the sweep variable(s) whereever they are used
 * 
 * Only one sweep variable is supported at the moment 
 *
 */
void v_updateVariables(double current_value, usageId_Pt usage_P)
{
  usageId_Pt entry_P;
  st_Entry_Pt st_P;
  generic_t gval;
  int	  type;
  gr_Id_t   id;
  st_Table_Pt symT_P;
  
  gval.d = current_value;
  type = GEN_TYPE_DOUBLE;
  for(entry_P = usage_P; entry_P; entry_P = entry_P->next_P) 
    { 
      switch(entry_P->type)
	{
	case USAGE_ID_OPTIONS_TYPE: 
	  St_ReplSym(capOptionsT_P, entry_P->name, type, gval);  
	  st_P = LookupSym(capOptionsT_P, entry_P->name);
	  if (st_P->usage_P) v_updateVariables(current_value, st_P->usage_P);
	  break;
	case USAGE_ID_ELEMENT_TYPE:
	  id.val = entry_P->val;
	  id.type = entry_P->type;
	  replaceElemParam(id, type, gval, entry_P->name);
	  break;
	case USAGE_ID_ANALYSIS_TYPE:
	  replaceAnalysisParam(entry_P->tableName, gval, entry_P->name);
	  break;
	case USAGE_ID_SYMBOL_TABLE:
	  symT_P = St_GetTable(entry_P->tableName);
	  if(symT_P == NULL) {
	    printf(
		   "ERROR: v_updateVariables: Internal error in handling variables.\n");
	    exit(987);
	  }
	  St_ReplSym(symT_P, entry_P->name, type, gval);  
	  st_P = LookupSym(symT_P, entry_P->name);
	  if (st_P->usage_P) v_updateVariables(current_value, st_P->usage_P);
	  break;
	case USAGE_ID_TERM_TYPE:
	default:
	  fprintf(stderr, 
		  "Usage of sweep, expression, changeable variables not \
supported on Terminals");
	}
    }
  return;  
}
