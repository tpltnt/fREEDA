/*
 * variables.h
 * 
 * Header file for the frequency sweep functions 
 */

extern int variableSequence;  /* This is the sequence number for sweeps and
                               * expressions */

void v_addSweep(double initial, double final, double step, int idVal,
		int idType, char *name);
double v_addExpression(char *expression, int idVal, int idType,
                     char *name);
void v_updateExpressionIds(int idVal, int idType);
void v_evaluateExpressions();
int v_evaluateOneExpression(expression_Pt entry_P);
void v_dumpSweeps();
void v_dumpExpressions();
void v_updateVariables(double current_value, usageId_Pt usage_P);
int v_setupNextSweep();
void v_resetSweeps();
void v_freeExpressions();
void v_freeSweeps();

void dumpUsage(FILE *output_F, usageId_Pt usage_P);
void v_updateSweepIds(int idVal, int idType);
