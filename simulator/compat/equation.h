/* equation.h */
#define		STACK_SIZE	50
extern st_Table_Pt PARAMSymT_P;
extern st_Table_Pt capOutputT_P, capOptionsT_P;
extern int expr(char *string, double *result, int *charsRead);
extern void init_equation(void);
extern void free_equation(void);
extern void Eq_DumpSymbolTable(FILE *output_F, st_Table_Pt st_P);
