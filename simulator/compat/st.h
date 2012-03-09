/*
 * st.h
 *	
 * definitions for the symbol table module
 *
 */

#define ST_MAX_VECTOR_LENGTH    1024

/* The definitions below are used to record information required to handle
 * the variable quantities:
 * sweeps   variables   and   expression handling
 *
 * A linked list of type  usageId_t is used to record every usage of the
 * variable quantities.
 * Special notes
 *
 * variables
 * A variable can be used in many places and so usageId_t is a linked list.
 * All variables are in the capOptionsT_P table and a variable is extracted
 * using the symbol table handling routines on capOptionsT_P
 *  capOptionsT_P   points into usageId_t
 *            which points into usageId_t
 *            which points into usageId_t
 *                  terminated by a null ?.next_P
 **
 * sweeps
 * Sweeps are recorded in a linked list of type sweep_t The top of the
 * linked list is pointed by sweepList_head.
 * The linked list is terminated by a null ?.next_P
 *  sweepList_head points into sweep_t
 *           which points into sweep_t
 *           which points into sweep_t
 *           terminated by a null ?.next_P
 * Each entry sweep_t points into usageId_t but there is only one entry here.
 *
 *
 * expressions
 * Expressions are recorded in a linked list of type expression_t. The top of
 * the linked list is pointed by expressionList_head.
 * The linked list is terminated by a null ?.next_P
 *  expressionList_head points into expression_t
 *           which points into expression_t
 *           which points into expression_t
 *           terminated by a null ?.next_P
 * Each entry expression_t points into usageId_t but there is only one entry
 * here.
 *
 *
 * The appropriate procedure for dereferencing is 
 * Go through the list of variables evaluating expressions as they are
 * encountered and updating the locations indicated in usageId_t.
 * Then evaluate the expressions that have not been evaluated are handled
 * and the locations indicated in usageId_t updated.
 */


/* The usageId_t data structure is used to record where information isued.
 * It is part of a linked list 
 */
typedef struct usageIdStruct {
  int     val,   /* Index for location where information is used.
		  * For  ?.type == USAGE_ID_TERM_TYPE
		  || USAGE_ID_ELEMENT_TYPE
		  * this is just id.val.  In this case
		  * gr_id_t == (usage_Id_t.val, usage_Id_t.type);
		  * Otherwise ?.type = USAGE_ID_OPTIONS_TYPE 
		  * and this is not defined.
		  */
    type;  /* Type of place where information is used.
	    * gr_id_t == (usage_Id_t.val, usage_Id_t.type);
	    */
  char	*tableName;  /* Name of the symbol table where it is used
		      * Only used when type = USAGE_ID_SYMBOL_TABLE 
		      * 
		      * If the type is USAGE_ID_ANALYSIS_TYPE, this
		      * is the name of the analysis type. */
  char	*name;  /* Name of a variable if ?.type =
		 * USAGE_ID_OPTIONS_TYPE (defined in a .OPTIONS
		 * statement) or of an element parameter otherwise. */
  struct usageIdStruct *next_P; /* Pointer to the next structure in
				 * chain.  This has meaning only for variables and not
				 * for sweeps and expressions. */
} usageId_t, *usageId_Pt;


/*The following defines are used to indicate where the information is used.  */
#define USAGE_ID_OPTIONS_TYPE	-1
/* The following defines must not be changed. */
#define USAGE_ID_NO_TYPE	GR_ID_NO_TYPE
#define USAGE_ID_TERM_TYPE	GR_ID_TERM_TYPE
#define USAGE_ID_ELEMENT_TYPE	GR_ID_NODE_TYPE
#define USAGE_ID_SYMBOL_TABLE	1000
#define USAGE_ID_ANALYSIS_TYPE	10000



/*
 * Define storage for sweep information.
 */
typedef struct sweep_str{
  double	initial,  
    final,
    step,
    current;	/* Current value during sweep */
  int	used;		/* 0 if never set, 1 if set */
  int	sequence;	/* Sequence index shared with
			 * expressions */
  usageId_t id;		/* This is the id of the variable/parameter
			 * where the sweep information is first
			 * defined. ?.next_P will always be NULL
			 */
  struct sweep_str	*next_sweep; /* next structure in linked list */
} sweep_t, *sweep_Pt;




/*
 * Define storage for expressions.
 * All in-line expressions are stored.  The place where the expression is
 * used is stored in ?.id
 */
typedef struct expression_str{
  char *expression;  	/* The expression */
  double	current;	/* The current value of the expression */
  int	used;		/* 0 if expression has not been evaluated
			 * 1 if expression has been evaluated */
  int	sequence;	/* Sequence index shared with
			 * expressions */
  usageId_t id;		/* This is the id of the variable/parameter
			 * where the expression is used
			 */
  struct expression_str
  *next_expression; /* next structure in linked list */
} expression_t, *expression_Pt;

extern expression_Pt expressionList_head; /* Defined in transim.c */


typedef struct sweep_expression_str {
  union  {
    expression_t regular_expression;
    sweep_t swp_expression;
  } expression_to_be_evaluated;
  struct sweep_expression_str *next_expression_to_be_evaluated;
} sweep_exp_t, *sweep_exp_Pt;

/*
 * A single type of symbol table is used to store all sorts of information.
 * It is used by the parser to store evrything about a line of input before
 * it is processed.
 * It is used to store option parameters.
 */

/*
 * Struct for individual entries of the symbol tables.
 */

typedef	struct st_EntryStruct {

  char	name[MAX_NAME_LEN];	/* symbol string		    */
  int		type,			/* symbol type			    */
    allocHere;		/* TRUE if symbol data alloc'ed	    */
  /* by symtab			    */
  struct st_EntryStruct *next_P;	/* next entry in chain		    */
  generic_t	gval;			/* symbol value			    */
  usageId_Pt  usage_P;		/* If this is a variable keep track
				 * of where it is used */
  int         marker;			/* marks if symbol has been read    */
} st_Entry_t, *st_Entry_Pt;

/*
 * Struct for table header
 */

typedef struct {

  char	name[MAX_NAME_LEN];	
  /* table name			    */
  int		size,		    /* size of table to be allocated	    */
    entryCt;	    /* # of entries in table currently	    */
  st_Entry_Pt	*tab_A;		    /* vector of ptrs to entries	    */

} st_Table_t, *st_Table_Pt;

#define ST_COULDNT_COERCE -3
#define ST_SYM_READ_FMT_ERR -2
#define ST_AT_EOF	    -1
#define	ST_SYM_FOUND	    0
#define ST_OK		    0
#define ST_SYM_NOT_FOUND    1
#define ST_NULL_TABLE	    2

#define ST_MAX_INPUT_LEN    MAX_STRING_LEN

#define ST_ALLOW_STRING_TO_INT 0x0001
#define ST_ALLOW_STRING_TO_FLOAT 0x0002
#define ST_ALLOW_STRING_TO_BOOLEAN 0x0004
#define ST_ALLOW_STRING_TO_CHAR 0x0008

#define ST_ALLOW_INT_TO_STRING 0x0010
#define ST_ALLOW_FLOAT_TO_STRING 0x0020
#define ST_ALLOW_CHAR_TO_STRING 0x0040

#define ST_ALLOW_FLOAT_TO_INT 0x0080


/*
 * GLOBAL functions defined in st_symtab.c
 */

extern st_Table_Pt St_NewTable(const char *name, int size);
extern void St_DelTable(st_Table_Pt st_P);
extern int St_GetSym(st_Table_Pt st_P, const char *name, int *type_P, 
  generic_Pt gval_P);
extern int St_DefSym(st_Table_Pt st_P, const char *name, int type, 
		     generic_t gval);
extern int St_ReplSym(st_Table_Pt st_P, const char *name, int type,
		      generic_t gval);
extern int St_DelSym(st_Table_Pt st_P, const char *name);
extern int St_MarkSym(st_Table_Pt st_P, const char *name);
extern st_Table_Pt St_GetTable(char *name);
extern void St_CleanUpTable(st_Table_Pt st_P);
extern void St_ListTable(st_Table_Pt st_P, char ***list_A, char **chars_P);
extern void St_SListTable(st_Table_Pt st_P, char ***list_A, char **chars_P);

extern void St_SListTableFree(char ***list_A, char **chars_P);
extern void St_ListTableFree(char ***list_A, char **chars_P);

extern void St_SetCoerceMask(int mask);
extern int St_GetCoerceMask(void);

extern int St_GetSymAsInt(st_Table_Pt st_P, const char *name, int *i_P);
extern int St_GetSymAsBoolean(st_Table_Pt st_P, const char *name, int *i_P);
extern int St_GetSymAsChar(st_Table_Pt st_P, const char *name, int *i_P);
extern int St_GetSymAsLong(st_Table_Pt st_P, const char *name, long *l_P);
extern int St_GetSymAsFloat(st_Table_Pt st_P, const char *name, float *f_P);
extern int St_GetSymAsDouble(st_Table_Pt st_P, const char *name, double *d_P);
extern int St_GetSymAsString(st_Table_Pt st_P, const char *name, char *s);
extern int St_GetSymAsVector(st_Table_Pt st_P, const char *name,
			     double **dP_P, int *length);
extern st_Entry_Pt LookupSym(st_Table_Pt st_P, const char *name);
extern int St_CheckEmpty(FILE *output_F, st_Table_Pt st_P);
extern int St_CheckAllMarked(FILE *output_F, st_Table_Pt st_P);
extern int St_GetCount(st_Table_Pt st_P);
extern void St_DelAllSym(st_Table_Pt st_P);

/* Case insensitive string operations */
extern int strcmpCI(char *s1, char *s2);

/* Copy string but first allocate space and then check memory */
void strcpyCareFul(char **s1, const char *s2);

extern void ml_CleanUp();

