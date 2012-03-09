/**************************************************************************
 *
 * File Name:
 *       st.c
 *
 * Description:
 *       routines for managing symbol tables 
 *       plus string manipulation routines
 *
 * Author:
 *       Joseph Nathan Hall
 *       Michael B. Steer
 *
 **************************************************************************/

#include <stdio.h>
#ifndef ultrix
#include <stdlib.h>
#endif
#include <string.h>
#include <ctype.h>

#include "standard.h"
#include "mltypes.h"
#include "generic.h"
#include "st.h"
/* #include "cap.h" */

int coercionBits = 0x0000;

/*
 * There's a symbol table containing all the names of the symbol tables,
 * and of course it contains its own name ...
 *
 * "tnt" stands for "Table Name Table" ...
 */ 
static st_Table_Pt  tnt_P = NULL;

expression_Pt expressionList_head = 0;

/*
 * If there are a lot of symbol tables around you might want to increase
 * TNT_SIZE.  For best results use a prime number.
 */

#define	TNT_SIZE    11

/*
 * Names of basic generic types.
 */

#define NAMES 7

static struct {
  char    *name;
  int	    num;
} types[NAMES] =

{{"int",	GEN_TYPE_INT},
 {"long",	GEN_TYPE_LONG},
 {"float",	GEN_TYPE_FLOAT},
 {"double",	GEN_TYPE_DOUBLE},
 {"char",	GEN_TYPE_CHAR},
 {"string",	GEN_TYPE_STRING},
 {"stringarray", GEN_TYPE_STRING_A}};

/*
 * Hash--
 *	Simple, stupid hash function.  Why be great, you know, when you can
 * be adequate so easily?
 */
int Hash(const char *s, int mod)
{

  register int h = 0;

  while (*s)
    h += (unsigned) *s++;

  return(h % mod);

}


/*
 * LookupSym--
 *	Search for a name in a table.  Returns NULL if not found.
 */
st_Entry_Pt LookupSym(st_Table_Pt st_P, const char *name)
{

  st_Entry_Pt	entry_P;
  int		h;

  h = Hash(name, st_P->size);

  for (entry_P = st_P->tab_A[h]; entry_P; entry_P = entry_P->next_P)
    if (!strncmp(entry_P->name, name, MAX_NAME_LEN))
      break;

  return(entry_P);

}



/*
 * AddSym--
 *	Add a name to a table and return a pointer to the new entry.
 *
 * If the name already exists then there will be two (or more) entries of the
 * same name.  The new entry will be the first to be accessed or deleted.
 *
 * So you need to be careful if what you reply want to do is to replace the
 * symbol.  You should use ReplSym instead. 
 */
st_Entry_Pt AddSym(st_Table_Pt st_P, const char *name)
{

  int		h;
  st_Entry_Pt	new_P;

  h = Hash(name, st_P->size);

  new_P = (st_Entry_Pt) malloc(sizeof(st_Entry_t));
  new_P->usage_P = NULL;
  strncpy(new_P->name, name, MAX_NAME_LEN - 1);
  new_P->name[MAX_NAME_LEN - 1] = 0;
  new_P->next_P = st_P->tab_A[h];
  st_P->tab_A[h] = new_P;

  st_P->entryCt++;

  return(new_P);

}



/*
 * St_NewTable--
 *	Create a new symbol table header.  Returns NULL if name isn't
 * unique with respect to the symbol tables currently in existence.
 */
st_Table_Pt St_NewTable(const char *name, int size)
{

  st_Table_Pt	st_P;
  generic_t	gval;

  /*
   * Create the name table if doesn't already exist.  Obviously we
   * can't use St_NewTable for this ...
   */

  if (!tnt_P) {
    /* MEMTODO : tnt_P and tnt_P->tab_A need to be freed */
    tnt_P = (st_Table_Pt) malloc(sizeof(st_Table_t));

    strncpy(tnt_P->name, "_TNT_", MAX_NAME_LEN - 1);
    tnt_P->name[MAX_NAME_LEN - 1] = 0;
    tnt_P->size = TNT_SIZE;
    tnt_P->tab_A = (st_Entry_Pt *) calloc(TNT_SIZE, sizeof(st_Entry_Pt));
    tnt_P->entryCt = 0;

    gval.v = (VOIDPTR) tnt_P;
    St_DefSym(tnt_P, "_TNT_", GEN_TYPE_VOID_PTR, gval);

  }

  /*
   * See if the new table name is unique.
   */

  if (LookupSym(tnt_P, name))
    return(NULL);

  /*
   * Create the new table.
   */

  st_P = (st_Table_Pt) malloc(sizeof(st_Table_t));

  strncpy(st_P->name, name, MAX_NAME_LEN - 1);
  st_P->name[MAX_NAME_LEN - 1] = 0;
  st_P->size = size;
  st_P->tab_A = (st_Entry_Pt *) calloc(size, sizeof(st_Entry_Pt));
  st_P->entryCt = 0;

  /*
   * Add the name of the new table to the "table name table" now.
   * gval.v is a pointer to the new table.
   */

  gval.v = (VOIDPTR) st_P;
  St_DefSym(tnt_P, name, GEN_TYPE_VOID_PTR, gval);

  return(st_P);

}



/*
 * St_DelTable--
 *	Delete a symbol table and associated storage.
 */
void St_DelTable(st_Table_Pt st_P)
{

  st_Entry_Pt entry_P, entry1_P;
  int		i;

  if (!st_P) return;

  St_CleanUpTable(st_P);

  for (i = 0; i < st_P->size; i++)
    for (entry_P = st_P->tab_A[i]; entry_P; entry_P = entry1_P) {
      entry1_P = entry_P->next_P;
      free(entry_P);
    }
  /*  	for (entry_P = st_P->tab_A[i]; entry_P; entry_P = entry_P->next_P) { */
  /*  	    free(entry_P); */
  /*  	} */

  if (strncmp(st_P->name, "_TNT_", MAX_NAME_LEN))
    St_DelSym(tnt_P, st_P->name);

  free(st_P->tab_A);
  free(st_P);

  return;

}



/*
 * St_CleanUpTable--
 * Clean out symbol table. That is set it to its prestine state.
 * The following objects are freed
 *         strings pointed to by symbols in the symbol table
 *         symbols in the symbol table
 *         vectors
 *         vectors of vectors
 */
void St_CleanUpTable(st_Table_Pt st_P)
{
  char **sym_A, **originalSym_A;
  char *chars_P;
  int i, type;
  generic_t gval;

  St_ListTable(st_P, &originalSym_A, &chars_P);
  if(originalSym_A) {
    sym_A = originalSym_A;
    while (*sym_A) /* Delete strings in the symbol table */
      {
	St_GetSym(st_P, *sym_A, &type, &gval);
	switch(type)
	  {
	  case GEN_TYPE_STRING:
	    free(gval.s);  /* free string */
	    break;
    
	  case GEN_TYPE_VOID_PTR:
	    free(gval.v); /* free vector */
	    break;
    
	  case GEN_TYPE_DOUBLE_VECTOR:
	    free(gval.dv.v);
	    break;

	  case GEN_TYPE_DOUBLE_VECTOR_OF_VECTORS:
	    for(i = gval.dvv.size-1; i; i--) {
	      free(gval.dvv.dv_P[i].v);
	    }
	    free(gval.dvv.dv_P);
	    break;

	  default:
	    break;
    
	  }
    
	/* if (type == GEN_TYPE_STRING) free(gval.s);*/ /* Free Strings*/
	sym_A++;
      }
    St_ListTableFree( &originalSym_A, &chars_P); ; /* Free storage */
  }
  return;
}

/*
 * St_ListTable--
 *	Returns an unsorted list of symbols in the table.  The list will be
 * terminated with a NULL pointer.  This routine frees the storage used by
 * the last call to St_ListTable or St_SListTable; a call with a NULL 
 * argument is a convenient way to free storage allocated by a previous call.
 */
void St_ListTable(st_Table_Pt st_P, char ***list_A, char **chars_P)
{

  st_Entry_Pt entry_P;
  int		i, j;

  if (!st_P) return;
  if(st_P->entryCt) {

    *list_A = (char **) malloc(sizeof(char *) * (st_P->entryCt + 1));
    *chars_P = (char *) malloc(sizeof(char) * MAX_NAME_LEN * st_P->entryCt);

    for (i = 0; i < st_P->entryCt; i++)
      (*list_A)[i] = *chars_P + MAX_NAME_LEN * i;
    (*list_A)[st_P->entryCt] = NULL;

    j = 0;
    for (i = 0; i < st_P->size; i++)
      for (entry_P = st_P->tab_A[i]; entry_P; entry_P = entry_P->next_P)
	strcpy((*list_A)[j++], entry_P->name);

    (*list_A)[st_P->entryCt] = NULL;
  }
  else {
    *list_A = NULL;
    *chars_P = NULL;
  }


  return;

}


/*
 * St_ListTableFree--
 *	Free temporary Storage used by previous call to St_ListTable.
 */
void St_ListTableFree(char ***list_A, char **chars_P)
{
  if (list_A && *list_A) {
    free(*list_A);
    *list_A = NULL;
  }
  if (chars_P && *chars_P) {
    free(*chars_P);
    *chars_P = NULL;
  }
}

/*
 * SCompare--
 *	String compare function for qsort below.
 */
/*static int SCompare(VOIDPTR s1, VOIDPTR s2) */
static int SCompare(CONST_VOIDPTR s1, CONST_VOIDPTR s2)
{
  return(strcmp(s1, s2));
}


/*
 * St_SListTable--
 *	Returns a sorted list of symbols in a table.  Otherwise is exactly
 * like St_ListTable, above.
 */
void St_SListTable(st_Table_Pt st_P, char ***list_A, char **chars_P)
{
  if(st_P->entryCt) {
    St_ListTable(st_P, list_A, chars_P);
    qsort(**list_A, st_P->entryCt, sizeof(char) * MAX_NAME_LEN, SCompare);
  }
  else {
    *list_A = NULL;
  }

  return;
}



/*
 * St_SListTableFree--
 *	Free temporary Storage used by previous call to St_SListTable.
 */
void St_SListTableFree(char ***list_A, char **chars_P)
{
  St_ListTableFree(list_A, chars_P); ; /* Free storage */
}
    

/* 
 * St_GetSym--
 *	Look for a symbol in a table.  Return type and ptr to val if found.
 */
int St_GetSym(st_Table_Pt st_P, const char *name, int *type_P, generic_Pt gval_P)
{

  st_Entry_Pt	entry_P;

  if (!st_P)
    return(ST_NULL_TABLE);

  if (!(entry_P = LookupSym(st_P, name)))
    return(ST_SYM_NOT_FOUND);

  *type_P = entry_P->type;
  *gval_P = entry_P->gval;

  return(ST_SYM_FOUND);

}


/* 
 * St_DefSym--
 *	Add a symbol to a table.  Returns ST_SYM_FOUND and does nothing if
 * name is already in table.
 */
int St_DefSym(st_Table_Pt st_P, const char *name, int type, generic_t gval)
{

  st_Entry_Pt entry_P;

  if (!st_P)
    return(ST_NULL_TABLE);

  if (LookupSym(st_P, name))
    return(ST_SYM_FOUND);

  entry_P = AddSym(st_P, name);

  /*
   * Assign data.  gval, if it's a pointer, points to an object that wasn't
   * allocated by symtab, so allocHere is set to FALSE.
   */

  entry_P->type = type;
  entry_P->gval = gval;
  entry_P->allocHere = FALSE;

  return(ST_SYM_NOT_FOUND);

}


/*
 * FreeGval--
 *	If gval is a pointer, frees space used by the object(s) it points to.
 */
void FreeGval(int type, generic_t gval)
{

  char **s_A;

  switch (type) {
  case GEN_TYPE_STRING:
  case GEN_TYPE_VOID_PTR:
    if (gval.v)
      free(gval.v);
    break;

  case GEN_TYPE_STRING_A:
    if (gval.v) {
      for (s_A = gval.s_A; *s_A; s_A++)
	free((VOIDPTR ) (*s_A++));
      free(gval.v);
    }
    break;
  }

  return;

}


/* 
 * St_ReplSym--
 *	Add or supersede a symbol in a table.
 */
int St_ReplSym(st_Table_Pt st_P, const char *name, int type, generic_t gval)
{

  st_Entry_Pt entry_P;
  int		status;

  if (!st_P)
    return(ST_NULL_TABLE);

  if (!(entry_P = LookupSym(st_P, name))) {
    entry_P = AddSym(st_P, name);
    status = ST_SYM_NOT_FOUND;
  } else {
    status = ST_SYM_FOUND;
    if (entry_P->allocHere)
      FreeGval(entry_P->type, entry_P->gval);
  }

  /*
   * Assign data.  gval, if it's a pointer, points to an object that wasn't
   * allocated by symtab, so allocHere is set to FALSE.
   */

  entry_P->type = type;
  entry_P->gval = gval;
  entry_P->allocHere = FALSE;

  return(status);

}



/*
 * St_DelSym--
 *	Delete a symbol from the table.
 */
int St_DelSym(st_Table_Pt st_P, const char *name)
{

  st_Entry_Pt	entry_P, last_P;
  int		h;

  if (!st_P)
    return(ST_NULL_TABLE);

  h = Hash(name, st_P->size);

  for (last_P = NULL, entry_P = st_P->tab_A[h]; entry_P; 
       last_P = entry_P, entry_P = entry_P->next_P)
    if (!strncmp(entry_P->name, name, MAX_NAME_LEN))
      break;

  if (!entry_P)
    return(ST_SYM_NOT_FOUND);

  if (last_P)
    last_P->next_P = entry_P->next_P;
  else
    st_P->tab_A[h] = NULL;

  /*
   * If gval is a pointer and its data was allocated by symtab, then
   * free the data.
   */

  if (entry_P->allocHere)
    FreeGval(entry_P->type, entry_P->gval);

  free(entry_P);
  st_P->entryCt--;

  return(ST_SYM_FOUND);

}


/*
 * St_MarkSym--
 *	Mark a symbol as read in the table.
 */
int St_MarkSym(st_Table_Pt st_P, const char *name)
{

  st_Entry_Pt	entry_P, last_P;
  int		h;

  if (!st_P)
    return(ST_NULL_TABLE);

  h = Hash(name, st_P->size);

  for (last_P = NULL, entry_P = st_P->tab_A[h]; entry_P; 
       last_P = entry_P, entry_P = entry_P->next_P)
    if (!strncmp(entry_P->name, name, MAX_NAME_LEN))
      break;

  if (!entry_P)
    return(ST_SYM_NOT_FOUND);
  else {
    entry_P->marker = 1;
    return(ST_SYM_FOUND);
  }

}

    
/*
 * St_GetTable--
 *	Get a table by name
 */
st_Table_Pt St_GetTable(char *name)
{

  /*      st_Table_Pt	st_P; */
  int		type;
  generic_t	gval;

  if (!tnt_P)
    return(NULL);

  if (St_GetSym(tnt_P, name, &type, &gval) != ST_SYM_FOUND)
    return(NULL);

  return((st_Table_Pt) gval.v);

}



/*
 * St_ReadSym--
 *	Read a symbol from a line of the stream specified by input_F into
 * the specified table.  May read several lines of the stream if the
 * data type is a string array.  Symbols read replace existing ones of
 * the same names.
 */
int St_ReadSym(FILE *input_F, st_Table_Pt st_P)
{

  char    s[ST_MAX_INPUT_LEN], name[ST_MAX_INPUT_LEN], s1[ST_MAX_INPUT_LEN],
    *p, *stmp, c;
  int	    i, done, type, status;
  generic_t g;
  st_Entry_Pt entry_P;

  if (!st_P)
    return(ST_NULL_TABLE);

  /*
   * Read two "words" into name and s1.  s1's case is folded.  The pointer
   * p is left at the start of the third word.  name may be enclosed
   * in quotes.
   */

  do {
    if (!fgets(s, ST_MAX_INPUT_LEN, input_F)) {
      return(ST_AT_EOF);
    }
  } while (sscanf(s, " %c", &c) == EOF);

  for (p = s; *p && isspace(*p); p++)
    ;

  if (*p == '\"') {
    p++;
    for (i = 0; *p && (*p != '\"'); p++, i++)
      name[i] = *p;
    if (*p == '\"')
      p++;
  } else {
    for (i = 0; *p && !isspace(*p); p++, i++)
      name[i] = *p;
  }
  name[i] = 0;

  for (; *p && isspace(*p); p++)
    ;

  for (i = 0; *p && !isspace(*p); p++, i++)
    s1[i] = tolower(*p);
  s1[i] = 0;

  for (; *p && isspace(*p); p++)
    ;

  /*
   * Match type name
   */

  type = 0;
  for (i = 0; i < NAMES; i++) {
    if (!strcmp(s1, types[i].name)) {
      type = types[i].num;
      break;
    }
  }

  /*
   * Exit early if type wasn't recognized.
   */

  if (!type) {
    return(ST_SYM_READ_FMT_ERR);
  }

  /*
   * Now, read in the value.
   */

  status = ST_OK;

  switch (type) {
  case GEN_TYPE_INT:
    g.i = 0;
    sscanf(p, "%d", &g.i);
    break;
  case GEN_TYPE_LONG:
    g.l = 0;
    sscanf(p, "%ld", &g.l);
    break;
  case GEN_TYPE_FLOAT:
    g.f = 0;
    sscanf(p, "%f", &g.f);
    break;
  case GEN_TYPE_DOUBLE:
    g.d = 0;
    sscanf(p, "%lf", &g.d);
    break;
  case GEN_TYPE_CHAR:
    i = 0;
    sscanf(p, "%d", &i);
    g.c = i;
    break;
  case GEN_TYPE_STRING:
    *s1 = 0;
    sscanf(p, "%*[^\"]%*c%[^\"]", s1);
    g.s = (char *) malloc(strlen(s1) + 1);
    strcpy(g.s, s1);
    break;
  case GEN_TYPE_STRING_A:
    for (done = FALSE, i = 0, stmp = NULL, g.s_A = NULL; !done;) {
      c = 0;
      *s1 = 0;
      if (sscanf(p, " %*[^\"]%*c%[^\"]%*c %c", s1, &c) == EOF) {
	char *p1;
	for (p1 = p; *p1 && isspace(*p1); p1++)
	  ;
	if (strncmp(p1, "\"\"", 2))
	  return(ST_SYM_READ_FMT_ERR);
      }
      if (!stmp) {
	stmp = (char *) malloc(strlen(s1) + 1);
	strcpy(stmp, s1);
      } else {
	stmp = (char *) realloc(stmp, strlen(stmp) + strlen(s1) + 1);
	strcat(stmp, s1);
      }
      if (c == ',') {
	if (!g.s_A) {
	  g.s_A = (char **) malloc(sizeof(char *));
	} else {
	  g.s_A = (char **) realloc(g.s_A, sizeof(char *) * (i + 1));
	}
	g.s_A[i++] = stmp;
      }
      if (*stmp) {
	if (c == ',')
	  stmp = NULL;
	do {
	  if (!fgets(s, ST_MAX_INPUT_LEN, input_F)) {
	    if (!g.s_A) {
	      g.s_A = (char **) malloc(sizeof(char *));
	    } else {
	      g.s_A = (char **) 
		realloc(g.s_A, sizeof(char *) * (i + 1));
	    }
	    g.s_A[i] = NULL;
	    status = ST_AT_EOF;
	    done = TRUE;
	  }
	} while (sscanf(s, "%c", &c) == EOF);
	p = s;
      }
      else
	done = TRUE;
    }
    break;
  }

  /*
   * Add the symbol to the table
   */

  if (!(entry_P = LookupSym(st_P, name)))
    entry_P = AddSym(st_P, name);

  /*
   * If gval is a pointer and its data was allocated by symtab, then
   * free the data before proceeding.
   */

  if (entry_P->allocHere)
    FreeGval(entry_P->type, entry_P->gval);

  /*
   * Assign data.  We allocated it if it was a string or string array,
   * so set allocHere to TRUE.
   */

  entry_P->type = type;
  entry_P->gval = g;
  if ((type == GEN_TYPE_STRING) || (type == GEN_TYPE_STRING_A))
    entry_P->allocHere = TRUE;
  else
    entry_P->allocHere = FALSE;

  return(status);

}



/*
 * St_ReadSymFile--
 *	Read symbols from the stream until EOF.
 */
void St_ReadSymFile(FILE *input_F, st_Table_Pt st_P)
{

  int status;

  while ((status = St_ReadSym(input_F, st_P)) != ST_AT_EOF) 
    {
      if (status == ST_SYM_READ_FMT_ERR) 
	{
	  fprintf(stderr, "Warning: error parsing symbol in file.");
	}
    }

  return;

}




/*
 * St_CheckEmpty--
 *	Checks to see if a table is empty if not write out a list of entrys
 *
 * Return codes:
 *   0  if the table was empty
 *   1  if there is at least one entry in the table
 */
int St_CheckEmpty(FILE *output_F, st_Table_Pt st_P)
{
  st_Entry_Pt	entry_P;
  int		bucket;

  if (!st_P) return(0);
  if(!st_P->entryCt) return(0);
  fprintf(output_F, "\n\nUndefined symbols:\n");
  for (bucket = 0; bucket < st_P->size; bucket++)
    {
      entry_P = st_P->tab_A[bucket];
      while (entry_P)
	{
	  fprintf(output_F, "\t%s\n", entry_P->name);
	  entry_P = entry_P->next_P;
	}
    }
  return(1);
}



/*
 * St_CheckAllMarked--
 *	Checks to see if all entries in a table have been used.
 *
 * Return codes:
 *   0  if one or more entries were not used
 *   1  if all entries were used
 */
int St_CheckAllMarked(FILE *output_F, st_Table_Pt st_P)
{
  int             flag = 0, unused = 0, size, h;
  st_Entry_Pt	entry_P;

  if (!st_P) return(0);

  /* Go through entire hash table */
  size = st_P->size;
  for(h = 0; h < size; h++)
    {
      for (entry_P = st_P->tab_A[h]; entry_P; entry_P = entry_P->next_P)
	{
	  flag = entry_P->marker;
	  if(!flag)
	    {
	      fprintf(output_F,
		      "Warning: unused parameter: %s\n", entry_P->name); 
	      unused =1;
	    }
	}
    }

  if(unused)
    return(0);
  else
    return(1);
}

/*
 * St_WriteSymFile--
 *	Write the contents of a symbol table to the given stream.
 */
void St_WriteSymFile(FILE *output_F, st_Table_Pt st_P)
{

  st_Entry_Pt	entry_P;
  /*      int		bucket; */
  char	name[256], **symName_A, *chars_P;

  if (!st_P) {
    return;
  }

  St_SListTable(st_P, &symName_A, &chars_P);
  for (; *symName_A; symName_A++) {

    entry_P = LookupSym(st_P, *symName_A);

    if (!entry_P)
      continue;

    *name = 0;
    if (strchr(entry_P->name, ' ')) {
      strcpy(name, "\"");
      strcat(name, entry_P->name);
      strcat(name, "\"");
    } else {
      strcpy(name, entry_P->name);
    }

    switch(entry_P->type) {

    case GEN_TYPE_INT:
      fprintf(output_F, "%s int %d\n", name,
	      entry_P->gval.i);
      break;

    case GEN_TYPE_LONG:
      fprintf(output_F, "%s long %ld\n", name,
	      entry_P->gval.l);
      break;

    case GEN_TYPE_FLOAT:
      fprintf(output_F, "%s float %.6e\n", name,
	      entry_P->gval.f);
      break;

    case GEN_TYPE_DOUBLE:
      fprintf(output_F, "%s double %.12e\n", name,
	      entry_P->gval.d);
      break;

    case GEN_TYPE_CHAR:
      fprintf(output_F, "%s char %d\n", name,
	      entry_P->gval.c);
      break;

    case GEN_TYPE_STRING:
      if (entry_P->gval.s)
	fprintf(output_F, "%s string \"%s\"\n", 
		name, entry_P->gval.s);
      else
	fprintf(output_F, "%s string \"\"\n", name);
      break;

    case GEN_TYPE_STRING_A:

      if (!entry_P->gval.s_A) {
	fprintf(output_F, "%s stringarray \"\"\n", name);
      } else {
	char	**s_A, *s;

	s_A = entry_P->gval.s_A;

	if (!s_A) {
	  fprintf(output_F, "%s stringarray \"\"\n", name);
	  break;
	}

	fprintf(output_F, "%s stringarray ", name);
	for (; *s_A; s_A++) {
	  s = *s_A;
	  while (strlen(s) > 60) {
	    fprintf(output_F, "\"%60s\"\n", s);
	    s += 60;
	    if (!*s)
	      break;
	  }
	  fprintf(output_F, "\"%s\",\n", s);
	}
	fprintf(output_F, "\"\"\n");
      }
      break;

    case GEN_TYPE_VOID_PTR:
    default:
      break;

    }

  }

  St_SListTableFree( &symName_A, &chars_P);

  return;
}


/*
 * Coercion functions -- retrieve a symbol and coerce it to a given type
 * if possible.
 */
int St_GetSymAsInt(st_Table_Pt st_P, const char *name, int *i_P)
{

  int type;
  generic_t gval;

  if (St_GetSym(st_P, name, &type, &gval) == ST_SYM_NOT_FOUND)
    return(ST_SYM_NOT_FOUND);

  switch(type) {
  case GEN_TYPE_FLOAT:
    *i_P = gval.f;
    return(ST_SYM_FOUND);
  case GEN_TYPE_DOUBLE:
    *i_P = gval.d;
    return(ST_SYM_FOUND);
  case GEN_TYPE_INT: 
    *i_P = gval.i; 
    return(ST_SYM_FOUND);
  case GEN_TYPE_LONG:
    *i_P = gval.l;
    return(ST_SYM_FOUND);
  case GEN_TYPE_CHAR:
    *i_P = gval.c;
    return(ST_SYM_FOUND);
  case GEN_TYPE_STRING:
    *i_P = atoi(gval.s);
  default:
    return(ST_COULDNT_COERCE);

  }

}



int St_GetSymAsLong(st_Table_Pt st_P, const char *name, long *l_P)
{

  int type;
  generic_t gval;

  if (St_GetSym(st_P, name, &type, &gval) == ST_SYM_NOT_FOUND)
    return(ST_SYM_NOT_FOUND);

  switch(type) {
  case GEN_TYPE_FLOAT:
    if (coercionBits & ST_ALLOW_FLOAT_TO_INT) {
      *l_P = gval.f;
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  case GEN_TYPE_DOUBLE:
    if (coercionBits & ST_ALLOW_FLOAT_TO_INT) {
      *l_P = gval.d;
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  case GEN_TYPE_INT: 
    *l_P = gval.i; 
    return(ST_SYM_FOUND);
  case GEN_TYPE_LONG:
    *l_P = gval.l;
    return(ST_SYM_FOUND);
  case GEN_TYPE_CHAR:
    *l_P = gval.c;
    return(ST_SYM_FOUND);
  case GEN_TYPE_STRING:
    if (coercionBits & ST_ALLOW_STRING_TO_INT) {
      *l_P = atoi(gval.s);
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  default:
    return(ST_COULDNT_COERCE);

  }

}



int St_GetSymAsBoolean(st_Table_Pt st_P, const char *name, int *i_P)
{

  int type;
  generic_t gval;

  if (St_GetSym(st_P, name, &type, &gval) == ST_SYM_NOT_FOUND)
    return(ST_SYM_NOT_FOUND);

  switch(type) {
  case GEN_TYPE_FLOAT:
    if (coercionBits & ST_ALLOW_FLOAT_TO_INT) {
      *i_P = (gval.f != 0);
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  case GEN_TYPE_DOUBLE:
    if (coercionBits & ST_ALLOW_FLOAT_TO_INT) {
      *i_P = (gval.d != 0);
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  case GEN_TYPE_INT: 
    *i_P = (gval.i != 0); 
    return(ST_SYM_FOUND);
  case GEN_TYPE_LONG:
    *i_P = (gval.l != 0);
    return(ST_SYM_FOUND);
  case GEN_TYPE_CHAR:
    *i_P = (gval.c != 0);
    return(ST_SYM_FOUND);
  case GEN_TYPE_STRING:
    if (coercionBits & ST_ALLOW_STRING_TO_BOOLEAN) {
      if (atoi(gval.s) || *gval.s == 'T' || *gval.s == 't' ||
	  *gval.s == '+')
	*i_P = TRUE;
      else
	*i_P = FALSE;
    } else {
      return(ST_COULDNT_COERCE);
    }
  default:
    return(ST_COULDNT_COERCE);

  }

}



int St_GetSymAsChar(st_Table_Pt st_P, const char *name, int *i_P)
{

  int type;
  generic_t gval;

  if (St_GetSym(st_P, name, &type, &gval) == ST_SYM_NOT_FOUND)
    return(ST_SYM_NOT_FOUND);

  switch(type) {
  case GEN_TYPE_FLOAT:
    if (coercionBits & ST_ALLOW_FLOAT_TO_INT) {
      *i_P = (int) gval.f & 0xFF;
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  case GEN_TYPE_DOUBLE:
    if (coercionBits & ST_ALLOW_FLOAT_TO_INT) {
      *i_P = (int) gval.d & 0xFF;
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  case GEN_TYPE_INT: 
    *i_P = gval.i & 0xFF; 
    return(ST_SYM_FOUND);
  case GEN_TYPE_LONG:
    *i_P = gval.l & 0xFF;
    return(ST_SYM_FOUND);
  case GEN_TYPE_CHAR:
    *i_P = gval.c;
    return(ST_SYM_FOUND);
  case GEN_TYPE_STRING:
    if (coercionBits & ST_ALLOW_STRING_TO_CHAR) {
      *i_P = *gval.s;
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  default:
    return(ST_COULDNT_COERCE);

  }

}



int St_GetSymAsFloat(st_Table_Pt st_P, const char *name, float *f_P)
{

  int type;
  generic_t gval;

  if (St_GetSym(st_P, name, &type, &gval) == ST_SYM_NOT_FOUND)
    return(ST_SYM_NOT_FOUND);

  switch(type) {
  case GEN_TYPE_FLOAT:
    *f_P = gval.f;
    return(ST_SYM_FOUND);
  case GEN_TYPE_DOUBLE:
    *f_P = gval.d;
    return(ST_SYM_FOUND);
  case GEN_TYPE_INT: 
    *f_P = gval.i; 
    return(ST_SYM_FOUND);
  case GEN_TYPE_LONG:
    *f_P = gval.l;
    return(ST_SYM_FOUND);
  case GEN_TYPE_CHAR:
    *f_P = gval.c;
    return(ST_SYM_FOUND);
  case GEN_TYPE_STRING:
    if (coercionBits & ST_ALLOW_STRING_TO_FLOAT) {
      *f_P = atof(gval.s);
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  default:
    return(ST_COULDNT_COERCE);

  }

}



int St_GetSymAsDouble(st_Table_Pt st_P, const char *name, double *d_P)
{

  int type;
  generic_t gval;

  if (St_GetSym(st_P, name, &type, &gval) == ST_SYM_NOT_FOUND)
    return(ST_SYM_NOT_FOUND);

  switch(type) {
  case GEN_TYPE_FLOAT:
    *d_P = gval.f;
    return(ST_SYM_FOUND);
  case GEN_TYPE_DOUBLE:
    *d_P = gval.d;
    return(ST_SYM_FOUND);
  case GEN_TYPE_INT: 
    *d_P = gval.i; 
    return(ST_SYM_FOUND);
  case GEN_TYPE_LONG:
    *d_P = gval.l;
    return(ST_SYM_FOUND);
  case GEN_TYPE_CHAR:
    *d_P = gval.c;
    return(ST_SYM_FOUND);
  case GEN_TYPE_STRING:
    *d_P = atof(gval.s);
    return(ST_SYM_FOUND);
  default:
    return(ST_COULDNT_COERCE);

  }

}



int St_GetSymAsString(st_Table_Pt st_P, const char *name, char *s)
{

  int type;
  generic_t gval;

  if (St_GetSym(st_P, name, &type, &gval) == ST_SYM_NOT_FOUND)
    return(ST_SYM_NOT_FOUND);

  switch(type) {
  case GEN_TYPE_FLOAT:
    if (coercionBits & ST_ALLOW_FLOAT_TO_STRING) {
      sprintf(s, "%.6g", gval.f);
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  case GEN_TYPE_DOUBLE:
    if (coercionBits & ST_ALLOW_FLOAT_TO_STRING) {
      sprintf(s, "%.12g", gval.d);
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  case GEN_TYPE_INT: 
    if (coercionBits & ST_ALLOW_INT_TO_STRING) {
      sprintf(s, "%1d", gval.i);
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  case GEN_TYPE_LONG:
    if (coercionBits & ST_ALLOW_INT_TO_STRING) {
      sprintf(s, "%1d", (int)gval.l);
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  case GEN_TYPE_CHAR:
    if (coercionBits & ST_ALLOW_CHAR_TO_STRING) {
      sprintf(s, "%c", gval.c);
      return(ST_SYM_FOUND);
    } else {
      return(ST_COULDNT_COERCE);
    }
  case GEN_TYPE_STRING:
    strcpy(s, gval.s);
    return(ST_SYM_FOUND);
  default:
    return(ST_COULDNT_COERCE);

  }

}


/*
 *
 * St_GetSymAsVector
 * 
 * Returns the pointer to a double vector.
 *
 */
int St_GetSymAsVector(st_Table_Pt st_P, const char *name, double **v_P, int *length)
{
  int type;
  generic_t gval;

  if (St_GetSym(st_P, name, &type, &gval) == ST_SYM_NOT_FOUND)
    return(ST_SYM_NOT_FOUND);

  switch(type) {
  case GEN_TYPE_FLOAT:
    *v_P = (double *) &(gval.f);
    *length = 1;
    return(ST_SYM_FOUND);
  case GEN_TYPE_DOUBLE:
    *v_P = (double *) &(gval.d);
    *length = 1;
    return(ST_SYM_FOUND);
  case GEN_TYPE_INT: 
    *v_P = (double *) &(gval.i);
    *length = 1;
    return(ST_SYM_FOUND);
  case GEN_TYPE_LONG:
    *v_P = (double *) &(gval.l);
    *length = 1;
    return(ST_SYM_FOUND);
  case GEN_TYPE_DOUBLE_VECTOR:
    *v_P = (double *) gval.dv.v;
    *length = gval.dv.size;
    return(ST_SYM_FOUND);
  default:
    return(ST_COULDNT_COERCE);
  }
}

/*
 *
 * St_GetCount
 * 
 * Returns the number of entries in a table.
 *
 */
int St_GetCount(st_Table_Pt st_P)
{
  return(st_P->entryCt);
}


/*
 * St_DelAllSym--
 *	Delete all the entries in a symbol table but not the table itself.
 */
void St_DelAllSym(st_Table_Pt st_P)
{
  st_Entry_Pt entry_P, next_P;
  int		i;

  if (!st_P) return;

  for (i = 0; i < st_P->size; i++)
    {
      for (entry_P = st_P->tab_A[i]; entry_P; entry_P = next_P)
	{
	  next_P = entry_P->next_P;
	  free(entry_P);
	}
      st_P->tab_A[i] = NULL;
    }
  return;
}

/*
 * strcmpCI
 *	Case insensitive strcmp
 * returns 0 if strings are identical
 * returns 1 if strings are different
 */
int strcmpCI(char *s1, char *s2)
{
  int c, d;
  while( (c=tolower(*s1++)) )
    {
      d = tolower(*s2++);
      if(c != d) return(1);
    }
  if(tolower(*s2)) return(1);
  return(0);
}

/*
 * strcpyCareFul
 *
 * Careful version of strcpy
 * Allocates space for s1 (must supply pointer for destination string.
 * copies s2 to *s1
 */
void strcpyCareFul(char **s1, const char *s2)
{
  *s1 = (char *) malloc(strlen(s2)+1);
  if(*s1 == NULL)
    {
      fprintf(stderr,"Error: Out of memory ... exiting");
      exit(1);
    }
  strcpy(*s1,s2);
  return;
}

