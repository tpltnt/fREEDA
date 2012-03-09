#ifndef parserGlobals_h
#define parserGlobals_h

#include "../compat/graphID.h"

/* termCount  = current # of terminals
 *Use to track number of terminals of an element
 *Of use in initAnalysistializing elements with variable number of
 *terminals. */
extern int thisTermCount;   /* current # of terminals */

/* thisTermList[i].val is the index of the i th terminal. i = 0, 1, 2 ... 
 * Can only use when constructing an element as this is set in the parser. */
extern gr_Id_t thisTermList[];  /* current list of terminals    */

/***************************
 * Function prototypes 
 ***************************/

/* defined in parser_functions.cc */

/* Use to report parser errors. Routine locates source of error */
void ErrMsg(const char *s);

#endif

