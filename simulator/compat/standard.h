/* standard.h */

typedef void *VOIDPTR;

typedef const void *CONST_VOIDPTR;

#ifndef PI
#define PI		3.14159265359
#endif

#define TWOPI		6.28318530718
#define RAD2DEG		57.2957795131
#define TRUE		1
#define FALSE		0
#define REALLYSMALL	1.0e-15
#define MAX_NAME_LEN	1024
#define MAX_STRING_LEN	1024
#define HASH_SIZE 	29      /* Size of various hash tables.
                                 * Best if it is a prime number and the       
                                 * larger it is the better for speed.  */

extern char   outputFilename[MAX_STRING_LEN];  /* name of output file */
/*
 * ARGUMENT LABELS
 *
 * These are used to label function arguments.
 *
 * IN indicates that the argument data is read only by the function, and 
 * not written.
 *
 * INOUT indicates that the argument data is read and written by the
 * function.
 *
 * OUT indicates that the argument data is written by the function,
 * without reference to its previous value.
 *
 * INOUT and OUT apply only to pointer arguments, since there is no
 * direct pass-by-reference mechanism in C.
 *
 * A fourth keyword, SHARED, is also defined.  It indicates that the
 * argument is a pointer to shared data (i. e., shared between PROCESSES
 * or simultaneously-running coroutines).  Data can be read or written,
 * but not meaningfully copied, since the value of the data pointed to
 * may CHANGE at any time.  A good example would be an address containing
 * a time.
 *
 * Example:
 *
 * char CharFunc(    IN readOnlyArg,
 *		    OUT writeOnlyArg, writeOnlyArg2,
 *		  INOUT readWriteArg,
 *		 SHARED theTime)
 */

#define	IN
#define OUT
#define INOUT
#define SHARED
