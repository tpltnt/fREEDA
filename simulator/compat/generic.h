/*
 * generic.h
 *
 * useful generic union type, used here and there
 *
 * Author:
 * Joseph Nathan Hall
 * Michael Steer
 *
 */


typedef struct {
  int	size;
  int	*v;     /* v is ant vector of length size */
} genericIntV_t, *genericIntV_Pt;

typedef struct {
  int	size;
  double	*v;     /* v is a double vector of length size */
} genericDoubleV_t, *genericDoubleV_Pt;

typedef struct {
  int	size;
  /* dv is a vector of pointers of length size */
  genericDoubleV_Pt dv_P;  
} genericDoubleVectors_t, *genericDoubleVectors_Pt;

typedef struct {
  int	size;
  char	**v;     /* v is a vector of length size of strings */
} genericStringV_t, *genericStringV_Pt;

typedef union {
  int	i;
  long	l;
  float	f;
  double	d;
  char	c;
  char	*s;
  char	**s_A;
  dcx_t	x;
  VOIDPTR	v;
  genericIntV_t iv;
  genericDoubleV_t dv;
  genericDoubleVectors_t dvv;
  genericStringV_t sv;
} generic_t, *generic_Pt;

#define GEN_TYPE_NONE 0
#define GEN_TYPE_INT 1
#define GEN_TYPE_LONG 2
#define GEN_TYPE_FLOAT 3
#define GEN_TYPE_DOUBLE 4
#define GEN_TYPE_CHAR 5
#define GEN_TYPE_STRING 6
#define GEN_TYPE_STRING_A 7
#define GEN_TYPE_COMPLEX 14
#define GEN_TYPE_VOID_PTR 8
#define GEN_TYPE_INT_VECTOR 9	    /* type for an genericInt_t vector */
#define GEN_TYPE_DOUBLE_VECTOR 10   /* type for an genericDoubleV_t vector*/
#define GEN_TYPE_STRING_VECTOR 11   /* type for an genericStringV_t vector*/
#define GEN_TYPE_BOOLEAN 12	    /* Noty a generic type, included here
				     * for parameter type purposes.       */
#define GEN_TYPE_ENUM 13	    /* Noty a generic type, included here
				     * for parameter type purposes.       */
#define GEN_TYPE_VARIABLE_USAGE 14  /* This is a variable usage. The
				     * VOIDPTR is used to point to
				     * st_Entry_t for the information on
				     * this variable. */
#define GEN_TYPE_SYMBOL_TABLE 15    /* Symbol table pointer */
#define GEN_TYPE_DOUBLE_VECTOR_OF_VECTORS 16 /* type for an 
					      * genericDoubleVectors_t vector*/


