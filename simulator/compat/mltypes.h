/*
 * mltypes.h--
 *	type definitions for ml module
 */

#ifndef mltypes_h
#define mltypes_h 1

typedef struct {
  float	re, im;
} cx_t, *cx_Pt;

typedef struct {
  double	re, im;
} dcx_t, *dcx_Pt;

typedef float *floatv_t, **floatm_t, ***floatmv_t, ****floatmm_t;
typedef double *doublev_t, **doublem_t, ***doublemv_t, ****doublemm_t;
typedef cx_t *cxv_t, **cxm_t, ***cxmv_t;
typedef dcx_t *dcxv_t, **dcxm_t, ***dcxmv_t,
  ***dcxsmv_t, ****dcxmm_t;   /* sprase matrix type */
typedef int *integerv_t, **integerm_t, ***integermv_t, ****integermm_t;
typedef char *charv_t, **charm_t, ***charmv_t, ****charmm_t;


/*
 *	Double precision and Complex matrix linked list structures.
 */
typedef struct gurnt temp_dstruct;
typedef struct gurnt {
  int	size;
  short	matInUse,
    vecInUse;
  double 	 **mat;
  double   *vec;
  temp_dstruct *next,
    *prior;
}delement;

typedef struct zweedot temp_cstruct;
typedef struct zweedot {
  int	size;
  short	matInUse,
    vecInUse;
  dcx_t	**mat;
  dcx_t    *vec;
  temp_cstruct *next,
    *prior;
}celement;

typedef struct thingy temp_istruct;
typedef struct thingy {
  int	size;
  int	matInUse,
    vecInUse;
  int		**mat;
  int		*vec;
  temp_istruct *next,
    *prior;
}ielement;


/*
 *	Integer linked list structure.
 */
typedef struct int_link_list_struct {
  int	x;
  struct int_link_list_struct *next,
    *prior;
}int_link_list_t, *int_link_list_Pt;

#endif
