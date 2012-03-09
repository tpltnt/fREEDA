/****************************************************************************
 * Miscellaneous declarations and definitions, many from old cap.h
  ****************************************************************************/

#ifndef graphID_h
#define graphID_h 1

/* IDENTIFIER TYPE */
typedef struct {
  unsigned val;
  short type;
} gr_Id_t, *gr_Id_Pt;

enum GR_Type {GR_ID_NO_VAL=0, GR_ID_NO_TYPE=0, GR_ID_TERM_TYPE,
  GR_ID_NODE_TYPE, GR_ID_CKTID_TYPE, 
  GR_ID_SUBCKT_TYPE, GR_ID_TOP_TYPE
};

#endif

