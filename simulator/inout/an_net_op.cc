// Network-related output commands

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

extern "C" 
{
	#include "../compat/stuff.h"
	#include "oper.h"
	#include "../compat/dsp.h"
	#include "../compat/ml.h"
	#include "oper_net.h"
	#include "ftvec.h"
}

#include "../network/CircuitManager.h"
#include <assert.h>

// Defined in freeda.cc
// This variable is required by the output routines.
extern Circuit* main_cir;

// An_Op_VT--
//	return the transient voltage vector at a terminal in the network
int An_Op_VT(capOutReq_t *arg, capOutReq_t *result)
{
  assert(TimeV_P);
  
  if (arg->type != CAP_OBJ_TERM) 
	{
    An_CalcError("Operator VT takes a terminal as an argument");
    return -1;
  }
	
  unsigned id = arg->obj.nid.val;
  Terminal* term = main_cir->getTerminal(id);
	
  const DenseDoubleVector& v_v = term->getTermData()->getRealV();
	
  if (!v_v.length()) 
	{
    An_CalcError("Data not available.");
    return -1;
  }
	
  assert(v_v.length() == unsigned(NoTimePoints));
	
  // Alloc memory
  doublev_t dv, xv;
  dv = Mlib_DNewVec(NoTimePoints);
  xv = Mlib_DNewVec(NoTimePoints);
	
  // Fill output vectors
  for (int i = 0; i < NoTimePoints; i++) 
	{
    dv[i] = v_v[i];
    xv[i] = TimeV_P[i];
  }
	
  result->obj.dv = dv;
  result->size = NoTimePoints;
  result->x = xv;
  result->type = CAP_OBJ_DOUBLEV;
  
  return 0;
}

/*
* An_Op_IT--
* 	return the transient current vector of a particular terminal
*      of an element in the network
*/
int An_Op_IT(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  assert(TimeV_P);
  
  if (arg1->type != CAP_OBJ_NODE) 
	{
    An_CalcError("Operator IT takes a element as argument 1");
    return -1;
  }
  
  if (arg2->type != CAP_OBJ_INT) 
	{
    An_CalcError("Operator IT takes an integer as second argument");
    return(-1);
  }
	
  unsigned id = arg1->obj.nid.val;
  Element* elem = main_cir->getElement(id);
  unsigned terminalIndex = arg2->obj.i;
	
  if (terminalIndex > elem->getNumTerms()) 
	{
    An_CalcError("Index greater than number of terminals");
    return(-1);
  }
  
  const DenseDoubleVector& i_v = elem->getElemData()->getRealI(terminalIndex);
	
  if (!i_v.length()) 
	{
    An_CalcError("Data not available.");
    return -1;
  }
	
  assert(i_v.length() == unsigned(NoTimePoints));
	
  // Alloc memory
  doublev_t dv, xv;
  dv = Mlib_DNewVec(NoTimePoints);
  xv = Mlib_DNewVec(NoTimePoints);
	
  // Fill output vectors
  for (int i = 0; i < NoTimePoints; i++) 
	{
    dv[i] = i_v[i];
    xv[i] = TimeV_P[i];
  }
	
  result->obj.dv = dv;
  result->size = NoTimePoints;
  result->x = xv;
  result->type = CAP_OBJ_DOUBLEV;
  
  return 0;
}

/*
* 	return the double transient voltage vector of a particular
*      port of an element in the network
*/
int An_Op_UT(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  assert(TimeV_P);
	
  if (arg1->type != CAP_OBJ_NODE)
    An_CalcError("Operator UT takes an element as argument 1");
  
  if (arg2->type != CAP_OBJ_INT) 
	{
    An_CalcError("Operator UT takes an integer as second argument");
    return(-1);
  }
	
  unsigned id = arg1->obj.nid.val;
  Element* elem = main_cir->getElement(id);
  unsigned terminalIndex = arg2->obj.i;
  
  if (terminalIndex > elem->getNumberOfStates()) 
	{
    An_CalcError("Index greater than number of states");
    return(-1);
  }
  
  const DenseDoubleVector& u_v = elem->getElemData()->getRealU(terminalIndex);
	
  if (!u_v.length()) 
	{
    An_CalcError("Data not available.");
    return -1;
  }
	
  assert(u_v.length() == unsigned(NoTimePoints));
	
  // Alloc memory
  doublev_t dv, xv;
  dv = Mlib_DNewVec(NoTimePoints);
  xv = Mlib_DNewVec(NoTimePoints);
	
  // Fill output vectors
  for (int i = 0; i < NoTimePoints; i++) 
	{
    dv[i] = u_v[i];
    xv[i] = TimeV_P[i];
  }
	
  result->obj.dv = dv;
  result->size = NoTimePoints;
  result->x = xv;
  result->type = CAP_OBJ_DOUBLEV;
	
  return 0;
}

/*
* 	return the double state variable vector of a particular
*      element in the network
*/
int An_Op_XT(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  assert(TimeV_P);
	
  if (arg1->type != CAP_OBJ_NODE)
    An_CalcError("Operator UT takes an element as argument 1");
  
  if (arg2->type != CAP_OBJ_INT) 
	{
    An_CalcError("Operator UT takes an integer as second argument");
    return(-1);
  }
	
  unsigned id = arg1->obj.nid.val;
  Element* elem = main_cir->getElement(id);
  unsigned terminalIndex = arg2->obj.i;
  
  if (terminalIndex > elem->getNumberOfStates()) 
	{
    An_CalcError("Index greater than number of states");
    return(-1);
  }
  
  const DenseDoubleVector& x_v = elem->getElemData()->getRealX(terminalIndex);
	
  if (!x_v.length()) 
	{
    An_CalcError("Data not available.");
    return -1;
  }
	
  assert(x_v.length() == unsigned(NoTimePoints));
	
  // Alloc memory
  doublev_t dv, xv;
  dv = Mlib_DNewVec(NoTimePoints);
  xv = Mlib_DNewVec(NoTimePoints);
	
  // Fill output vectors
  for (int i = 0; i < NoTimePoints; i++) 
	{
    dv[i] = x_v[i];
    xv[i] = TimeV_P[i];
  }
	
  result->obj.dv = dv;
  result->size = NoTimePoints;
  result->x = xv;
  result->type = CAP_OBJ_DOUBLEV;
	
  return 0;
}

/*
* 
* Return the complex voltage vector at a terminal in the network
*
*/
int An_Op_VF(capOutReq_t *arg, capOutReq_t *result)
{
  assert(FreqV_P);
  
  if (arg->type != CAP_OBJ_TERM) 
	{
    An_CalcError("Operator VF takes a terminal as an argument");
    return -1;
  }
	
  unsigned id = arg->obj.nid.val;
  Terminal* term = main_cir->getTerminal(id);
	
  const DenseComplexVector& v_v = term->getTermData()->getPhasorV();
	
  if (!v_v.length()) 
	{
    An_CalcError("Data not available.");
    return -1;
  }
	
  assert(v_v.length() == unsigned(NoFreqPoints));
	
  // Alloc memory
  doublev_t xv;
  dcxv_t dcxv;
	
  dcxv = Mlib_CNewVec(NoFreqPoints);
  xv = Mlib_DNewVec(NoFreqPoints);
	
  // Fill output vectors
  for (int i = 0; i < NoFreqPoints; i++) 
	{
    dcxv[i].re = v_v[i].real();
    dcxv[i].im = v_v[i].imag();
    xv[i] = FreqV_P[i];
  }
	
  result->obj.dcxv = dcxv;
  result->size = NoFreqPoints;
  result->x = xv;
  result->type = CAP_OBJ_DCXV;
  
  return 0;
}

/*
* 
*	return the complex current vector of a particular terminal
*      of an element in the network
*/
int An_Op_IF(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  assert(FreqV_P);
  if (arg1->type != CAP_OBJ_NODE) 
	{
    An_CalcError("Operator IF takes a element as argument 1");
    return -1;
  }
  
  if (arg2->type != CAP_OBJ_INT) 
	{
    An_CalcError("Operator IF takes an integer as second argument");
    return(-1);
  }
	
  unsigned id = arg1->obj.nid.val;
  Element* elem = main_cir->getElement(id);
  unsigned terminalIndex = arg2->obj.i;
	
  if (terminalIndex > unsigned(elem->getNumTerms())) 
	{
    An_CalcError("Index greater than number of terminals");
    return(-1);
  }
  
  const DenseComplexVector& i_v = elem->getElemData()->getPhasorI(terminalIndex);
	
  if (!i_v.length()) 
	{
    An_CalcError("Data not available.");
    return -1;
  }
	
  assert(i_v.length() == unsigned(NoFreqPoints));
	
  // Alloc memory
  doublev_t xv;
  dcxv_t dcxv;
	
  dcxv = Mlib_CNewVec(NoFreqPoints);
  xv = Mlib_DNewVec(NoFreqPoints);
	
  // Fill output vectors
  for (int i = 0; i < NoFreqPoints; i++) 
	{
    dcxv[i].re = i_v[i].real();
    dcxv[i].im = i_v[i].imag();
    xv[i] = FreqV_P[i];
  }
	
  result->obj.dcxv = dcxv;
  result->size = NoFreqPoints;
  result->x = xv;
  result->type = CAP_OBJ_DCXV;
	
  return 0;
}


/*
* 
*	return the complex state variable vector of a particular terminal
*      of an element in the network
*/
int An_Op_XF(capOutReq_t *arg1, capOutReq_t *arg2, capOutReq_t *result)
{
  assert(FreqV_P);
	
  if (arg1->type != CAP_OBJ_NODE)
    An_CalcError("Operator XF takes an element as argument 1");
  
  if (arg2->type != CAP_OBJ_INT) 
	{
    An_CalcError("Operator XF takes an integer as second argument");
    return(-1);
  }
	
  unsigned id = arg1->obj.nid.val;
  Element* elem = main_cir->getElement(id);
  unsigned terminalIndex = arg2->obj.i;
  
  if (terminalIndex > elem->getNumberOfStates()) 
	{
    An_CalcError("Index greater than number of states");
    return(-1);
  }
  
  const DenseComplexVector& x_v = elem->getElemData()->getPhasorX(terminalIndex);
	
  if (!x_v.length()) 
	{
    An_CalcError("Data not available.");
    return -1;
  }
	
  assert(x_v.length() == unsigned(NoFreqPoints));
	
  // Alloc memory
  doublev_t xv;
  dcxv_t dcxv;
	
  dcxv = Mlib_CNewVec(NoFreqPoints);
  xv = Mlib_DNewVec(NoFreqPoints);
	
  // Fill output vectors
  for (int i = 0; i < NoFreqPoints; i++) 
	{
    dcxv[i].re = x_v[i].real();
    dcxv[i].im = x_v[i].imag();
    xv[i] = FreqV_P[i];
  }
	
  result->obj.dcxv = dcxv;
  result->size = NoFreqPoints;
  result->x = xv;
  result->type = CAP_OBJ_DCXV;
	
  return 0;
}

