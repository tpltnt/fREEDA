// This may look like C code, but it is really -*- C++ -*-
//
// This is an ideal bi-directional circuit delay element.
//
//                 1                        2
//         I1 --->  o--------.---.     .---.----------o  <--- I2
//                           |   |     |   |
//           	   +  	     \  /^\   /^\  \         +
//           	  V1	     /   |     |   /         V2
//           	   -	     \  \ /   \ /  \         -
//                           |   |     |   |
//           	    o--------'---'     '---'----------o
//                 3                        4
//
// by Chris Saunders

#ifndef Delay_h
#define Delay_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class Delay : public ADInterface
{
public:

  Delay(const string& iname);

  ~Delay() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // This is a multi-reference element.
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list);
  // fill MNAM
  //  virtual void fillMNAM(FreqMNAM* mnam);
  //  virtual void fillMNAM(TimeMNAM* mnam);
  //  virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const ;
  //  virtual unsigned getExtraRC(const unsigned& eqn_number, const MNAMType& type) ;
  virtual void eval(AD * x, AD * effort, AD * flow);
  virtual void init() throw(string&) ;

  private:

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double delay;
  double res;

  unsigned temp_eval_counter;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
