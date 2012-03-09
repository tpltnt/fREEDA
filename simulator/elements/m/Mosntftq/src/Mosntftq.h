// MOSFET LEVEL 16 MODEL
//
//          Drain 2
//                  o
//                  |
//                  |
//              |---+
//              |
// Gate 1 o-----|------o 4 Bulk
//              |
//              |---+
//                  |
//                  |
//                  o
//         Source 3
//
//
//	  Author:
//		718 Project

#ifndef Mosntftq_h
#define Mosntftq_h 1

#include "../../../../network/ADInterface.h"
#include "../../../../analysis/TimeDomainSV.h"

class Mosntftq : public ADInterface
{
	public:

  Mosntftq(const string& iname);

  ~Mosntftq() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:
  virtual void eval(AD * x, AD * effort, AD * flow);

  // Some constants
  double  vt,A,D;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double l,w,u,vto,tts,ctr,rd,vthm,beta,roff,cgso,cgdo,cl,lambda,rs,tnom;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
