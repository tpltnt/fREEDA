// PHILIPS MOS9 --CapacitorJn MODEL

//	      |\ |

//   0 o------| >|------o 1

//	      |/ |

// anode                 cathode

//	  Authors:
//		Yogesh Ramdoss, Kuldip Gothi, Xuemin Yang,
//    Ajit Rajagopalan, Dapeng ding, Xin Cai

#ifndef CapacitorJn_h
#define CapacitorJn_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class CapacitorJn : public ADInterface

{
	public:

  CapacitorJn(const string& iname);
  ~CapacitorJn() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  // Generic state variable evaluation routine
  virtual void eval(AD * x, AD * effort, AD * flow);

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double level;
  double ab;
  double ls;
  double lg;
  double dta;
  double tr;
  double vr;
  double jsgbr;
  double jsdbr;
  double jsgsr;
  double jsdsr;
  double jsggr;
  double jsdgr;
  double nb;
  double ns;
  double ng;
  double vb;
  double cjbr;
  double cjsr;
  double cjgr;
  double vdbr;
  double vdsr;
  double vdgr;
  double pb;
  double ps;
  double pg;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
