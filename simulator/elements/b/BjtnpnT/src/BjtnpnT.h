//		  Spice BJT model with charge conservation model
//
//
//
//                  0    NCollector
//                  |
//                  |
//	         | /
//               |/
//   NBase  0----|
//               |
//               |
//                  |
//                  |
//                  0    NEmitter
//
//	  Author:
//		 Senthil  Velu

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

#ifndef BjtnpnT_h
#define BjtnpnT_h 1

class BjtnpnT : public ADInterface
{
	public:

  BjtnpnT(const string& iname);

  ~BjtnpnT() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list);

	private:

  // Generic state variable evaluation routine
	virtual void eval(AD * x, AD * effort, AD * flow);

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double is, bf, nf, vaf, ikf, ise, ne, br, nr, var, ikr, isc,
	nc, re, rb, rbm, irb, rc, eg, cje, vje, mje, cjc, vjc, mjc,
	xcjc, fc, tf, xtf, vtf, itf, tr, xtb, xti, tre1, tre2, trb1, trb2,
	trm1, trm2, trc1, trc2, tnom, b, ts, cjs, mjs, vjs, area, ns, iss;
  double bm1;
  bool kirchhoff;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif

