//
//
//   UCSD HBT model
//
//
//
//                  0 Collector
//                  |
//                  |
//	               |
//               ||
//    Base  0----|
//               |    --------0 Power dissipation
//               ||   --------0 Thermal ground
//                 |
//                  |
//                  |
//                  0 Emitter
//
//
//
//	Author: Jian Ding, Sonali Luniya
//


#ifndef HBTnpnxT_h
#define HBTnpnxT_h 1

#include "../../../../network/ADInterface.h"
#include "../../../../analysis/TimeDomainSV.h"

class HBTnpnxT : public ADInterface
{
  public:

  HBTnpnxT(const string& iname);

  ~HBTnpnxT() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
      TerminalVector& term_list);

  private:
  // Generic state variable evaluation routine
  virtual void eval(AD * x, AD * effort, AD * flow);

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double type, afn, bf, bfn, bkdn, br, bvc, ccmin, cemin, cjc, cjcx, cje;
  double cjs, cth, cxmin, dtmax, eaa, eab, eac, eae, eax, eg, fa, fc, fce;
  double fex, icrit0, ics, ik, ikrk, is, isa, isb, isc, iscx, ise, isex, itc;
  double itc2, kfn, mjc, mjcx, mje, mjs, na, nb, nbc, nc, ncs, ncx, ne, nex;
  double nf, nr, rbi, rbx, rci, rcx, re, rex, rth, tbcxs, tbexs, tfb, tfc0;
  double tkrk, tnc, tne, tnex, tnom, tr, trx, tre, tvjc, tvjci, tvjcx, tvje;
  double tvjs, vaf, var, vjc, vjci, vjcx, vje, vjs, vkrk, vtc, xcjc, xrb;
  double xrc, xre, xrex, xrt, xtb, xti, xtikrk, xtitc, xtitc2, xttf, xttkrk;
  double xtvkrk, b, ts, bm1;
  bool kirchhoff;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
