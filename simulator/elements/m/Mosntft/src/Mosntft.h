// MOSFET LEVEL 16 MODEL
//
//          Drain 2
//                  o
//                  |
//                  |
//              |---+
//              |
// Gate 1 o-----|
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
//
// - Model based on parameters extracted by HP.

#ifndef Mosntft_h
#define Mosntft_h 1

#include "../../../../network/CircuitManager.h"
#include "../../../../network/ADInterface.h"

class Mosntft : public ADInterface
{
	public:

  Mosntft(const string& iname);

  ~Mosntft() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:
  virtual void eval(AD * x, AD * effort, AD * flow);

  // Some constants
  double leff, weff, vt, cfm, tfm;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double L, LD, W, WD;
  double U0, VT0, PHI, NFS, NSS, T1, T2, E1, E2, THETA, ETA;
  double VMAX, G0, DEFF, VX, TVST, PSI, GAMMA, VTIME, TREF, T, RD, RS;
  double CGS0, CGD0, CSC, FREQ, FEFF, TAU,NU, CHI,K2;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
