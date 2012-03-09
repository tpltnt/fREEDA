// MOSFET EKV MODEL
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
//	  Authors:
//		William Cornwell
//		Mohamed Guled
//		Jordan Novgrod
//		Mac Perry
//		Ryon Stewart

#ifndef Mosnekv2_h
#define Mosnekv2_h 1
#define NMOS 1  //
#define PMOS -1 //

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class Mosnekv2 : public ADInterface
{
	public:

  Mosnekv2(const string& iname);

  ~Mosnekv2() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Constants
  double q, k, tref,ca;

  // Basic definitions
  double epsilonsi, epsilonox, tnom, temp;

  //Parameters
  // Device input variables
  double l, w, np, ns;
  // EKV intrinsic model process related parameters
  double cox, xj, dw, dl;
  // EKV basic intrinsic model parameters
  double vto, gamma, phi, kp, e0, ucrit;
  // Optional EKV intrinsic model parameters
  double tox, nsub, vfb, uo, vmax, theta;
  // Channel length modulation and charge sharing parameters
  double lambda, weta, leta;
  // Reverse short-channel effect parameters
  double q0,lk;
  // Impact ionization related parameters
  double iba, ibb, ibn;
  // Intrinsic model temperature parameters
  double tcv, bex, ucex, ibbt;
  // Matching parameters
  double avto, akp, agamma;
  // Flicker noise parameters
  double kf, af;
  // Setubp parameters
  double nqs, satlim, xqc;

  //starting here
  double scale, weff, leff, vtoa, kpa, gammaa, ce, xi, deltavrsce;

  // Parameter information
  static ParmInfo pinfo[];

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;
};

#endif
