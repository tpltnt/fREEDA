// Physical resistor model, subtype p
//
//                                     __
//                                      /|
//          Terminal 0                 /             Terminal 1
//                  o--------------/\/\/\--------------o
//                         |         /           |
//                         |        /            |
//                       -----                 -----
//                         ^                     ^
//                        ---                   ---
//                         |                     |
//                         |                     |
//                         -----------------------
//                                    |
//                                    |
//                                    o
//                          Terminal 2 (substrate)
//
//
//	Author:  NCSU ECE718 student
//
//	Model is based on description of Cadence physical resistor model,
//	  "Resistor models in Cadence Spectre", found at
//	  http://www.uta.edu/ronc/cadence/Resistor Models.pdf
//
//	Code for the diodes was taken from the SPDiode model by Carlos E. Christoffersen
//

#ifndef ResistorPhyP_h
#define ResistorPhyP_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class ResistorPhyP : public ADInterface
{
	public:

  ResistorPhyP(const string& iname);

  ~ResistorPhyP() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double r, coeff0, coeff1, coeff2, coeff3, coeff4, coeff5;
  bool polyarg;
  double tc1, tc2, tnom, tdev;
  double is, n, ibv, bv, fc, cj0, vj, m, tt, area, rs;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
