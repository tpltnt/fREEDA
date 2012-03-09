// Physical resistor model, subtype poly
//
//                                     __
//                                      /|
//          Terminal 0                 /             Terminal 1
//                  o--------------/\/\/\--------------o
//                         |         /           |
//                         |        /            |
//                       -----                 -----
//                       -----                 -----
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

#ifndef ResistorPhyPoly_h
#define ResistorPhyPoly_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class ResistorPhyPoly : public ADInterface
{
	public:

  ResistorPhyPoly(const string& iname);

  ~ResistorPhyPoly() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Some constants -- not using this one now
  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double r, coeff0, coeff1, coeff2, coeff3, coeff4, coeff5;
  bool polyarg;
  double tc1, tc2, tc1c, tc2c, tnom, tdev, c;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
