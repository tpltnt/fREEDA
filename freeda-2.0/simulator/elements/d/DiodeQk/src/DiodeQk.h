//
//	      |\ |
//   0 o------| >|------o 1
//	      |/ |
// anode                 cathode
//
//
//	  Author:
//		James Murray
//

#ifndef DiodeQk_h
#define DiodeQk_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class DiodeQk : public ADInterface
{
	public:

 DiodeQk(const string& iname);

  ~DiodeQk() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

	private:

  virtual void eval(AD * x, AD * effort, AD * flow);

  // Some constants
  double Vo;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double js;
  double alfa;
  double r0;
  double area;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
