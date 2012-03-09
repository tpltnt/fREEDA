//  Ideal mixer
//
//  Author:
//  Mark Buff


#ifndef AbmMixer_h
#define AbmMixer_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class AbmMixer : public ADInterface
{
	public:

  AbmMixer(const string& iname);

  ~AbmMixer() {}

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
  int op;
  int inputs;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
